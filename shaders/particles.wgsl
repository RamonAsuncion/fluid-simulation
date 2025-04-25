/* 
 * Copyright (c) 2025 SingChun LEE @ Bucknell University. CC BY-NC 4.0.
 * 
 * This code is provided mainly for educational purposes at Bucknell University.
 *
 * This code is licensed under the Creative Commons Attribution-NonCommerical 4.0
 * International License. To view a copy of the license, visit 
 *   https://creativecommons.org/licenses/by-nc/4.0/
 * or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
 *
 * You are free to:
 *  - Share: copy and redistribute the material in any medium or format.
 *  - Adapt: remix, transform, and build upon the material.
 *
 * Under the following terms:
 *  - Attribution: You must give appropriate credit, provide a link to the license,
 *                 and indicate if changes where made.
 *  - NonCommerical: You may not use the material for commerical purposes.
 *  - No additional restrictions: You may not apply legal terms or technological 
 *                                measures that legally restrict others from doing
 *                                anything the license permits.
 */
struct Particle {
  pos : vec2f,
  initPos : vec2f,
  vel : vec2f,
  initVel : vec2f,
  lifeTime : vec2f,
};
//boundary box for the fluid
const left_bound = -1;
const right_bound = 1;
const top_bound = 1;
const bot_bound = -1;
//grid parameters
const cell_size = .08;
const grid_length = i32(ceil((right_bound - left_bound)/cell_size));
const grid_height = i32(ceil((top_bound - bot_bound)/cell_size));
const max_per_cell = 512;
const use_acceleration = 1; //non-zero to use acceleration structures
//fluid simulation parameters
const mass = 1;
const smoothing_radius = cell_size; //this way we only need to check adjacent grid cells, MUST change code to change this value
const smoothing_rate = 4;
const pi = 3.14159265;
const e = 2.71828;
const func_volume = (smoothing_rate + 1)/(2 * pi * smoothing_radius);
const target_density = .015;
const pressure_multiplier= 10;
const gravity_multiplier = 0.2;
const viscosity_multiplier = 1;
const max_vel = .2;
const max_force = .02;
const velocity_damping = .8; //multiplies final computed velocity, between 0.5 and 1 prolly
const steps_per_update = .4; //lower value = slower/more stable simulation, somewhere between .1 and 1


// TODO 4: Write the bind group spells here using array<Particle>
@group(0) @binding(0) var<storage> particlesIn: array<Particle>;
@group(0) @binding(1) var<storage, read_write> particlesOut: array<Particle>;
@group(0) @binding(2) var<storage, read_write> timeBuffer: array<u32>;
@group(0) @binding(3) var<storage> gridIn: array<i32>;
@group(0) @binding(4) var<storage, read_write> gridOut: array<atomic<i32>>;
@group(0) @binding(5) var<storage, read_write> densityBuffer: array<f32>;

// currently, in 2D, for each particles (instances), we draw multiple vertices
// and for each vertex, we draw the sample on the circle using line-strips

// in 3D, for each particles (instances), we draw multiple vertices
// and for each vertex, we draw the sample on the sphere using triangle-lists

@vertex
fn vertexMain(@builtin(instance_index) idx: u32, @builtin(vertex_index) vIdx: u32) -> @builtin(position) vec4f {
  let particle = particlesIn[idx];
  let size = 0.0250 / 2f;
  let theta = 2. * pi / 8 * f32(vIdx);
  let x = cos(theta) * size;
  let y = sin(theta) * size;
  return vec4f(vec2f(x + particle.pos.x, y + particle.pos.y), 0, 1);
}

@fragment
fn fragmentMain() -> @location(0) vec4f {
  // return vec4f(238.f/255, 118.f/255, 35.f/255, 1); // (R, G, B, A)
  return vec4f(1, 0, 0, 1); // red is cool
}

@compute @workgroup_size(256)
fn computeMain(@builtin(global_invocation_id) global_id: vec3u) {
  let idx = global_id.x;
  
  if (idx < arrayLength(&particlesIn)) {
    particlesOut[idx] = particlesIn[idx];
    
    let particle = particlesIn[idx];
    //f = ma
    let forces = calculateForces(idx);
    let accel = forces / mass;
    var newVel = particle.vel + accel * steps_per_update;
    
    //cap the velocity
    if (abs(newVel.x) > max_vel) {
      newVel.x = max_vel * (abs(newVel.x) / newVel.x);
    }
    if (abs(newVel.y) > max_vel) {
      newVel.y = max_vel * (abs(newVel.y) / newVel.y);
    }
    
    //damp the velocity
    newVel *= velocity_damping;

    // update particle position
    var newPos = particle.pos + newVel * steps_per_update;
    

    //keep the particles in the bounding box, damp a bit from bouncing off the sides
    if (newPos.x < left_bound){
      newPos.x = left_bound;
      newVel.x *= -.9;
    }
    else if (newPos.x > right_bound){
      newPos.x = right_bound;
      newVel.x *= -.9;
    }
    if (newPos.y < bot_bound){
      newPos.y = bot_bound;
      newVel.y *= -.9;
    }
    else if (newPos.y > top_bound){
      newPos.y = top_bound;
      newVel.y *= -.9;
    }
    particlesOut[idx].pos = newPos;
    particlesOut[idx].vel = newVel;
  }
}

//find all the forces that should be applied to a given particle, and the net direction of these forces
fn calculateForces(idx: u32) -> vec2f{
  
  var particle = particlesIn[idx];
  var forces = vec2f(0, 0);
  //apply gravity
  forces += vec2f(0, -9.81) * gravity_multiplier;

  // calculate pressure force
  forces -= pressureAproximation(particle.pos, idx) * pressure_multiplier; //negative for some reason...
  
  //calculate viscosity force
  forces += viscocityApproximation(idx) * viscosity_multiplier;

  forces = max(min(forces, vec2f(max_force, max_force)), vec2f(-max_force, -max_force));
  return forces;
}

fn densityApproximation(position: vec2f) -> f32{
  if (use_acceleration == 0){
    return rawDensityCalculation(position);
  } 
  return acceleratedDensityCalculation(position);
}

fn rawDensityCalculation(position : vec2f) -> f32{
  var density = 0.0;
  var particle2 : Particle;
  for (var i = 0; i < i32(arrayLength(&particlesIn)); i++){
    particle2 = particlesIn[i];
    var distance = length(position - particle2.pos);
    density += mass * smoothingFunction(smoothing_radius, distance);
  }
  return density;
}

fn acceleratedDensityCalculation(position: vec2f) -> f32{
  var density = 0.0;
  var particle2 : Particle;
  //find the cell that the particle is in
  var x = i32(floor((position.x - left_bound) / cell_size));
  var y = i32(floor((position.y - bot_bound) / cell_size));
  for(var adjacent_x = x-1; adjacent_x < x+2; adjacent_x++){
    for (var adjacent_y = y-1; adjacent_y < y+2; adjacent_y++){ //iterate over adjacent cells from the cell_pos
        //boundary checking
        if (adjacent_x < 0 || adjacent_x > grid_length){ continue; }
        if (adjacent_y < 0 || adjacent_y > grid_height){ continue; }

        var cell_start = (adjacent_y * grid_length + adjacent_x) * max_per_cell; //recall that the first index of each cell is used as an index
        for (var i = cell_start + 1; i < cell_start + max_per_cell; i++){
          var pIdx = gridIn[i]; // the particle index
          if (pIdx == -1) {  //check if no more particles in this cell
            break;
          }
        particle2 = particlesIn[pIdx];
        var distance = length(position - particle2.pos);
        density += mass * smoothingFunction(smoothing_radius, distance);   
      }
    }
  }
  // for(var i = 0; i < i32(arrayLength(&gridIn)); i++){
  //  var pIdx = gridIn[i]; // the particle index
  // // for (var pIdx = 0; pIdx < i32(arrayLength(&particlesIn)); pIdx++) {
  //   if (pIdx == -1) {  //check if no more particles in this cell
  //     continue;
  //   }
  //   particle2 = particlesIn[pIdx];
  //   var distance = length(position - particle2.pos);
  //   density += mass * smoothingFunction(smoothing_radius, distance);
  // }
  return density;
}

//given a certain density, how hard should the fluid push(pressure)
//further away from target density, the higher the pressure
fn densityToPressure(density : f32) -> f32{
  return density - target_density;
}

fn pressureAproximation(position: vec2f, particle_idx: u32) -> vec2f{
  if (use_acceleration == 0){
    return rawPressureCalculation(position, particle_idx);
  }
  return acceleratedPressureCalculation(position, particle_idx);
}

fn rawPressureCalculation(position: vec2f, particle_idx: u32) -> vec2f{
  var pressure_force = vec2f(0, 0);
  var current_density = densityBuffer[particle_idx];
  var particle2 : Particle;
  
  for (var i = 0; i < i32(arrayLength(&particlesIn)); i++){
    pressure_force += pressureFromPoint(particle_idx, i);
  }
  return pressure_force;
}

fn acceleratedPressureCalculation(position: vec2f, particle_idx: u32) -> vec2f{
  var pressure_force = vec2f(0, 0);
  
  //find the cell that the particle is in
  var x = i32(floor((position.x - left_bound) / cell_size));
  var y = i32(floor((position.y - bot_bound) / cell_size));
  for(var adjacent_x = x-1; adjacent_x < x+2; adjacent_x++){
    for (var adjacent_y = y-1; adjacent_y < y+2; adjacent_y++){ //iterate over adjacent cells from the cell_pos
        //boundary checking
        if (adjacent_x < 0 || adjacent_x > grid_length){ continue; }
        if (adjacent_y < 0 || adjacent_y > grid_height){ continue; }

        var cell_start = (adjacent_y * grid_length + adjacent_x) * max_per_cell; //recall that the first index of each cell is used as an index
        for (var i = cell_start + 1; i < cell_start + max_per_cell; i++){
          var pIdx = gridIn[i]; // the particle index
          if (pIdx == -1) {  //check if no more particles in this cell
            break;
          }
          pressure_force += pressureFromPoint(particle_idx, pIdx);
      }
    }
  }
  // for(var i = 0; i < i32(arrayLength(&gridIn)); i++){
  //  var pIdx = gridIn[i]; // the particle index
  // // for (var pIdx = 0; pIdx < i32(arrayLength(&particlesIn)); pIdx++) {
  //   if (pIdx == -1) {  //check if no more particles in this cell
  //     continue;
  //   }
  //   pressure_force += pressureFromPoint(particle_idx, pIdx);
  // }
  return pressure_force;
}

fn pressureFromPoint(particle1Idx: u32, particle2Idx: i32) -> vec2f{
  var temp = particlesIn[particle1Idx].pos - particlesIn[particle2Idx].pos;
  var distance = length(temp);
  var pressure_force = vec2f(0, 0);
  if (distance > 0.0000001) {
    var direction = temp / distance;
    var current_density = densityBuffer[particle1Idx];
    var other_density = densityBuffer[particle2Idx];
    var slope = smoothingDerivative(smoothing_radius, distance);
    var shared_pressure = getSharedPressure(current_density, other_density);
    pressure_force = direction * slope * mass * shared_pressure;
    if (length(pressure_force) > 1000.f) {
      pressure_force *= (1000.f / length(pressure_force));
    }
  }
  return pressure_force;
}

fn getSharedPressure(density1: f32, density2: f32) -> f32{
  var pressure1 = densityToPressure(density1);
  var pressure2 = densityToPressure(density2);
  return (pressure1 + pressure2) / 2;
}


fn viscocityApproximation(idx: u32) -> vec2f{
  if (use_acceleration == 0){
    return rawViscosityCalculation(idx);
  }
  return acceleratedViscocityCalculation(idx);
}
fn rawViscosityCalculation(idx : u32) -> vec2f{
  var viscosity_force = vec2f(0,0);
  var particle = particlesIn[idx];
  for (var i = 0; i < i32(arrayLength(&particlesIn)); i++){
    var particle2 = particlesIn[i];
    var distance = length(particle.pos - particle2.pos);
    var influence = viscocitySmoothFunction(smoothing_radius, distance);
    viscosity_force += (particle2.vel - particle.vel) * influence;
  }
  return viscosity_force;
}

fn acceleratedViscocityCalculation(idx: u32) -> vec2f{
  var viscosity_force = vec2f(0,0);
  var particle = particlesIn[idx];
  var particle2 : Particle;
  //find the cell that the particle is in
  var x = i32(floor((particle.pos.x - left_bound) / cell_size));
  var y = i32(floor((particle.pos.y - bot_bound) / cell_size));
  for(var adjacent_x = x-1; adjacent_x < x+2; adjacent_x++){
    for (var adjacent_y = y-1; adjacent_y < y+2; adjacent_y++){ //iterate over adjacent cells from the cell_pos
      //boundary checking
      if (adjacent_x < 0 || adjacent_x > grid_length){ continue; }
      if (adjacent_y < 0 || adjacent_y > grid_height){ continue; }

      var cell_start = (adjacent_y * grid_length + adjacent_x) * max_per_cell; //recall that the first index of each cell is used as an index
      for (var i = cell_start + 1; i < cell_start + max_per_cell; i++){
        var pIdx = gridIn[i]; // the particle index
        if (pIdx == -1) {  //check if no more particles in this cell
          break;
        }
        var particle2 = particlesIn[pIdx];
        var distance = length(particle.pos - particle2.pos);
        var influence = viscocitySmoothFunction(smoothing_radius, distance);
        viscosity_force += (particle2.vel - particle.vel) * influence; 
      }
    }
  }
  return viscosity_force;
}

//function to calculate how the influence a particle has on the density decays as the distance away from a particle grows
fn smoothingFunction(smoothing_radius : f32, distance : f32) -> f32{
  //for calc reasons, this is the volume of the function here. Need to keep that constant so we divide by it (?) just watch sebas' video
  var influence = pow(1 - (min(abs(distance), smoothing_radius))/smoothing_radius, smoothing_rate);
  return influence / func_volume;
}

fn smoothingDerivative(smoothing_radius : f32, distance: f32) -> f32{
  if (distance > smoothing_radius){
    return 0.0;
  }
  var deriv = -(smoothing_rate / smoothing_radius) * pow(1 - abs(distance) / smoothing_radius, smoothing_rate - 1) - 1;
  return deriv / func_volume;
}

fn viscocitySmoothFunction(smoothing_radius: f32, distance: f32) -> f32{
  //use a bell curve to get a smooth shape centered around 0. mean of 0, variance of smoothing radius
  if (distance < .000001){
    return 1;
  }
  else if (distance > smoothing_radius){
    return 0;
  }
  return pow(e, -distance * distance / (2 * smoothing_radius));
}

@compute @workgroup_size(256)
fn clearGridStructure(@builtin(global_invocation_id) global_id: vec3u){
  //empty out then write the new grid positions
  if (global_id.x < arrayLength(&gridOut)){
    for (var i = 0; i < max_per_cell; i++){
      atomicStore(&gridOut[global_id.x * max_per_cell + u32(i)], i32(-1));
    }
  }
}

@compute @workgroup_size(256)
fn computeGridStructure(@builtin(global_invocation_id) global_id: vec3u){
  //write the new grid positions
  //global_id is the index of the particle, so there must be a dispatch for each of the particles
  if (global_id.x < arrayLength(&particlesOut)){
    var x = i32(floor((particlesOut[global_id.x].pos.x - left_bound) / cell_size));
    var y = i32(floor((particlesOut[global_id.x].pos.y - bot_bound) / cell_size));
    var cell_start = (y * grid_length + x) * max_per_cell;
    
    var idx = atomicAdd(&gridOut[cell_start], 1); // the first position of each cell stores the next open position in the given cell, retrieve and then increment this value
    atomicStore(&gridOut[cell_start + idx + 1], i32(global_id.x)); //store the particle idx at the position pointed to by the first position of the cell, +1 to skip the index spot
  }
}

@compute @workgroup_size(256)
fn computeDensity(@builtin(global_invocation_id) global_id: vec3u){
  if (global_id.x < arrayLength(&densityBuffer)){
    densityBuffer[global_id.x] = densityApproximation(particlesIn[global_id.x].pos);
  }
}