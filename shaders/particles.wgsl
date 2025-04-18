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

// TODO 3: Define a struct to store a particle
struct Particle {
  pos : vec2f,
  initPos : vec2f,
  vel : vec2f,
  initVel : vec2f,
  lifeTime : vec2f,
};
const values_per_particle = 10;
//boundary box for the fluid
const left_bound = -.5;
const right_bound = .5;
const top_bound = .3;
const bot_bound = -.3;
//grid parameters
const grid_size = .07;
const grid_length = i32((abs(left_bound) + abs(right_bound))/grid_size);
const grid_height = i32((abs(bot_bound) + abs(top_bound))/grid_size);
const max_per_cell = 64;
const use_acceleration = 1; //non-zero to use acceleration structures
//fluid simulation parameters
const mass = 1.0;
const smoothingRadius = grid_size; //this way we only need to check adjacent grid cells
const smoothingRate = 4;
const pi = 3.14159265;
const e = 2.71828;
const func_volume = (smoothingRate + 1)/(2 * pi * smoothingRadius);
const num_particles = 256;
const targetDensity = 1;
const pressureMultiplier= 1;
const gravityMultiplier = 0.05;
const viscosityMultiplier = 2;
const maxVel = .1;
const steps_per_update = .5; //lower value = slower/more stable simulation, somewhere between .1 and 1


// TODO 4: Write the bind group spells here using array<Particle>
@group(0) @binding(0) var<storage> particlesIn: array<Particle>;
@group(0) @binding(1) var<storage, read_write> particlesOut: array<Particle>;
@group(0) @binding(2) var<storage, read_write> timeBuffer: array<u32>;
@group(0) @binding(3) var<storage> gridIn: array<i32>;
@group(0) @binding(4) var<storage, read_write> gridOut: array<atomic<i32>>;


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
    if (abs(newVel.x) > maxVel) {
      newVel.x = min(abs(newVel.x), maxVel) * (abs(newVel.x) / newVel.x);
    }
    if (abs(newVel.y) > maxVel) {
      newVel.y = min(abs(newVel.y), maxVel) * (abs(newVel.y) / newVel.y);
    }
    
    //damp the velocity
    newVel *= .80;

    // update particle position
    var newPos = particle.pos + newVel * steps_per_update;
    

    //keep the particles on the screen
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
  forces += vec2f(0, -9.81) * gravityMultiplier;

  //precompute densities to make things faster
  var densities : array<f32, num_particles>;
  for (var i = 0; i < i32(arrayLength(&particlesIn)); i++){
    densities[i] = pointDensity(particlesIn[i].pos);
  }

  // calculate pressure force
  forces -= pressureAproximation(particle.pos, idx, densities) * pressureMultiplier;
  
  //calculate viscosity force
  forces += rawViscosityCalculation(idx) * viscosityMultiplier;


  forces = max(min(forces, vec2f(0.01, 0.01)), vec2f(-0.01, -0.01));
  return forces;
}



//returns the density at a given point
fn pointDensity(position : vec2f) -> f32{
  // if (use_acceleration == 0){
    return rawDensityCalculation(position);
  // } 
  // return acceleratedDensityCalculation(position);
}

fn rawDensityCalculation(position : vec2f) -> f32{
  var density = 0.0;
  var particle2 : Particle;
  for (var i = 0; i < i32(arrayLength(&particlesIn)); i++){
    particle2 = particlesIn[i];
    var temp = position - particle2.pos;
    var distance = sqrt(dot(temp, temp));
    density += mass * smoothingFunction(smoothingRadius, distance);
  }
  return density;
}

fn acceleratedDensityCalculation(position: vec2f) -> f32{
  var density = 0.0;
  var particle2 : Particle;
  //find the cell that the particle is in
  var x = floor((position.x + 1) / grid_size); //how many grid cells does it take to go from -1 to position.x TODO: these lines dont take in to considerating if the botleft is moved
  var y = floor((position.y + 1) / grid_size);
  var cell_pos = vec2f(x,y);
  var radius_in_cells = i32(ceil(smoothingRadius / grid_size));
  //the cells that are within the smoothing distance of the position must be checked, iterate over the cells, check if they are within the smoothing radius, then compute those particles
  for (var dx = i32(-radius_in_cells); dx <= radius_in_cells; dx++) {
    for (var dy = i32(-radius_in_cells); dy <= radius_in_cells; dy++) {
      var cell_x = i32(x) + dx;
      var cell_y = i32(y) + dy;

      // Boundary check
      if (cell_x < 0 || cell_x >= grid_length || cell_y < 0 || cell_y >= grid_length) {
        continue;
      }
      //the particles contained in cell (cell_x, cell_y) are close enough to probably contribute to the function, and thus must be checked
      var cell_offset = (cell_y * grid_length + cell_x) * max_per_cell;
      for(var particle_index = cell_offset; particle_index < cell_offset + max_per_cell; particle_index++){
        var i = gridIn[particle_index];
        if (i == -1){ break; }
        //here is where we can access only adjacent particles within smoothing radius
        particle2 = particlesIn[i];
        var temp = position - particle2.pos;
        var distance = sqrt(dot(temp, temp));
        density += mass * smoothingFunction(smoothingRadius, distance);
      }
      }
    }
  return density;
}

//given a certain density, how hard should the fluid push(pressure)
//further away from target density, the higher the pressure
fn densityToPressure(density : f32) -> f32{
  return density - targetDensity;
}

fn pressureAproximation(position: vec2f, particle_idx: u32, densities : array<f32, num_particles>) -> vec2f{
  if (use_acceleration == 0){
    return rawPressureCalculation(position, particle_idx, densities);
  }
  return acceleratedPressureCalculation(position, particle_idx, densities);
}

fn rawPressureCalculation(position: vec2f, particle_idx: u32, densities : array<f32, num_particles>) -> vec2f{
  var pressure_force = vec2f(0, 0);
  var current_density = densities[particle_idx];
  var particle2 : Particle;
  
  for (var i = 0; i < i32(arrayLength(&particlesIn)); i++){
    particle2 = particlesIn[i];
    var temp = position - particle2.pos;
    var distance = length(position - particle2.pos);
    if (distance > 0.0000001) {
      var direction = temp / distance;
      var density = densities[i];
      var slope = smoothingDerivative(smoothingRadius, distance);
      var shared_pressure = (current_density + density) / 2;
      pressure_force += direction * slope * mass * shared_pressure;
      if (length(pressure_force) > 1000.f) {
        pressure_force *= (1000.f / length(pressure_force));
      }
    }
  }
  return pressure_force;
}

fn acceleratedPressureCalculation(position: vec2f, particle_idx: u32, densities : array<f32, num_particles>) -> vec2f{
  var pressure_force = vec2f(0, 0);
  var current_density = densities[particle_idx];
  var particle2 : Particle;

  //find the cell that the particle is in
  var x = i32(floor((position.x - left_bound) / grid_size));
  var y = i32(floor((position.y - bot_bound) / grid_size));
  
  for(var x_offset = -1; x_offset < 2; x_offset++){
    for (var y_offset = -1; y_offset < 2; y_offset++){ //iterate over adjacent cells from the cell_pos
        var adjacent_x = x + x_offset;
        var adjacent_y = y + y_offset;
        //boundary checking
        if (adjacent_x < 0 || adjacent_x < grid_length){ continue; }
        if (adjacent_y < 0 || adjacent_y < grid_height){ continue; }

        var cell_start = (adjacent_y * grid_length + adjacent_x) * max_per_cell;
        for (var i = 0; i < max_per_cell; i++){
          var index = cell_start + i;
          if (index == -1) { break; }
          var particle2 = particlesIn[index];
          var temp = position - particle2.pos;
          var distance = sqrt(dot(temp, temp));
          if (distance > 0.0000001) {
            var direction = temp / distance;
            var density = densities[i];
            var slope = smoothingDerivative(smoothingRadius, distance);
            var shared_pressure = (current_density + density) / 2;
            pressure_force += direction * slope * mass * shared_pressure;
            if (length(pressure_force) > 1000.f) {
              pressure_force *= (1000.f / length(pressure_force));
            }
        }
      }
    }
  }
  return pressure_force;
}

fn getSharedPressure(density1: f32, density2: f32) -> f32{
  var pressure1 = densityToPressure(density1);
  var pressure2 = densityToPressure(density2);
  return (pressure1 + pressure2) / 2;
}

//function to calculate how the influence a particle has on the density decays as the distance away from a particle grows
fn smoothingFunction(smoothingRadius : f32, distance : f32) -> f32{
  //for calc reasons, this is the volume of the function here. Need to keep that constant so we divide by it (?) just watch sebas' video
  var influence = pow(1 - (min(abs(distance), smoothingRadius))/smoothingRadius, smoothingRate);
  return influence / func_volume;
}

fn smoothingDerivative(smoothingRadius : f32, distance: f32) -> f32{
  if (distance > smoothingRadius){
    return 0.0;
  }
  var deriv = -(smoothingRate / smoothingRadius) * pow(1 - abs(distance) / smoothingRadius, smoothingRate - 1) - 1;
  return deriv / func_volume;
}

fn rawViscosityCalculation(idx : u32) -> vec2f{
  var viscosity_force = vec2f(0,0);
  var particle = particlesIn[idx];
  for (var i = 0; i < i32(arrayLength(&particlesIn)); i++){
    var particle2 = particlesIn[i];
    var distance = length(particle.pos - particle2.pos);
    var influence = viscocitySmoothFunction(smoothingRadius, distance);
    viscosity_force += (particle2.vel - particle.vel) * influence;
  }
  return viscosity_force;
}

fn viscocitySmoothFunction(smoothingRadius: f32, distance: f32) -> f32{
  //use a bell curve to get a smooth shape centered around 0. mean of 0, variance of smoothing radius
  if (distance < .000001){
    return 1;
  }
  else if (distance > smoothingRadius){
    return 0;
  }
  var value = pow(e, -distance * distance / (2 * smoothingRadius));
  return value;
}



fn intersectForce(particle1: Particle, particle2: Particle) -> vec2f{
  //check if the distance from the first particle is 2 radius away from the particles
  //we are calculating the force that should be applied onto the idx1 particle
  var radius = 0.0125 / 2f;
  var distance = sqrt(pow(particle1.pos.x - particle2.pos.x, 2) + pow(particle1.pos.y - particle2.pos.y, 2));
  //if the two particle centers are within 2 radius of eachother, then they are colliding
  if (distance < 2 * radius){
    //the repulsion force is along the line between the two centers of the particles, pointing towards idx1
    var force = particle1.pos - particle2.pos;
    var scaled_force = force * .1;
    return force;
  }
  return vec2f(0, 0);
}

@compute @workgroup_size(256)
fn computeGridStructure(@builtin(global_invocation_id) global_id: vec3u){
  //zero out then write the new grid positions
  for (var i = 0; i < i32(arrayLength(&gridOut)); i++){
    atomicStore(&gridOut[i], -1);
  }
  for(var i = 0; i < i32(arrayLength(&particlesOut)); i++){
    var x = floor((particlesOut[i].pos.x + 1) / grid_size);
    var y = floor((particlesOut[i + 1].pos.y + 1) / grid_size);
    var index = i32(y * max_per_cell + x);
    for (var offset = 0; offset < max_per_cell; offset++){
      if (atomicLoad(&gridOut[index + offset]) == -1){
        atomicStore(&gridOut[index + offset], i);
        break;
      }
    }
  }
}