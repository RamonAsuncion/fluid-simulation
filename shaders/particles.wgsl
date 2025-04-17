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

//fluid simulation parameters
const mass = 1.0;
const smoothingRadius = 0.07;
const smoothingRate = 4;
const pi = 3.14159265;
const func_volume = (smoothingRate + 1)/(2 * pi * smoothingRadius);
const num_particles = 256;
const targetDensity = 1;
const pressureMultiplier= 1;
const gravityMultiplier = 0.01;
//grid parameters
const grid_size = .05;
const grid_length = i32(2/grid_size);
const max_per_cell = 64;
//boundary box for the fluid
const left_bound = -.3;
const right_bound = .3;
const top_bound = .3;
const bot_bound = -.3;
const use_acceleration = 0; //non-zero to use acceleration structures
const maxVel = .1;


// TODO 4: Write the bind group spells here using array<Particle>
@group(0) @binding(0) var<storage> particlesIn: array<Particle>;
@group(0) @binding(1) var<storage, read_write> particlesOut: array<Particle>;
@group(0) @binding(2) var<storage, read_write> timeBuffer: array<u32>;
@group(0) @binding(3) var<storage> gridIn: array<i32>;
@group(0) @binding(4) var<storage, read_write> gridOut: array<i32>;


@vertex
fn vertexMain(@builtin(instance_index) idx: u32, @builtin(vertex_index) vIdx: u32) -> @builtin(position) vec4f {
  // TODO 5: Revise the vertex shader to draw circle to visualize the particles
  let particle = particlesIn[idx];
  let size = 0.0125 / 2f;
  let theta = 2. * pi / 8 * f32(vIdx);
  let x = cos(theta) * size;
  let y = sin(theta) * size;
  return vec4f(vec2f(x + particle.pos.x, y + particle.pos.y), 0, 1);
}

@fragment
fn fragmentMain() -> @location(0) vec4f {
  // return vec4f(238.f/255, 118.f/255, 35.f/255, 1); // (R, G, B, A)
  return vec4f(1, 0, 0, 1);
}

@compute @workgroup_size(256)
fn computeMain(@builtin(global_invocation_id) global_id: vec3u) {
  let idx = global_id.x;
  
  if (idx < arrayLength(&particlesIn)) {
    particlesOut[idx] = particlesIn[idx];
    
    let particle = particlesIn[idx];
    //f = ma
    let forces = calculateForces(idx);
    //let forces = vec2f(0, -0.1);
    let accel = forces / mass;
    var newVel = particle.vel + accel;
    
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
    var newPos = particle.pos + newVel;
    

    //keep the particles on the screen
    if (newPos.x < left_bound){
      newPos.x = left_bound;
      newVel.x *= -1;
    }
    else if (newPos.x > right_bound){
      newPos.x = right_bound;
      newVel.x *= -1;
    }
    if (newPos.y < bot_bound){
      newPos.y = bot_bound;
      newVel.y *= -1;
    }
    else if (newPos.y > top_bound){
      newPos.y = top_bound;
      newVel.y *= -1;
    }
    particlesOut[idx].pos = newPos;
    particlesOut[idx].vel = newVel;

    //zero out then write the new grid positions
    // for (var i = 0; i < i32(arrayLength(&gridOut)); i++){
    //   gridOut[i] = -1;
    // }
    // for(var i = 0; i < i32(arrayLength(&particlesOut)); i++){
    //   var x = floor((particlesOut[i].pos.x + 1) / grid_size);
    //   var y = floor((particlesOut[i + 1].pos.y + 1) / grid_size);
    //   var index = i32(y * max_per_cell + x);
    //   for (var offset = 0; offset < max_per_cell; offset++){
    //     if (gridOut[index + offset] == -1){
    //       gridOut[index + offset] = i;
    //       break;
    //     }
    //   }
    // }
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
  forces = max(min(forces, vec2f(0.01, 0.01)), vec2f(-0.01, -0.01));
  return forces;
}



//returns the density at a given point
fn pointDensity(position : vec2f) -> f32{
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
  return pressure_force;
}

fn acceleratedPressureCalculation(position: vec2f, particle_idx: u32, densities : array<f32, num_particles>) -> vec2f{
  var pressureForce = vec2f(0, 0);
  var current_density = densities[particle_idx];
  var particle2 : Particle;

  //find the cell that the particle is in
  var x = floor((position.x + 1) / grid_size);
  var y = floor((position.y + 1) / grid_size);
  var cell_pos = vec2f(x,y); //world coord for the bot left of the cell containing the position
  
  //the cells that are within the smoothing distance of the position must be checked, iterate over the cells, check if they are within the smoothing radius, then compute those particles
  for(var cell_x = 0; cell_x < grid_length; cell_x++){
    for (var cell_y = 0; cell_y < grid_length; cell_y++){
      var temp = cell_pos - vec2f(f32(cell_x), f32(cell_y)) * grid_size;
      var distance = sqrt(dot(temp, temp)); //distance from the position's cell to the current cell, comparing bot left corners
      if (distance < smoothingRadius){
        //the particles contained in cell (cell_x, cell_y) are close enough to probably contribute to the function, and thus must be checked
        var cell_offset = (cell_y * grid_length + cell_x) * max_per_cell;
        for(var particle_index = cell_offset; particle_index < cell_offset + max_per_cell; particle_index++){
          var i = gridIn[particle_index];
          if (i == -1){ break; }
          //here is where we can access only adjacent particles within smoothing radius
          particle2 = particlesIn[i];
          var temp = position - particle2.pos;
          var distance = sqrt(dot(temp, temp));
          var direction = (position - particle2.pos) / distance;
          var density = densities[i];
          var slope = smoothingDerivative(smoothingRadius, distance);
          var pressure = densityToPressure(density);
          var sharedPressure = (current_density + density) / 2;

          pressureForce += direction * slope * mass * sharedPressure / current_density;
        }
      }
    }
  }
  return pressureForce;
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