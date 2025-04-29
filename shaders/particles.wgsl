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

// struct to store 3D PGA multivector
struct MultiVector {
  s: f32, 
  exey: f32, 
  exez: f32, 
  eyez: f32, 
  eoex: f32, 
  eoey: f32, 
  eoez: f32, 
  exeyez: f32, 
  eoexey: f32, 
  eoexez: f32, 
  eoeyez: f32,
  ex: f32, 
  ey: f32, 
  ez: f32, 
  eo: f32,
  eoexeyez: f32
}

// struct to store 3D PGA pose
struct Camera {
  motor: MultiVector,
  focal: vec2f,
  res: vec2f,
}

struct MouseInteraction {
  position: vec2f,
  isDown: f32,
  radius: f32
};

struct BoundaryBox {
  left: f32,
  right: f32,
  top: f32,
  bottom: f32,
  front: f32,
  back: f32
};

/*
The paper talks about "color field" to identity surface particles.
I want to create a buffer to identify surface particles and use the 
surface information to target mouse interface.

*/

// const surfaceThreshold = 0.4; // detecting surface particles

// struct SurfaceInfo {
//   isOnSurface: i32, // 1 if particle is on surface, 0 otherwise
//   normal: vec2f,
//   curvature: f32,
// };

// define a constant
const EPSILON : f32 = 0.00000001;

// Fluid simulation parameters
const mass = 1.0;
const smoothingRadius = .1;
const pi = 3.14159265;
const num_particles = 64;
const targetDensity = 1;
const pressureMultiplier= 0.000001;
const gravityMultiplier = 0.000001;

// Grid parameters
const grid_size = .05;
const grid_length = i32(2/grid_size);
const max_per_cell = 64;
const use_acceleration = 0; //non-zero to use acceleration structures
const maxVel = .1;

@group(0) @binding(0) var<storage> particlesIn: array<Particle>;
@group(0) @binding(1) var<storage, read_write> particlesOut: array<Particle>;
// @group(0) @binding(2) var<storage, read_write> timeBuffer: array<u32>;
@group(0) @binding(2) var<storage, read_write> timeBuffer: array<f32>;
@group(0) @binding(3) var<storage> gridIn: array<i32>;
@group(0) @binding(4) var<storage, read_write> gridOut: array<i32>;
@group(0) @binding(5) var<uniform> boundaryBox: BoundaryBox;
@group(0) @binding(6) var<uniform> cameraPose: Camera ;
@group(0) @binding(7) var<uniform> mouse: MouseInteraction;

// @vertex
// fn vertexMain(@builtin(instance_index) idx: u32, @builtin(vertex_index) vIdx: u32) -> @builtin(position) vec4f {
//   let particle = particlesIn[idx];
//   let size = 0.0125 / 2f;
//   let theta = 2. * pi / 8 * f32(vIdx);
//   let x = cos(theta) * size;
//   let y = sin(theta) * size;
  // let worldPos = vec3f(x + particle.pos.x, y + particle.pos.y, 0.0);
  // let transformedPos = applyMotorToPoint(worldPos, cameraPose.motor);
//   return vec4f(vec2f(x + particle.pos.x, y + particle.pos.y), 0, 1);
// }


@vertex
fn vertexMain(@builtin(vertex_index) vIdx: u32) -> @builtin(position) vec4f {
  let corners = array<vec3f, 8>(
    vec3f(boundaryBox.left, boundaryBox.bottom, boundaryBox.back),
    vec3f(boundaryBox.right, boundaryBox.bottom, boundaryBox.back),
    vec3f(boundaryBox.right, boundaryBox.top, boundaryBox.back),
    vec3f(boundaryBox.left, boundaryBox.top, boundaryBox.back),
    vec3f(boundaryBox.left, boundaryBox.bottom, boundaryBox.front),
    vec3f(boundaryBox.right, boundaryBox.bottom, boundaryBox.front),
    vec3f(boundaryBox.right, boundaryBox.top, boundaryBox.front),
    vec3f(boundaryBox.left, boundaryBox.top, boundaryBox.front)
  );

  // 2 triangles per face
  let triangleIndices = array<u32, 36>(
    0, 1, 2,  0, 2, 3,
    4, 6, 5,  4, 7, 6,
    0, 3, 7,  0, 7, 4,
    1, 5, 6,  1, 6, 2,
    0, 4, 5,  0, 5, 1,
    3, 2, 6,  3, 6, 7
  );
  
  let vertex = corners[triangleIndices[vIdx]];
  let transformedPos = applyMotorToPoint(vertex, cameraPose.motor);
  
  let depth = 3.0;
  let perspectiveFactor = 1.0 / (1.0 - transformedPos.z/depth);
  
  let x = transformedPos.x * perspectiveFactor;
  let y = transformedPos.y * perspectiveFactor;
  
  return vec4f(x, y, transformedPos.z, 1.0);
}

@fragment
fn fragmentMain() -> @location(0) vec4f {
  return vec4f(238.f/255, 118.f/255, 35.f/255, 1); // (R, G, B, A)
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
    var newVel = particle.vel + accel;
    //damp the velocity
    newVel *= .90;

    //cap the velocity
    if (newVel.x > maxVel || newVel.x < -maxVel){
      newVel.x = maxVel;
    }
    if (newVel.y > maxVel || newVel.y < -maxVel){
      newVel.y = maxVel;
    }

    // var newPos = particle.pos + particle.vel; 
    var newPos = particle.pos + particle.vel * timeBuffer[1]; 
    
    //keep the particles on the screen
    if (newPos.x < boundaryBox.left) {
      newPos.x = boundaryBox.left;
      newVel.x *= -1;
    } else if (newPos.x > boundaryBox.right) {
      newPos.x = boundaryBox.right;
      newVel.x *= -1;
    }
    if (newPos.y < boundaryBox.bottom) {
      newPos.y = boundaryBox.bottom;
      newVel.y *= -1;
    } else if (newPos.y > boundaryBox.top) {
      newPos.y = boundaryBox.top;
      newVel.y *= -1;
    }
    particlesOut[idx].pos = newPos;
    particlesOut[idx].vel = newVel;

    //zero out then write the new grid positions
    for (var i = 0; i < i32(arrayLength(&gridOut)); i++){
      gridOut[i] = -1;
    }
    for(var i = 0; i < i32(arrayLength(&particlesOut)); i++){
      var x = floor((particlesOut[i].pos.x + 1) / grid_size);
      var y = floor((particlesOut[i + 1].pos.y + 1) / grid_size);
      var index = i32(y * max_per_cell + x);
      for (var offset = 0; offset < max_per_cell; offset++){
        if (gridOut[index + offset] == -1){
          gridOut[index + offset] = i;
          break;
        }
      }
    }
  }
}

// TODO:
// bring mouse position to simulation
// calculate distance from mouse to particles
// apply a force that decreases with distance (using a radial falloff function) - Recommendation.
// scale the force based on mouse movement speed or button state
fn calculateMouseForce(position: vec2f) -> vec2f {
  if (mouse.isDown < 0.5) {
    return vec2f(0.0, 0.0); // not pressed
  }
  
  // get distance
  let distVec = position - mouse.position;
  let distSq = dot(distVec, distVec);
  let radius = mouse.radius;
  
  // apply repulsion force
  if (distSq < radius * radius) {
    let dist = sqrt(distSq);
    
    // normalize
    var dir = distVec;
    if (dist > 0.0001) {
      dir = distVec / dist;
    }
    
    // apply stronger force to closer particles (inverse linear falloff)
    let strength = 0.01 * (1.0 - dist / radius);
    
    return dir * strength;
  }
  
  return vec2f(0.0, 0.0);
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
  forces += pressureApproximation(particle.pos, idx, densities);

  // forces += calculateMouseForce(particle.pos);

  return forces;
}

//returns the density at a given point
fn pointDensity(position : vec2f) -> f32{
  if (use_acceleration == 0) {
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
  // var x = floor((particlesOut[i].pos.x + 1) / grid_size);
  // var y = floor((particlesOut[i].pos.y + 1) / grid_size);
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
  var difference = density - targetDensity;
  var pressure = difference * pressureMultiplier;
  return pressure;
}

fn pressureApproximation(position: vec2f, particle_idx: u32, densities : array<f32, num_particles>) -> vec2f{
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
    var direction = (position - particle2.pos) / distance;
    var density = densities[i];
    var slope = smoothingDerivative(smoothingRadius, distance);
    var pressure = densityToPressure(density);
    var shared_pressure = (current_density + density) / 2;
    pressure_force += direction * slope * mass * shared_pressure / current_density;
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
  var func_volume = pi * pow(smoothingRadius, 4) / 6;
  var value = max(0, smoothingRadius * smoothingRadius - distance * distance * distance);
  return (smoothingRadius - distance) * (smoothingRadius - distance) / func_volume;
}

fn smoothingDerivative(smoothingRadius : f32, distance: f32) -> f32{
  if (distance > smoothingRadius){
    return 0.0;
  }
  var f = smoothingRadius * smoothingRadius - distance * distance;
  var scale = 12 / (pi * pow(smoothingRadius, 4));
  return scale * (distance - smoothingRadius); 
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


// the geometric product 
fn geometricProduct(a: MultiVector, b: MultiVector) -> MultiVector { 
  // The geometric product rules are:
  //   1. eoeo = 0, exex = 1 and eyey = 1, ezez = 1
  //   2. eoex + exeo = 0, eoey + eyeo = 0, eoez + ezeo = 0
  //   3. exey + eyex = 0, exez + ezex = 0, eyez + ezey = 0
  // This results in the following product table:
  var r: MultiVector;
  r.s = a.s * b.s - a.exey * b.exey - a.exez * b.exez - a.eyez * b.eyez - a.exeyez * b.exeyez + a.ex * b.ex + a.ey * b.ey + a.ez * b.ez; // scalar
  r.exey = a.s * b.exey + a.exey * b.s - a.exez * b.eyez + a.eyez * b.exez + a.exeyez * b.ez + a.ex * b.ey - a.ey * b.ex + a.ez * b.exeyez; // exey
  r.exez = a.s * b.exez + a.exey * b.eyez + a.exez * b.s - a.eyez * b.exey - a.exeyez * b.ey + a.ex * b.ez - a.ey * b.exeyez - a.ez * b.ex; // exez
  r.eyez = a.s * b.eyez - a.exey * b.exez + a.exez * b.exey + a.eyez * b.s + a.exeyez * b.ex + a.ex * b.exeyez + a.ey * b.ez - a.ez * b.ey; // eyez
  r.eoex = a.s * b.eoex + a.exey * b.eoey + a.exez * b.eoez - a.eyez * b.eoexeyez + a.eoex * b.s - a.eoey * b.exey - a.eoez * b.exez + a.exeyez * b.eoeyez + a.eoexey * b.ey + a.eoexez * b.ez - a.eoeyez * b.exeyez - a.ex * b.eo + a.ey * b.eoexey + a.ez * b.eoexez + a.eo * b.ex - a.eoexeyez * b.eyez; // eoex
  r.eoey = a.s * b.eoey - a.exey * b.eoex + a.exez * b.eoexeyez + a.eyez * b.eoez + a.eoex * b.exey + a.eoey * b.s - a.eoez * b.eyez - a.exeyez * b.eoexez - a.eoexey * b.ex + a.eoexez * b.exeyez + a.eoeyez * b.ey - a.ex * b.eoexey - a.ey * b.eo + a.ez * b.eoeyez + a.eo * b.ey + a.eoexeyez * b.exez; // eoey
  r.eoez = a.s * b.eoez - a.exey * b.eoexeyez - a.exez * b.eoex - a.eyez * b.eoey + a.eoex * b.exez + a.eoey * b.eyez + a.eoez * b.s + a.exeyez * b.eoexey - a.eoexey * b.exeyez - a.eoexez * b.ex - a.eoeyez * b.ey - a.ex * b.eoexez - a.ey * b.eoeyez - a.ez * b.eo + a.eo * b.ez - a.eoexeyez * b.exey; // eoez
  r.exeyez = a.s * b.exeyez + a.exey * b.ez - a.exez * b.ey + a.eyez * b.ex + a.exeyez * b.s + a.ex * b.eyez - a.ey * b.exez + a.ez * b.exey; // exeyez
  r.eoexey = a.s * b.eoexey + a.exey * b.eo - a.exez * b.eoeyez + a.eyez * b.eoexez + a.eoex * b.ey - a.eoey * b.ex + a.eoez * b.exeyez - a.exeyez * b.eoez + a.eoexey * b.s - a.eoexez * b.eyez + a.eoeyez * b.exez - a.ex * b.eoey + a.ey * b.eoex - a.ez * b.eoexeyez + a.eo * b.exey + a.eoexeyez * b.ez; // eoexey
  r.eoexez = a.s * b.eoexez + a.exey * b.eoeyez + a.exez * b.eo - a.eyez * b.eoexey + a.eoex * b.ez - a.eoey * b.exeyez - a.eoez * b.ex + a.exeyez * b.eoey + a.eoexey * b.eyez + a.eoexez * b.s - a.eoeyez * b.exey - a.ex * b.eoez + a.ey * b.eoexeyez + a.ez * b.eoex + a.eo * b.exez - a.eoexeyez * b.ey; // eoexez
  r.eoeyez = a.s * b.eoeyez - a.exey * b.eoexez + a.exez * b.eoexey + a.eyez * b.eo + a.eoex * b.exeyez + a.eoey * b.ez - a.eoez * b.ey - a.exeyez * b.eoex - a.eoexey * b.exez + a.eoexez * b.exey + a.eoeyez * b.s - a.ex * b.eoexeyez - a.ey * b.eoez + a.ez * b.eoey + a.eo * b.eyez + a.eoexeyez * b.ex; // eoeyez
  r.ex = a.s * b.ex + a.exey * b.ey + a.exez * b.ez - a.eyez * b.exeyez - a.exeyez * b.eyez + a.ex * b.s - a.ey * b.exey - a.ez * b.exez; // ex
  r.ey = a.s * b.ey - a.exey * b.ex + a.exez * b.exeyez + a.eyez * b.ez + a.exeyez * b.exez + a.ex * b.exey + a.ey * b.s - a.ez * b.eyez; // ey
  r.ez = a.s * b.ez - a.exey * b.exeyez - a.exez * b.ex - a.eyez * b.ey - a.exeyez * b.exey + a.ex * b.exez + a.ey * b.eyez + a.ez * b.s; // ez
  r.eo = a.s * b.eo - a.exey * b.eoexey - a.exez * b.eoexez - a.eyez * b.eoeyez + a.eoex * b.ex + a.eoey * b.ey + a.eoez * b.ez + a.exeyez * b.eoexeyez - a.eoexey * b.exey - a.eoexez * b.exez - a.eoeyez * b.eyez - a.ex * b.eoex - a.ey * b.eoey - a.ez * b.eoez + a.eo * b.s - a.eoexeyez * b.exeyez; // eo
  r.eoexeyez = a.s * b.eoexeyez + a.exey * b.eoez - a.exez * b.eoey + a.eyez * b.eoex + a.eoex * b.eyez - a.eoey * b.exez + a.eoez * b.exey - a.exeyez * b.eo + a.eoexey * b.ez - a.eoexez * b.ey + a.eoeyez * b.ex - a.ex * b.eoeyez + a.ey * b.eoexez - a.ez * b.eoexey + a.eo * b.exeyez + a.eoexeyez * b.s; // eoexeyez
  return r;
}

// the reverse of a Multivector
fn reverse(a: MultiVector) -> MultiVector {
  // The reverse is the reverse order of the basis elements
  //  the reverse of a scalar is the scalar
  //  the reverse of exey is eyex = -exey
  //  the reverse of exez is ezex = -exez
  //  the reverse of eyez is ezey = -eyez
  //  the reverse of eoex is exeo = -eoex
  //  the reverse of eoey is eyeo = -eoey
  //  the reverse of eoez is ezeo = -eoez
  //  the reverse of exeyez is ezeyex = exezey = -exeyez
  //  the reverse of eoexey is eyexeo = eoeyex = -eoexey
  //  the reverse of eoexez is ezexeo = eoezex = -eoexez
  //  the reverse of eoeyez is ezeyeo = eoezey = -eoeyez
  //  the reverse of ex, ey, ez, eo are ex, ey, ez, eo
  //  the reverse of eoexeyez is ezeyexeo = -eoezeyex = -eoexezey = eoexeyez
  // So, for [s, exey, exez, eyez, eoex, eoey, eoez, exeyez, eoexey, eoexez, eoeyez, ex, ey, ez, eo, eoexeyez],
  // Its reverse is [s, -exey, -exez, eyez, -eoex, -eoey, -eoez, -exeyez, -eoexey, -eoexez, -eoeyez, ex, ey, ez, eo, eoexeyez].
  return MultiVector(a.s, -a.exey, -a.exez, -a.eyez, -a.eoex, -a.eoey, -a.eoez, -a.exeyez, -a.eoexey, -a.eoexez, -a.eoeyez, a.ex, a.ey, a.ez, a.eo, a.eoexeyez);
}

fn applyMotor(p: MultiVector, m: MultiVector) -> MultiVector {
  // To apply a motor to a point, we use the sandwich operation
  // The formula is m * p * reverse of m
  // Here * is the geometric product
  return geometricProduct(m, geometricProduct(p, reverse(m)));
}

fn motorNorm(m: MultiVector) -> f32 {
  // The norm of a motor is square root of the sum of square of the terms:
  // we have
  var sum = 0.;
  sum += m.s * m.s;
  sum += m.exey * m.exey;
  sum += m.exez * m.exez;
  sum += m.eyez * m.eyez;
  sum += m.eoex * m.eoex;
  sum += m.eoey * m.eoey;
  sum += m.eoez * m.eoez;
  sum += m.exeyez * m.exeyez;
  sum += m.eoexey * m.eoexey;
  sum += m.eoexez * m.eoexez;
  sum += m.eoeyez * m.eoeyez;
  sum += m.ex * m.ex;
  sum += m.ey * m.ey;
  sum += m.ez * m.ez;
  sum += m.eo * m.eo;
  sum += m.eoexeyez * m.eoexeyez;
  return sqrt(sum);
}

fn createTranslator(d: vec3f) -> MultiVector {
  // Given dx and dy describing the moveming in the x and y directions,
  // the translator is given by 1 + dx/2 exeo + dy/2 eyeo + dz/2 ezeo
  // In code, we always store the coefficents of
  //    scalar, exey, exez, eyez, eoex, eoey, eoez, exeyez, eoexey, eoexez, eoeyez, ex, ey, ez, eo, eoexeyez
  // Hence the implementation is as below
  return MultiVector(1, 0, 0, 0, -d.x / 2, -d.y / 2, -d.z / 2, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

fn extractTranslator(m: MultiVector) -> MultiVector {
  // Given a general motor, we can extract the translator part
  return MultiVector(1, 0, 0, 0, m.eoex, m.eoey, m.eoez, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

fn createDir(d: vec3f) -> MultiVector {
  // A direction is given by dx eyez + dy ezex + dz exey
  //    scalar, exey, exez, eyez, eoex, eoey, eoez, exeyez, eoexey, eoexez, eoeyez, ex, ey, ez, eo, eoexeyez
  return MultiVector(0, d.z, -d.y, d.x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

fn createLine(s: vec3f, d: vec3f) -> MultiVector {
  // A line is given by a starting point (sx, sy, sz) and a direction (dx, dy, dz)
  //  in this form: dx eyez + dy ezex + dz exey + (dy * sz - dz * sy) exeo + (dz * sx - dx * sz) eyeo + (dx * sy - dy * sx) ezeo
  let n = createDir(d); // represent the input direction in PGA
  let dir = normalizeMotor(n); // normalize the direction to make sure it is a unit direction
  // Note dir.exey = dz, dir.exez = -dy, dir.eyez = dx
  return MultiVector(0, dir.exey, dir.exez, dir.eyez, -(-dir.exez * s.z - dir.exey * s.y), -(dir.exey * s.x - dir.eyez * s.z), -(dir.eyez * s.y + dir.exez * s.x), 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

fn createRotor(angle: f32, d: vec3f, spt: vec3f) -> MultiVector {
  // Given an angle and a rotation axis direction (dx, dy, dz) and a start point of the rotation axis,
  // the rotor is given by cos(angle / 2 ) + sin(angle / 2 ) L
  //  where L is the line in 3D PGA formed by the direction and the start point
  let c = cos(angle / 2);
  let s = sin(angle / 2);
  let L = createLine(spt, d);
  return MultiVector(c, s * L.exey, s * L.exez, s * L.eyez, s * L.eoex, s * L.eoey, s * L.eoez, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

fn extractRotor(m: MultiVector) -> MultiVector {
  // Given a general motor, we can extract the rotor part
  return MultiVector(m.s, m.exey, m.exez, m.eyez, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

fn createPoint(p: vec3f) -> MultiVector {
  // Given a point in 3D with coordinates (x, y, z)
  // A point in PGA is given by exeyez + x eoezey + y eoexez + z eoeyex
  // In code, we always store the coefficents of 
  //    scalar, exey, exez, eyez, eoex, eoey, eoez, exeyez, eoexey, eoexez, eoeyez, ex, ey, ez, eo, eoexeyez
  return MultiVector(0, 0, 0, 0, 0, 0, 0, 1, -p.z, p.y, -p.x, 0, 0, 0, 0, 0);
}

fn extractPoint(p: MultiVector) -> vec3f {
  // to extract the 3d point from a exeyez + b eoezey + c eoexez + d eoeyex
  // we have x = -b/a and y = c/a and z = -d/a
  return vec3f(-p.eoeyez / p.exeyez, p.eoexez / p.exeyez, -p.eoexey / p.exeyez);
}

fn createPlane(n: vec3f, d: f32) -> MultiVector {
  // Given a plane in 3D with normal (nx, ny, nz) and distance from the origin d
  // A plane in PGA is given by nx ex + ny ey + nz ez - deo
  // In code, we always store the coefficents of 
  //    scalar, exey, exez, eyez, eoex, eoey, eoez, exeyez, eoexey, eoexez, eoeyez, ex, ey, ez, eo, eoexeyez
  return MultiVector(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, n.x, n.y, n.z, -d, 0);
}

fn createPlaneFromPoints(p1: vec3f, p2: vec3f, p3: vec3f) -> MultiVector {
  // Given three poitns (x1, y1, z1), (x2, y2, z2), (x3, y3, z3)
  // A plane in PGA is given by 
  //        ((y2 * z3 - y3 * z2) -      (y1 * z3 - y3 * z1) +      (y1 * z2 - y2 * z1)) ex 
  // -      ((x2 * z3 - x3 * z2) -      (x1 * z3 - x3 * z1) +      (x1 * z2 - x2 * z1)) ey 
  // +      ((x2 * y3 - x3 * y2) -      (x1 * y3 - x3 * y1) +      (x1 * y2 - x2 * y1)) ez 
  // + (x1 * (y2 * z3 - y3 * z2) - x2 * (y1 * z3 - y3 * z1) + x3 * (y1 * z2 - y2 * z1)) eo
  let nx =          (p2[1] * p3[2] - p3[1] * p2[2]) -         (p1[1] * p3[2] - p3[1] * p1[2]) +         (p1[1] * p2[2] - p2[1] * p1[2]);
  let ny =          (p2[0] * p3[2] - p3[0] * p2[2]) -         (p1[0] * p3[2] - p3[0] * p1[2]) +         (p1[0] * p2[2] - p2[0] * p1[2]);
  let nz =          (p2[0] * p3[1] - p3[0] * p2[1]) -         (p1[0] * p3[1] - p3[0] * p1[1]) +         (p1[0] * p2[1] - p2[0] * p1[1]);
  let d = (p1[0] * (p2[1] * p3[2] - p3[1] * p2[2]) - p2[0] * (p1[1] * p3[2] - p3[1] * p1[2]) + p3[0] * (p1[1] * p2[2] - p2[1] * p1[2]));
  return createPlane(vec3f(nx, -ny, nz), d);
}


fn normalizeMotor(m: MultiVector) -> MultiVector {
  // To normalize a motor, we divide each coefficient by its norm
  let mnorm = motorNorm(m);
  if (mnorm == 0.0) {
    return MultiVector(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  }
  return MultiVector(m.s / mnorm, m.exey / mnorm, m.exez / mnorm, m.eyez / mnorm, m.eoex / mnorm, m.eoey / mnorm, m.eoez / mnorm, m.exeyez / mnorm, m.eoexey / mnorm, m.eoexez / mnorm, m.eoeyez / mnorm, m.ex / mnorm, m.ey / mnorm, m.ez / mnorm, m.eo / mnorm, m.eoexeyez / mnorm);
}

fn applyMotorToPoint(p: vec3f, m: MultiVector) -> vec3f {
  // apply the motor m to transform the point p
  // this code covert the 3d point p into PGA and apply the motor to transform it
  // then extra the result from PGA
  let new_p = applyMotor(createPoint(p), m);
  return extractPoint(new_p);
};

fn applyMotorToDir(d: vec3f, m: MultiVector) -> vec3f {
  // apply the motor m to transform the direction d
  // this code convert the 3d direction d into PGA, then extract the rotor from the motor
  // and transform the direction using the rotor
  // last, extra the result from PGA
  let r = extractRotor(m);
  let new_d = applyMotor(createPoint(d), r);
  return extractPoint(new_d);
}

// a function to transform the direction to the model coordiantes
// fn transformDir(d: vec3f) -> vec3f {
//   // transform the direction using the camera pose
//   var out = applyMotorToDir(d, cameraPose.motor);
//   // transform it using the object pose
//   out = applyMotorToDir(out, reverse(box.motor));
//   out /= box.scale.xyz;
//   return out;
// }

// a function to transform the start pt to the model coordiantes
// fn transformPt(pt: vec3f) -> vec3f {
//   // transform the point using the camera pose
//   var out = applyMotorToPoint(pt, cameraPose.motor);
//   // transform it using the object pose
//   out = applyMotorToPoint(out, reverse(box.motor));
//   out /= box.scale.xyz;
//   return out;
// }