struct Particle {
  pos : vec4f,
  initPos : vec4f,
  vel : vec4f,
  initVel : vec4f,
  lifeTime : vec2f,
  dummy : vec2f,
};

struct MouseInteraction {
  position: vec2f,
  isDown: f32,
  radius: f32,
  attractMode: f32 // 1 attract, 0 repel (todo: multiply with 1 or -1)
};


struct BoundaryBox {
  left: f32,
  right: f32,
  top: f32,
  bottom: f32,
  front: f32,
  back: f32
};

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

struct Camera {
  motor: MultiVector,
  focal: vec2f,
  res: vec2f,
}

struct SimulationConstants {
  pressureMultiplier: f32,
  gravityMultiplier: f32,
  max_per_cell: f32,
  reserved2: f32,
}

// TODO: use the mouse pose?

const use_binding_box = 1;
const fallback_left = -1.0;
const fallback_right = 1.0;
const fallback_top = 1.0;
const fallback_bottom = -1.0;
const fallback_front = 1.0;
const fallback_back = 2.0;

//grid parameters
const cell_size = .04;
//const max_per_cell = 256;

//fluid simulation parameters
const mass = 1;
const smoothing_radius = cell_size; //this way we only need to check adjacent grid cells, MUST change code to change this value
const smoothing_rate = 1;
const pi = 3.14159265;
const e = 2.71828;
const func_volume = (smoothing_rate + 1)/(2 * pi * smoothing_radius);
const target_density = .1;
const pressure_multiplier= 2;
const gravity_multiplier = .5;
const viscosity_multiplier = 5;
const max_vel = .1;
const max_force = .02;
const velocity_damping = .8; //multiplies final computed velocity, between 0.5 and 1 prolly
// const steps_per_update = .4; //lower value = slower/more stable simulation, somewhere between .1 and 1

@group(0) @binding(0) var<storage> particlesIn: array<Particle>;
@group(0) @binding(1) var<storage, read_write> particlesOut: array<Particle>;
@group(0) @binding(2) var<storage, read_write> timeBuffer: array<f32>;
@group(0) @binding(3) var<storage> gridIn: array<i32>;
@group(0) @binding(4) var<storage, read_write> gridOut: array<atomic<i32>>;
@group(0) @binding(5) var<uniform> boundaryBox: BoundaryBox;
@group(0) @binding(6) var<uniform> cameraPose: Camera;
@group(0) @binding(7) var<uniform> mouse: MouseInteraction;
@group(0) @binding(8) var<storage, read_write> densityBuffer: array<f32>;
@group(0) @binding(9) var<uniform> constants: SimulationConstants;


fn getSimulationSpeed() -> f32 {
  return timeBuffer[1] * 0.9 + 0.1; // 0->0.1 1->1.0
}

fn getLeftBound() -> f32 {
  if (use_binding_box != 0) {
    return boundaryBox.left;
  }
  return fallback_left;
}

fn getRightBound() -> f32 {
  if (use_binding_box != 0) {
    return boundaryBox.right;
  }
  return fallback_right;
}

fn getTopBound() -> f32 {
  if (use_binding_box != 0) {
    return boundaryBox.top;
  }
  return fallback_top;
}

fn getBottomBound() -> f32 {
  if (use_binding_box != 0) {
    return boundaryBox.bottom;
  }
  return fallback_bottom;
}

fn getFrontBound() -> f32 {
  if (use_binding_box != 0) {
    return boundaryBox.front;
  }
  return fallback_front;
}

fn getBackBound() -> f32 {
  if (use_binding_box != 0) {
    return boundaryBox.back;
  }
  return fallback_back;
}

fn getGridLength() -> i32 {
  return i32(ceil((getRightBound() - getLeftBound())/cell_size)) + 1;
}

fn getGridHeight() -> i32 {
  return i32(ceil((getTopBound() - getBottomBound())/cell_size)) + 1;
}

fn getGridWidth() -> i32 {
  return i32(ceil((getBackBound() - getFrontBound())/cell_size)) + 1;
}

fn transformDir(d: vec3f) -> vec3f {
  var out = applyMotorToDir(d, cameraPose.motor);
  return out;
}

fn transformPt(pt: vec3f) -> vec3f {
  var out = applyMotorToPoint(pt, cameraPose.motor);
  return out;
}

fn getBoundingBoxCenter() -> vec3f {
  return vec3f(
    (getLeftBound() + getRightBound()) * 0.5,
    (getBottomBound() + getTopBound()) * 0.5,
    (getFrontBound() + getBackBound()) * 0.5
  );
}

fn transformPointFromCenter(pt: vec3f) -> vec3f {
  let center = getBoundingBoxCenter();
  let centered = pt - center;
  let transformed = transformPt(centered);
  return transformed;
}

@vertex
fn ballBasedVertex(@builtin(instance_index) idx: u32, @builtin(vertex_index) vIdx: u32) -> @builtin(position) vec4f {
  let particle = particlesIn[idx];
  let r = .009;
  let pi = 3.14159265;
  let k = 6u;
  
  // a sphere in three dimensions can be parameterized by angles t and p, 0 <= t <= 2 pi. -pi/2 <= p <= pi/2
  
  //sweep a little bit of area of a sphere, some dt dp, will create a square. This square will be made into two triangles that share 2 vertices
  //idx refers to the particle index that we are computing. There are 6 vertices to compute per vIdx
  // vIdx / 6 refers to the angle offset we use. vIdx % 6 refers to which specific vertex we are drawing
  let angle_offset = u32(vIdx / 6);
  let vertex_number = i32(vIdx % 6);
  //I chose that each sphere will be broken into a 2k by k grid with (i,j) indices. let k = 3
  let i = (angle_offset) % (2 * k);
  let j = u32((angle_offset) / (2 * k));
  //there are 2k total divisions on theta. to get out angle out we do theta per division * number divisons
  let delta_theta = (2 * pi) / f32(2 * k - 1);
  var theta =  delta_theta * f32(i);
  //same for phi, k total divisions, j is the number of divisions we traveled
  let delta_phi = (pi) / f32(k - 1);
  var phi = -pi/2 + delta_phi * f32(j); 
  // 
  switch (vertex_number) {
    case 0: {
      //bot left for 1st triangle
      break;
    }
    case 1: {
      //bot right, need to move in theta direction
      theta += delta_theta;
      break;
    }
    case 2: {
      //top right for 1st triangle, both directions need to move
      theta += delta_theta;
      phi += delta_phi;
      break;
    }
    case 3: {
      //no offsets
      break;
    }
    case 4: {
      //top right for 2nd triangle, both need to move
      theta += delta_theta;
      phi += delta_phi;
      break;
    }
    case 5: {
      //top left, need phi change
      phi += delta_phi;
      break;
    }
    default: {
      break;
    }
  }

  let x = r * sin(phi) * sin(theta) + particle.pos.x;
  let y = r * sin(phi) * cos(theta) + particle.pos.y;
  let z = r * cos(phi) + particle.pos.z;
  
  // errr why why why why ljkfhsdalfjasdljflasdjfl;kjasd
  let center = getBoundingBoxCenter();
  let worldPos = vec3f(x, y, z);
  
  let relativePos = worldPos; // - center;
  let transformedPos = transformPt(relativePos);
  let finalPos = transformedPos + vec3f(0, 0, 2.0);
  
  var near = .01;
  var far = 100.;
  var projectedPt = vec2f(finalPos.x / finalPos.z, finalPos.y / finalPos.z);
  var depth = (finalPos.z - near) / (far - near);

  return vec4f(projectedPt, depth, 1);
}

@vertex
fn vertexMain(@builtin(instance_index) idx: u32, @builtin(vertex_index) vIdx: u32) -> @builtin(position) vec4f {
  let particle = particlesIn[idx];
  let size = 0.0250 / 2f;
  let theta = 2. * pi / 8 * f32(vIdx);
  let x = cos(theta) * size;
  let y = sin(theta) * size;
  let z = cos(theta) * size;
  return vec4f(x + particle.pos.x, y + particle.pos.y, z + particle.pos.z, 1);
}

@fragment
fn fragmentMain() -> @location(0) vec4f {
  // let particle = particlesIn[idx];
  // const maxVel = vec3f(.4f, .4f, .4f);
  // const maxMagnitude = sqrt(dot(maxVel, maxVel));
  // var particleMagnitude = sqrt(dot(particle.vel, particle.vel));
  // if (particleMagnitude > maxMagnitude) {
  //   particleMagnitude = maxMagnitude;
  // }
  // var r = (255.f*(particleMagnitude/maxMagnitude))/255.f;
  // var g = 0/255.f;
  // var b = 1 - (255.f*(particleMagnitude/maxMagnitude))/255.f;
  // var a = 1;
  // var color = vec4f(r, g, b, a);

  // return color; // (R, G, B, A)
  return vec4f(1, 0, 0, 1);
}

// https://youtu.be/rSKMYc1CQHE?feature=shared&t=1846
fn applyMouseForce(position: vec3f, velocity: vec3f) -> vec3f {
  var newVelocity = velocity;

  // button is pressed
  if (mouse.isDown > 0.5) {
    let mousePos = vec3f(mouse.position.x, mouse.position.y, 1); // x,y, z is fucking annoying
    let mouseRadius = mouse.radius; // interaction radius

    // dist for particle and mouse pos
    let distance = length(position - mousePos);

    if (distance < mouseRadius) {
      var direction: vec3f;
      var forceMagnitude: f32;

      if (mouse.attractMode > 0.5) {
        // attract: pull toward mouse
        direction = normalize(mousePos - position);

        // trying to fix clumping issues
        forceMagnitude = 0.05 * (1.0 - (distance / mouseRadius) * (distance / mouseRadius));

        if (distance < 0.02) {
          forceMagnitude = 0.0;
        }

      } else {
        // repel: push away from mouse
        direction = normalize(position - mousePos);

        // base force with linear falloff
        // increase the base force to create a wider
        // area of influence fot the mouse
        forceMagnitude = 0.08 * (1.0 - distance / mouseRadius);
      }

      // apply force by adding to velocity
      newVelocity += direction * forceMagnitude;
    }
  }

  return newVelocity;
}

@compute @workgroup_size(256)
fn computeMain(@builtin(global_invocation_id) global_id: vec3u) {
  let idx = global_id.x;
  
  if (idx < arrayLength(&particlesIn)) {
    // particlesOut[idx] = particlesIn[idx];
    
    let particle = particlesIn[idx];
    //f = ma
    let forces = calculateForces(idx);
    //let forces = vec3f(0, -0.1, 0);
    let accel = forces / mass;
    var newVel = particle.vel.xyz + accel * getSimulationSpeed();
    // * steps_per_update;

    
    
    //cap the velocity
    if (abs(newVel.x) > max_vel) {
      newVel.x = max_vel * (abs(newVel.x) / newVel.x);
    }
    if (abs(newVel.y) > max_vel) {
      newVel.y = max_vel * (abs(newVel.y) / newVel.y);
    }
    if (abs(newVel.z) > max_vel) {
      newVel.z = max_vel * (abs(newVel.z) / newVel.z);
    }
    
    //damp the velocity
    newVel *= velocity_damping;

    // update particle position
    var newPos = particle.pos.xyz + newVel * getSimulationSpeed();
    // * steps_per_update;
    
    //keep the particles in the bounding box, damp a bit from bouncing off the sides
    if (newPos.x < getLeftBound()){
      newPos.x = getLeftBound();
      newVel.x *= -.9;
    }
    else if (newPos.x > getRightBound()){
      newPos.x = getRightBound();
      newVel.x *= -.9;
    }
    if (newPos.y < getBottomBound()){
      newPos.y = getBottomBound();
      newVel.y *= -.9;
    }
    else if (newPos.y > getTopBound()){
      newPos.y = getTopBound();
      newVel.y *= -.9;
    }
    if (newPos.z < getFrontBound()){
      newPos.z = getBottomBound();
      newVel.z *= -.9;
    }
    else if (newPos.z > getBackBound()){
      newPos.z = getTopBound();
      newVel.z *= -.9;
    }
    
    particlesOut[idx].pos = vec4f(newPos, 0);
    particlesOut[idx].vel = vec4f(newVel, 0);
  }
}

//find all the forces that should be applied to a given particle, and the net direction of these forces
fn calculateForces(idx: u32) -> vec3f{
  
  var particle = particlesIn[idx];
  var forces = vec3f(0, 0, 0);
  //apply gravity
  forces += vec3f(0, 9.81, 0) * constants.gravityMultiplier;

  // // calculate pressure force
  forces -= pressureApproximation(particle.pos.xyz, idx) * constants.pressureMultiplier;
  
  // //calculate viscosity force
  // forces += viscosityApproximation(idx) * viscosity_multiplier;

  // forces += applyMouseForce(particle.pos.xyz, particle.vel.xyz);

  forces = max(min(forces, vec3f(max_force, max_force, max_force)), vec3f(-max_force, -max_force, -max_force));
  return forces;
}

fn densityApproximation(position: vec3f) -> f32{
  let use_acceleration = timeBuffer[2];
  if (use_acceleration == 0){
    return rawDensityCalculation(position);
  }
  return acceleratedDensityCalculation(position);
}

fn rawDensityCalculation(position : vec3f) -> f32{
  var density = 0.0;
  for (var i = 0; i < i32(arrayLength(&particlesIn)); i++){
    var distance = length(position - particlesIn[i].pos.xyz);
    //density = max(density + mass * smoothingFunction(smoothing_radius, distance), 1);
    density += mass * smoothingFunction(smoothing_radius, distance);
  }
  return density;
}

fn acceleratedDensityCalculation(position: vec3f) -> f32{
  var density = 0.0;
  var particle2 : Particle;
  //find the cell that the particle is in
  var x = i32(floor((position.x - getLeftBound()) / cell_size));
  var y = i32(floor((position.y - getBottomBound()) / cell_size));
  var z = i32(floor((position.z - getFrontBound()) / cell_size));

  for(var adjacent_x = x-1; adjacent_x < x+2; adjacent_x++){
    for (var adjacent_y = y-1; adjacent_y < y+2; adjacent_y++){ //iterate over adjacent cells from the cell_pos
      for (var adjacent_z = z-1; adjacent_z < z+2; adjacent_z++){

          //boundary checking
          if (adjacent_x < 0 || adjacent_x > getGridLength()){ continue; }
          if (adjacent_y < 0 || adjacent_y > getGridHeight()){ continue; }
          if (adjacent_z < 0 || adjacent_z > getGridWidth()){ continue; }

          //var cell_start = ((z * this._grid_height * this._grid_length) + (y * this._grid_length) + x) * this._max_per_cell;

          var cell_start = getGridIndex(x,y,z);
          var len = gridIn[cell_start] + 2; //recall that the first index of each cell is used as an index
          for (var i = cell_start + 1; i < cell_start + len; i++){
            var pIdx = gridIn[i]; // the particle index
            if (pIdx == -1) {  //check if no more particles in this cell
              break;
            }
          particle2 = particlesIn[pIdx];
          var distance = length(position - particle2.pos.xyz);
          density += mass * smoothingFunction(smoothing_radius, distance);
        }
      }
    }
  }
  return density;
}

//given a certain density, how hard should the fluid push(pressure)
//further away from target density, the higher the pressure
fn densityToPressure(density : f32) -> f32{
  return density - target_density;
}

fn pressureApproximation(position: vec3f, particle_idx: u32) -> vec3f{
  let use_acceleration = timeBuffer[2];
  //if (use_acceleration == 0){
    return rawPressureCalculation(position, particle_idx);
  //}
  //return acceleratedPressureCalculation(position, particle_idx);
}

fn rawPressureCalculation(position: vec3f, particle_idx: u32) -> vec3f{
  var pressure_force = vec3f(0, 0, 0);
  //var current_density = densityBuffer[particle_idx];
  var particle2 : Particle;
  for (var i = 0; i < i32(arrayLength(&particlesIn)); i++){
      pressure_force += pressureFromPoint(particle_idx, i);
  }
  return pressure_force;
  //return vec3f(0, 0.1, 0);
}

fn acceleratedPressureCalculation(position: vec3f, particle_idx: u32) -> vec3f{
  var pressure_force = vec3f(0, 0, 0);
  
  //find the cell that the particle is in
  var x = i32(floor((position.x - getLeftBound()) / cell_size));
  var y = i32(floor((position.y - getBottomBound()) / cell_size));
  var z = i32(floor((position.z - getFrontBound()) / cell_size));

  for(var adjacent_x = x-1; adjacent_x < x+2; adjacent_x++){
    if (adjacent_x < 0 || adjacent_x > getGridLength()){ continue; }
    for (var adjacent_y = y-1; adjacent_y < y+2; adjacent_y++){ //iterate over adjacent cells from the cell_pos
      //boundary checking
      if (adjacent_y < 0 || adjacent_y > getGridHeight()){ continue; }
      for (var adjacent_z = z-1; adjacent_z < z+2; adjacent_z++){
        
        if (adjacent_z < 0 || adjacent_z > getGridWidth()){ continue; }

        var cell_start = getGridIndex(x,y,z); //recall that the first index of each cell is used as an index
        var len = gridIn[cell_start] + 2;
        for (var i = cell_start + 1; i < cell_start + len; i++){
          var pIdx = gridIn[i]; // the particle index
          if (pIdx == -1) {  //check if no more particles in this cell
            break;
          }
          if (pIdx != i32(particle_idx)) {
            pressure_force += pressureFromPoint(particle_idx, pIdx);
          }
        }
      }
    }
  }
  return pressure_force;
  //return vec3f(0, 0.1, 0);
}

fn pressureFromPoint(particle1Idx: u32, particle2Idx: i32) -> vec3f{
  var temp = particlesIn[particle1Idx].pos.xyz - particlesIn[particle2Idx].pos.xyz;
  var distance = length(temp);
  var pressure_force = vec3f(0, 0, 0);
  if (distance > 0.00000001) {
    var direction = temp / distance;
    var current_density = densityBuffer[particle1Idx];
    var other_density = densityBuffer[particle2Idx];
    var slope = smoothingDerivative(smoothing_radius, distance);
    var shared_pressure = getSharedPressure(current_density, other_density);
    pressure_force = direction * slope * mass * shared_pressure;
    
    if (length(pressure_force) > 1000.f) {
      pressure_force *= (1000.f / length(pressure_force));
    }
    if (distance < smoothing_radius/2){ // if the particles are too close, push them harder apart
      pressure_force *= 10;
    }
    if (distance < smoothing_radius/4){
      pressure_force *= 10;
    }
    if (distance < smoothing_radius/6){
      pressure_force *= 10;
    }
  }
  return pressure_force;
}

fn getSharedPressure(density1: f32, density2: f32) -> f32{
  var pressure1 = densityToPressure(density1);
  var pressure2 = densityToPressure(density2);
  return (pressure1 + pressure2) / 2;
}

fn viscosityApproximation(idx: u32) -> vec3f{
  let use_acceleration = timeBuffer[2];
  if (use_acceleration == 0){
    return rawViscosityCalculation(idx);
  }
  return acceleratedViscosityCalculation(idx);
}

fn rawViscosityCalculation(idx : u32) -> vec3f{
  var viscosity_force = vec3f(0,0,0);
  var particle = particlesIn[idx];
  for (var i = 0; i < i32(arrayLength(&particlesIn)); i++){
    var particle2 = particlesIn[i];
    var distance = length(particle.pos.xyz - particle2.pos.xyz);
    var influence = viscositySmoothFunction(smoothing_radius, distance);
    viscosity_force += (particle2.vel.xyz - particle.vel.xyz) * influence;
  }
  return viscosity_force;
}

fn acceleratedViscosityCalculation(idx: u32) -> vec3f{
  var viscosity_force = vec3f(0,0,0);
  var particle = particlesIn[idx];
  var particle2 : Particle;
  //find the cell that the particle is in
  var x = i32(floor((particle.pos.x - getLeftBound()) / cell_size));
  var y = i32(floor((particle.pos.y - getBottomBound()) / cell_size));
  var z = i32(floor((particle.pos.z - getFrontBound()) / cell_size));
  for(var adjacent_x = x-1; adjacent_x < x+2; adjacent_x++){
    for (var adjacent_y = y-1; adjacent_y < y+2; adjacent_y++){ //iterate over adjacent cells from the cell_pos
      for (var adjacent_z = z-1; adjacent_z < z+2; adjacent_z++){
        
        
        //boundary checking
        if (adjacent_x < 0 || adjacent_x > getGridLength()){ continue; }
        if (adjacent_y < 0 || adjacent_y > getGridHeight()){ continue; }
        if (adjacent_z < 0 || adjacent_z > getGridWidth()){ continue; }

        var cell_start = (adjacent_y * getGridLength() + adjacent_x) * i32(constants.max_per_cell); //recall that the first index of each cell is used as an index
        for (var i = cell_start + 1; i < cell_start + i32(constants.max_per_cell); i++){
          var pIdx = gridIn[i]; // the particle index
          if (pIdx == -1) {  //check if no more particles in this cell
            break;
          }
          var particle2 = particlesIn[pIdx];
          var distance = length(particle.pos.xyz - particle2.pos.xyz);
          var influence = viscositySmoothFunction(smoothing_radius, distance);
          viscosity_force += (particle2.vel.xyz - particle.vel.xyz) * influence; 
        }
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

fn viscositySmoothFunction(smoothing_radius: f32, distance: f32) -> f32{
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
    for (var i = 0; i < i32(constants.max_per_cell); i++){
      atomicStore(&gridOut[global_id.x * u32(constants.max_per_cell) + u32(i)], i32(-1));
    }
  }
}

@compute @workgroup_size(256)
fn computeGridStructure(@builtin(global_invocation_id) global_id: vec3u){
  //write the new grid positions
  if (global_id.x < arrayLength(&particlesOut)){
    var x = i32(floor((particlesOut[global_id.x].pos.x - getLeftBound()) / cell_size));
    var y = i32(floor((particlesOut[global_id.x].pos.y - getBottomBound()) / cell_size));
    var z = i32(floor((particlesOut[global_id.x].pos.z - getFrontBound()) / cell_size));
    if (x < 0 || x >= getGridLength() || y < 0 || y >= getGridHeight() || z < 0 || z >= getGridWidth()) {
      return;
    }
    
    var cell_start = (z * getGridLength() * getGridHeight() + y * getGridLength() + x) * i32(constants.max_per_cell);
    
    var idx = atomicAdd(&gridOut[cell_start], 1);
    if (idx >= i32(constants.max_per_cell) - 2) {
      return;
    }
    
    atomicStore(&gridOut[cell_start + idx + 2], i32(global_id.x));
  }
}

@compute @workgroup_size(256)
fn computeDensity(@builtin(global_invocation_id) global_id: vec3u){
  if (global_id.x < arrayLength(&densityBuffer)){
    densityBuffer[global_id.x] = densityApproximation(particlesIn[global_id.x].pos.xyz);
  }
}

fn getGridIndex(x: i32, y:i32, z:i32) -> i32{
  return (z * getGridLength() * getGridHeight() + y * getGridLength() + x) * i32(constants.max_per_cell); 
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