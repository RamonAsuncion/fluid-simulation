struct Particle {
  pos : vec3f,
  initPos : vec3f,
  vel : vec3f,
  initVel : vec3f,
  lifeTime : vec2f,
  dummy : vec2f,
};

struct MouseInteraction {
  position: vec2f,
  isDown: f32,
  radius: f32
};

struct Camera {
  motor: MultiVector,
  focal: vec2f,
  res: vec2f,
}

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

const use_binding_box = 1;
const fallback_left = -1.0;
const fallback_right = 1.0;
const fallback_top = 1.0;
const fallback_bottom = -1.0;
const fallback_front = 1.0;
const fallback_back = 2.0;

//grid parameters
const cell_size = .04;
const max_per_cell = 512;

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

@vertex
fn ballBasedVertex(@builtin(instance_index) idx: u32, @builtin(vertex_index) vIdx: u32) -> @builtin(position) vec4f {
  let particle = particlesIn[idx];
  let r = 0.01;
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

  let x = r * sin(phi) * sin(theta);
  let y = r * sin(phi) * cos(theta);
  let z = r * cos(phi);

  return vec4f(x + particle.pos.x, y + particle.pos.y, z + particle.pos.z, 1);
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

@compute @workgroup_size(256)
fn computeMain(@builtin(global_invocation_id) global_id: vec3u) {
  let idx = global_id.x;
  
  if (idx < arrayLength(&particlesIn)) {
    // particlesOut[idx] = particlesIn[idx];
    
    let particle = particlesIn[idx];
    //f = ma
    let forces = calculateForces(idx);
    let accel = forces / mass;
    var newVel = particle.vel + accel * getSimulationSpeed();
    // * steps_per_update;
    
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
    var newPos = particle.pos + newVel * getSimulationSpeed();
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
    particlesOut[idx].pos = newPos;
    particlesOut[idx].vel = newVel;
  }
}

//find all the forces that should be applied to a given particle, and the net direction of these forces
fn calculateForces(idx: u32) -> vec3f{
  
  var particle = particlesIn[idx];
  var forces = vec3f(0, 0, 0);
  //apply gravity
  forces += vec3f(0, -9.81, 0) * gravity_multiplier;

  // calculate pressure force
  forces -= pressureApproximation(particle.pos, idx) * pressure_multiplier; //negative for some reason...
  
  //calculate viscosity force
  forces += viscosityApproximation(idx) * viscosity_multiplier;

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
    var distance = length(position - particlesIn[i].pos);
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
          var distance = length(position - particle2.pos);
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
  if (use_acceleration == 0){
    return rawPressureCalculation(position, particle_idx);
  }
  return acceleratedPressureCalculation(position, particle_idx);
}

fn rawPressureCalculation(position: vec3f, particle_idx: u32) -> vec3f{
  var pressure_force = vec3f(0, 0, 0);
  var current_density = densityBuffer[particle_idx];
  var particle2 : Particle;
  for (var i = 0; i < i32(arrayLength(&particlesIn)); i++){
      pressure_force += pressureFromPoint(particle_idx, i);
  }
  return pressure_force;
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
}

fn pressureFromPoint(particle1Idx: u32, particle2Idx: i32) -> vec3f{
  var temp = particlesIn[particle1Idx].pos - particlesIn[particle2Idx].pos;
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
    var distance = length(particle.pos - particle2.pos);
    var influence = viscositySmoothFunction(smoothing_radius, distance);
    viscosity_force += (particle2.vel - particle.vel) * influence;
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

        var cell_start = (adjacent_y * getGridLength() + adjacent_x) * max_per_cell; //recall that the first index of each cell is used as an index
        for (var i = cell_start + 1; i < cell_start + max_per_cell; i++){
          var pIdx = gridIn[i]; // the particle index
          if (pIdx == -1) {  //check if no more particles in this cell
            break;
          }
          var particle2 = particlesIn[pIdx];
          var distance = length(particle.pos - particle2.pos);
          var influence = viscositySmoothFunction(smoothing_radius, distance);
          viscosity_force += (particle2.vel - particle.vel) * influence; 
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
    for (var i = 0; i < max_per_cell; i++){
      atomicStore(&gridOut[global_id.x * max_per_cell + u32(i)], i32(-1));
    }
  }
}

@compute @workgroup_size(256)
fn computeGridStructure(@builtin(global_invocation_id) global_id: vec3u){
  //write the new grid positions
  if (global_id.x < arrayLength(&particlesOut)){
    var x = i32(floor((particlesOut[global_id.x].pos.x - getLeftBound()) / cell_size));
    var y = i32(floor((particlesOut[global_id.x].pos.y - getBottomBound()) / cell_size));
    if (x < 0 || x >= getGridLength() || y < 0 || y >= getGridHeight()) {
      return;
    }
    
    var cell_start = (y * getGridLength() + x) * max_per_cell;
    
    var idx = atomicAdd(&gridOut[cell_start], 1);
    if (idx >= max_per_cell - 2) {
      return;
    }
    
    atomicStore(&gridOut[cell_start + idx + 2], i32(global_id.x));
  }
}

@compute @workgroup_size(256)
fn computeDensity(@builtin(global_invocation_id) global_id: vec3u){
  if (global_id.x < arrayLength(&densityBuffer)){
    densityBuffer[global_id.x] = densityApproximation(particlesIn[global_id.x].pos);
  }
}

fn getGridIndex(x: i32, y:i32, z:i32) -> i32{
  return (z * getGridLength() * getGridHeight() + y * getGridLength() + x) * max_per_cell; 
}