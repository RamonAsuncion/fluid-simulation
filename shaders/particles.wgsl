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

// Updated camera struct to use matrices instead of MultiVector
struct Camera {
  viewMatrix: mat4x4f,
  invViewMatrix: mat4x4f,
  cameraPosition: vec4f,
}

struct SimulationConstants {
  pressureMultiplier: f32,
  gravityMultiplier: f32,
  max_per_cell: f32,
  reserved2: f32,
}

struct Point {
  pos: vec3f,
  scalar: f32
}

struct GridCell {
  vertices: array<Point, 8>
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
@group(0) @binding(10)var<storage> trianglesIn: array<vec3f>;
@group(0) @binding(11)var<storage, read_write> trianglesOut: array<vec3f>;

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

fn transformPt(pt: vec3f) -> vec3f {
  var point4 = vec4f(pt, 1.0);
  // apply the view matrix transformation:
  return (cameraPose.viewMatrix * point4).xyz;
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

  let worldPos = vec3f(x, y, z);
  let viewPos = transformPt(worldPos);

  // perspective  projection
  var near = .01;
  var far = 100.;
  var aspect = 1.33;
  var fov = 0.8; // ~45 degrees
  var f = 1.0 / tan(fov / 2.0);
  var nf = 1.0 / (near - far);

  var projectedX = viewPos.x * f / aspect;
  var projectedY = viewPos.y * f;
  var depth = (viewPos.z * nf + near * far * nf);

  return vec4f(projectedX, projectedY, depth, -viewPos.z);
}

@vertex
fn surfaceBasedVertex(@builtin(instance_index) idx: u32, @builtin(vertex_index) vIdx: u32) -> @builtin(position) vec4f {
  var start_index = idx * 16;
  var specific_index = start_index + vIdx;

  var worldPos = trianglesIn[specific_index];

  // Check for sentinel value (-1, -1, -1)
  if (all(worldPos == vec3f(-1.0, -1.0, -1.0))) {
    // Push off-screen or behind camera
    return vec4f(0.0, 0.0, -10.0, 1.0);
  }

  let viewPos = transformPt(worldPos);

  // perspective projection
  var near = 0.01;
  var far = 100.0;
  var aspect = 1.33;
  var fov = 0.8;
  var f = 1.0 / tan(fov / 2.0);
  var nf = 1.0 / (near - far);

  var projectedX = viewPos.x * f / aspect;
  var projectedY = viewPos.y * f;
  var depth = (viewPos.z * nf + near * far * nf);

  return vec4f(projectedX, projectedY, depth, -viewPos.z);
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
  let r = timeBuffer[3];
  let g = timeBuffer[4];
  let b = timeBuffer[5];
  // return vec4f(1, 0, 0, 1);
  return vec4f(r,g,b,1);
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
    let particle = particlesIn[idx];
    let forces = calculateForces(idx);
    let accel = forces / mass;
    var newVel = particle.vel.xyz + accel * getSimulationSpeed();
    
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
      newPos.z = getFrontBound();
      newVel.z *= -.9;
    }
    else if (newPos.z > getBackBound()){
      newPos.z = getBackBound();
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
  return rawPressureCalculation(position, particle_idx);
}

fn rawPressureCalculation(position: vec3f, particle_idx: u32) -> vec3f{
  var pressure_force = vec3f(0, 0, 0);
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
  if (global_id.x < arrayLength(&gridOut)){
    for (var i = 0; i < i32(constants.max_per_cell); i++){
      atomicStore(&gridOut[global_id.x * u32(constants.max_per_cell) + u32(i)], i32(-1));
    }
  }
}

@compute @workgroup_size(256)
fn computeGridStructure(@builtin(global_invocation_id) global_id: vec3u){
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

@compute @workgroup_size(256)
fn computeCubeMarch(@builtin(global_invocation_id) global_id: vec3u){
  // scalar function for grid cubes is based on DENSITY of particles around each vertex
  // Isolevel determines how many particles need to be near a vertex to count it for the rendering of the isosurface
  // Need to find a way to get the density of particles for each vertex - distance from all particles to each corner
  // Also need to get vertices for each grid cell
  // Compare distances from all particles to each corner to isolevel
  if (global_id.x < arrayLength(&gridOut)){
    var start_index = global_id.x * u32(constants.max_per_cell);
    var num_particles = atomicLoad(&gridOut[start_index]);
    // get grid indices
    var grid_x = i32(global_id.x) % getGridLength();
    var grid_y = (i32(global_id.x) / getGridLength()) % getGridHeight();
    var grid_z = i32(global_id.x) / (getGridHeight() * getGridWidth());

    // get grid world positions for all 8 vertices - then scalar values
    var current_grid_cell: GridCell;
    var world_x = f32(grid_x) * cell_size - f32(getLeftBound());
    var world_y = f32(grid_y) * cell_size - f32(getBottomBound());
    var world_z = f32(grid_z) * cell_size - f32(getFrontBound());
    for (var i = 0; i < 8; i++) {
      // world position
      current_grid_cell.vertices[i].pos = vec3f(world_x, world_y, world_z);
      var remainder = i % 4;
      if (remainder == 1 || remainder == 2) {
        current_grid_cell.vertices[i].pos.x += cell_size;
      }
      if (remainder == 2 || remainder == 3) {
        current_grid_cell.vertices[i].pos.z += cell_size;
      }
      if (i < 4) {
        current_grid_cell.vertices[i].pos.y += cell_size;
      }

      // scalar value
      current_grid_cell.vertices[i].scalar = 0.0;
      for (var j = 1; j < num_particles; j++) {
        var particle_idx = atomicLoad(&gridOut[j]);
        var position = particlesIn[particle_idx].pos.xyz;
        var distance = length(current_grid_cell.vertices[i].pos - position);
        current_grid_cell.vertices[i].scalar += 1.0 / (distance * distance + 0.01);
      }
    }

    // create triangle information from grid cell
    polygonise(current_grid_cell, 2.0, global_id.x);
  }
}

fn getGridIndex(x: i32, y:i32, z:i32) -> i32{
  return (z * getGridLength() * getGridHeight() + y * getGridLength() + x) * i32(constants.max_per_cell); 
}

// Linearly interpolate the position where an isosurface cuts
// an edge between two vertices, each with their own scalar value
fn vertexInterp(isolevel: f32, p1: vec3f, p2: vec3f, valp1: f32, valp2: f32) -> vec3f {
    var p = vec3f(0, 0, 0);
    if (abs(isolevel - valp1) < 0.00001) {
        return p1;
    }
    if (abs(isolevel - valp2) < 0.00001) {
        return p2;
    }
    if (abs(valp1 - valp2) < 0.00001) {
        return p1;
    }
    let mu = (isolevel - valp1) / (valp2 - valp1);
    p.x = p1.x + mu * (p2.x - p1.x);
    p.y = p1.y + mu * (p2.y - p1.y);
    p.z = p1.z + mu * (p2.z - p1.z);

    return p;
}

fn polygonise(grid: GridCell, isolevel: f32, grid_index: u32) {
  var cubeindex = 0;
  var vertex_list = array<vec3f, 12>(vec3f(0, 0, 0), vec3f(0, 0, 0), vec3f(0, 0, 0), vec3f(0, 0, 0), vec3f(0, 0, 0), vec3f(0, 0, 0), vec3f(0, 0, 0), vec3f(0, 0, 0), vec3f(0, 0, 0), vec3f(0, 0, 0), vec3f(0, 0, 0), vec3f(0, 0, 0));

  // Determine the index into the edge table which
  // tells us which vertices are inside of the surface
  if (grid.vertices[0].scalar < isolevel) {
    cubeindex |= 1;
  }
  if (grid.vertices[1].scalar < isolevel) {
    cubeindex |= 2;
  }
  if (grid.vertices[2].scalar < isolevel) {
    cubeindex |= 4;
  }
  if (grid.vertices[3].scalar < isolevel) {
    cubeindex |= 8;
  }
  if (grid.vertices[4].scalar < isolevel) {
    cubeindex |= 16;
  }
  if (grid.vertices[5].scalar < isolevel) {
    cubeindex |= 32;
  }
  if (grid.vertices[6].scalar < isolevel) {
    cubeindex |= 64;
  }
  if (grid.vertices[7].scalar < isolevel) {
    cubeindex |= 128;
  }
  
  // Cube is entirely in/out of the surface
  if (edgeTable[cubeindex] == 0) {
    return;
  }

  // Find the vertices where the surface intersects the cube
  if ((edgeTable[cubeindex] & 1) != 0)
    {vertex_list[0] = vertexInterp(isolevel, grid.vertices[0].pos, grid.vertices[1].pos, grid.vertices[0].scalar, grid.vertices[1].scalar);}
  if ((edgeTable[cubeindex] & 2) != 0)
    {vertex_list[1] = vertexInterp(isolevel, grid.vertices[1].pos, grid.vertices[2].pos, grid.vertices[1].scalar, grid.vertices[2].scalar);}
  if ((edgeTable[cubeindex] & 4) != 0)
    {vertex_list[2] = vertexInterp(isolevel, grid.vertices[2].pos, grid.vertices[3].pos, grid.vertices[2].scalar, grid.vertices[3].scalar);}
  if ((edgeTable[cubeindex] & 8) != 0)
    {vertex_list[3] = vertexInterp(isolevel, grid.vertices[3].pos, grid.vertices[0].pos, grid.vertices[3].scalar, grid.vertices[0].scalar);}
  if ((edgeTable[cubeindex] & 16) != 0)
    {vertex_list[4] = vertexInterp(isolevel, grid.vertices[4].pos, grid.vertices[5].pos, grid.vertices[4].scalar, grid.vertices[5].scalar);}
  if ((edgeTable[cubeindex] & 32) != 0)
    {vertex_list[5] = vertexInterp(isolevel, grid.vertices[5].pos, grid.vertices[6].pos, grid.vertices[5].scalar, grid.vertices[6].scalar);}
  if ((edgeTable[cubeindex] & 64) != 0)
    {vertex_list[6] = vertexInterp(isolevel, grid.vertices[6].pos, grid.vertices[7].pos, grid.vertices[6].scalar, grid.vertices[7].scalar);}
  if ((edgeTable[cubeindex] & 128) != 0)
    {vertex_list[7] = vertexInterp(isolevel, grid.vertices[7].pos, grid.vertices[4].pos, grid.vertices[7].scalar, grid.vertices[4].scalar);}
  if ((edgeTable[cubeindex] & 256) != 0)
    {vertex_list[8] = vertexInterp(isolevel, grid.vertices[0].pos, grid.vertices[4].pos, grid.vertices[0].scalar, grid.vertices[4].scalar);}
  if ((edgeTable[cubeindex] & 512) != 0)
    {vertex_list[9] = vertexInterp(isolevel, grid.vertices[1].pos, grid.vertices[5].pos, grid.vertices[1].scalar, grid.vertices[5].scalar);}
  if ((edgeTable[cubeindex] & 1024) != 0)
    {vertex_list[10] = vertexInterp(isolevel, grid.vertices[2].pos, grid.vertices[6].pos, grid.vertices[2].scalar, grid.vertices[6].scalar);}
  if ((edgeTable[cubeindex] & 2048) != 0)
    {vertex_list[11] = vertexInterp(isolevel, grid.vertices[3].pos, grid.vertices[7].pos, grid.vertices[3].scalar, grid.vertices[7].scalar);}
  // Create the triangle
  var buffer_idx = grid_index * 16;
  for (var i = 0u; triTable[cubeindex][i] != -1 && i < 16; i++) {
    trianglesOut[buffer_idx + i] = vertex_list[triTable[cubeindex][i]];
  }
}

const edgeTable = array<u32, 256>(
  0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
  0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
  0x190, 0x99, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
  0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
  0x230, 0x339, 0x33, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
  0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
  0x3a0, 0x2a9, 0x1a3, 0xaa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
  0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
  0x460, 0x569, 0x663, 0x76a, 0x66, 0x16f, 0x265, 0x36c,
  0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
  0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff, 0x3f5, 0x2fc,
  0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
  0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55, 0x15c,
  0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
  0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc,
  0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
  0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
  0xcc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
  0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
  0x15c, 0x55, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
  0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
  0x2fc, 0x3f5, 0xff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
  0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
  0x36c, 0x265, 0x16f, 0x66, 0x76a, 0x663, 0x569, 0x460,
  0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
  0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa, 0x1a3, 0x2a9, 0x3a0,
  0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
  0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33, 0x339, 0x230,
  0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
  0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99, 0x190,
  0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
  0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
);

const triTable = array<array<i32, 16>, 256>(
  array<i32, 16>(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1),
  array<i32, 16>(8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1),
  array<i32, 16>(3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1),
  array<i32, 16>(4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1),
  array<i32, 16>(4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1),
  array<i32, 16>(9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1),
  array<i32, 16>(10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1),
  array<i32, 16>(5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1),
  array<i32, 16>(5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1),
  array<i32, 16>(8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1),
  array<i32, 16>(2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1),
  array<i32, 16>(2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1),
  array<i32, 16>(11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1),
  array<i32, 16>(5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1),
  array<i32, 16>(11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1),
  array<i32, 16>(11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1),
  array<i32, 16>(2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1),
  array<i32, 16>(6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1),
  array<i32, 16>(3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1),
  array<i32, 16>(6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1),
  array<i32, 16>(6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1),
  array<i32, 16>(8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1),
  array<i32, 16>(7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1),
  array<i32, 16>(3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1),
  array<i32, 16>(0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1),
  array<i32, 16>(9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1),
  array<i32, 16>(8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1),
  array<i32, 16>(5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1),
  array<i32, 16>(0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1),
  array<i32, 16>(6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1),
  array<i32, 16>(10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1),
  array<i32, 16>(1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1),
  array<i32, 16>(0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1),
  array<i32, 16>(3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1),
  array<i32, 16>(6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1),
  array<i32, 16>(9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1),
  array<i32, 16>(8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1),
  array<i32, 16>(3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1),
  array<i32, 16>(10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1),
  array<i32, 16>(10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1),
  array<i32, 16>(2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1),
  array<i32, 16>(7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1),
  array<i32, 16>(2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1),
  array<i32, 16>(1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1),
  array<i32, 16>(11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1),
  array<i32, 16>(8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1),
  array<i32, 16>(0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1),
  array<i32, 16>(7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1),
  array<i32, 16>(7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1),
  array<i32, 16>(10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1),
  array<i32, 16>(0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1),
  array<i32, 16>(7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1),
  array<i32, 16>(6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1),
  array<i32, 16>(4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1),
  array<i32, 16>(10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1),
  array<i32, 16>(8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1),
  array<i32, 16>(1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1),
  array<i32, 16>(10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1),
  array<i32, 16>(10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1),
  array<i32, 16>(9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1),
  array<i32, 16>(7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1),
  array<i32, 16>(3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1),
  array<i32, 16>(7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1),
  array<i32, 16>(3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1),
  array<i32, 16>(6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1),
  array<i32, 16>(9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1),
  array<i32, 16>(1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1),
  array<i32, 16>(4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1),
  array<i32, 16>(7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1),
  array<i32, 16>(6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1),
  array<i32, 16>(0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1),
  array<i32, 16>(6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1),
  array<i32, 16>(0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1),
  array<i32, 16>(11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1),
  array<i32, 16>(6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1),
  array<i32, 16>(5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1),
  array<i32, 16>(9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1),
  array<i32, 16>(1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1),
  array<i32, 16>(10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1),
  array<i32, 16>(0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1),
  array<i32, 16>(11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1),
  array<i32, 16>(9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1),
  array<i32, 16>(7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1),
  array<i32, 16>(2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1),
  array<i32, 16>(9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1),
  array<i32, 16>(9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1),
  array<i32, 16>(1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1),
  array<i32, 16>(0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1),
  array<i32, 16>(10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1),
  array<i32, 16>(2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1),
  array<i32, 16>(0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1),
  array<i32, 16>(0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1),
  array<i32, 16>(9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1),
  array<i32, 16>(5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1),
  array<i32, 16>(5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1),
  array<i32, 16>(8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1),
  array<i32, 16>(9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1),
  array<i32, 16>(1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1),
  array<i32, 16>(3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1),
  array<i32, 16>(4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1),
  array<i32, 16>(9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1),
  array<i32, 16>(11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1),
  array<i32, 16>(2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1),
  array<i32, 16>(9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1),
  array<i32, 16>(3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1),
  array<i32, 16>(1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1),
  array<i32, 16>(4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1),
  array<i32, 16>(0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1),
  array<i32, 16>(1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
  array<i32, 16>(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)
);