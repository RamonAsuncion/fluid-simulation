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

// TODO 4: Write the bind group spells here using array<Particle>
@group(0) @binding(0) var<storage> particlesIn: array<Particle>;
@group(0) @binding(1) var<storage, read_write> particlesOut: array<Particle>;
@group(0) @binding(2) var<storage, read_write> timeBuffer: array<u32>;
// name the binded variables particlesIn and particlesOut


@vertex
fn vertexMain(@builtin(instance_index) idx: u32, @builtin(vertex_index) vIdx: u32) -> @builtin(position) vec4f {
  // TODO 5: Revise the vertex shader to draw circle to visualize the particles
  let particle = particlesIn[idx];
  let size = 0.0125 / 2f;
  let pi = 3.14159265;
  let theta = 2. * pi / 8 * f32(vIdx);
  let x = cos(theta) * size;
  let y = sin(theta) * size;
  return vec4f(vec2f(x + particle.pos.x, y + particle.pos.y), 0, 1);
}

@fragment
fn fragmentMain() -> @location(0) vec4f {
  return vec4f(238.f/255, 118.f/255, 35.f/255, 1); // (R, G, B, A)
}

@compute @workgroup_size(256)
fn computeMain(@builtin(global_invocation_id) global_id: vec3u) {
  // TODO 6: Revise the compute shader to update the particles using the velocity
  let idx = global_id.x;
  const maxVel = .1;
  if (idx < arrayLength(&particlesIn)) {
    particlesOut[idx] = particlesIn[idx];
    
    // TOOD 7: Add boundary checking and respawn the particle when it is offscreen
    let particle = particlesIn[idx];

    let g = vec2f(0, -.000001);
    let wind = generateWind(f32(particle.pos.y), 1.5, 0.00003); 
    let accel = g + wind;
    var newVel = particle.vel + accel;
    if (newVel.x > maxVel || newVel.x < -maxVel){
      newVel.x = maxVel;
    }
    if (newVel.y > maxVel || newVel.y < -maxVel){
      newVel.y = maxVel;
    }

    var newPos = particle.pos + particle.vel; 
    
    
    //wrap around the screen if they go off. They will start slightly off screen but velocity should carry them back in.
    if (newPos.x < -1 || newPos.x > 1){
      newPos.x *= -1;
    }
    if (newPos.y < -1 || newPos.y > 1){
      newPos.y *= -1;
    }
    particlesOut[idx].lifeTime[0] = particle.lifeTime[0] - 1;
    if (particle.lifeTime[0] <= 0){
      newPos = particle.initPos;
      newVel = particle.initVel;
      particlesOut[idx].lifeTime[0] = particlesOut[idx].lifeTime[1];
    }
    //cap velocity
    
    particlesOut[idx].pos = newPos;
    particlesOut[idx].vel = newVel;
  }
}


fn generateWind(y: f32, frequency: f32, strength: f32) -> vec2f {
  let angle = sin(y * frequency) * 3.14159265 * f32(timeBuffer[0]);
  return vec2f(cos(angle), sin(angle)) * strength;
}