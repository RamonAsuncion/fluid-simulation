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
  vel: vec2f,
  rg : vec2f,
  ba : vec2f,
  lifeTime : vec2f,
};

struct VertexOutput {
  @builtin(position) pos: vec4f,
  @location(0) color: vec4f // pass the color
};

@group(0) @binding(0) var<storage> particlesIn: array<Particle>;
@group(0) @binding(1) var<storage, read_write> particlesOut: array<Particle>;


@vertex
fn vertexMain(@builtin(instance_index) idx: u32, @builtin(vertex_index) vIdx: u32) -> VertexOutput {
  let particle = particlesIn[idx];
  let size = 0.0125 / 2f;
  let pi = 3.14159265;
  let theta = 2. * pi / 8 * f32(vIdx);
  let x = cos(theta) * size;
  let y = sin(theta) * size;
  var out: VertexOutput;
  out.pos = vec4f(vec2f(x + particle.pos.x, y + particle.pos.y), 0, 1);
  
  out.color = vec4f(particle.rg.x, particle.rg.y, particle.ba.x, particle.ba.y);
  return out;
}

@fragment
fn fragmentMain(@location(0) color: vec4f) -> @location(0) vec4f {
  return vec4f(color[0]/255, color[1]/255, color[2]/255, color[3]/255); // (R, G, B, A)
}


@compute @workgroup_size(256)
fn computeMain(@builtin(global_invocation_id) global_id: vec3u) {
  let idx = global_id.x;
  if (idx < arrayLength(&particlesIn)) {
    particlesOut[idx] = particlesIn[idx];
    let particle = particlesIn[idx];

    let g = vec2f(0, -.0002);
    var newVel = particle.vel + g;
    var newPos = particle.pos + particle.vel; 
    
    particlesOut[idx].lifeTime[0] = particle.lifeTime[0] - 1;

    if (particle.lifeTime[0] <= 0){
      //how do I kill this thing?
    }
    
    particlesOut[idx].pos = newPos;
    particlesOut[idx].vel = newVel;
  }
}