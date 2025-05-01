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
  pos : vec3f,
  initPos : vec3f,
  vel : vec3f,
  initVel : vec3f,
  lifeTime : vec2f,
  dummy : vec2f,
};

// TODO 4: Write the bind group spells here using array<Particle>
@group(0) @binding(0) var<storage> particlesIn: array<Particle>;
@group(0) @binding(1) var<storage, read_write> particlesOut: array<Particle>;
@group(0) @binding(2) var<storage, read_write> timeBuffer: array<u32>;
// name the binded variables particlesIn and particlesOut


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
  let size = 0.0125 / 2f;
  let pi = 3.14159265;
  let theta = 2. * pi / 8 * f32(vIdx);
  let x = cos(theta) * size;
  let y = sin(theta) * size;
  return vec4f(x + particle.pos.x, y + particle.pos.y, particle.pos.z, 1);
}

@fragment
fn fragmentMain(@builtin(position) p: vec4f) -> @location(0) vec4f {
  return vec4f(p.z, 0, 0, 1); // (R, G, B, A)
  // return p;
}

@compute @workgroup_size(256)
fn computeMain(@builtin(global_invocation_id) global_id: vec3u) {
  // TODO 6: Revise the compute shader to update the particles using the velocity
  let idx = global_id.x;
  const maxVel = .1;
  if (idx < arrayLength(&particlesIn)) {
    // particlesOut[idx] = particlesIn[idx];
    
    // // // TOOD 7: Add boundary checking and respawn the particle when it is offscreen
    // let particle = particlesIn[idx];

    // let g = vec3f(0, -.000001, 0);
    // let accel = g;
    // var newVel = particle.vel + accel;

    // var newPos = particle.pos + particle.vel; 
    
    
    // //wrap around the screen if they go off. They will start slightly off screen but velocity should carry them back in.
    // if (newPos.x < -1 || newPos.x > 1){
    //   newPos.x *= -1;
    // }
    // if (newPos.y < -1 || newPos.y > 1){
    //   newPos.y *= -1;
    // }
    // particlesOut[idx].lifeTime[0] = particle.lifeTime[0] - 1;
    // if (particle.lifeTime[0] <= 0){
    //   newPos = particle.initPos;
    //   newVel = particle.initVel;
    //   particlesOut[idx].lifeTime[0] = particlesOut[idx].lifeTime[1];
    // }
    // particlesOut[idx].pos = newPos;
    // particlesOut[idx].vel = newVel;

    particlesOut[idx].pos = particlesIn[idx].pos;
    particlesOut[idx].vel = vec3f(0, 0, 0);
  }
}


fn generateWind(y: f32, frequency: f32, strength: f32) -> vec2f {
  let angle = sin(y * frequency) * 3.14159265 * f32(timeBuffer[0]);
  return vec2f(cos(angle), sin(angle)) * strength;
}