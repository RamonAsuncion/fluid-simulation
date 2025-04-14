/*!
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

import SceneObject from '/lib/DSViz/SceneObject.js'

export default class ParticleSystemObject extends SceneObject {
  constructor(device, canvasFormat, numParticles = 4096 * 4) {
    super(device, canvasFormat);
    this._numParticles = numParticles;
    this._step = 0;
  }
  
  async createGeometry() { 
    await this.createParticleGeometry();
  }
  
  async createParticleGeometry() {
    // Create particles
    this._particles = new Float32Array(this._numParticles * 10); // [x, y, ix, iy, vx, vy, ivx, ivy, l1, l2]
    // TODO 1 - create ping-pong buffers to store and update the particles in GPU
    // name the ping-pong buffers _particleBuffers
    
    this._particleBuffers = [
      this._device.createBuffer({
        label: "particle buffer 1 " + this.getName(),
        size: this._particles.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      }),
      this._device.createBuffer({
        label: "particle buffer 2" + this.getName(),
        size: this._particles.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      })
    ];
    // calling the resetParticles function to reset the particle buffers
    this.resetParticles();

    this._time = new Uint32Array([0]);
    this._timeBuffer = this._device.createBuffer({
      label: "particle buffer 1 " + this.getName(),
        size: this._time.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      })
    this._device.queue.writeBuffer(this._timeBuffer, 0, this._time);
  }
    
  resetParticles() {
    for (let i = 0; i < this._numParticles; ++i) {
      // random position between [-1, 1] x [-1, 1]
      this._particles[10 * i + 0] = (Math.random() * 2 - 1); // [-1, 1] 
      this._particles[10 * i + 1] = (Math.random() * 2 - 1);
      // store the initial positions
      this._particles[10 * i + 2] = this._particles[10 * i + 0];
      this._particles[10 * i + 3] = this._particles[10 * i + 1];
      // TODO 10: update the velocity
      this._particles[10 * i + 4] = (Math.random() * 2 - 1) * .001; //[-.001, .001]
      this._particles[10 * i + 5] = (Math.random() * 2 - 1) * .001;
      //store initial velocity
      this._particles[10 * i + 6] = this._particles[10 * i + 4];
      this._particles[10 * i + 7] = this._particles[10 * i + 5];
      //store lifetime
      let rand = (Math.random()) * 300 + 300;
      this._particles[10 * i + 8] = rand;
      this._particles[10 * i + 9] = rand;
    }
    // Copy from CPU to GPU
    this._step = 0;
    this._device.queue.writeBuffer(this._particleBuffers[this._step % 2], 0, this._particles);
  }

  updateTimeBuffer(){
    const size = 500;
    if (this._step < size){
      this._device.queue.writeBuffer(this._timeBuffer, 0, new Uint32Array([1]));
    }
    else if (this._step < size * 2){
      this._device.queue.writeBuffer(this._timeBuffer, 0, new Uint32Array([0]));
    }
    else if (this._step < size * 3){
      this._device.queue.writeBuffer(this._timeBuffer, 0, new Uint32Array([-1]));
    }
    else{
      this._step = 0;
    }
  }
  
  updateGeometry() { }
  
  async createShaders() {
    let shaderCode = await this.loadShader("/shaders/particles.wgsl");
    this._shaderModule = this._device.createShaderModule({
      label: "Particles Shader " + this.getName(),
      code: shaderCode,
    });
    // TODO 2 - Create the bind group layout for using the ping-pong buffers in the GPU
    // name the bind group layout _bindGroupLayout
    this._bindGroupLayout = this._device.createBindGroupLayout({
      label: "Grid Bind Group Layout " + this.getName(),
      entries: [{
        binding: 0,
        visibility: GPUShaderStage.VERTEX | GPUShaderStage.COMPUTE,
        buffer: { type: "read-only-storage"} // particle input buffer
      }, {
        binding: 1,
        visibility: GPUShaderStage.COMPUTE,
        buffer: { type: "storage"} // particle output buffer
      }, {
        binding: 2,
        visibility: GPUShaderStage.COMPUTE,
        buffer: {type: "storage"}
      }]
    });
    
    
    
    // create the pipeline layout using the bind group layout
    this._pipelineLayout = this._device.createPipelineLayout({
      label: "Particles Pipeline Layout",
      bindGroupLayouts: [ this._bindGroupLayout ],
    });
  }
  
  async createRenderPipeline() { 
    await this.createParticlePipeline();
  }
  
  async createParticlePipeline() {
    this._particlePipeline = this._device.createRenderPipeline({
      label: "Particles Render Pipeline " + this.getName(),
      layout: this._pipelineLayout,
      vertex: {
        module: this._shaderModule, 
        entryPoint: "vertexMain",
      },
      fragment: {
        module: this._shaderModule,
        entryPoint: "fragmentMain",
        targets: [{
          format: this._canvasFormat
        }]
      },
      primitives: {
        typology: 'line-strip'
      }
    }); 
    // Create bind group to bind the particle buffers
    this._bindGroups = [
      this._device.createBindGroup({
        layout: this._particlePipeline.getBindGroupLayout(0),
        entries: [
          {
            binding: 0,
            resource: { buffer: this._particleBuffers[0] }
          },
          {
            binding: 1,
            resource: { buffer: this._particleBuffers[1] }
          },
          {
            binding: 2,
            resource: { buffer: this._timeBuffer }
          }
        ],
      }),
      this._device.createBindGroup({
        layout: this._particlePipeline.getBindGroupLayout(0),
        entries: [
          {
            binding: 0,
            resource: { buffer: this._particleBuffers[1] }
          },
          {
            binding: 1,
            resource: { buffer: this._particleBuffers[0] }
          },
          {
            binding: 2,
            resource: { buffer: this._timeBuffer }
          }
        ],
      })
    ];
  }
  
  render(pass) { 
    pass.setPipeline(this._particlePipeline); 
    pass.setBindGroup(0, this._bindGroups[this._step % 2]);
    pass.draw(128, this._numParticles);
  }
  
  async createComputePipeline() { 
    this._computePipeline = this._device.createComputePipeline({
      label: "Particles Compute Pipeline " + this.getName(),
      layout: this._pipelineLayout,
      compute: {
        module: this._shaderModule,
        entryPoint: "computeMain",
      }
    });
  }
  
  compute(pass) { 
    pass.setPipeline(this._computePipeline);
    pass.setBindGroup(0, this._bindGroups[this._step % 2]);
    pass.dispatchWorkgroups(Math.ceil(this._numParticles / 256));
    ++this._step
    this.updateTimeBuffer();
    console.log("updating time: " + this._step);
  }

  mouseInteraction(x, y, attract){
    const radius = .2;
    const force = .007 * attract;
    //brute force way, find every single particle that is within the radius from the x,y coords
    var dist = 0;
    for(var i = 0; i < this._particles.length; i++){
     let ox = this._particles[10 * i] - x;
     let oy = this._particles[10 * i + 1] - y;
     dist = Math.sqrt(ox *  ox + oy * oy);

     if(dist < radius){
      //calcuate px,py in mouse click's coordinate as origin
      let angle = Math.atan(oy / ox);
      if (ox < 0){ //arctan only does left side of graph, need to correct if on the other side
        angle += Math.PI;
      }
      
      //add the effect to the particle's velocity
      this._particles[10 * i + 4] += force * Math.cos(angle) / dist;
      this._particles[10 * i + 5] += force * Math.sin(angle) / dist;
      //write the new values to the buffers
      this._device.queue.writeBuffer(this._particleBuffers[0], (10 * i + 4) * Float32Array.BYTES_PER_ELEMENT, new Float32Array([force * Math.cos(angle) * radius / dist, force * Math.sin(angle) * radius / dist]));
     }
    }
  }
}