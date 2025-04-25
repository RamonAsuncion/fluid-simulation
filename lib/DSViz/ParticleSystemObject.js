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

import SceneObject from './SceneObject.js'
let VALUES_PER_PARTICLE = 10;

export default class ParticleSystemObject extends SceneObject {
  constructor(device, canvasFormat, numParticles = 4096 * 2) {
    super(device, canvasFormat);
    this._num_particles = numParticles;
    this._cell_size = .08;
    this._max_per_cell = 512;
    this._step = 0;
    this._left_bound = -1;
    this._right_bound = 1;
    this._top_bound = 1;
    this._bot_bound = -1;
    this._grid_length = Math.ceil((this._right_bound - this._left_bound) / this._cell_size);
    this._grid_height = Math.ceil((this._top_bound - this._bot_bound) / this._cell_size);
    this._num_grid_cells = this._grid_length * this._grid_height;
    this._initialized = false;
  }
  
  async createGeometry() { 
    await this.createParticleGeometry();
  }
  
  async createParticleGeometry() {
    // Create particles
    this._particles = new Float32Array(this._num_particles * VALUES_PER_PARTICLE); // [x, y, ix, iy, vx, vy, ivx, ivy, l1, l2]
    
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
    this._grid = new Int32Array(Array(this._num_grid_cells * this._max_per_cell).fill(-1));
    // calling the resetParticles function to reset the particle adn grid buffers
    this.resetParticles();
    await this.fill_grid();
    this._initialized = true;
    this._gridBuffers = [
      this._device.createBuffer({
        label: "grid buffer 1 " + this.getName(),
        size: this._grid.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      }),
      this._device.createBuffer({
        label: "grid buffer 2" + this.getName(),
        size: this._grid.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      })
    ];
    this._device.queue.writeBuffer(this._gridBuffers[0], 0, this._grid);
    
    this._time = new Uint32Array([0]);
    this._timeBuffer = this._device.createBuffer({
      label: "time buffer " + this.getName(),
        size: this._time.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      });
    this._device.queue.writeBuffer(this._timeBuffer, 0, this._time);
    
    this._density = new Float32Array(this._num_particles);
    this._densityBuffer = this._device.createBuffer({
      label: "density buffer " + this.getName(),
        size: this._density.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });
  }
    
  resetParticles() {
    for (let i = 0; i < this._num_particles; ++i) {
      // random position between [-1, 1] x [-1, 1]
      this._particles[VALUES_PER_PARTICLE * i + 0] = (Math.random() * 2 - 1); // [-1, 1] x coord
      this._particles[VALUES_PER_PARTICLE * i + 1] = (Math.random() * 2 - 1); // y coord
      // store the initial positions
      this._particles[VALUES_PER_PARTICLE * i + 2] = this._particles[VALUES_PER_PARTICLE * i + 0];
      this._particles[VALUES_PER_PARTICLE * i + 3] = this._particles[VALUES_PER_PARTICLE * i + 1];
      // TODO VALUES_PER_PARTICLE: update the velocity
      this._particles[VALUES_PER_PARTICLE * i + 4] = 0; //no initial velocity
      this._particles[VALUES_PER_PARTICLE * i + 5] = 0;
      //store initial velocity
      this._particles[VALUES_PER_PARTICLE * i + 6] = this._particles[VALUES_PER_PARTICLE * i + 4];
      this._particles[VALUES_PER_PARTICLE * i + 7] = this._particles[VALUES_PER_PARTICLE * i + 5];
      //store lifetime
      let rand = (Math.random()) * 300 + 300;
      this._particles[VALUES_PER_PARTICLE * i + 8] = rand;
      this._particles[VALUES_PER_PARTICLE * i + 9] = rand;
    }
    // Copy from CPU to GPU
    this._step = 0;
    this._device.queue.writeBuffer(this._particleBuffers[this._step % 2], 0, this._particles);
  }

  async fill_grid() {
    for (var i = 0; i < this._num_particles; i++){
      var x = Math.floor((this._particles[i * VALUES_PER_PARTICLE] - this._left_bound) / this._cell_size);
      var y = Math.floor((this._particles[i * VALUES_PER_PARTICLE + 1] - this._bot_bound) / this._cell_size);
      var cell_start = (y * this._grid_length + x) * this._max_per_cell;
      var stored = false;
      for (var index = cell_start; index < cell_start + this._max_per_cell; index++){
        if (this._grid[index] == -1){
          this._grid[index] = i;
          stored = true;
          console.log("stored at coord: ",x,y);
          break;
        }
      }
      if (!stored){
        console.log("failed to store particle " + (i) + ", maybe increase max_per_cell - start at ", cell_start, "x, y = ", x, y, this._cell_size);
      }
    }
    console.log("done storing values");
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
    let shaderCode = await this.loadShader("../../shaders/particles.wgsl");
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
        buffer: {type: "storage"} // time buffer
      }, {
        binding: 3,
        visibility: GPUShaderStage.VERTEX | GPUShaderStage.COMPUTE,
        buffer: { type: "read-only-storage"} // grid input buffer
      }, {
        binding: 4,
        visibility: GPUShaderStage.COMPUTE,
        buffer: { type: "storage"} // grid output buffer
      }, 
      {
        binding: 5,
        visibility: GPUShaderStage.COMPUTE,
        buffer: {type: "storage"} //density buffer
      }
    
    ]
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
          },
          {
            binding: 3,
            resource: { buffer: this._gridBuffers[0] }
          },
          {
            binding: 4,
            resource: { buffer: this._gridBuffers[1] }
          },
          {
            binding: 5,
            resource: { buffer: this._densityBuffer }
          }
        ]
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
          },
          {
            binding: 3,
            resource: { buffer: this._gridBuffers[1] }
          },
          {
            binding: 4,
            resource: { buffer: this._gridBuffers[0] }
          },
          {
            binding: 5,
            resource: { buffer: this._densityBuffer }
          },
        ],
      })
    ];
  }
  
  render(pass) { 
    if (this._initialized) {
      pass.setPipeline(this._particlePipeline); 
      pass.setBindGroup(0, this._bindGroups[this._step % 2]);
      pass.draw(128, this._num_particles);
    }
  }
  
  async createComputePipeline() { 
    this._densityPipeline = this._device.createComputePipeline({
      label: "Density Pipeline " + this.getName(),
      layout: this._pipelineLayout,
      compute: {
        module: this._shaderModule,
        entryPoint: "computeDensity",
      }
    });
    this._computePipeline = this._device.createComputePipeline({
      label: "Particles Compute Pipeline " + this.getName(),
      layout: this._pipelineLayout,
      compute: {
        module: this._shaderModule,
        entryPoint: "computeMain",
      }
    });
    this._clearGridPipeline = this._device.createComputePipeline({
      label: "Grid Compute Pipeline " + this.getName(),
      layout: this._pipelineLayout,
      compute: {
        module: this._shaderModule,
        entryPoint: "clearGridStructure",
      }
    });
    this._gridPipeline = this._device.createComputePipeline({
      label: "Grid Compute Pipeline " + this.getName(),
      layout: this._pipelineLayout,
      compute: {
        module: this._shaderModule,
        entryPoint: "computeGridStructure",
      }
    });
  }
  
  compute(pass) { 
    if (this._initialized) {
      pass.setBindGroup(0, this._bindGroups[this._step % 2]);

      pass.setPipeline(this._densityPipeline); //precompute densities
      pass.dispatchWorkgroups(Math.ceil(this._num_particles / 256));

      pass.setPipeline(this._computePipeline); //perform particle computation
      pass.dispatchWorkgroups(Math.ceil(this._num_particles / 256));

      pass.setPipeline(this._clearGridPipeline); //clear grid structure
      pass.dispatchWorkgroups(Math.ceil(this._num_grid_cells / 256));

      pass.setPipeline(this._gridPipeline); //write new particle values to grid structure for next time
      pass.dispatchWorkgroups(Math.ceil(this._num_particles / 256));

      ++this._step
      this.updateTimeBuffer();
    }
  }

  mouseInteraction(x, y, attract){
    const radius = .08;
    const force = 10 * attract;
    //brute force way, find every single particle that is within the radius from the x,y coords
    var dist = 0;
    for(var i = 0; i < this._particles.length; i++){
     let ox = this._particles[VALUES_PER_PARTICLE * i] - x;
     let oy = this._particles[VALUES_PER_PARTICLE * i + 1] - y;
     dist = Math.sqrt(ox *  ox + oy * oy);

     if(dist < radius){
      //calcuate px,py in mouse click's coordinate as origin
      let angle = Math.atan(oy / ox);
      if (ox < 0){ //arctan only does left side of graph, need to correct if on the other side
        angle += Math.PI;
      }
      
      //add the effect to the particle's velocity
      this._particles[VALUES_PER_PARTICLE * i + 4] += force * Math.cos(angle) / dist;
      this._particles[VALUES_PER_PARTICLE * i + 5] += force * Math.sin(angle) / dist;
      //write the new values to the buffers
      this._device.queue.writeBuffer(this._particleBuffers[this._step % 2], (VALUES_PER_PARTICLE * i + 4) * Float32Array.BYTES_PER_ELEMENT, new Float32Array([force * Math.cos(angle) * radius / dist, force * Math.sin(angle) * radius / dist]));
     }
    }
  }
}