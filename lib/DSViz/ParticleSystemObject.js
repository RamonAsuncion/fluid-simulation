import SceneObject from "./SceneObject.js";
let VALUES_PER_PARTICLE = 10;

export default class ParticleSystemObject extends SceneObject {
  constructor(device, canvasFormat, camera = null, numParticles = 10000) {
    super(device, canvasFormat);
    this._num_particles = numParticles;
    this._cell_size = 0.04;
    this._max_per_cell = 512;
    this._step = 0;
    this._camera = camera;
    this._simulationSpeed = 1.0;
    this._use_acceleration = false;

    this._left_bound = -1;
    this._right_bound = 1;
    this._top_bound = 1;
    this._bot_bound = -1;

    this._grid_length = Math.ceil(
      (this._right_bound - this._left_bound) / this._cell_size
    );
    this._grid_height = Math.ceil(
      (this._top_bound - this._bot_bound) / this._cell_size
    );
    this._num_grid_cells = this._grid_length * this._grid_height;
    this._initialized = false;
    this._mousePosition = new Float32Array([0, 0, 0, 0.1]); // [x, y, isDown, radius]
  }

  async createGeometry() {
    this._boundaryBoxBuffer = this._device.createBuffer({
      label: "boundary box buffer",
      size: 24,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });

    const initialBoundaryBox = new Float32Array([
      this._left_bound, // left
      this._right_bound, // right
      this._top_bound, // top
      this._bot_bound, // bottom
      0.5, // front (for 3D later)
      -0.5, // back (for 3D later)
    ]);

    this._device.queue.writeBuffer(
      this._boundaryBoxBuffer,
      0,
      initialBoundaryBox
    );

    if (this._camera) {
      this._cameraBuffer = this._device.createBuffer({
        label: "Camera " + this.getName(),
        size: 80,
        usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
      });

      const cameraIdentity = new Float32Array(20);
      cameraIdentity[0] = 1;
      this._device.queue.writeBuffer(this._cameraBuffer, 0, cameraIdentity);
    }

    this._mouseBuffer = this._device.createBuffer({
      label: "Mouse interaction buffer",
      size: this._mousePosition.byteLength,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });

    this._device.queue.writeBuffer(this._mouseBuffer, 0, this._mousePosition);

    await this.createParticleGeometry();
  }

  setMousePosition(x, y) {
    this._mousePosition[0] = x;
    this._mousePosition[1] = y;
    this._device.queue.writeBuffer(this._mouseBuffer, 0, this._mousePosition);
  }

  setMouseDown(isDown) {
    this._mousePosition[2] = isDown ? 1.0 : 0.0;
    this._device.queue.writeBuffer(this._mouseBuffer, 0, this._mousePosition);
  }

  setMouseRadius(radius = 0.1) {
    this._mousePosition[3] = radius;
    this._device.queue.writeBuffer(this._mouseBuffer, 0, this._mousePosition);
  }

  async createParticleGeometry() {
    this._particles = new Float32Array(
      this._num_particles * VALUES_PER_PARTICLE
    );

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
      }),
    ];
    this._grid = new Int32Array(
      Array(this._num_grid_cells * this._max_per_cell).fill(-1)
    );
    // calling the resetParticles function to reset the particle and grid buffers
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
      }),
    ];
    this._device.queue.writeBuffer(this._gridBuffers[0], 0, this._grid);

    this._time = new Float32Array([
      0,
      this._simulationSpeed,
      this._use_acceleration,
    ]);
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

  async fill_grid() {
    for (var i = 0; i < this._num_particles; i++) {
      var x = Math.floor(
        (this._particles[i * VALUES_PER_PARTICLE] - this._left_bound) /
          this._cell_size
      );
      var y = Math.floor(
        (this._particles[i * VALUES_PER_PARTICLE + 1] - this._bot_bound) /
          this._cell_size
      );

      //bound checking
      if (x < 0 || x >= this._grid_length || y < 0 || y >= this._grid_height) {
        console.log(`Particle ${i} outside grid bounds: (${x}, ${y})`);
        continue;
      }

      var cell_start = (y * this._grid_length + x) * this._max_per_cell;
      var stored = false;
      for (
        var index = cell_start;
        index < cell_start + this._max_per_cell;
        index++
      ) {
        if (this._grid[index] == -1) {
          this._grid[index] = i;
          stored = true;
          break;
        }
      }
      if (!stored) {
        console.log(
          "failed to store particle " +
            i +
            ", maybe increase max_per_cell - start at ",
          cell_start,
          "x, y = ",
          x,
          y,
          this._cell_size
        );
      }
    }
    console.log("Grid filled with particles");
  }

  updateAccelerationMode() {
    this._use_acceleration = !this._use_acceleration;
    console.log(`Using acceleration ${this._use_acceleration}`);
  }

  resetParticles() {
    // reset the grid
    this._grid = new Int32Array(
      Array(this._num_grid_cells * this._max_per_cell).fill(-1)
    );

    for (let i = 0; i < this._num_particles; ++i) {
      this._particles[VALUES_PER_PARTICLE * i + 0] =
        this._left_bound +
        Math.random() * (this._right_bound - this._left_bound);

      // y pos
      // START THIS BITCH AT THE TOP!!!!!!!!!!!
      this._particles[VALUES_PER_PARTICLE * i + 1] =
        this._top_bound - Math.random() * 0.05;

      // initial stuff
      this._particles[VALUES_PER_PARTICLE * i + 2] =
        this._particles[VALUES_PER_PARTICLE * i + 0];
      this._particles[VALUES_PER_PARTICLE * i + 3] =
        this._particles[VALUES_PER_PARTICLE * i + 1];

      // no velocity
      this._particles[VALUES_PER_PARTICLE * i + 4] = 0;
      this._particles[VALUES_PER_PARTICLE * i + 5] = 0;

      // initial velocity
      this._particles[VALUES_PER_PARTICLE * i + 6] =
        this._particles[VALUES_PER_PARTICLE * i + 4];
      this._particles[VALUES_PER_PARTICLE * i + 7] =
        this._particles[VALUES_PER_PARTICLE * i + 5];

      // lifetime
      let rand = Math.random() * 300 + 300;
      this._particles[VALUES_PER_PARTICLE * i + 8] = rand;
      this._particles[VALUES_PER_PARTICLE * i + 9] = rand;
    }

    this._step = 0;
    this._device.queue.writeBuffer(
      this._particleBuffers[this._step % 2],
      0,
      this._particles
    );

    // reset density
    if (this._densityBuffer) {
      this._density = new Float32Array(this._num_particles);
      this._device.queue.writeBuffer(this._densityBuffer, 0, this._density);
    }

    console.log("Particles reset at top of boundary box");
  }

  updateTimeBuffer() {
    const size = 500;
    if (this._step < size) {
      this._device.queue.writeBuffer(this._timeBuffer, 0, new Uint32Array([1]));
    } else if (this._step < size * 2) {
      this._device.queue.writeBuffer(this._timeBuffer, 0, new Uint32Array([0]));
    } else if (this._step < size * 3) {
      this._device.queue.writeBuffer(
        this._timeBuffer,
        0,
        new Uint32Array([-1])
      );
    } else {
      this._step = 0;
    }

    if (this._simulationSpeed !== undefined) {
      const speedBuffer = new Float32Array([this._simulationSpeed]);
      this._device.queue.writeBuffer(this._timeBuffer, 4, speedBuffer);
    }

    if (this._use_acceleration !== undefined) {
      const accelerationModeBuffer = new Float32Array([
        this._use_acceleration ? 1 : 0,
      ]);
      this._device.queue.writeBuffer(
        this._timeBuffer,
        8,
        accelerationModeBuffer
      );
    }
  }

  updateGeometry() {}

  async createShaders() {
    let shaderCode = await this.loadShader("../../shaders/particles.wgsl");
    this._shaderModule = this._device.createShaderModule({
      label: "Particles Shader " + this.getName(),
      code: shaderCode,
    });
    let bindGroupEntries = [
      {
        binding: 0,
        visibility: GPUShaderStage.VERTEX | GPUShaderStage.COMPUTE,
        buffer: { type: "read-only-storage" }, // particle input buffer
      },
      {
        binding: 1,
        visibility: GPUShaderStage.COMPUTE,
        buffer: { type: "storage" }, // particle output buffer
      },
      {
        binding: 2,
        visibility: GPUShaderStage.COMPUTE,
        buffer: { type: "storage" }, // time buffer
      },
      {
        binding: 3,
        visibility: GPUShaderStage.VERTEX | GPUShaderStage.COMPUTE,
        buffer: { type: "read-only-storage" }, // grid input buffer
      },
      {
        binding: 4,
        visibility: GPUShaderStage.COMPUTE,
        buffer: { type: "storage" }, // grid output buffer
      },
    ];
    if (this._boundaryBoxBuffer) {
      bindGroupEntries.push({
        binding: 5,
        visibility: GPUShaderStage.VERTEX | GPUShaderStage.COMPUTE,
        buffer: { type: "uniform" }, // boundary box buffer
      });
    }
    if (this._cameraBuffer) {
      bindGroupEntries.push({
        binding: 6,
        visibility: GPUShaderStage.VERTEX | GPUShaderStage.COMPUTE,
        buffer: { type: "uniform" }, // camera buffer
      });
    }
    if (this._mouseBuffer) {
      bindGroupEntries.push({
        binding: 7,
        visibility: GPUShaderStage.VERTEX | GPUShaderStage.COMPUTE,
        buffer: { type: "uniform" }, // mouse buffer
      });
    }
    bindGroupEntries.push({
      binding: 8,
      visibility: GPUShaderStage.COMPUTE,
      buffer: { type: "storage" }, // density buffer
    });

    this._bindGroupLayout = this._device.createBindGroupLayout({
      label: "Grid Bind Group Layout " + this.getName(),
      entries: bindGroupEntries,
    });

    this._pipelineLayout = this._device.createPipelineLayout({
      label: "Particles Pipeline Layout",
      bindGroupLayouts: [this._bindGroupLayout],
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
        targets: [
          {
            format: this._canvasFormat,
          },
        ],
      },
      primitives: {
        topology: "line-strip",
      },
    });

    let bindGroupEntries1 = [
      {
        binding: 0,
        resource: { buffer: this._particleBuffers[0] },
      },
      {
        binding: 1,
        resource: { buffer: this._particleBuffers[1] },
      },
      {
        binding: 2,
        resource: { buffer: this._timeBuffer },
      },
      {
        binding: 3,
        resource: { buffer: this._gridBuffers[0] },
      },
      {
        binding: 4,
        resource: { buffer: this._gridBuffers[1] },
      },
    ];

    let bindGroupEntries2 = [
      {
        binding: 0,
        resource: { buffer: this._particleBuffers[1] },
      },
      {
        binding: 1,
        resource: { buffer: this._particleBuffers[0] },
      },
      {
        binding: 2,
        resource: { buffer: this._timeBuffer },
      },
      {
        binding: 3,
        resource: { buffer: this._gridBuffers[1] },
      },
      {
        binding: 4,
        resource: { buffer: this._gridBuffers[0] },
      },
    ];

    if (this._boundaryBoxBuffer) {
      bindGroupEntries1.push({
        binding: 5,
        resource: { buffer: this._boundaryBoxBuffer },
      });
      bindGroupEntries2.push({
        binding: 5,
        resource: { buffer: this._boundaryBoxBuffer },
      });
    }

    if (this._cameraBuffer) {
      bindGroupEntries1.push({
        binding: 6,
        resource: { buffer: this._cameraBuffer },
      });
      bindGroupEntries2.push({
        binding: 6,
        resource: { buffer: this._cameraBuffer },
      });
    }

    if (this._mouseBuffer) {
      bindGroupEntries1.push({
        binding: 7,
        resource: { buffer: this._mouseBuffer },
      });
      bindGroupEntries2.push({
        binding: 7,
        resource: { buffer: this._mouseBuffer },
      });
    }

    bindGroupEntries1.push({
      binding: 8,
      resource: { buffer: this._densityBuffer },
    });
    bindGroupEntries2.push({
      binding: 8,
      resource: { buffer: this._densityBuffer },
    });

    this._bindGroups = [
      this._device.createBindGroup({
        layout: this._particlePipeline.getBindGroupLayout(0),
        entries: bindGroupEntries1,
      }),
      this._device.createBindGroup({
        layout: this._particlePipeline.getBindGroupLayout(0),
        entries: bindGroupEntries2,
      }),
    ];
  }

  updateCameraPose(camera) {
    if (this._cameraBuffer && camera) {
      this._device.queue.writeBuffer(
        this._cameraBuffer,
        0,
        camera._pose ||
          new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
      );
      this._device.queue.writeBuffer(
        this._cameraBuffer,
        64,
        camera._focal || new Float32Array([1, 1])
      );
      this._device.queue.writeBuffer(
        this._cameraBuffer,
        72,
        camera._resolutions || new Float32Array([1, 1])
      );
    }
  }

  updateBoundaryBox(width) {
    const top = 1;
    const bottom = -1;
    const front = 0.5;
    const back = -0.5;

    const boundaryBox = new Float32Array([
      -width / 2, // left
      width / 2, // right
      top,
      bottom,
      front,
      back,
    ]);

    this._left_bound = -width / 2;
    this._right_bound = width / 2;
    this._top_bound = top;
    this._bot_bound = bottom;

    this._grid_length = Math.ceil(width / this._cell_size);
    this._grid_height = Math.ceil((top - bottom) / this._cell_size);
    this._num_grid_cells = this._grid_length * this._grid_height;
    this._grid = new Int32Array(
      Array(this._num_grid_cells * this._max_per_cell).fill(-1)
    );

    console.log(
      `Updating box to: ${width} (${this._left_bound} to ${this._right_bound})`
    );

    this._device.queue.writeBuffer(this._boundaryBoxBuffer, 0, boundaryBox);
    // this.resetParticles();
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
      },
    });
    this._computePipeline = this._device.createComputePipeline({
      label: "Particles Compute Pipeline " + this.getName(),
      layout: this._pipelineLayout,
      compute: {
        module: this._shaderModule,
        entryPoint: "computeMain",
      },
    });
    this._clearGridPipeline = this._device.createComputePipeline({
      label: "Grid Compute Pipeline " + this.getName(),
      layout: this._pipelineLayout,
      compute: {
        module: this._shaderModule,
        entryPoint: "clearGridStructure",
      },
    });
    this._gridPipeline = this._device.createComputePipeline({
      label: "Grid Compute Pipeline " + this.getName(),
      layout: this._pipelineLayout,
      compute: {
        module: this._shaderModule,
        entryPoint: "computeGridStructure",
      },
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

      ++this._step;
      this.updateTimeBuffer();
    }
  }

  mouseInteraction(x, y, attract) {
    const radius = 0.08;
    const force = 10 * attract;
    //brute force way, find every single particle that is within the radius from the x,y coords
    var dist = 0;
    for (var i = 0; i < this._particles.length; i++) {
      let ox = this._particles[VALUES_PER_PARTICLE * i] - x;
      let oy = this._particles[VALUES_PER_PARTICLE * i + 1] - y;
      dist = Math.sqrt(ox * ox + oy * oy);

      if (dist < radius) {
        //calcuate px,py in mouse click's coordinate as origin
        let angle = Math.atan(oy / ox);
        if (ox < 0) {
          //arctan only does left side of graph, need to correct if on the other side
          angle += Math.PI;
        }

        //add the effect to the particle's velocity
        this._particles[VALUES_PER_PARTICLE * i + 4] +=
          (force * Math.cos(angle)) / dist;
        this._particles[VALUES_PER_PARTICLE * i + 5] +=
          (force * Math.sin(angle)) / dist;
        //write the new values to the buffers
        this._device.queue.writeBuffer(
          this._particleBuffers[this._step % 2],
          (VALUES_PER_PARTICLE * i + 4) * Float32Array.BYTES_PER_ELEMENT,
          new Float32Array([
            (force * Math.cos(angle) * radius) / dist,
            (force * Math.sin(angle) * radius) / dist,
          ])
        );
      }
    }
  }
}
