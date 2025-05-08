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

import SceneObject from "./SceneObject.js";

let VALUES_PER_PARTICLE = 16;
let DEBUG = false;

export default class ParticleSystemObject extends SceneObject {
  constructor(device, canvasFormat, camera = null, numParticles = 4096) {
    super(device, canvasFormat);
    this._num_particles = numParticles;
    this._cell_size = 0.04;
    this._max_per_cell = 256;
    this._step = 0;
    this._camera = camera;
    this._simulationSpeed = 0.3;
    this._use_acceleration = true;
    this._camera = camera;
    this._left_bound = -1;
    this._right_bound = 1;
    this._top_bound = 1;
    this._bot_bound = -1;
    this._front_bound = 1.5;
    this._back_bound = 2.5;

    this._r_value = 1;
    this._g_value = 0;
    this._b_value = 0;

    this._grid_length = Math.ceil(
      (this._right_bound - this._left_bound) / this._cell_size
    );
    this._grid_height = Math.ceil(
      (this._top_bound - this._bot_bound) / this._cell_size
    );
    this._grid_width = Math.ceil(
      (this._back_bound - this._front_bound) / this._cell_size
    );
    this._num_grid_cells =
      this._grid_length * this._grid_height * this._grid_width;
    this._initialized = false;
    this._particle_mouse_radius = 0.25;
    this._pressureMultiplier = 2.0;
    this._gravityMultiplier = 0.5;
    this._isoValue = 0.0001;
    this._showMarchingCubes = false;
    this._mousePosition = new Float32Array([
      0,
      0,
      0,
      this._particle_mouse_radius,
      0,
      0,
      0,
      0,
      0, // padding
      0, // padding
    ]); // [x, y, isDown, radius, attractMode, lastMoveTime, velocityX, velocityY, isMoving]
  }

  async createGeometry() {
    this._cameraBuffer = this._device.createBuffer({
      label: "Camera " + this.getName(),
      size: 144, // 16 floats for view matrix, 16 floats for projection matrix, 4 floats for camera position
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });

    if (this._camera) {
      const cameraData = new Float32Array(36);
      cameraData[0] = 1;
      this._device.queue.writeBuffer(this._cameraBuffer, 0, cameraData);
    }

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
      this._front_bound, // front (for 3D later)
      this._back_bound, // back (for 3D later)
    ]);

    this._device.queue.writeBuffer(
      this._boundaryBoxBuffer,
      0,
      initialBoundaryBox
    );

    this._mouseBuffer = this._device.createBuffer({
      label: "Mouse interaction buffer",
      size: this._mousePosition.byteLength,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });

    this._device.queue.writeBuffer(this._mouseBuffer, 0, this._mousePosition);

    this._pressureMultiplier = 2.0;
    this._gravityMultiplier = 0.5;
    this.updateShaderConstants();

    await this.createParticleGeometry();

    this.start();
  }

  start() {
    if (!this._initialized) return;
    this.setMousePosition(0, 0);
    this.setMouseDown(true);
    setTimeout(() => {
      this.setMouseDown(false);
    }, 50);
    if (this._camera) {
      this.updateCameraPose(this._camera);
    }
    if (DEBUG) console.log("Simulation started automatically");
  }

  setMousePosition(x, y) {
    this._mousePosition[0] = x;
    this._mousePosition[1] = y;
    this._device.queue.writeBuffer(this._mouseBuffer, 0, this._mousePosition);
  }

  setAttractMode(isAttract) {
    this._mousePosition[4] = isAttract ? 1.0 : 0.0;
    this._device.queue.writeBuffer(this._mouseBuffer, 0, this._mousePosition);
  }

  setMouseDown(isDown) {
    // console.log(`Is down ${isDown}?`);
    this._mousePosition[2] = isDown ? 1.0 : 0.0;
    this._device.queue.writeBuffer(this._mouseBuffer, 0, this._mousePosition);
  }

  setMouseRadius(radius = 0.1) {
    this._particle_mouse_radius = radius;
    this._mousePosition[3] = radius;
    this._device.queue.writeBuffer(this._mouseBuffer, 0, this._mousePosition);
  }

  setPressureMultiplier(value) {
    this._pressureMultiplier = value;
    this.updateShaderConstants();
  }

  setGravityMultiplier(value) {
    this._gravityMultiplier = value;
    this.updateShaderConstants();
  }

  setIsoValue(value) {
    this._isoValue = value;
    this.updateShaderConstants();
  }

  updateMouseVelocity(velocityX, velocityY, isMoving) {
    this._mousePosition[5] = performance.now(); // lastMove timestamp
    this._mousePosition[6] = velocityX; // velX
    this._mousePosition[7] = velocityY; // velY
    this._mousePosition[8] = isMoving ? 1.0 : 0.0; // isMoving flag
    this._device.queue.writeBuffer(this._mouseBuffer, 0, this._mousePosition);
  }

  updateShaderConstants() {
    if (!this._constantsBuffer) {
      this._constantsBuffer = this._device.createBuffer({
        label: "Simulation constants buffer",
        size: 16,
        usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
      });
    }

    const constants = new Float32Array([
      this._pressureMultiplier,
      this._gravityMultiplier,
      this._max_per_cell,
      this._isoValue,
    ]);

    this._device.queue.writeBuffer(this._constantsBuffer, 0, constants);
  }

  async createParticleGeometry() {
    this._particles = new Float32Array(
      this._num_particles * VALUES_PER_PARTICLE
    );

    this._particleBuffers = [
      this._device.createBuffer({
        label: "particle buffer 1" + this.getName(),
        size: this._particles.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      }),
      this._device.createBuffer({
        label: "particle buffer 2" + this.getName(),
        size: this._particles.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      }),
    ];
    if (DEBUG)
      console.log(
        "new grid length: " + this._num_grid_cells,
        " ",
        this._max_per_cell
      );
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
      this._r_value,
      this._g_value,
      this._b_value,
    ]);
    this._timeBuffer = this._device.createBuffer({
      label: "time buffer " + this.getName(),
      size: this._time.byteLength,
      usage:
        GPUBufferUsage.STORAGE |
        GPUBufferUsage.COPY_DST |
        GPUBufferUsage.UNIFORM,
    });
    this._device.queue.writeBuffer(this._timeBuffer, 0, this._time);

    this._density = new Float32Array(this._num_particles);
    this._densityBuffer = this._device.createBuffer({
      label: "density buffer " + this.getName(),
      size: this._density.byteLength,
      usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    });

    this._triangles = new Float32Array(
      Array(this._num_grid_cells * 16 * 3).fill(-1)
    );
    this._triangleBuffers = [
      this._device.createBuffer({
        label: "triangle buffer 1" + this.getName(),
        size: this._triangles.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      }),
      this._device.createBuffer({
        label: "triangle buffer 2" + this.getName(),
        size: this._triangles.byteLength,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
      }),
    ];

    this._device.queue.writeBuffer(
      this._triangleBuffers[0],
      0,
      this._triangles.buffer
    );
    this._device.queue.writeBuffer(
      this._triangleBuffers[1],
      0,
      this._triangles.buffer
    );

    // what about the bind ground? which buffer do they point to in the bind group?
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
      var z = Math.floor(
        (this._particles[i * VALUES_PER_PARTICLE + 2] - this._front_bound) /
          this._cell_size
      );

      //bound checking
      if (
        x < 0 ||
        x > this._grid_length ||
        y < 0 ||
        y > this._grid_height ||
        z < 0 ||
        z > this._grid_width
      ) {
        if (DEBUG)
          console.log(`Particle ${i} outside grid bounds: (${x}, ${y}, ${z})`);
        continue;
      }
      var cell_start =
        (z * this._grid_height * this._grid_length +
          y * this._grid_length +
          x) *
        this._max_per_cell;
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
        if (DEBUG)
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
    if (DEBUG) console.log("Grid filled with particles");
    if (DEBUG)
      console.log(
        "grid dimensions:",
        this._grid_length,
        this._grid_height,
        this._grid_width
      );
  }

  resetParticles() {
    // reset the grid
    this._grid = new Int32Array(
      Array(this._num_grid_cells * this._max_per_cell).fill(-1)
    );

    for (let i = 0; i < this._num_particles; ++i) {
      // x pos
      this._particles[VALUES_PER_PARTICLE * i + 0] =
        this._left_bound +
        Math.random() * (this._right_bound - this._left_bound);

      // y pos
      this._particles[VALUES_PER_PARTICLE * i + 1] =
        this._top_bound - Math.random() * 0.05;

      // z pos
      this._particles[VALUES_PER_PARTICLE * i + 2] =
        this._front_bound +
        Math.random() * (this._back_bound - this._front_bound);

      this._particles[VALUES_PER_PARTICLE * i + 3] = 0; // dummy variable

      // initial position
      this._particles[VALUES_PER_PARTICLE * i + 4] =
        this._particles[VALUES_PER_PARTICLE * i + 0];
      this._particles[VALUES_PER_PARTICLE * i + 5] =
        this._particles[VALUES_PER_PARTICLE * i + 1];
      this._particles[VALUES_PER_PARTICLE * i + 6] =
        this._particles[VALUES_PER_PARTICLE * i + 2];

      this._particles[VALUES_PER_PARTICLE * i + 7] = 0; // dummy variable

      // no velocity
      this._particles[VALUES_PER_PARTICLE * i + 8] = 0;
      this._particles[VALUES_PER_PARTICLE * i + 9] = 0;
      this._particles[VALUES_PER_PARTICLE * i + 10] = 0;

      this._particles[VALUES_PER_PARTICLE * i + 11] = 0; // dummy dummy dummy dummy

      // initial velocity
      this._particles[VALUES_PER_PARTICLE * i + 12] =
        this._particles[VALUES_PER_PARTICLE * i + 8];
      this._particles[VALUES_PER_PARTICLE * i + 13] =
        this._particles[VALUES_PER_PARTICLE * i + 9];
      this._particles[VALUES_PER_PARTICLE * i + 14] =
        this._particles[VALUES_PER_PARTICLE * i + 10];

      this._particles[VALUES_PER_PARTICLE * i + 15] = 0; //dummy
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

    // reset the triangles
    if (this._triangleBuffers) {
      this._device.queue.writeBuffer(
        this._triangleBuffers[0],
        0,
        this._triangles.buffer
      );
      this._device.queue.writeBuffer(
        this._triangleBuffers[1],
        0,
        this._triangles.buffer
      );
    }

    if (DEBUG) console.log("Particles reset at top of boundary box");
  }

  async changeParticleCount(new_value) {
    this._num_particles = new_value;
    this._initialized = false;
    await this.createGeometry();
    await this.createShaders();
    await this.createRenderPipeline();
    await this.createComputePipeline();
    if (DEBUG) console.log("particle length:", this._particles.length);
  }

  updateTimeBuffer() {
    if (DEBUG)
      console.log(
        "array: ",
        0,
        this._simulationSpeed,
        this._use_acceleration ? 1 : 0,
        this._r_value,
        this._g_value,
        this._b_value
      );
    const time = new Float32Array([
      0,
      this._simulationSpeed,
      this._use_acceleration ? 1 : 0,
      this._r_value,
      this._g_value,
      this._b_value,
    ]);
    this._device.queue.writeBuffer(this._timeBuffer, 0, time);
  }

  setSimulationSpeed(speed) {
    this._simulationSpeed = speed;
    this.updateTimeBuffer();
    if (DEBUG) console.log(`Simulation speed set to: ${speed}`);
  }

  setRGB(r, g, b) {
    this._r_value = r;
    this._g_value = g;
    this._b_value = b;
  }

  updateGeometry() {}

  async createShaders() {
    const baseUrl = window.location.origin;
    let loadUrl = "../../shaders/particles.wgsl";
    if (baseUrl.includes("ramonasuncion.github.io")) {
      console.log("Github page", baseUrl);
      loadUrl = "/fluid-simulation/shaders/particles.wgsl";
    }
    let shaderCode = await this.loadShader(loadUrl);
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
        visibility: GPUShaderStage.COMPUTE | GPUShaderStage.FRAGMENT,
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
      {
        binding: 10,
        visibility: GPUShaderStage.VERTEX | GPUShaderStage.COMPUTE,
        buffer: { type: "read-only-storage" }, // triangle input buffer
      },
      {
        binding: 11,
        visibility: GPUShaderStage.COMPUTE,
        buffer: { type: "storage" }, // triangle output buffer
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

    if (this._constantsBuffer) {
      bindGroupEntries.push({
        binding: 9,
        visibility: GPUShaderStage.COMPUTE,
        buffer: { type: "uniform" }, // simulation constants buffer
      });
    }

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
        typology: "line-strip",
      },
    });
    this._ballPipeline = this._device.createRenderPipeline({
      label: "Particles Render Pipeline " + this.getName(),
      layout: this._pipelineLayout,
      vertex: {
        module: this._shaderModule,
        entryPoint: "ballBasedVertex",
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
        typology: "triangle-list",
      },
      depthStencil: {
        format: "depth24plus",
        depthWriteEnabled: true, // enable z-buffer - depth test
        depthCompare: "less", // closer pixels overwrite farther ones
      },
    });

    this._surfacePipeline = this._device.createRenderPipeline({
      label: "Surface Render Pipeline " + this.getName(),
      layout: this._pipelineLayout,
      vertex: {
        module: this._shaderModule,
        entryPoint: "surfaceBasedVertex",
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
        typology: "triangle-list",
      },
      depthStencil: {
        format: "depth24plus",
        depthWriteEnabled: true, // enable z-buffer - depth test
        depthCompare: "less", // closer pixels overwrite farther ones
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
      {
        binding: 10,
        resource: { buffer: this._triangleBuffers[0] },
      },
      {
        binding: 11,
        resource: { buffer: this._triangleBuffers[1] },
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
      {
        binding: 10,
        resource: { buffer: this._triangleBuffers[1] },
      },
      {
        binding: 11,
        resource: { buffer: this._triangleBuffers[0] },
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
    if (this._constantsBuffer) {
      bindGroupEntries1.push({
        binding: 9,
        resource: { buffer: this._constantsBuffer },
      });
      bindGroupEntries2.push({
        binding: 9,
        resource: { buffer: this._constantsBuffer },
      });
    }

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
    if (DEBUG) console.log("Camera pose updated");
    if (this._cameraBuffer && camera) {
      const cameraData = new Float32Array(36);
      // view matrix
      if (camera.viewMatrix) {
        cameraData.set(camera.viewMatrix, 0);
      }
      // inverse view matrix
      if (camera.invViewMatrix) {
        cameraData.set(camera.invViewMatrix, 16);
      }
      // camera position
      const position = [0, 0, camera.currentDistance, 1];
      cameraData.set(position, 32);

      this._device.queue.writeBuffer(this._cameraBuffer, 0, cameraData);
    }
  }

  getCenterBoundaryBox() {
    const centerX = (this._left_bound + this._right_bound) / 2;
    const centerY = this._top_bound; // use top value "the bottom"
    const centerZ = (this._front_bound + this._back_bound) / 2;
    return [centerX, centerY, centerZ];
  }

  updateBoundaryBox(width) {
    const top = 1;
    const bottom = -1;
    const front = 1.0;
    const back = 2.0;

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

    const new_grid_length = Math.ceil(width / this._cell_size);
    const new_grid_height = Math.ceil((top - bottom) / this._cell_size);
    const new_grid_width = Math.ceil((back - front) / this._cell_size);
    const new_num_grid_cells =
      new_grid_height * new_grid_length * new_grid_width;

    const MAX_CELLS = 18000000;
    if (new_num_grid_cells * this._max_per_cell > MAX_CELLS) {
      console.warn(
        `Grid too large ${new_num_grid_cells} * ${this._max_per_cell}`
      );
      return;
    }

    this._grid_length = new_grid_length;
    this._grid_height = new_grid_height;
    this._num_grid_cells = new_num_grid_cells;

    this._grid = new Int32Array(
      Array(this._num_grid_cells * this._max_per_cell).fill(-1)
    );

    if (DEBUG)
      console.log(
        `Updating box to: ${width} (${this._left_bound} to ${this._right_bound})`
      );

    this._device.queue.writeBuffer(this._boundaryBoxBuffer, 0, boundaryBox);
  }

  render(pass) {
    pass.setBindGroup(0, this._bindGroups[this._step % 2]);
    if (this._initialized) {
      if (this._showMarchingCubes) {
        pass.setPipeline(this._surfacePipeline);
        pass.draw(15, this._num_grid_cells);
      }
      else {
        pass.setPipeline(this._ballPipeline);
        var k = 6;
        pass.draw(6 * 2 * k * k, this._num_particles);
      }
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
    this._trianglePipeline = this._device.createComputePipeline({
      label: "Triangle Compute Pipeline " + this.getName(),
      layout: this._pipelineLayout,
      compute: {
        module: this._shaderModule,
        entryPoint: "computeCubeMarch",
      },
    });
  }
  setShowMarchingCubes(enabled) {
    this._showMarchingCubes = enabled;
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

      if (this._showMarchingCubes) {
        // reset the triangle
        this._device.queue.writeBuffer(
          this._triangleBuffers[(this._step + 1) % 2],
          0,
          this._triangles.buffer
        );

        pass.setPipeline(this._trianglePipeline); //set buffer for storing triangles from marching cubes algorithm
        pass.dispatchWorkgroups(Math.ceil(this._num_grid_cells / 256));
      }

      ++this._step;
      this.updateTimeBuffer();
    }
  }
}
