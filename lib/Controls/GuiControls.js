export default class GuiControls {
  constructor(particles) {
    this._loadLilGui().then(() => {
      this._init(particles);
    });
  }

  async _loadLilGui() {
    return new Promise((resolve, reject) => {
      if (window.lil && window.lil.GUI) {
        resolve();
        return;
      }

      const script = document.createElement("script");
      script.src = "https://cdn.jsdelivr.net/npm/lil-gui@0.20";
      script.onload = () => resolve();
      script.onerror = () =>
        reject(new Error("Failed to load lil-gui library"));
      document.head.appendChild(script);
    });
  }

  setInitialValues(values) {
    if (this.params) {
      Object.assign(this.params, values);
    } else {
      setTimeout(() => {
        this.setInitialValues(values);
      }, 50);
    }
  }

  _init(particles) {
    this.particles = particles;
    this.params = {
      particleCount: 10000,
      diffuseColor: { r: 0.5, g: 0.5, b: 1.0 },
      surfaceTension: 0.5,
      speed: 1,
      debugMode: false,
      boxWidth: 1.9,
      showBoundingBox: true,

      resetSimulation: () => {
        this.particles.resetParticles();
        console.log("Reset simulation requested");
      },
    };

    const GUI = window.lil.GUI;
    this.gui = new GUI({ title: "Controls" });
    const particlesFolder = this.gui.addFolder("Particles");

    particlesFolder
      .add(this.params, "particleCount", {
        "10,000": 10000,
        "25,000": 25000,
        "50,000": 50000,
        "100,000": 100000,
      })
      .name("Number of particles")
      .onChange((value) => {
        // TODO: Update particle count
        console.log(`Particle count changed to: ${value}`);
      });

    particlesFolder
      .add(this.params, "speed", 0, 1)
      .name("Simulation Speed")
      .onChange((value) => {
        this.particles.simulationSpeed = value;
        this.particles.setSimulationSpeed(value);
        console.log(`Speed changed to: ${value}`);
      });

    particlesFolder
      .add(this.params, "debugMode")
      .name("Particle Mode (Debug)")
      .onChange((value) => {
        // TODO: Toggle between debug mode
        console.log(`Debug mode ${value ? "enabled" : "disabled"}`);
      });

    const colorFolder = this.gui.addFolder("Diffuse Color");
    colorFolder
      .addColor(this.params, "diffuseColor")
      .name("Diffuse Color")
      .onChange((value) => {
        // TODO: Update diffuse color
        console.log(
          `Color changed to: R=${value.r}, G=${value.g}, B=${value.b}`
        );
      });

    const physicsFolder = this.gui.addFolder("Physics");
    physicsFolder
      .add(this.params, "surfaceTension", 0, 1)
      .name("Surface Tension")
      .onChange((value) => {
        // TODO: Update surface tension
        console.log(`Surface tension changed to: ${value}`);
      });

    physicsFolder
      .add(this.params, "boxWidth", 0.05, 1.9)
      .name("Box Width")
      .onChange((value) => {
        this.particles.updateBoundaryBox(value);
        console.log(`Box width changed to: ${value}`);
      });

    const actionsFolder = this.gui.addFolder("Actions");
    actionsFolder.add(this.params, "resetSimulation").name("Reset Simulation");
  }
}
