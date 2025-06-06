import Renderer from "./lib/Viz/3DRenderer.js";
import ParticleSystemObject from "./lib/DSViz/ParticleSystemObject.js";
import StandardTextObject from "./lib/DSViz/StandardTextObject.js";
import GuiControls from "./lib/Controls/GuiControls.js";
import { Camera } from "./camera.js";

let DEBUG = false;

async function init() {
  const urlParams = new URLSearchParams(window.location.search);
  const mode = urlParams.get("mode") || "3D";

  if (mode !== "3D") {
    window.location.href = "/?mode=2D";
    return;
  }

  var frameCnt = 0;
  var tgtFPS = 60;
  var secPerFrame = 1 / tgtFPS;
  var frameInterval = secPerFrame * 1000;
  var lastCalled;
  var playing = true;
  var isDragging = false;
  var attract = -1;

  let mouseDown = false;
  let lastMouseX = 0;
  let lastMouseY = 0;
  let mouseX = 0;
  let mouseY = 0;

  let lastMouseMoveTime = 0;
  let mouseVelocityX = 0;
  let mouseVelocityY = 0;
  let mouseMoving = 0;
  let mouseMovementTimeout = null;

  const canvasTag = document.createElement("canvas");
  canvasTag.id = "renderCanvas";
  document.body.appendChild(canvasTag);
  const renderer = new Renderer(canvasTag);
  await renderer.init();

  const camera = new Camera(canvasTag);
  /**
   * camera distance from target,
   * initial target point,
   * field of view (60 degrees)
   * zoom rate
   */
  camera.reset(3, [0, 0, 0], Math.PI / 3, 0.2);

  const particles = new ParticleSystemObject(
    renderer._device,
    renderer._canvasFormat,
    camera
  );

  const boxCenter = particles.getCenterBoundaryBox();
  camera.target = boxCenter;
  camera.reset(3, boxCenter, Math.PI / 3, 0.2);

  await renderer.appendSceneObject(particles);

  let controls;
  try {
    controls = new GuiControls(particles);
    controls.setInitialValues({ boxWidth: 1.9, simulationMode: "3D" });
  } catch (e) {
    console.error("GUI controls not available: ", e);
  }

  let fps = "??";
  var fpsText = new StandardTextObject("fps: " + fps);

  var instructionsText = new StandardTextObject(
    "Controls:\n" +
      "R - Reset simulation\n" +
      "P - Pause menu\n" +
      "Drag - Move particles\n" +
      "Drag + Click - Rotate camera"
  );
  instructionsText._textCanvas.style.top = "60px";

  let renderFrame = () => {
    let elapsed = Date.now() - lastCalled;
    if (elapsed > frameInterval) {
      ++frameCnt;
      lastCalled = Date.now() - (elapsed % frameInterval);
      if (playing) {
        renderer.render();
      }
    }
    requestAnimationFrame(renderFrame);
  };

  /**
   * todo: shift sets the attraction mode (fix later/add interaction)
   */

  window.addEventListener("keydown", (e) => {
    if (e.key === "Shift") {
      particles.setAttractMode(true);
      particles.setMouseDown(true);
    }

    switch (e.key) {
      case "q":
      case "Q":
        attract *= -1;
        break;
      case "p":
      case "P":
        playing = !playing;
        break;
      case "r":
      case "R":
        particles.resetParticles();
        break;
    }
  });

  window.addEventListener("keyup", (e) => {
    if (e.key === "Shift") {
      particles.setAttractMode(false);
      if (!isDragging) {
        particles.setMouseDown(false);
      }
    }
  });

  canvasTag.addEventListener("mousedown", (e) => {
    const rect = canvasTag.getBoundingClientRect();
    mouseX = ((e.clientX - rect.left) / rect.width) * 2 - 1;
    mouseY = (1 - (e.clientY - rect.top) / rect.height) * 2 - 1;

    lastMouseX = e.clientX;
    lastMouseY = e.clientY;

    mouseDown = true;
    isDragging = true;

    particles.setMousePosition(mouseX, mouseY);
    particles.setMouseDown(true);
    if (DEBUG)
      console.log(`Mouse down at (${mouseX.toFixed(2)}, ${mouseY.toFixed(2)})`);
  });

  particles.setMouseRadius(0.5);

  canvasTag.addEventListener("mousemove", (e) => {
    const rect = canvasTag.getBoundingClientRect();
    const currentTime = performance.now();
    const timeDelta = currentTime - lastMouseMoveTime;
    mouseX = ((e.clientX - rect.left) / rect.width) * 2 - 1;
    mouseY = (1 - (e.clientY - rect.top) / rect.height) * 2 - 1;
    if (timeDelta > 5) {
      mouseVelocityX = (mouseX - lastMouseX) / (timeDelta / 1000);
      mouseVelocityY = (mouseY - lastMouseY) / (timeDelta / 1000);
      lastMouseMoveTime = currentTime;
      mouseMoving = 1;
      clearTimeout(mouseMovementTimeout);
      mouseMovementTimeout = setTimeout(() => {
        mouseMoving = 0;

        mouseVelocityX = 0;
        mouseVelocityY = 0;
        particles.updateMouseVelocity(
          mouseVelocityX,
          mouseVelocityY,
          mouseMoving
        );
      }, 100);
    }
    particles.setMousePosition(mouseX, mouseY);
    particles.updateMouseVelocity(mouseVelocityX, mouseVelocityY, mouseMoving);
    particles.updateCameraPose(camera);
    lastMouseX = e.clientX;
    lastMouseY = e.clientY;
    if (DEBUG)
      console.log(`Mouse move at (${mouseX.toFixed(2)}, ${mouseY.toFixed(2)})`);
  });

  canvasTag.addEventListener("mouseup", (e) => {
    console.log("Mouse down:", particles._mousePosition[2]);
    console.log(
      "Mouse position:",
      particles._mousePosition[0],
      particles._mousePosition[1]
    );
    console.log("Mouse radius:", particles._mousePosition[3]);
    console.log("Attract mode:", particles._mousePosition[4]);
    mouseDown = false;
    isDragging = false;
    particles.setMouseDown(false);
  });

  // https://developer.mozilla.org/en-US/docs/Web/API/WheelEvent/deltaY
  canvasTag.addEventListener(
    "wheel",
    (e) => {
      e.preventDefault();
      particles.updateCameraPose(camera);
    },
    { passive: false }
  );

  lastCalled = Date.now();
  renderFrame();
  setInterval(() => {
    fpsText.updateText("fps: " + frameCnt);
    frameCnt = 0;
  }, 1000); // call every 1000 ms
  return renderer;
}

init()
  .then((ret) => {
    console.log(ret);
  })
  .catch((error) => {
    const pTag = document.createElement("p");
    pTag.innerHTML = navigator.userAgent + "</br>" + error.message;
    document.body.appendChild(pTag);
    document.getElementById("renderCanvas").remove();
  });
