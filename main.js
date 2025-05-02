import Renderer from "./lib/Viz/3DRenderer.js";
import ParticleSystemObject from "./lib/DSViz/ParticleSystemObject.js";
import StandardTextObject from "./lib/DSViz/StandardTextObject.js";
import GuiControls from "./lib/Controls/GuiControls.js";
import Camera from "/lib/Viz/3DCamera.js";

/**
 *
 * TODO: Ctrl + click rotates the camera.
 * TODO: move around mouse effects the particles
 */

async function init() {
  var frameCnt = 0;
  var tgtFPS = 60;
  var secPerFrame = 1 / tgtFPS;
  var frameInterval = secPerFrame * 1000;
  var lastCalled;
  var playing = true;
  let lastMouseX = 0;
  let lastMouseY = 0;
  let mouseDown = 0;
  let rotateSpeed = 0.01;
  var isDragging = false;
  var attract = -1;

  const canvasTag = document.createElement("canvas");
  canvasTag.id = "renderCanvas";
  document.body.appendChild(canvasTag);
  const renderer = new Renderer(canvasTag);
  await renderer.init();

  var camera = new Camera();

  const particles = new ParticleSystemObject(
    renderer._device,
    renderer._canvasFormat,
    camera
  );

  await renderer.appendSceneObject(particles);

  let controls;
  try {
    controls = new GuiControls(particles);
    controls.setInitialValues({ boxWidth: 1.9 });
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
      "Drag + Ctrl - Rotate camera"
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
    const mouseX = ((e.clientX - rect.left) / rect.width) * 2 - 1;
    const mouseY = (1 - (e.clientY - rect.top) / rect.height) * 2 - 1;
    lastMouseX = mouseX;
    lastMouseY = mouseY;

    particles.setMousePosition(mouseX, mouseY);
    particles.setMouseDown(true);
    isDragging = true;
    console.log("currently dragging");
  });

  canvasTag.addEventListener("mousemove", (e) => {
    const rect = canvasTag.getBoundingClientRect();
    const mouseX = ((e.clientX - rect.left) / rect.width) * 2 - 1;
    const mouseY = (1 - (e.clientY - rect.top) / rect.height) * 2 - 1;
    lastMouseX = mouseX;
    lastMouseY = mouseY;

    if (isCtrlPressed) {
      camera.rotateY(deltaX * rotateSpeed);
      camera.rotateX(deltaY * rotateSpeed);
      particles.updateCameraPose(camera);
    }

    particles.setMousePosition(mouseX, mouseY);

    if (isDragging) {
      particles.setMousePosition(mouseX, mouseY);
    }
  });

  canvasTag.addEventListener("mouseup", (e) => {
    particles.setMouseDown(false);
    isDragging = false;
  });

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
