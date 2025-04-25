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

// Check your browser supports: https://github.com/gpuweb/gpuweb/wiki/Implementation-Status#implementation-status
// Need to enable experimental flags chrome://flags/
// Chrome & Edge 113+ : Enable Vulkan, Default ANGLE Vulkan, Vulkan from ANGLE, Unsafe WebGPU Support, and WebGPU Developer Features (if exsits)
// Firefox Nightly: sudo snap install firefox --channel=latext/edge or download from https://www.mozilla.org/en-US/firefox/channel/desktop/

import Renderer from "./lib/Viz/2DRenderer.js";
import ParticleSystemObject from "./lib/DSViz/ParticleSystemObject.js";
import StandardTextObject from "./lib/DSViz/StandardTextObject.js";
import GuiControls from "./lib/Controls/GuiControls.js";
import Camera from "/lib/Viz/3DCamera.js";

async function init() {
  // Create a canvas tag
  const canvasTag = document.createElement("canvas");
  canvasTag.id = "renderCanvas";
  document.body.appendChild(canvasTag);
  // Create a 2d animated renderer
  const renderer = new Renderer(canvasTag);
  await renderer.init();
  // Create a 3D Camera
  var camera = new Camera();
  // Camera mode (false orthogonal, true projective)
  var lastMouseX = 0;
  var lastMouseY = 0;
  var mouseDown = false;
  var rotationSpeed = 0.01;
  const particles = new ParticleSystemObject(
    renderer._device,
    renderer._canvasFormat,
    camera
  );

  await renderer.appendSceneObject(particles);
  let fps = "??";
  var fpsText = new StandardTextObject("fps: " + fps);

  // Add instructions text with commands
  var instructionsText = new StandardTextObject(
    "Controls:\n" +
      "R - Reset simulation\n" +
      "P - Pause menu\n" +
      "Drag - Move camera\n" +
      "Scroll - Zoom in/out"
  );
  // Position the instructions text below FPS counter
  instructionsText._textCanvas.style.top = "60px";

  // Initialize GUI controls
  const controls = new GuiControls(particles);
  controls.setInitialValues({ boxWidth: 0.1 });

  // run animation at 60 fps
  var frameCnt = 0;
  var tgtFPS = 60;
  var secPerFrame = 1 / tgtFPS;
  var frameInterval = secPerFrame * 1000;
  var lastCalled;
  let renderFrame = () => {
    let elapsed = Date.now() - lastCalled;
    if (elapsed > frameInterval) {
      ++frameCnt;
      lastCalled = Date.now() - (elapsed % frameInterval);
      renderer.render();
    }
    requestAnimationFrame(renderFrame);
  };

  var isDragging = false;
  var mouseX = 0;
  var mouseY = 0;
  let isCtrlPressed = false;

  // Set up keyboard interaction
  window.addEventListener("keydown", (e) => {
    if (e.key === "Control") {
      isCtrlPressed = true;
    }
    switch (e.key) {
      case "r":
      case "R":
        // TODO: Reset simulation
        console.log("Reset requested via keyboard");
        break;
      case "p":
      case "P":
        // TODO: Toggle pause menu
        console.log("Pause menu toggled");
        break;
      // TODO: Add more keyboard interactions
    }
  });

  window.addEventListener("keyup", (e) => {
    if (e.key === "Control") {
      isCtrlPressed = false;
    }
  });

  canvasTag.addEventListener("mousedown", (e) => {
    mouseX = (e.clientX / window.innerWidth) * 2 - 1;
    mouseY = (-e.clientY / window.innerHeight) * 2 + 1;
    lastMouseX = e.clientX;
    lastMouseY = e.clientY;
    mouseDown = true;
    isDragging = true;
  });

  canvasTag.addEventListener("mousedown", (e) => {
    mouseX = (e.clientX / window.innerWidth) * 2 - 1;
    mouseY = (-e.clientY / window.innerHeight) * 2 + 1;
    console.log(`X: ${mouseX} Y: ${mouseY}`);
    lastMouseX = e.clientX;
    lastMouseY = e.clientY;
    mouseDown = true;
    isDragging = true;
  });

  // Update your mousemove event handler
  canvasTag.addEventListener("mousemove", (e) => {
    if (!mouseDown) return;

    const deltaX = e.clientX - lastMouseX;
    const deltaY = e.clientY - lastMouseY;

    if (isCtrlPressed) {
      camera.moveX(-deltaX * 0.001);
      camera.moveY(deltaY * 0.001);
    } else {
      camera.rotateY(deltaX * rotationSpeed);
      camera.rotateX(-deltaY * rotationSpeed);
    }

    particles.updateCameraPose(camera);

    lastMouseX = e.clientX;
    lastMouseY = e.clientY;
  });

  canvasTag.addEventListener("mouseup", (e) => {
    mouseDown = false;
    isDragging = false;
  });

  // https://developer.mozilla.org/en-US/docs/Web/API/WheelEvent/deltaY
  canvasTag.addEventListener(
    "wheel",
    (e) => {
      e.preventDefault();
      const zoomAmount = e.deltaY * 0.005;
      camera.moveZ(zoomAmount);
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
