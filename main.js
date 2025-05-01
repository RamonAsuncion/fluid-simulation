import Renderer from "./lib/Viz/3DRenderer.js";
import ParticleSystemObject from "./lib/DSViz/ParticleSystemObject.js";
import StandardTextObject from "./lib/DSViz/StandardTextObject.js";
import GuiControls from "./lib/Controls/GuiControls.js";

async function init() {
  const canvasTag = document.createElement("canvas");
  canvasTag.id = "renderCanvas";
  document.body.appendChild(canvasTag);
  const renderer = new Renderer(canvasTag);
  await renderer.init();

  const particles = new ParticleSystemObject(
    renderer._device,
    renderer._canvasFormat
  );

  await renderer.appendSceneObject(particles);

  let controls;
  try {
    controls = new GuiControls(particles);
    controls.setInitialValues({ boxWidth: 1.9 });
  } catch (e) {
    console.log("GUI controls not available: ", e);
  }

  let fps = "??";
  var fpsText = new StandardTextObject(
    "fps: " +
      fps +
      "\nClick and drag to interact!\nq: push/pull\np: pause\nr: reset"
  );

  // run animation at 60 fps
  var frameCnt = 0;
  var tgtFPS = 60;
  var secPerFrame = 1 / tgtFPS;
  var frameInterval = secPerFrame * 1000;
  var lastCalled;
  var playing = true;
  let renderFrame = () => {
    //console.log("rendering");
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

  var isDragging = false;
  var attract = -1;

  window.addEventListener("keydown", (e) => {
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
        console.log("Simulation reset");
        break;
    }
  });

  canvasTag.addEventListener("mousedown", (e) => {
    const rect = canvasTag.getBoundingClientRect();
    const mouseX = ((e.clientX - rect.left) / rect.width) * 2 - 1;
    const mouseY = (1 - (e.clientY - rect.top) / rect.height) * 2 - 1;
    if (particles.setMouseDown) {
      particles.setMouseDown(true);
      particles.setMousePosition(mouseX, mouseY);
    }

    particles.mouseInteraction(mouseX, mouseY, attract);
    isDragging = true;
  });

  canvasTag.addEventListener("mousemove", (e) => {
    const rect = canvasTag.getBoundingClientRect();
    const mouseX = ((e.clientX - rect.left) / rect.width) * 2 - 1;
    const mouseY = (1 - (e.clientY - rect.top) / rect.height) * 2 - 1;
    if (particles.setMousePosition) {
      particles.setMousePosition(mouseX, mouseY);
    }
    if (isDragging) {
      particles.mouseInteraction(mouseX, mouseY, attract);
    }
  });

  canvasTag.addEventListener("mouseup", (e) => {
    if (particles.setMouseDown) {
      particles.setMouseDown(false);
    }

    isDragging = false;
  });

  // make a wheel

  lastCalled = Date.now();
  renderFrame();
  setInterval(() => {
    fpsText.updateText(
      "fps: " +
        frameCnt +
        "\nClick and drag to interact!\nq: push/pull\np: pause\nr: reset"
    );
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
