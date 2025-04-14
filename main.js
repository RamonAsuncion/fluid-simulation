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

import Renderer from './lib/Viz/2DRenderer.js'
import FireWorkObject from './lib/DSViz/FireWorkObject.js'
import ParticleSystemObject from './lib/DSViz/ParticleSystemObject.js'
import StandardTextObject from './lib/DSViz/StandardTextObject.js'

async function init() {
  // Create a canvas tag
  const canvasTag = document.createElement('canvas');
  canvasTag.id = "renderCanvas";
  document.body.appendChild(canvasTag);
  // Create a 2d animated renderer
  const renderer = new Renderer(canvasTag);
  await renderer.init();
  const particles = new ParticleSystemObject(renderer._device, renderer._canvasFormat);
  await renderer.appendSceneObject(particles);
  let fps = '??';
  var fpsText = new StandardTextObject('fps: ' + fps + '\nClick and drag to interact!\nq: push/pull\nw: firework');
  
  // run animation at 60 fps
  var frameCnt = 0;
  var tgtFPS = 60;
  var secPerFrame = 1. / tgtFPS;
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
  var activeFireworks = [];
  var attract = -1;

  window.addEventListener("keydown", (e) => {
    switch (e.key) {
      case 'q': case 'Q':
        attract *= -1;
        break;
      case 'w': case 'W':
        var firework = createFirework(renderer);
        activeFireworks.push(firework);
        console.log("fireworks!");
    }
  });
  canvasTag.addEventListener('mousedown', (e) => {
    var mouseX = (e.clientX / window.innerWidth) * 2 - 1;
    var mouseY = (-e.clientY / window.innerHeight) * 2 + 1;
    particles.mouseInteraction(mouseX, mouseY, attract);
    isDragging = true;
  });
  canvasTag.addEventListener('mousemove', (e) => {
      var mouseX = (e.clientX / window.innerWidth) * 2 - 1;
      var mouseY = (-e.clientY / window.innerHeight) * 2 + 1;
      if (isDragging){
        particles.mouseInteraction(mouseX, mouseY, attract);
      }
    });
    canvasTag.addEventListener('mouseup', (e) => {
      isDragging = false;
    });

  lastCalled = Date.now();
  renderFrame();
  setInterval(() => { 
    fpsText.updateText('fps: ' + frameCnt + '\nClick and drag to interact!\nq: push/pull\nw: firework');
    frameCnt = 0;
  }, 1000); // call every 1000 ms
  return renderer;
}

init().then( ret => {
  console.log(ret);
}).catch( error => {
  const pTag = document.createElement('p');
  pTag.innerHTML = navigator.userAgent + "</br>" + error.message;
  document.body.appendChild(pTag);
  document.getElementById("renderCanvas").remove();
});

async function createFirework(renderer){
  const randX = Math.random() * 2 - 1;
  const randY = Math.random() * 2 - 1;
  const fireworks = new FireWorkObject(renderer._device, renderer._canvasFormat, randX, randY);
  await renderer.appendSceneObject(fireworks);
  return fireworks;
}