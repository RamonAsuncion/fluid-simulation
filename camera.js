import { mat4 } from "https://wgpu-matrix.org/dist/3.x/wgpu-matrix.module.js";

// from https://github.com/matsuoka-601/Splash/blob/3df5d621e6f153d7269037cbd8606bd83796f3fd/camera.ts

export const renderUniformsValues = new ArrayBuffer(272);
export const renderUniformsViews = {
  texelSize: new Float32Array(renderUniformsValues, 0, 2),
  sphereSize: new Float32Array(renderUniformsValues, 8, 2),
  invProjectionMatrix: new Float32Array(renderUniformsValues, 16, 16),
  projectionMatrix: new Float32Array(renderUniformsValues, 80, 16),
  viewMatrix: new Float32Array(renderUniformsValues, 144, 16),
  invViewMatrix: new Float32Array(renderUniformsValues, 208, 16),
};

export class Camera {
  constructor(canvas) {
    this.canvas = canvas;
    this.isDragging = false;
    this.prevX = 0;
    this.prevY = 0;
    this.prevHoverX = 0;
    this.prevHoverY = 0;
    this.currentHoverX = 0;
    this.currentHoverY = 0;
    this.currentXtheta = -Math.PI / 2;
    this.currentYtheta = (-Math.PI / 12) * 0.8;
    this.maxYTheta = (-Math.PI / 12) * 0.8;
    this.minYTheta = (-0.99 * Math.PI) / 2;
    this.sensitivity = 0.005;
    this.currentDistance = 3;
    this.maxDistance = 5;
    this.minDistance = 1;
    this.target = [0, 0, 0];
    this.fov = Math.PI / 3;
    this.zoomRate = 0.2;

    // Expose view matrices for shader access
    this.viewMatrix = new Float32Array(16);
    this.invViewMatrix = new Float32Array(16);
    mat4.identity(this.viewMatrix);
    mat4.identity(this.invViewMatrix);

    // Event listeners for camera control
    if (this.canvas) {
      this.canvas.addEventListener("mousedown", (event) => {
        this.isDragging = true;
        this.prevX = event.clientX;
        this.prevY = event.clientY;
      });

      this.canvas.addEventListener("wheel", (event) => {
        event.preventDefault();
        var scrollDelta = event.deltaY;
        this.currentDistance += (scrollDelta > 0 ? 1 : -1) * this.zoomRate;
        if (this.currentDistance < this.minDistance)
          this.currentDistance = this.minDistance;
        if (this.currentDistance > this.maxDistance)
          this.currentDistance = this.maxDistance;
        this.recalculateView();
      });

      this.canvas.addEventListener("mousemove", (event) => {
        this.currentHoverX = event.clientX;
        this.currentHoverY = event.clientY;
        if (this.isDragging) {
          const deltaX = this.prevX - event.clientX;
          const deltaY = this.prevY - event.clientY;
          this.currentXtheta += this.sensitivity * deltaX;
          this.currentYtheta += this.sensitivity * deltaY;
          if (this.currentYtheta > this.maxYTheta)
            this.currentYtheta = this.maxYTheta;
          if (this.currentYtheta < this.minYTheta)
            this.currentYtheta = this.minYTheta;
          this.prevX = event.clientX;
          this.prevY = event.clientY;
          this.recalculateView();
        }
      });

      this.canvas.addEventListener("mouseup", () => {
        if (this.isDragging) this.isDragging = false;
      });

      this.recalculateView();
    }
  }

  reset(initDistance, target, fov, zoomRate) {
    this.isDragging = false;
    this.prevX = 0;
    this.prevY = 0;
    this.currentXtheta = -Math.PI / 2;
    this.currentYtheta = (-Math.PI / 12) * 0.8;
    this.maxYTheta = (-Math.PI / 12) * 0.8;
    this.minYTheta = (-0.99 * Math.PI) / 2;
    this.sensitivity = 0.005;
    this.currentDistance = initDistance;
    this.maxDistance = 1.3 * this.currentDistance;
    this.minDistance = 0.8 * this.currentDistance;
    this.target = target;
    this.fov = fov;
    this.zoomRate = zoomRate;

    const aspect = this.canvas.clientWidth / this.canvas.clientHeight;
    const projection = mat4.perspective(fov, aspect, 0.1, 300);
    renderUniformsViews.projectionMatrix.set(projection);
    renderUniformsViews.invProjectionMatrix.set(mat4.inverse(projection));

    this.recalculateView();
  }

  recalculateView() {
    var mat = mat4.identity();
    mat4.translate(mat, this.target, mat);
    mat4.rotateY(mat, this.currentXtheta, mat);
    mat4.rotateX(mat, this.currentYtheta, mat);
    mat4.translate(mat, [0, 0, this.currentDistance], mat);
    var position = mat4.multiply(mat, [0, 0, 0, 1]);

    const view = mat4.lookAt(
      [position[0], position[1], position[2]], // position
      this.target, // target
      [0, 1, 0] // up
    );

    renderUniformsViews.viewMatrix.set(view);
    renderUniformsViews.invViewMatrix.set(mat4.inverse(view));

    // Update local copies for direct access by ParticleSystemObject
    this.viewMatrix.set(view);
    this.invViewMatrix.set(mat4.inverse(view));
  }

  calcMouseVelocity() {
    if (this.isDragging) {
      return [0, 0];
    }

    let [mousePlaneX, mousePlaneY] = this.calcPlaneCoord(
      this.currentHoverX,
      this.currentHoverY
    );
    let [prevMousePlaneX, prevMousePlaneY] = this.calcPlaneCoord(
      this.prevHoverX,
      this.prevHoverY
    );

    let mouseVelocityX = mousePlaneX - prevMousePlaneX;
    let mouseVelocityY = mousePlaneY - prevMousePlaneY;

    let clampValue = 4.0;
    if (mouseVelocityX > clampValue) mouseVelocityX = clampValue;
    if (mouseVelocityX < -clampValue) mouseVelocityX = -clampValue;
    if (mouseVelocityY > clampValue) mouseVelocityY = clampValue;
    if (mouseVelocityY < -clampValue) mouseVelocityY = -clampValue;

    return [mouseVelocityX, mouseVelocityY, 0, 0];
  }

  calcPlaneCoord(x, y) {
    let normalizedX = x / this.canvas.width;
    let normalizedY = y / this.canvas.height;
    let ndcX = 2.0 * normalizedX - 1.0;
    let ndcY = (1.0 - normalizedY) * 2.0 - 1.0;

    let viewSpaceMouseRay = [
      ndcX *
        Math.tan(this.fov / 2.0) *
        (this.canvas.width / this.canvas.height),
      ndcY * Math.tan(this.fov / 2.0),
      -1.0,
    ];

    return [
      viewSpaceMouseRay[0] * this.currentDistance,
      viewSpaceMouseRay[1] * this.currentDistance,
    ];
  }

  setNewPrevMouseCoord() {
    this.prevHoverX = this.currentHoverX;
    this.prevHoverY = this.currentHoverY;
  }

  stepAngle() {
    this.currentXtheta += 0.012;
    this.recalculateView();
  }
}
