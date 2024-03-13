#include "camera.h"
#include "rtweekend.h"

int main(void) {

  spherelist world = {0};

  spherelistadd(&world, sp(v3(0, 0, -1), 0.5));
  spherelistadd(&world, sp(v3(0, -100.5, -1), 100));

  camera cam = {1.0, 100, 10};

  cam.aspectratio = 16.0 / 9.0;
  cam.imagewidth = 400;
  cam.samplesperpixel = 100;

  render(&cam, &world);
  return 0;
}
