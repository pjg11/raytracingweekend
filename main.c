#include "rtweekend.h"

int main(void) {

  camera cam = {1.0, 100, 10, 10};
  double r = cos(M_PI/4);

  spherelist world = {0};
  material left = lambertian(v3(0, 0, 1)),
           right = lambertian(v3(1, 0, 0));

  spherelistadd(&world, sp(v3(-r, 0, -1), r, left));
  spherelistadd(&world, sp(v3(r, 0, -1), r, right));

  cam.aspectratio = 16.0 / 9.0;
  cam.imagewidth = 400;
  cam.samplesperpixel = 100;
  cam.maxdepth = 50;

  render(&cam, &world);
  return 0;
}
