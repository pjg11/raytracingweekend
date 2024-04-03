#include "rtweekend.h"

int main(void) {

  camera cam = {1.0, 100, 10, 10};
  spherelist world = {0};
  material ground = lambertian(v3(0.8, 0.8, 0)),
           center = lambertian(v3(0.1, 0.2, 0.5)),
           left = dielectric(1.5),
           right = metal(v3(0.8, 0.6, 0.2), 0);

  spherelistadd(&world, sp(v3(0, -100.5, -1), 100, ground));
  spherelistadd(&world, sp(v3(0, 0.0, -1), 0.5, center));
  spherelistadd(&world, sp(v3(-1, 0.0, -1), 0.5, left));
  spherelistadd(&world, sp(v3(-1, 0.0, -1), -0.4, left));
  spherelistadd(&world, sp(v3(1, 0.0, -1), 0.5, right));

  cam.aspectratio = 16.0 / 9.0;
  cam.imagewidth = 400;
  cam.samplesperpixel = 100;
  cam.maxdepth = 50;

  render(&cam, &world);
  return 0;
}
