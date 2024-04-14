#include "rtweekend.h"

int main(void) {
  spherelist world = {0};
  camera cam = {1, 100, 10, 10, 90, {0, 0, -1}, {0, 0, 0}, {0, 1, 0}, 0, 10};

  material ground = lambertian(v3(0.5, 0.5, 0.5));
  spherelistadd(&world, sp(v3(0, -1000, 0), 1000, ground));

  for (int a = -11; a < 11; a++) {
    for (int b = -11; b < 11; b++) {
      double choosemat = randomdouble();
      vec3 center = v3(a + 0.9 * randomdouble(), 0.2, b + 0.9 * randomdouble());

      if (v3length(v3sub(center, v3(4, 0.2, 0))) > 0.9) {
        material spherematerial;
        if (choosemat < 0.8) {
          // diffuse
          spherematerial = lambertian(v3mul(v3random(), v3random()));
          spherelistadd(&world, sp(center, 0.2, spherematerial));
        } else if (choosemat < 0.95) {
          // metal
          spherematerial =
              metal(v3randominterval(0.5, 1), 0.5 * randomdouble());
          spherelistadd(&world, sp(center, 0.2, spherematerial));
        } else {
          // glass
          spherematerial = dielectric(1.5);
          spherelistadd(&world, sp(center, 0.2, spherematerial));
        }
      }
    }
  }

  material material1 = dielectric(1.5);
  spherelistadd(&world, sp(v3(0, 1, 0), 1.0, material1));

  material material2 = lambertian(v3(0.4, 0.2, 0.1));
  spherelistadd(&world, sp(v3(-4, 1, 0), 1.0, material2));

  material material3 = metal(v3(0.7, 0.6, 0.5), 0.0);
  spherelistadd(&world, sp(v3(4, 1, 0), 1.0, material3));

  cam.aspectratio = 16.0 / 9.0;
  cam.imagewidth = 1200;
  cam.samplesperpixel = 500;
  cam.maxdepth = 50;

  cam.vfov = 20;
  cam.lookfrom = v3(13, 2, 3);
  cam.lookat = v3(0, 0, 0);
  cam.vup = v3(0, 1, 0);

  cam.defocusangle = 0.6;
  cam.focusdist = 10.0;

  render(&cam, &world);

  return 0;
}
