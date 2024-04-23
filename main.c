#include "rtweekend.h"

int main(int argc, char *argv[]) {
  spherelist world = {0};
  camera cam = {1, 100, 10, 10, 90, {0, 0, -1}, {0, 0, 0}, {0, 1, 0}, 0, 10};
  int a = -11, b = -11, small = 0;

  fastsrand(1713755957);

  if (argc == 2 && !strcmp(argv[1], "-small"))
    small = 1;
  else if (argc > 1) {
    fprintf(stderr, "Usage: main [-small]\n");
    return 1;
  }

  material ground = lambertian(v3(0.5, 0.5, 0.5));
  spherelistadd(&world, sp(v3(0, -1000, 0), 1000, ground));

  for (a = -11; a < 11; a++) {
    for (b = -11; b < 11; b++) {
      float choosemat = randomfloat();
      float radius = 0.2;
      vec3 center =
          v3(a + 0.9 * randomfloat(), radius, b + 0.9 * randomfloat());

      if (v3length(v3sub(center, v3(4, 0.2, 0))) > 0.9) {
        material spherematerial;
        if (choosemat < 0.8) {
          // diffuse
          spherematerial = lambertian(v3mul(v3random(), v3random()));
          spherelistadd(&world, sp(center, radius, spherematerial));
        } else if (choosemat < 0.95) {
          // metal
          spherematerial = metal(v3randominterval(0.5, 1), 0.5 * randomfloat());
          spherelistadd(&world, sp(center, radius, spherematerial));
        } else {
          // glass
          spherematerial = dielectric(1.5);
          spherelistadd(&world, sp(center, radius, spherematerial));
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
  cam.imagewidth = small ? 400 : 1200;
  cam.samplesperpixel = small ? 100 : 500;
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
