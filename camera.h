#ifndef CAMERA_H
#define CAMERA_H

#include "rtweekend.h"

typedef struct {
  double aspectratio;
  int imageheight, samplesperpixel, maxdepth, imagewidth;
  vec3 center, pixel100loc, pixeldelu, pixeldelv;
} camera;

void render(camera *c, spherelist *world);

#endif // CAMERA_H
