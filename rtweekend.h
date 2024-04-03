#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
  double x, y, z;
} vec3;

vec3 v3(double x, double y, double z);

typedef struct {
  vec3 orig, dir;
} ray;

typedef struct hitrecord hitrecord;

typedef struct {
  int (*scatter)(ray in, hitrecord *rec, vec3 *attenuation, ray *scattered);
  vec3 albedo;
  double fuzz;
} material;

material lambertian(vec3 albedo);
material metal(vec3 albedo, double fuzz);

struct hitrecord {
  vec3 point, normal;
  double t;
  int frontface;
  material mat;
};

typedef struct {
  vec3 center;
  double radius;
  material mat;
} sphere;

sphere sp(vec3 center, double radius, material mat);

typedef struct {
  sphere *spheres;
  int n, max;
} spherelist;

void spherelistadd(spherelist *l, sphere s);

typedef struct {
  double aspectratio;
  int imageheight, samplesperpixel, maxdepth, imagewidth;
  vec3 center, pixel100loc, pixeldelu, pixeldelv;
} camera;

void render(camera *c, spherelist *world);

#endif // RTWEEKEND_H
