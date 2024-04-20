#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NTHREADS 4

double randomdouble(void);

typedef struct {
  double x, y, z;
} vec3;

vec3 v3(double x, double y, double z);
vec3 v3add(vec3 v, vec3 w);
vec3 v3sub(vec3 v, vec3 w);
vec3 v3neg(vec3 v);
vec3 v3mul(vec3 v, vec3 w);
vec3 v3scale(vec3 v, double c);
double v3dot(vec3 v, vec3 w);
double v3length(vec3 v);
vec3 v3unit(vec3 v);
vec3 v3cross(vec3 v, vec3 w);
vec3 v3random(void);
vec3 v3randominterval(double min, double max);
vec3 v3randomunit(void);
vec3 v3randomunitdisk(void);

typedef struct {
  vec3 orig, dir;
} ray;

enum material { LAMBERTIAN, METAL, DIELECTRIC };

typedef struct {
  enum material type;
  union {
    struct lambertian {
      vec3 albedo;
    } lambertian;
    struct metal {
      vec3 albedo;
      double fuzz;
    } metal;
    struct dielectric {
      double ir;
    } dielectric;
  } data;
} material;

material lambertian(vec3 albedo);
material metal(vec3 albedo, double fuzz);
material dielectric(double ir);

typedef struct {
  vec3 point, normal;
  double t;
  int frontface;
  material mat;
} hitrecord;

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
  double aspectratio, vfov;
  int imagewidth, samplesperpixel, maxdepth;
  vec3 lookfrom, lookat, vup;
  double defocusangle, focusdist;

  int imageheight;
  vec3 center, pixel100loc, pixeldelu, pixeldelv, u, w, v, defdisku, defdiskv;
} camera;

void render(camera *c, spherelist *world);

typedef struct {
  camera *c;
  spherelist *world;
  int num, step;
  vec3 *pixels;
} threaddata;

#endif // RTWEEKEND_H
