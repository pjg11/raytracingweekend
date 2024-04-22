#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef NTHREADS
#define NTHREADS 4
#endif

extern const float pi, MAX_RAND;
unsigned int g_seed;

// From https://stackoverflow.com/a/26237777

inline void fast_srand(int seed) {
  g_seed = seed;
}

inline int fast_rand(void) {
  g_seed = (214013*g_seed+2531011);
  return (g_seed>>16)&(int)MAX_RAND;
}

inline float randomfloat(void) { return fast_rand() / (MAX_RAND + 1.0); }

#if defined(__ARM_NEON)

#include <arm_neon.h>
typedef float32x4_t vec3;

#else

typedef struct {
  float x, y, z;
} vec3;

#endif

float v3x(vec3 v);
float v3y(vec3 v);
float v3z(vec3 v);

vec3 v3(float x, float y, float z);
vec3 v3add(vec3 v, vec3 w);
vec3 v3sub(vec3 v, vec3 w);
vec3 v3neg(vec3 v);
vec3 v3mul(vec3 v, vec3 w);
vec3 v3scale(vec3 v, float c);
float v3dot(vec3 v, vec3 w);
float v3length(vec3 v);
vec3 v3unit(vec3 v);
vec3 v3cross(vec3 v, vec3 w);
vec3 v3random(void);
vec3 v3randominterval(float min, float max);
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
      float fuzz;
    } metal;
    struct dielectric {
      float ir;
    } dielectric;
  } data;
} material;

material lambertian(vec3 albedo);
material metal(vec3 albedo, float fuzz);
material dielectric(float ir);

typedef struct {
  vec3 point, normal;
  float t;
  int frontface;
  material mat;
} hitrecord;

typedef struct {
  vec3 center;
  float radius, rsquare;
  material mat;
} sphere;

sphere sp(vec3 center, float radius, material mat);

typedef struct {
  sphere *spheres;
  int n, max;
} spherelist;

void spherelistadd(spherelist *l, sphere s);

typedef struct {
  float aspectratio, vfov;
  int imagewidth, samplesperpixel, maxdepth;
  vec3 lookfrom, lookat, vup;
  float defocusangle, focusdist;

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
