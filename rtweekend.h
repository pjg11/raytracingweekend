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

#ifndef MAXRAND
#define MAXRAND 0x7FFF
#endif

// fastsrand() and fastrand() from https://stackoverflow.com/a/26237777
unsigned int gseed;

inline void fastsrand(unsigned int seed) { gseed = seed; }

inline int fastrand(void) {
  gseed = (214013 * gseed + 2531011);
  return (gseed >> 16) & (int)MAXRAND;
}

inline float randomfloat(void) { return fastrand() / (MAXRAND + 1.0); }

inline float clamp(float x) {
  float tmin = 0.000, tmax = 0.999;
  return x < tmin ? tmin : x > tmax ? tmax : x;
}

inline float degtorad(float deg) { return M_PI * deg / 180.0; }

#if defined(__ARM_NEON)

#include <arm_neon.h>
typedef float32x4_t vec3;

#else

typedef struct {
  float x, y, z;
} vec3;

#endif

vec3 v3(float x, float y, float z);
vec3 v3sub(vec3 v, vec3 w);
vec3 v3mul(vec3 v, vec3 w);
vec3 v3random(void);
vec3 v3randominterval(float min, float max);
float v3length(vec3 v);

typedef struct {
  vec3 orig, dir;
} ray;

ray r(vec3 from, vec3 to);

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
