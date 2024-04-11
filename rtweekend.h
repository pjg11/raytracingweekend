#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define NTHREADS 8

double randomdouble(void);
double randomintervaldouble(double min, double max);

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

typedef struct hitrecord hitrecord;

typedef struct {
  int (*scatter)(ray in, hitrecord *rec, vec3 *attenuation, ray *scattered);
  vec3 albedo;
  double fuzz;
  double ir;
} material;

material lambertian(vec3 albedo);
material metal(vec3 albedo, double fuzz);
material dielectric(double ir);

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
  double aspectratio, vfov;
  int imagewidth, samplesperpixel, maxdepth;
  vec3 lookfrom, lookat, vup;
  double defocusangle, focusdist;

  int imageheight;
  vec3 center, pixel100loc, pixeldelu, pixeldelv, 
       u, w, v, defdisku, defdiskv;
} camera;

void render(camera *c, spherelist *world);

typedef struct {
  camera *c;
  spherelist *world;
  int i, j;
  vec3 *pixelcolor;
  
} threaddata;

#endif // RTWEEKEND_H
