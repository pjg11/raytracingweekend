#ifndef RTWEEKEND_H
#define RTWEEKEND_H

#include <math.h>
#include <stdlib.h>

typedef struct {
  double x, y, z;
} vec3;

vec3 v3(double x, double y, double z);
vec3 v3add(vec3 v, vec3 w);
vec3 v3sub(vec3 v, vec3 w);
vec3 v3scale(vec3 v, double c);
vec3 v3mul(vec3 v, vec3 w);
double v3dot(vec3 v, vec3 w);
double v3length(vec3 v);
vec3 v3unit(vec3 v);
double randomdouble(void);
vec3 v3random(void);
vec3 v3randominterval(double min, double max);
int v3nearzero(vec3 v);
vec3 reflect(vec3 v, vec3 w);

typedef struct {
  vec3 orig, dir;
} ray;

typedef struct hitrecord hitrecord;

typedef struct {
  int (*scatter)(ray in, hitrecord *rec, vec3 *attenuation, ray *scattered,
                 vec3 albedo);
  vec3 albedo;
  double fuzz;
} material;

material lambertian(vec3 albedo);
material metal(vec3 albedo, double fuzz);

struct hitrecord {
  vec3 point, normal;
  material mat;
  double t;
  int frontface;
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

vec3 randomunitvector(void);

#endif // RTWEEKEND_H
