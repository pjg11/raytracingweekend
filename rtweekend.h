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
double v3dot(vec3 v, vec3 w);
double v3length(vec3 v);
vec3 v3unit(vec3 v);

typedef struct {
  vec3 orig, dir;
} ray;

typedef struct {
  vec3 center;
  double radius;
} sphere;

sphere sp(vec3 center, double radius);

typedef struct {
  sphere *spheres;
  int n, max;
} spherelist;

void spherelistadd(spherelist *l, sphere s);

#endif // RTWEEKEND_H
