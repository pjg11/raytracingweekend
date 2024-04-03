#include "rtweekend.h"

vec3 v3(double x, double y, double z) {
  vec3 v;
  v.x = x;
  v.y = y;
  v.z = z;
  return v;
}

vec3 v3add(vec3 v, vec3 w) {
  v.x += w.x;
  v.y += w.y;
  v.z += w.z;
  return v;
}

vec3 v3sub(vec3 v, vec3 w) {
  v.x -= w.x;
  v.y -= w.y;
  v.z -= w.z;
  return v;
}

vec3 v3scale(vec3 v, double c) {
  v.x *= c;
  v.y *= c;
  v.z *= c;
  return v;
}

double v3dot(vec3 v, vec3 w) { return v.x * w.x + v.y * w.y + v.z * w.z; }

double v3length(vec3 v) { return sqrt(v.x * v.x + v.y * v.y + v.z * v.z); }

vec3 v3unit(vec3 v) { return v3scale(v, 1.0 / v3length(v)); }

double randomdouble(void) { return rand() / (RAND_MAX + 1.0); }

vec3 v3random(void) {
  return v3(randomdouble(), randomdouble(), randomdouble());
}

vec3 v3randominterval(double min, double max) {
  return v3add(v3(min, min, min), v3scale(v3random(), max - min));
}

sphere sp(vec3 center, double radius, material mat) {
  sphere s;
  s.center = center;
  s.radius = radius;
  s.mat = mat;
  return s;
}

void spherelistadd(spherelist *l, sphere s) {
  if (l->n == l->max) {
    l->max = l->max ? 2 * l->max : 1;
    l->spheres = realloc(l->spheres, l->max * sizeof(*l->spheres));
  }
  l->spheres[l->n++] = s;
}
