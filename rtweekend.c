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

int v3nearzero(vec3 v) {
	double s = 1e-8;
	return (fabs(v.x) < s) && (fabs(v.y) < s) && (fabs(v.z) < s);
}

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

vec3 randominunitsphere(void) {
  while (1) {
    vec3 p = v3randominterval(-1, 1);
    if (v3length(p) < 1) {
      return p;
    }
  }
}

vec3 randomunitvector(void) { return v3unit(randominunitsphere()); }

int lambertianscatter(ray in, hitrecord *rec, vec3 *attenuation, ray *scattered, vec3 albedo) {
	vec3 scatterdir = v3add(rec->normal, randomunitvector());
	if(v3nearzero(scatterdir)) {
		scatterdir = rec->normal;
	}
	scattered->orig = rec->point;
	scattered->dir = scatterdir;
	*attenuation = albedo;
	return 1;
}

material lambertian(vec3 albedo) {
	material mat;
	mat.scatter = &lambertianscatter;
	mat.albedo = albedo;
	return mat;
}
