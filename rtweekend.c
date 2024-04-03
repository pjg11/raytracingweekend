#include "rtweekend.h"

double randomdouble(void) { return rand() / (RAND_MAX + 1.0); }

double clamp(double x) {
  double tmin = 0.000, tmax = 0.999;
  if (x < tmin)
    return tmin;
  if (x > tmax)
    return tmax;
  return x;
}

double togamma(double linear) { return sqrt(linear); }

// Vector functions

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

vec3 v3neg(vec3 v) {
  vec3 zero = {0};
  return v3sub(zero, v);
}

vec3 v3mul(vec3 v, vec3 w) {
  v.x *= w.x;
  v.y *= w.y;
  v.z *= w.z;
  return v;
}

vec3 v3scale(vec3 v, double c) {
  v.x *= c;
  v.y *= c;
  v.z *= c;
  return v;
}

double v3length(vec3 v) { return sqrt(v.x * v.x + v.y * v.y + v.z * v.z); }
double v3dot(vec3 v, vec3 w) { return v.x * w.x + v.y * w.y + v.z * w.z; }
vec3 v3unit(vec3 v) { return v3scale(v, 1.0 / v3length(v)); }

vec3 v3random(void) {
  return v3(randomdouble(), randomdouble(), randomdouble());
}

vec3 v3randominterval(double min, double max) {
  return v3add(v3(min, min, min), v3scale(v3random(), max - min));
}

vec3 v3randomunit(void) {
  vec3 v;
  do
    v = v3randominterval(-1, 1);
  while (v3dot(v, v) >= 1);
  return v3unit(v);
}

int v3nearzero(vec3 v) {
  double s = 1e-8;
  return fabs(v.x) < s && fabs(v.y) < s && fabs(v.z) < s;
}

vec3 rayat(ray r, double t) { return v3add(r.orig, v3scale(r.dir, t)); }

// Sphere functions
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

void setfacenormal(ray r, vec3 outwardnormal, hitrecord *rec) {
  rec->frontface = v3dot(r.dir, outwardnormal) < 0;
  rec->normal = rec->frontface ? outwardnormal : v3neg(outwardnormal);
}

int spherehit(sphere s, ray r, double tmin, double tmax, hitrecord *rec) {
  vec3 oc = v3sub(r.orig, s.center);
  double a = v3dot(r.dir, r.dir);
  double halfb = v3dot(oc, r.dir);
  double c = v3dot(oc, oc) - s.radius * s.radius;
  double discriminant = halfb * halfb - a * c;

  if (discriminant < 0)
    return 0;

  double sqrtd = sqrt(discriminant);
  double root = (-halfb - sqrtd) / a;
  if (root <= tmin || root >= tmax) {
    root = (-halfb + sqrtd) / a;
    if (root <= tmin || root >= tmax)
      return 0;
  }

  rec->t = root;
  rec->point = rayat(r, rec->t);
  setfacenormal(r, v3scale(v3sub(rec->point, s.center), 1.0 / s.radius), rec);
  rec->mat = s.mat;
  return 1;
}

int spherelisthit(spherelist *l, ray r, double tmin, double tmax,
                  hitrecord *rec) {
  int hitanything = 0, i;
  double closest = tmax;

  for (i = 0; i < l->n; i++) {
    if (spherehit(l->spheres[i], r, tmin, closest, rec)) {
      hitanything++;
      closest = rec->t;
    }
  }
  return hitanything;
}

// Material functions
vec3 reflect(vec3 v, vec3 n) { return v3sub(v, v3scale(n, 2 * v3dot(v, n))); }

int lambertianscatter(ray in, hitrecord *rec, vec3 *attenuation,
                      ray *scattered) {
  material mat = rec->mat;
  vec3 scatterdir = v3add(rec->normal, v3randomunit());
  if (v3nearzero(scatterdir))
    scatterdir = rec->normal;

  scattered->orig = rec->point;
  scattered->dir = scatterdir;
  *attenuation = mat.albedo;
  return 1;
}

material lambertian(vec3 albedo) {
  material mat = {0};
  mat.scatter = &lambertianscatter;
  mat.albedo = albedo;
  return mat;
}

int metalscatter(ray in, hitrecord *rec, vec3 *attenuation, ray *scattered) {
  material mat = rec->mat;
  vec3 reflected = reflect(v3unit(in.dir), rec->normal);

  scattered->orig = rec->point;
  scattered->dir = v3add(reflected, v3scale(v3randomunit(), mat.fuzz));
  *attenuation = mat.albedo;
  return 1;
}

material metal(vec3 albedo, double fuzz) {
  material mat = {0};
  mat.scatter = &metalscatter;
  mat.albedo = albedo;
  mat.fuzz = fuzz > 1 ? 1 : fuzz;
  return mat;
}

// Camera functions
vec3 pixelsamplesquare(camera *c) {
  double px = -0.5 + randomdouble();
  double py = -0.5 + randomdouble();
  return v3add(v3scale(c->pixeldelu, px), v3scale(c->pixeldelv, py));
}

ray getray(camera *c, int i, int j) {
  vec3 pixelcenter = v3add(c->pixel100loc, v3add(v3scale(c->pixeldelu, i),
                                                 v3scale(c->pixeldelv, j)));
  vec3 pixelsample = v3add(pixelcenter, pixelsamplesquare(c));
  ray r = {c->center, v3sub(pixelsample, c->center)};
  return r;
}

vec3 raycolor(ray r, int depth, spherelist *world) {
  hitrecord rec;
  vec3 black = {0};

  if (depth <= 0)
    return black;

  if (spherelisthit(world, r, 0.001, INFINITY, &rec)) {
    ray scattered;
    vec3 attenuation;
    if (rec.mat.scatter(r, &rec, &attenuation, &scattered))
      return v3mul(attenuation, raycolor(scattered, depth - 1, world));
    return black;
  } else {
    vec3 dir = v3unit(r.dir);
    double a = 0.5 * (dir.y + 1);
    return v3add(v3scale(v3(1, 1, 1), 1 - a), v3scale(v3(0.5, 0.7, 1), a));
  }
}

void writecolor(FILE *out, vec3 color, int samplesperpixel) {
  int s = 256;
  double scale = 1.0 / samplesperpixel, r = togamma(color.x * scale),
         g = togamma(color.y * scale), b = togamma(color.z * scale);

  fprintf(out, "%d %d %d\n", (int)(s * clamp(r)), (int)(s * clamp(g)),
          (int)(s * clamp(b)));
}

void initialize(camera *c) {
  double focallength = 1, viewportheight = 2, viewportwidth;
  vec3 viewportu, viewportv, viewportupperleft;

  c->imageheight = c->imagewidth / c->aspectratio;

  if (c->imageheight < 1) {
    c->imageheight = 1;
  }

  viewportwidth = viewportheight * ((double)(c->imagewidth) / c->imageheight);
  viewportu = v3(viewportwidth, 0, 0);
  viewportv = v3(0, -viewportheight, 0);

  c->pixeldelu = v3scale(viewportu, 1.0 / c->imagewidth);
  c->pixeldelv = v3scale(viewportv, 1.0 / c->imageheight);

  viewportupperleft = v3sub(
      v3sub(v3sub(c->center, v3(0, 0, focallength)), v3scale(viewportu, 0.5)),
      v3scale(viewportv, 0.5));

  c->pixel100loc =
      v3add(viewportupperleft, v3scale(v3add(c->pixeldelu, c->pixeldelv), 0.5));
}

void render(camera *c, spherelist *world) {
  int i, j, sample;
  initialize(c);

  printf("P3\n%d %d\n255\n", c->imagewidth, c->imageheight);

  for (j = 0; j < c->imageheight; j++) {
    fprintf(stderr, "\rScanlines remaining: %d ", c->imageheight - j);
    for (i = 0; i < c->imagewidth; i++) {
      vec3 pixelcolor = v3(0, 0, 0);
      for (sample = 0; sample < c->samplesperpixel; ++sample) {
        ray r = getray(c, i, j);
        pixelcolor = v3add(pixelcolor, raycolor(r, c->maxdepth, world));
      }
      writecolor(stdout, pixelcolor, c->samplesperpixel);
    }
  }

  fprintf(stderr, "\rDone.                       \n");
}
