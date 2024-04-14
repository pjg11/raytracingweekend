#include "rtweekend.h"

double randomdouble(void) { return rand() / (RAND_MAX + 1.0); }

double clamp(double x) {
  double tmin = 0.000, tmax = 0.999;
  return x < tmin ? tmin : x > tmax ? tmax : x;
}
double degtorad(double deg) { return M_PI * deg / 180.0; }

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

double v3dot(vec3 v, vec3 w) { return v.x * w.x + v.y * w.y + v.z * w.z; }
double v3length(vec3 v) { return sqrt(v3dot(v, v)); }
vec3 v3unit(vec3 v) { return v3scale(v, 1.0 / v3length(v)); }

vec3 v3cross(vec3 v, vec3 w) {
  return v3(v.y * w.z - v.z * w.y, v.z * w.x - v.x * w.z,
            v.x * w.y - v.y * w.x);
}

vec3 v3random(void) {
  return v3(randomdouble(), randomdouble(), randomdouble());
}

vec3 v3randominterval(double min, double max) {
  return v3add(v3(min, min, min), v3scale(v3random(), max - min));
}

vec3 v3randomunit(void) {
  while (1) {
    vec3 v = v3randominterval(-1, 1);
    if (v3dot(v, v) < 1)
      return v3unit(v);
  }
}

vec3 v3randomunitdisk(void) {
  while (1) {
    vec3 v = v3(-1.0 + 2.0 * randomdouble(), -1.0 + 2.0 * randomdouble(), 0);
    if (v3dot(v, v) < 1)
      return v;
  }
}

int v3nearzero(vec3 v) {
  double s = 1e-8;
  return fabs(v.x) < s && fabs(v.y) < s && fabs(v.z) < s;
}

vec3 rayat(ray r, double t) { return v3add(r.orig, v3scale(r.dir, t)); }

ray r(vec3 from, vec3 to) {
  ray r;
  r.orig = from;
  r.dir = v3sub(to, from);
  return r;
}

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

vec3 refract(vec3 uv, vec3 n, double etaioveretat) {
  double costheta = fmin(v3dot(v3neg(uv), n), 1.0);
  vec3 routperp = v3scale(v3add(uv, v3scale(n, costheta)), etaioveretat),
       routparallel = v3scale(n, -sqrt(fabs(1 - v3dot(routperp, routperp))));
  return v3add(routperp, routparallel);
}

double reflectance(double cosine, double refidx) {
  double r0 = (1.0 - refidx) / (1.0 + refidx);
  r0 *= r0;
  return r0 + (1.0 - r0) * pow((1.0 - cosine), 5);
}

material lambertian(vec3 albedo) {
  material mat;
  mat.type = LAMBERTIAN;
  mat.data.lambertian.albedo = albedo;
  return mat;
}

material metal(vec3 albedo, double fuzz) {
  material mat;
  mat.type = METAL;
  mat.data.metal.albedo = albedo;
  mat.data.metal.fuzz = fuzz > 1 ? 1 : fuzz;
  return mat;
}

material dielectric(double ir) {
  material mat;
  mat.type = DIELECTRIC;
  mat.data.dielectric.ir = ir;
  return mat;
}

int scatter(ray in, hitrecord *rec, vec3 *attenuation, ray *scattered) {
  material mat = rec->mat;
  switch (mat.type) {
  case LAMBERTIAN: {
    struct lambertian data = mat.data.lambertian;
    vec3 dir = v3add(rec->normal, v3randomunit());
    if (v3nearzero(dir))
      dir = rec->normal;
    scattered->orig = rec->point;
    scattered->dir = dir;
    *attenuation = data.albedo;
    return 1;
  }

  case METAL: {
    struct metal data = mat.data.metal;
    vec3 reflected = reflect(v3unit(in.dir), rec->normal);
    scattered->orig = rec->point;
    scattered->dir = v3add(reflected, v3scale(v3randomunit(), data.fuzz));
    *attenuation = data.albedo;
    return 1;
  }

  case DIELECTRIC: {
    struct dielectric data = mat.data.dielectric;
    double refractionratio = rec->frontface ? 1.0 / data.ir : data.ir;
    *attenuation = v3(1, 1, 1);
    vec3 unitdir = v3unit(in.dir);
    double costheta = fmin(v3dot(v3neg(unitdir), rec->normal), 1.0),
           sintheta = sqrt(1.0 - costheta * costheta);
    scattered->orig = rec->point;
    scattered->dir =
        refractionratio * sintheta > 1.0 ||
                reflectance(costheta, refractionratio) > randomdouble()
            ? reflect(unitdir, rec->normal)
            : refract(unitdir, rec->normal, refractionratio);
    return 1;
  }

  default:
    return 0;
  }
}

// Camera functions
vec3 pixelsamplesquare(camera *c) {
  double px = -0.5 + randomdouble();
  double py = -0.5 + randomdouble();
  return v3add(v3scale(c->pixeldelu, px), v3scale(c->pixeldelv, py));
}

vec3 defocusdisksample(camera *c) {
  vec3 p = v3randomunitdisk();
  return v3add(c->center,
               v3add(v3scale(c->defdisku, p.x), v3scale(c->defdiskv, p.y)));
}

ray getray(camera *c, int i, int j) {
  vec3 pixelcenter = v3add(c->pixel100loc, v3add(v3scale(c->pixeldelu, i),
                                                 v3scale(c->pixeldelv, j)));
  vec3 pixelsample = v3add(pixelcenter, pixelsamplesquare(c));

  return r((c->defocusangle <= 0) ? c->center : defocusdisksample(c),
           pixelsample);
}

vec3 raycolor(ray r, int depth, spherelist *world) {
  hitrecord rec;
  vec3 black = {0};

  if (depth <= 0)
    return black;

  if (spherelisthit(world, r, 0.001, INFINITY, &rec)) {
    ray scattered;
    vec3 attenuation;
    if (scatter(r, &rec, &attenuation, &scattered))
      return v3mul(attenuation, raycolor(scattered, depth - 1, world));
    else
      return black;
  } else {
    vec3 dir = v3unit(r.dir);
    double a = 0.5 * (dir.y + 1);
    return v3add(v3scale(v3(1, 1, 1), 1 - a), v3scale(v3(0.5, 0.7, 1), a));
  }
}

void writecolor(FILE *out, vec3 color, int samplesperpixel) {
  int s = 256;
  color = v3scale(color, 1.0 / samplesperpixel);
  color = v3(clamp(sqrt(color.x)), clamp(sqrt(color.y)), clamp(sqrt(color.z)));

  fprintf(out, "%d %d %d\n", (int)(s * color.x), (int)(s * color.y),
          (int)(s * color.z));
}

void initialize(camera *c) {
  double viewportheight, viewportwidth, h, defocusradius;
  vec3 viewportu, viewportv, viewportupperleft;

  c->imageheight = c->imagewidth / c->aspectratio;
  if (c->imageheight < 1)
    c->imageheight = 1;

  c->center = c->lookfrom;

  h = tan(degtorad(c->vfov) / 2);
  viewportheight = 2 * h * c->focusdist;
  viewportwidth = viewportheight * ((double)c->imagewidth / c->imageheight);

  c->w = v3unit(v3sub(c->lookfrom, c->lookat));
  c->u = v3unit(v3cross(c->vup, c->w));
  c->v = v3cross(c->w, c->u);

  viewportu = v3scale(c->u, viewportwidth);
  viewportv = v3scale(v3neg(c->v), viewportheight);
  c->pixeldelu = v3scale(viewportu, 1.0 / c->imagewidth);
  c->pixeldelv = v3scale(viewportv, 1.0 / c->imageheight);

  viewportupperleft = v3sub(v3sub(v3sub(c->center, v3scale(c->w, c->focusdist)),
                                  v3scale(viewportu, 0.5)),
                            v3scale(viewportv, 0.5));

  c->pixel100loc =
      v3add(viewportupperleft, v3scale(v3add(c->pixeldelu, c->pixeldelv), 0.5));

  defocusradius = c->focusdist * tan(degtorad(c->defocusangle / 2));
  c->defdisku = v3scale(c->u, defocusradius);
  c->defdiskv = v3scale(c->v, defocusradius);
}

void *linesrender(void *args) {
  threaddata *td = (threaddata *)args;
  camera *c = td->c;
  for (int j = td->num; j < c->imageheight; j += td->step) {
    for (int i = 0; i < c->imagewidth; i++) {
      vec3 pixelcolor = {0};
      for (int sample = 0; sample < c->samplesperpixel; ++sample) {
        pixelcolor = v3add(
            pixelcolor, raycolor(getray(td->c, i, j), c->maxdepth, td->world));
      }

      *(td->pixels + j * c->imagewidth + i) = pixelcolor;
    }
  }
  return 0;
}

void render(camera *c, spherelist *world) {
  int i, k, npixels;
  vec3 *pixels;
  pthread_t threads[NTHREADS];
  threaddata threadargs[NTHREADS];

  initialize(c);
  npixels = c->imageheight * c->imagewidth;
  pixels = calloc(npixels, sizeof(*pixels));

  printf("P3\n%d %d\n255\n", c->imagewidth, c->imageheight);

  for (k = 0; k < NTHREADS; k++) {
    threaddata *args = threadargs + k;
    args->c = c;
    args->world = world;
    args->num = k;
    args->step = NTHREADS;
    args->pixels = pixels;
    int err = pthread_create(threads + k, NULL, linesrender, threadargs + k);
    if (err) {
      fprintf(stderr, "error: pthread_create\n");
      return;
    }
  }

  for (k = 0; k < NTHREADS; k++) {
    int err = pthread_join(threads[k], NULL);
    if (err) {
      fprintf(stderr, "error: pthread_join %d\n", err);
      return;
    }
  }

  for (i = 0; i < npixels; i++) {
    writecolor(stdout, pixels[i], c->samplesperpixel);
  }

  free(pixels);
}
