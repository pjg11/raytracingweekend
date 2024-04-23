#include "rtweekend.h"

/* --------------------------------------------------------------------------
 * Vector Functions
 * -------------------------------------------------------------------------- */

// ARM Neon version from
// https://blog.jacobvosmaer.nl/0022-ray-tracing-weekend
#if defined(__ARM_NEON)

float v3x(vec3 v) { return v[0]; }
float v3y(vec3 v) { return v[1]; }
float v3z(vec3 v) { return v[2]; }

vec3 v3(float x, float y, float z) {
  float32_t ar[4];
  ar[0] = x;
  ar[1] = y;
  ar[2] = z;
  ar[3] = 0;
  return vld1q_f32(ar);
}

vec3 v3add(vec3 v, vec3 w) { return vaddq_f32(v, w); }
vec3 v3sub(vec3 v, vec3 w) { return vsubq_f32(v, w); }
vec3 v3mul(vec3 v, vec3 w) { return vmulq_f32(v, w); }
vec3 v3scale(vec3 v, float c) { return vmulq_n_f32(v, c); }
vec3 v3neg(vec3 v) { return vnegq_f32(v); }
float v3dot(vec3 v, vec3 w) { return vaddvq_f32(vmulq_f32(v, w)); }

#else

static vec3 zero;

float v3x(vec3 v) { return v.x; }
float v3y(vec3 v) { return v.y; }
float v3z(vec3 v) { return v.z; }

vec3 v3(float x, float y, float z) {
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

vec3 v3mul(vec3 v, vec3 w) {
  v.x *= w.x;
  v.y *= w.y;
  v.z *= w.z;
  return v;
}

vec3 v3scale(vec3 v, float c) {
  v.x *= c;
  v.y *= c;
  v.z *= c;
  return v;
}

float v3dot(vec3 v, vec3 w) { return v.x * w.x + v.y * w.y + v.z * w.z; }
vec3 v3neg(vec3 v) { return v3sub(zero, v); }

#endif

float v3length(vec3 v) { return sqrtf(v3dot(v, v)); }

vec3 v3unit(vec3 v) { return v3scale(v, 1.0 / v3length(v)); }

vec3 v3cross(vec3 v, vec3 w) {
  return v3(v3y(v) * v3z(w) - v3z(v) * v3y(w),
            v3z(v) * v3x(w) - v3x(v) * v3z(w),
            v3x(v) * v3y(w) - v3y(v) * v3x(w));
}

vec3 v3random(void) { return v3(randomfloat(), randomfloat(), randomfloat()); }

vec3 v3randominterval(float min, float max) {
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
    vec3 v = v3(-1.0 + 2.0 * randomfloat(), -1.0 + 2.0 * randomfloat(), 0);
    if (v3dot(v, v) < 1)
      return v;
  }
}

int v3nearzero(vec3 v) {
  float s = 1e-8;
  return fabsf(v3x(v)) < s && fabsf(v3y(v)) < s && fabsf(v3z(v)) < s;
}

vec3 rayat(ray r, float t) { return v3add(r.orig, v3scale(r.dir, t)); }

ray r(vec3 from, vec3 to) {
  ray r;
  r.orig = from;
  r.dir = v3sub(to, from);
  return r;
}

/* --------------------------------------------------------------------------
 * Sphere Functions
 * -------------------------------------------------------------------------- */

sphere sp(vec3 center, float radius, material mat) {
  sphere s;
  s.center = center;
  s.radius = radius;
  s.rsquare = radius * radius;
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

// https://stackoverflow.com/a/76500658
float32x4x4_t zip(float32x4x4_t a) {
  float32x4x4_t b = {
      vzip1q_s32(a.val[0], a.val[2]), vzip2q_s32(a.val[0], a.val[2]),
      vzip1q_s32(a.val[1], a.val[3]), vzip2q_s32(a.val[1], a.val[3])};
  return b;
}

// ARM Neon version from
// https://www.jakef.science/posts/simd-parallelism/#four-spheres-at-a-time
int spherelisthit(spherelist *l, ray r, float tmin, float tmax,
                  hitrecord *rec) {
  int i, j;
  float closest = tmax, a = v3dot(r.dir, r.dir);
  sphere *spheres = l->spheres, *s;

#if defined(__ARM_NEON)
  vec3 rorigx = vdupq_n_f32(r.orig[0]), rorigy = vdupq_n_f32(r.orig[1]),
       rorigz = vdupq_n_f32(r.orig[2]);

  for (i = 0; i < l->n; i += 4) {
    float32x4x4_t xyzr = {spheres[i].center, spheres[i + 1].center,
                          spheres[i + 2].center, spheres[i + 3].center};
    xyzr = zip(zip(xyzr));

    xyzr.val[3][0] = spheres[i].rsquare;
    xyzr.val[3][1] = spheres[i + 1].rsquare;
    xyzr.val[3][2] = spheres[i + 2].rsquare;
    xyzr.val[3][3] = spheres[i + 3].rsquare;

    // vec3 oc = v3sub(r.orig, s.center);
    vec3 ocx = v3sub(rorigx, xyzr.val[0]), ocy = v3sub(rorigy, xyzr.val[1]),
         ocz = v3sub(rorigz, xyzr.val[2]);

    // float halfb = v3dot(oc, r.dir);
    vec3 halfbx = v3scale(ocx, v3x(r.dir)), halfby = v3scale(ocy, v3y(r.dir)),
         halfbz = v3scale(ocz, v3z(r.dir)),
         halfb = v3add(halfbz, v3add(halfbx, halfby));

    // float c = v3dot(oc, oc) - s.radius * s.radius;
    vec3 cx = v3mul(ocx, ocx), cy = v3mul(ocy, ocy), cz = v3mul(ocz, ocz),
         c = v3sub(v3add(cz, v3add(cx, cy)), xyzr.val[3]);

    // float discriminant = halfb * halfb - a * c;
    vec3 discrimiant = v3sub(v3mul(halfb, halfb), v3scale(c, a));

    if (vmaxvq_f32(discrimiant) < 0.0)
      continue;

    for (j = 0; j < 4; j++) {
      if (discrimiant[j] < 0)
        continue;

      float sqrtd = sqrtf(discrimiant[j]);
      float root = (-halfb[j] - sqrtd) / a;
      if (root <= tmin || root >= closest) {
        root = (-halfb[j] + sqrtd) / a;
        if (root <= tmin || root >= closest)
          continue;
      }

      rec->point = rayat(r, root);
      s = spheres + i + j;
      closest = root;
    }
    setfacenormal(r, v3scale(v3sub(rec->point, s->center), 1.0 / s->radius),
                  rec);
    rec->mat = s->mat;
  }

#else

  for (i = 0; i < l->n; i++) {
    sphere s = spheres[i];
    vec3 oc = v3sub(r.orig, s.center);
    float halfb = v3dot(oc, r.dir);
    float c = v3dot(oc, oc) - s.rsquare;
    float discriminant = halfb * halfb - a * c;

    if (discriminant < 0)
      continue;

    float sqrtd = sqrtf(discriminant);
    float root = (-halfb - sqrtd) / a;
    if (root <= tmin || root >= closest) {
      root = (-halfb + sqrtd) / a;
      if (root <= tmin || root >= closest)
        continue;
    }

    rec->t = root;
    rec->point = rayat(r, rec->t);
    setfacenormal(r, v3scale(v3sub(rec->point, s.center), 1.0 / s.radius), rec);
    rec->mat = s.mat;
    closest = rec->t;
  }

#endif
  return closest != tmax ? 1 : 0;
}

/* --------------------------------------------------------------------------
 * Material Functions
 * -------------------------------------------------------------------------- */

vec3 reflect(vec3 v, vec3 n) { return v3sub(v, v3scale(n, 2 * v3dot(v, n))); }

vec3 refract(vec3 uv, vec3 n, float etaioveretat) {
  float costheta = fminf(v3dot(v3neg(uv), n), 1.0);
  vec3 routperp = v3scale(v3add(uv, v3scale(n, costheta)), etaioveretat),
       routparallel = v3scale(n, -sqrtf(fabsf(1 - v3dot(routperp, routperp))));
  return v3add(routperp, routparallel);
}

float reflectance(float cosine, float refidx) {
  float r0 = (1.0 - refidx) / (1.0 + refidx);
  r0 *= r0;
  return r0 + (1.0 - r0) * powf((1.0 - cosine), 5);
}

material lambertian(vec3 albedo) {
  material mat;
  mat.type = LAMBERTIAN;
  mat.data.lambertian.albedo = albedo;
  return mat;
}

material metal(vec3 albedo, float fuzz) {
  material mat;
  mat.type = METAL;
  mat.data.metal.albedo = albedo;
  mat.data.metal.fuzz = fuzz > 1 ? 1 : fuzz;
  return mat;
}

material dielectric(float ir) {
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
    float refractionratio = rec->frontface ? 1.0 / data.ir : data.ir;
    *attenuation = v3(1, 1, 1);
    vec3 unitdir = v3unit(in.dir);
    float costheta = fminf(v3dot(v3neg(unitdir), rec->normal), 1.0),
          sintheta = sqrtf(1.0 - costheta * costheta);
    scattered->orig = rec->point;
    scattered->dir =
        refractionratio * sintheta > 1.0 ||
                reflectance(costheta, refractionratio) > randomfloat()
            ? reflect(unitdir, rec->normal)
            : refract(unitdir, rec->normal, refractionratio);
    return 1;
  }

  default:
    return 0;
  }
}

/* --------------------------------------------------------------------------
 * Camera Functions
 * -------------------------------------------------------------------------- */

vec3 pixelsamplesquare(camera *c) {
  float px = -0.5 + randomfloat();
  float py = -0.5 + randomfloat();
  return v3add(v3scale(c->pixeldelu, px), v3scale(c->pixeldelv, py));
}

vec3 defocusdisksample(camera *c) {
  vec3 p = v3randomunitdisk();
  return v3add(c->center, v3add(v3scale(c->defdisku, v3x(p)),
                                v3scale(c->defdiskv, v3y(p))));
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
  vec3 color, attenuation = {1, 1, 1};

  while (depth > 0) {
    if (spherelisthit(world, r, 0.001, INFINITY, &rec)) {
      if (scatter(r, &rec, &color, &r)) {
        attenuation = v3mul(attenuation, color);
        depth -= 1;
      }
    } else {
      vec3 dir = v3unit(r.dir);
      float a = 0.5 * (v3y(dir) + 1);
      return v3mul(attenuation, v3add(v3scale(v3(1, 1, 1), 1 - a),
                                      v3scale(v3(0.5, 0.7, 1), a)));
    }
  }
  return v3(0, 0, 0);
}

void writecolor(FILE *out, vec3 color, int samplesperpixel) {
  int s = 256;
  color = v3scale(color, 1.0 / samplesperpixel);
  color = v3(clamp(sqrtf(v3x(color))), clamp(sqrtf(v3y(color))),
             clamp(sqrtf(v3z(color))));

  fprintf(out, "%d %d %d\n", (int)(s * v3x(color)), (int)(s * v3y(color)),
          (int)(s * v3z(color)));
}

void initialize(camera *c) {
  float viewportheight, viewportwidth, h, defocusradius;
  vec3 viewportu, viewportv, viewportupperleft;

  c->imageheight = c->imagewidth / c->aspectratio;
  if (c->imageheight < 1)
    c->imageheight = 1;

  c->center = c->lookfrom;

  h = tanf(degtorad(c->vfov) / 2);
  viewportheight = 2 * h * c->focusdist;
  viewportwidth = viewportheight * ((float)c->imagewidth / c->imageheight);

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

  defocusradius = c->focusdist * tanf(degtorad(c->defocusangle / 2));
  c->defdisku = v3scale(c->u, defocusradius);
  c->defdiskv = v3scale(c->v, defocusradius);
}

/* --------------------------------------------------------------------------
 * Render Functions
 * -------------------------------------------------------------------------- */

void *linesrender(void *args) {
  int i, j, sample;
  threaddata *td = (threaddata *)args;
  camera *c = td->c;

  for (j = td->num; j < c->imageheight; j += td->step) {
    for (i = 0; i < c->imagewidth; i++) {
      vec3 pixelcolor = {0};
      for (sample = 0; sample < c->samplesperpixel; sample++) {
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
    int err = pthread_create(threads + k, NULL, linesrender, args);
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

  for (i = 0; i < c->imageheight; i++) {
    for (k = 0; k < c->imagewidth; k++)
      writecolor(stdout, *(pixels + i * c->imagewidth + k), c->samplesperpixel);
  }

  free(pixels);
  free(world->spheres);
}
