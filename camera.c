#include "camera.h"
#include <stdio.h>

double clamp(double x) {
  double tmin = 0.000, tmax = 0.999;
  if (x < tmin)
    return tmin;
  if (x > tmax)
    return tmax;
  return x;
}

double togamma(double linear) { return sqrt(linear); }

void writecolor(FILE *out, vec3 color, int samplesperpixel) {
  int s = 256;
  double scale = 1.0 / samplesperpixel, r = togamma(color.x * scale),
         g = togamma(color.y * scale), b = togamma(color.z * scale);

  fprintf(out, "%d %d %d\n", (int)(s * clamp(r)), (int)(s * clamp(g)),
          (int)(s * clamp(b)));
}

vec3 at(ray r, double t) { return v3add(r.orig, v3scale(r.dir, t)); }

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

void setfacenormal(ray r, vec3 outwardnormal, hitrecord *rec) {
  rec->frontface = v3dot(r.dir, outwardnormal) < 0;
  rec->normal = rec->frontface ? outwardnormal : v3scale(outwardnormal, -1);
}

int spherehit(sphere s, ray r, double tmin, double tmax, hitrecord *rec) {
  vec3 oc = v3sub(r.orig, s.center), outwardnormal;
  double a = v3dot(r.dir, r.dir);
  double halfb = v3dot(oc, r.dir);
  double c = v3dot(oc, oc) - s.radius * s.radius;

  double discriminant = halfb * halfb - a * c;
  if (discriminant < 0) {
    return 0;
  }
  double sqrtd = sqrt(discriminant);

  double root = (-halfb - sqrtd) / a;
  if (root <= tmin || tmax <= root) {
    root = (-halfb - sqrtd) / a;
    if (root <= tmin || tmax <= root) {
      return 0;
    }
  }

  rec->t = root;
  rec->point = at(r, rec->t);
  outwardnormal = v3scale(v3sub(rec->point, s.center), 1 / s.radius);
  setfacenormal(r, outwardnormal, rec);
  rec->mat = s.mat;
  return 1;
}

int spherelisthit(spherelist *l, ray r, double tmin, double tmax,
                  hitrecord *rec) {
  hitrecord tmp;
  int hitanything = 0, i;
  double closest = tmax;

  for (i = 0; i < l->n; i++) {
    if (spherehit(l->spheres[i], r, tmin, closest, &tmp)) {
      hitanything = 1;
      closest = tmp.t;
      *rec = tmp;
    }
  }
  return hitanything;
}

vec3 randomonhemisphere(vec3 normal) {
  vec3 onunitsphere = randomunitvector();
  if (v3dot(onunitsphere, normal) > 0.0)
    return onunitsphere;
  else
    return v3scale(onunitsphere, -1);
}

vec3 raycolor(ray r, int depth, spherelist *world) {
  hitrecord rec;
  vec3 dir;
  double a;
  vec3 black = {0};

  if (depth <= 0) {
    return black;
  }

  if (spherelisthit(world, r, 0.001, INFINITY, &rec)) {
    ray scattered;
    vec3 attenuation;
    if (rec.mat.scatter(r, &rec, &attenuation, &scattered, rec.mat.albedo))
      return v3mul(attenuation, raycolor(scattered, depth - 1, world));
    return black;
  }

  dir = v3unit(r.dir);
  a = 0.5 * (dir.y + 1.0);
  return v3add(v3scale(v3(1.0, 1.0, 1.0), 1.0 - a),
               v3scale(v3(0.5, 0.7, 1.0), a));
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
