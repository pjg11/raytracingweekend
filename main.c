#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
  double x, y, z;
} vec3;

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

void writecolor(FILE *out, vec3 color) {
  double s = 255.99;
  fprintf(out, "%d %d %d\n", (int)(s * color.x), (int)(s * color.y),
          (int)(s * color.z));
}

typedef struct {
  vec3 orig, dir;
} ray;

vec3 at(ray r, double t) { return v3add(r.orig, v3scale(r.dir, t)); }

typedef struct {
  vec3 center;
  double radius;
} sphere;

typedef struct {
  vec3 point, normal;
  double t;
  int frontface;
} hitrecord;

void setfacenormal(ray r, vec3 outwardnormal, hitrecord *rec) {
  rec->frontface = v3dot(r.dir, outwardnormal) < 0;
  rec->normal = rec->frontface ? outwardnormal : v3scale(outwardnormal, -1);
}

int spherehit(sphere s, ray r, double tmin, double tmax, hitrecord *rec) {
  vec3 oc = v3sub(r.orig, s.center);
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
  vec3 outwardnormal = v3scale(v3sub(rec->point, s.center), 1 / s.radius);
  setfacenormal(r, outwardnormal, rec);

  return 1;
}

sphere makesphere(vec3 center, double radius) {
  sphere s;
  s.center = center;
  s.radius = radius;
  return s;
}

typedef struct {
  sphere *spheres;
  int n, max;
} spherelist;

void spherelistadd(sphere s, spherelist *l) {
  if (l->n == l->max) {
    l->max = l->max ? 2 * l->max : 1;
    l->spheres = realloc(l->spheres, l->max * sizeof(*l->spheres));
  }
  l->spheres[l->n++] = s;
}

int spherelisthit(spherelist *l, ray r, double tmin, double tmax,
                  hitrecord *rec) {
  hitrecord tmp;
  int hitanything = 0;
  double closest = tmax;
  for (int i = 0; i < l->n; i++) {
    if (spherehit(l->spheres[i], r, tmin, closest, &tmp)) {
      hitanything = 1;
      closest = tmp.t;
      *rec = tmp;
    }
  }
  return hitanything;
}

vec3 ray_color(ray r, spherelist *world) {
  hitrecord rec;
  if (spherelisthit(world, r, 0, INFINITY, &rec)) {
    return v3scale(v3add(rec.normal, v3(1, 1, 1)), 0.5);
  }

  vec3 dir = v3unit(r.dir);
  double a = 0.5 * (dir.y + 1.0);
  return v3add(v3scale(v3(1.0, 1.0, 1.0), 1.0 - a),
               v3scale(v3(0.5, 0.7, 1.0), a));
}

int main(void) {

  double aspectratio = 16.0 / 9.0, focallength = 1.0, viewportheight = 2.0,
         viewportwidth;

  int imagewidth = 400, imageheight = imagewidth / aspectratio;

  if (imageheight < 1) {
    imageheight = 1;
  }

  spherelist world = {0};
  spherelistadd(makesphere(v3(0, 0, -1), 0.5), &world);
  spherelistadd(makesphere(v3(0, -100.5, -1), 100), &world);

  viewportwidth = viewportheight * ((double)(imagewidth) / imageheight);

  vec3 cameracenter = {0}, viewportu = v3(viewportwidth, 0, 0),
       viewportv = v3(0, -viewportheight, 0),
       pixeldelu = v3scale(viewportu, 1.0 / imagewidth),
       pixeldelv = v3scale(viewportv, 1.0 / imageheight),
       viewportupperleft =
           v3sub(v3sub(v3sub(cameracenter, v3(0, 0, focallength)),
                       v3scale(viewportu, 0.5)),
                 v3scale(viewportv, 0.5)),
       pixel00loc =
           v3add(viewportupperleft, v3scale(v3add(pixeldelu, pixeldelv), 0.5));

  printf("P3\n%d %d\n255\n", imagewidth, imageheight);

  for (int j = 0; j < imageheight; j++) {
    fprintf(stderr, "\rScanlines remaining: %d", imageheight - j);
    for (int i = 0; i < imagewidth; i++) {
      vec3 pixelcenter = v3add(
          pixel00loc, v3add(v3scale(pixeldelu, i), v3scale(pixeldelv, j)));
      vec3 raydirection = v3sub(pixelcenter, cameracenter);
      ray r = {cameracenter, raydirection};
      writecolor(stdout, ray_color(r, &world));
    }
  }

  fprintf(stderr, "\rDone.                       \n");
}
