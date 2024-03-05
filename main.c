#include <stdio.h>

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

void writecolor(FILE *out, vec3 color) {
  double s = 255.99;
  fprintf(out, "%d %d %d\n", (int)(s * color.x), (int)(s * color.y), (int)(s * color.z));
}

int main(void) {
  int image_width = 256;
  int image_height = 256;
  
  printf("P3\n%d %d\n255\n",image_width,image_height);
  for(int j = 0; j < image_height; j++) {
	fprintf(stderr, "\rScanlines remaining: %d", image_height - j);
	for (int i = 0; i < image_width; i++) {
	  vec3 color = v3((double) i / (image_width - 1),
					  (double) j / (image_height - 1), 0.0);
	  writecolor(stdout, color);
	}
  }
  fprintf(stderr, "\rDone.                       \n");
}
