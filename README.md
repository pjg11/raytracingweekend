# Ray Tracing in (More Than) One Weekend

A attempt at the [Ray Tracing in One Weekend](https://raytracing.github.io/books/RayTracingInOneWeekend.html) book, in C.

## Usage

```
$ make
$ ./main > section14.ppm
```

The resulting image is the cover image of the book: 1200x675px, 500 samples/pixel. Running `./main -small` generates a smaller image for testing purposes - 400x225px, 100 samples/pixel.

The original, single-threaded implementation ran in 18m 41s with optimization flags enabled. Further optimizations reduced the time to 3m 31s, a 5.3x speedup. A writeup is in the works, so here's a summary:

* Multi-threading (4 threads) - 62% faster
* Implementation changes - 18% faster
  * Changing from function pointers to enum + union for `scatter()`
  * Change `raycolor()` from recursive to iterative
  * Use custom random function instead of stdlib `rand()`
* SIMD - 37% faster
  * Add [ARM Neon](https://community.arm.com/arm-community-blogs/b/operating-systems-blog/posts/arm-neon-programming-quick-reference) instructions in `vec3` functions
  * Vectorize `spherelisthit()` further to calculate 4 spheres at once

## Resources
- [jacobvosmaer/raytracingweekend](https://github.com/jacobvosmaer/raytracingweekend/) + accompaniying [blog post](http://blog.jacobvosmaer.nl/0022-ray-tracing-weekend/)
- [jfeintzeig/ray_tracer](https://github.com/JFeintzeig/ray_tracer) + accompanying [blog post](https://www.jakef.science/posts/simd-parallelism/)
- [ARM Neon Instruction Set](https://developer.arm.com/architectures/instruction-sets/intrinsics/#f:@navigationhierarchiessimdisa=[Neon])
