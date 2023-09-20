# CPPMANDELBROT
- Mandelbrot set renderer accelerated with AVX and OpenMP in C++.

![sample](samples/mandelbrot.gif)

## Build and Execution
```
$ mkdir build
$ cd build
$ cmake ..
$ make
$ ../zoom.py cppmandelbrot # Generate images included in samples.
```

## Algorithms
### Naive Algorithm
- Coloring is determined by the number of iterations of the Mandelbrot function until the complex number assigned to each pixel in the region exceeds a threshold.
- This algorithm is called the escape time algorithm.
- `calc_mandelbrot_countmap` runs this algorith for every pixel in the region.
- This function is parallelized using AVX and OpenMP.

### Border Tracing / Edge Checking Algorithm
- Most of the computation cost of the naive algorithm is consumed in computing pixels belonging to the Mandelbrot set.
- Therefore, one possible appoach is to execute the escape time algorithm only at the edge pixels of escape times.
- Specifically, it is as follows
  1. Calculate the escape time of all pixels on perimeter of the region.
  2. Find the boundaries of the escape time on perimeter and push it to the QUEUE.
  3. If the QUEUE is not empty, repeat the following.
	 1. Pop a pixel from QUEUE.
	 2. Calculate the escape time of the popped pixel.
	 3. Check the boundaries between the pixel and calculated pixels of its 4-neighbour.
	 4. Push the unknown 8-neighbours of the pixel that touch the boundary to the QUEUE.
- The AVX version pops pixels from the QUEUE in groups of 4.
- Divide the area vertically and horizontally for parallelization.
  - This is the reason for the vertical stripes in the images included in samples.
- This algorithm is implemented as `calc_mandelbrot_boundary`, `calc_submandelbrot_boundary` and `pop4elms`.

## References
- https://en.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set
- https://geocities.restorativland.org/CapeCanaveral/5003/mandel.htm
- https://github.com/shapoco/accelbrot

## Author
- Mugi Noda (void-hoge)

## License
- GPLv3
