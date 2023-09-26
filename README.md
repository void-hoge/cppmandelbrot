# CPPMANDELBROT
- Mandelbrot set renderer accelerated with AVX and OpenMP in C++.

![mandelbrot](samples/mandelbrot.gif)
![buddhabrot](samples/buddhabrot.png)

## Build and Execution
```
$ mkdir build
$ cd build
$ cmake ..
$ make
$ ./mandel > countmap
$ ../saveimage.py countmap mandelbrot.png
$ ./buddha > countmap
$ ../saveimage.py countmap buddhabrot.png
```

## Options
- AVX2
  - macro: `ENABLE_AVX` in `mandelbrot/mandelbrot.hpp`
- GMP
  - macro: `ENABLE_GMP` in `mandelbrot/mandelbrot.hpp`

## Algorithms
### Naive Algorithm
- The color of pixels in the rendering region is determined by the number of iterations it takes for the absolute value of the output of the Mandelbrot function to exceed a threshold (2.0).
- The algorithm that repeatedly applies the function in this manner is called the "escape time algorithm."
- `calc_mandelbrot_countmap` executes this algorithm for all pixels in the rendering region.
- It calculates 4 pixels in parallel using AVX2 and parallelized row-wise computations with OpenMP.

### Border Tracing / Edge Checking Algorithm
- The majority of the computation time in the naive algorithm is consumed by the pixels that belong to the Mandelbrot set.
- Hence, an algorithm that computes only the edge regions of the escape time can be considered.
- Specifically, it is as follows
  1. Calculate the escape time for pixels on the outer perimeter of the rendering region.
  2. Search for areas on the perimeter where the escape time is changing and push neighboring pixels to QUEUE.
  3. If the QUEUE is not empty, repeat the following.
	 1. Pop a pixel from QUEUE.
	 2. Calculate the escape time of the popped pixel.
	 3. Check whether there are boundaries between the popped pixel and the pixels for which escape times have been calculated in its 4 neighbours.
	 4. Push into the queue those pixels among the 8 neighbours of the popped pixel that are adjacent to a boundary, haven't been added to the QUEUE, and haven't their escape times calculated yet.
- When processed with AVX2, 4 pixels are popped from the QUEUE at a time.
- To process with OpenMP, it's necessary to have a QUEUE for each thread.
  - To achieve this, the rendering region is divided vertically and horizontally, and the above algorithm is executed.
- This algorithm is implemented as `calc_mandelbrot_boundary`, `calc_submandelbrot_boundary` and `pop4elms`.

### Buddhabrot
- Buddhabrot is a probability distribution of the trajectories of points that do not belong to the Mandelbrot set.
- To save computation for points belonging to the set, the Mandelbrot set is precomputed using the algorithm mentioned above.
  - Specifically, generate complex numbers randomly, and if the complex numbers at corners of that pixel belong to the Mandelbrot set, skip the computation.

## References
- https://en.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set
- https://geocities.restorativland.org/CapeCanaveral/5003/mandel.htm
- https://github.com/shapoco/accelbrot
- https://en.wikipedia.org/wiki/Buddhabrot

## Author
- Mugi Noda (void-hoge)

## License
- GPLv3
