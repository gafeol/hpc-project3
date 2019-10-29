
#include <stdlib.h>
#include <stdio.h>

#include <unistd.h>
#include <time.h>
#include <sys/time.h>

#include "pngwriter.h"
#include "consts.h"


unsigned long get_time ()
{
    struct timeval tp;
    gettimeofday (&tp, NULL);
    return tp.tv_sec * 1000000 + tp.tv_usec;
}

int main (int argc, char** argv)
{
    printf("Mandel Parallel - runned with %d threads\n", omp_get_max_threads());
	png_data* pPng = png_create (IMAGE_WIDTH, IMAGE_HEIGHT);
	
	double x, y, x2, y2, cx, cy;
	
	double fDeltaX = (MAX_X - MIN_X) / (double) IMAGE_WIDTH;
	double fDeltaY = (MAX_Y - MIN_Y) / (double) IMAGE_HEIGHT;
	
	long nTotalIterationsCount = 0;
	unsigned long nTimeStart = get_time ();
	
	long i, j, n;
        
	// do the calculation
    #pragma omp parallel for schedule(dynamic) private(i, j, x, y, cx, cy, x2, y2, n) shared(nTotalIterationsCount)
	for (j = 0; j < IMAGE_HEIGHT; j++)
	{
        cx = MIN_X;
        cy = MIN_Y + fDeltaY*j;
        long partial_n = 0;
		for (i = 0; i < IMAGE_WIDTH; i++)
		{			
			// compute the orbit z, f(z), f²(z), f³(z), ...
			// count the iterations until the orbit leaves the circle |z|=2.
			// stop if the number of iterations exceeds the bound MAX_ITERS.

            double zx = 0, zy = 0;
            n = 0;

            double ox = zx, oy = zy;
            double oxx = zx*zx, oyy = zy*zy, oxy = zx*zy;
            while(oxx + oyy < 2. && n < MAX_ITERS){
                zx = oxx - oyy + cx;
                zy = oxy*2 + cy;
                n++;

                oxx = zx*zx; 
                oyy = zy*zy; 
                oxy = zx*zy;
            }
            partial_n += n;
            // n indicates if the point belongs to the mandelbrot set
			// plot the number of iterations at point (i, j)
			int c = ((long) n * 255) / MAX_ITERS;
			png_plot (pPng, i, j, c, c, c);
			
			cx += fDeltaX;
		}
		
        #pragma omp critical
        {nTotalIterationsCount += partial_n;}
	}
	
	unsigned long nTimeEnd = get_time ();
	
	// print benchmark data
	printf ("Total time:                 %g ms\n", (nTimeEnd - nTimeStart) / 1000.0);
	printf ("Image size:                 %ld x %ld = %ld Pixels\n",
		(long) IMAGE_WIDTH, (long) IMAGE_HEIGHT, (long) (IMAGE_WIDTH * IMAGE_HEIGHT));
	printf ("Total number of iterations: %ld\n", nTotalIterationsCount);
	printf ("Avg. time per pixel:        %g µs\n", (nTimeEnd - nTimeStart) / (double) (IMAGE_WIDTH * IMAGE_HEIGHT));
	printf ("Avg. time per iteration:    %g µs\n", (nTimeEnd - nTimeStart) / (double) nTotalIterationsCount);
	printf ("Iterations/second:          %g\n", nTotalIterationsCount / (double) (nTimeEnd - nTimeStart) * 1e6);
	// assume there are 8 floating point operations per iteration
	printf ("MFlop/s:                    %g\n", nTotalIterationsCount * 8.0 / (double) (nTimeEnd - nTimeStart));
	
	png_write (pPng, "mandel.png");
	return 0;
}
