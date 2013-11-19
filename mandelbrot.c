/*
 * This file contains the implementation of the interface. You can change
 * everything in this file (but it still has to contain the function
 * generateMandelbrot). For example, for a proper SSE implementation of
 * the color mapping, you may want to process multiple pixels at once.
 */
#include "mandelbrot.h"

/*
 * Calculates a color mapping for a given iteration number by exploiting the
 * YUV color space. Returns the color as 8-bit unsigned char per channel (RGB).
 *
 * Arguments:
 *	index - Number of iterations that resulted from iterating the Mandelbrot series.
 *	maxIterations - Parameter that was also used for series iteration.
 *
 * Returns:
 *	The associated color as an array of 3 8-bit unsigned char values.
 */
void
colorMapYUV(int index, int maxIterations, unsigned char* color)
{
    float y = 0.2;
    float u = -1.0+2.0*(index/maxIterations);
    float v = 0.5 - (index/maxIterations);

    float b = y + u/0.492;
    float r = y + v/0.877;
    float g = 1.704*y-0.509*r-0.194*b;

    //printf ("r: %f, g: %f b: %f \n", r,g,b);
    color[0] = (char) r;
    color[1] = (char) g;
    color[2] = (char) b;
}

/*
 * Executes the complex series for a given parameter c for up to maxIterations
 * and saves the last series component in last.
 *
 * Arguments:
 *  c - Additive component (complex number) for Mandelbrot series.
 *	maxIterations - Maximum number of iterations that are executed to determine a series' boundedness
 *	last - Pointer to a complex float number that can be used for storing the last component in a series - useful for color mapping
 *
 * Returns:
 *	The number of iterations that were executed before the series escaped our
 *	circle or - if the point is part of the Mandelbrot set a special
 *	(user-defined) value.
 */
int
testEscapeSeriesForPoint(complex float c, int maxIterations, complex float * last)
{
    complex float z = 0+0*I;
    int iteration = 0;

    while ((sqrt(pow(creal(z), 2) + pow(cimag(z), 2)) <= RADIUS) && (iteration < maxIterations)) {
        z = z*z+c;
        iteration = iteration + 1;
    }

    if (iteration < maxIterations) {
        int mu = log(log(sqrt(pow(creal(z), 2) + pow(cimag(z), 2))) / log(2));
        iteration = iteration + 1 - mu;
    }

    if (iteration == maxIterations) {
        iteration = -1;
    }

    last = &z;
    return iteration;
}



/*
 * Generates an image of a Mandelbrot set.
 */
unsigned char *
generateMandelbrot(
    complex float upperLeft, 
    complex float lowerRight, 
    int maxIterations, 
    int width, 
    int height)
{
    // Allocate image buffer, row-major order, 3 channels.
    unsigned char *image = malloc(height * width * 3);
    // TODO: Generate and color map the image.
    for(int y = 0; y < height; y++) {
        for(int x = 0; x < width; x++) {
            // TODO: Call testEscapeSeriesForPoint() here.
            // TODO: Color map the result.
            int iter = testEscapeSeriesForPoint(x+y*I, maxIterations, 0);
            printf("iter: %d \n", iter);

            int offset = (y * width + x) * 3;
            colorMapYUV(iter, maxIterations, image + offset);
        }
    }

    return image;
}

