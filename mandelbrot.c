/*
 * This file contains the implementation of the interface. You can change
 * everything in this file (but it still has to contain the function
 * generateMandelbrot). For example, for a proper SSE implementation of
 * the color mapping, you may want to process multiple pixels at once.
 */
#include "mandelbrot.h"
#include "stdio.h"

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
    float r,g,b;

    if (index == -1){
        r = 0.0;
        g = 0.0;
        b = 0.0;

    }
    else {
        float y = 0.2;
        float u = -1.0 + 2.0 * ((float)index / (float)maxIterations);
        float v = 0.5 - ((float)index / (float)maxIterations);

        b = y + u / 0.492;
        r = y + v / 0.877;
        g = 1.704 * y - 0.509 * r - 0.194 * b;

        b = b * 255.0;
        r = r * 255.0;
        g = g * 255.0;
    }

    //printf ("it: %d,r: %f, g: %f b: %f \n",index,r, g, b);
    //printf ("it: %d,r: %d, g: %d b: %d \n",index,(unsigned char)r, (unsigned char)g, (unsigned char)b);
    //printf("---\n");
    color[0] = (unsigned char) r;
    color[1] = (unsigned char) g;
    color[2] = (unsigned char) b;
}


/*
 *Addition von komplexen Zahlen
 */

complex float addComplex(complex float a, complex float b) {
    complex float result = (crealf(a) + crealf(b)) + (cimagf(a) + cimagf(b)) * I;
    return result;
}

/*
 *Multiplikation von komplexen Zahlen
 */

complex float mulComplex(complex float a, complex float b) {
    complex float result = ((crealf(a)*crealf(b)) - (cimagf(a)*cimagf(b))) + ((crealf(a)*cimagf(b)) + (cimagf(a)*crealf(b))) * I;
    return result;
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
    //printf("real: %f, imag: %f\n",crealf(c),cimagf(c));
    complex float z = 0+0*I;
    int iteration = 0;

    while ((sqrt(pow(crealf(z), 2.0) + pow(cimagf(z), 2.0)) <= RADIUS) && (iteration < maxIterations)) {
        z = addComplex(mulComplex(z,z),c);
        iteration = iteration + 1;
    }

    //printf("Real: %f, Imag: %f\n",crealf(z),cimag(z));

    if (iteration < maxIterations) {
        int mu = log(log(sqrt(pow(crealf(z), 2.0) + pow(cimagf(z), 2.0))) / log(2.0));
        iteration = iteration + 1 - mu;
    }

    if (iteration == maxIterations) {
        iteration = -1;
    }

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
    float widthPiece = (crealf(upperLeft) - crealf(lowerRight)) / width;
    float heightPiece = (cimagf(upperLeft) - cimagf(lowerRight)) / height;
    widthPiece = fabs(widthPiece);
    heightPiece = fabs(heightPiece);
    //printf("%f %f \n",widthPiece, heightPiece);

    for(int y = 0; y < height; y++) {
        for(int x = 0; x < width; x++) {
            double real = crealf(upperLeft) + (widthPiece * x);
            double imag = cimagf(upperLeft) - (heightPiece * y);
            //printf ("x: %d, y: %d\n",x,y);
            //printf("real: %f, imag: %f\n",real,imag);
            //printf("---\n");

            int value = testEscapeSeriesForPoint(real+imag*I, maxIterations, 0);
            //printf("Value: %d\n",value);
            int offset = (y * width + x) * 3;
            colorMapYUV(value, maxIterations, image + offset);
        }
    }

    return image;
}

