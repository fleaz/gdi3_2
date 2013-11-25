/*
 * This file contains the implementation of the interface. You can change
 * everything in this file (but it still has to contain the function
 * generateMandelbrot). For example, for a proper SSE implementation of
 * the color mapping, you may want to process multiple pixels at once.
 */
#include "mandelbrot.h"
#include "stdio.h"
#include <stdlib.h>
#include <xmmintrin.h>

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
colorMapYUV(int* index, int maxIterations, unsigned char* color)
{

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
 *Berechnet den Betrag einer komplexen Zahl
 */

float absComplex(complex float a) {
    float result = sqrt(crealf(a)*crealf(a) + cimagf(a)*cimagf(a));
    return result;
}

/*
 *
 */
__m128 sseABS(__m128 real, __m128 imag){
    real = _mm_mul_ps(real, real);
    imag = _mm_mul_ps(imag, imag);
    real = _mm_add_ps(real, imag);

    return _mm_sqrt_ps(real);
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
int*
testEscapeSeriesForPoint(float r1,float r2,float r3,float r4, float i1, int maxIterations, complex float * last)
{
    //printf("real: %f, imag: %f\n",crealf(c),cimagf(c));
    //printf("C1 r:%f, i:%f\n",r1,i1);
    //printf("C2 r:%f, i:%f\n",r2,i1);
    //printf("C3 r:%f, i:%f\n",r3,i1);
    //printf("C4 r:%f, i:%f\n",r4,i1);
    __m128 realC = _mm_set_ps(r1, r2, r3, r4);
    __m128 imagC = _mm_set_ps(i1, i1, i1, i1);

    __m128 realZ = _mm_set_ps(0.0f, 0.0f, 0.0f, 0.0f);
    __m128 imagZ = _mm_set_ps(0.0f, 0.0f, 0.0f, 0.0f);

    __m128 tempZ = _mm_set_ps(0.0f, 0.0f, 0.0f, 0.0f);

    __m128 comp = _mm_set_ps(0.0f, 0.0f, 0.0f, 0.0f);
    __m128 radius = _mm_set_ps(RADIUS,RADIUS,RADIUS,RADIUS);

    int it[] = {0,0,0,0};
    int loop = 0;

    while (loop <= maxIterations) {
        comp = _mm_cmple_ps (sseABS(realZ, imagZ), radius);
        int x = _mm_movemask_ps(comp);

        //printf("%d %d %d %d %d \n",x,it[0],it[1],it[2],it[3]);
        if (x == 0){
            break;
        }
        if (x % 2 == 1 && it[0] < maxIterations){
            it[0]++;
        }
        x = x >> 1;

        if (x % 2 == 1 && it[1] < maxIterations){
            it[1]++;
        }
        x = x >> 1;

        if (x % 2 == 1 && it[2] < maxIterations){
            it[2]++;
        }
        x = x >> 1;

        if (x % 2 == 1 && it[3] < maxIterations){
            it[3]++;
        }

        tempZ = _mm_sub_ps(_mm_mul_ps(realZ, realZ),_mm_mul_ps(imagZ, imagZ));
        imagZ = _mm_add_ps(_mm_mul_ps(realZ, imagZ),_mm_mul_ps(imagZ,realZ));
        realZ = tempZ;

        realZ = _mm_add_ps(realZ, realC);
        imagZ = _mm_add_ps(imagC, imagZ);

        //printf("R: %f %f %f %f \n",rZ[0],rZ[1],rZ[2],rZ[3]);
        //printf("I: %f %f %f %f \n\n",iZ[0],iZ[1],iZ[2],iZ[3]);

        loop++;
    }

    float rZ[4];
    _mm_store_ps (rZ, realZ);
    float iZ[4];
    _mm_store_ps (iZ, imagZ);

    for(int i=0; i < 4; i++){
        if( it[i] < maxIterations){
            int mu = log(log(absComplex(rZ[i]+iZ[i]*I)) / log(2.0)) / log(2);
            it[i] = it[i] + 1 - mu;
        }
        else{
            it[i] = -1;
        }
    }

    //printf("%d %d %d %d \n",it[0],it[1],it[2],it[3]);
    //mu berechnen
    // -1 wenn it = maxIt
    int *p = it[0];
    return p;
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
    float widthPiece = fabs((crealf(upperLeft) - crealf(lowerRight)) / width);
    float heightPiece = fabs((cimagf(upperLeft) - cimagf(lowerRight)) / height);
    //widthPiece = fabs(widthPiece);
    //heightPiece = fabs(heightPiece);
    //printf("%f %f \n",widthPiece, heightPiece);

    for(int y = 0; y < height; y++) {
        for(int x = 0; x < width; x++) {

            //printf ("x: %d, y: %d\n",x,y);
            //printf("real: %f, imag: %f\n",real,imag);
            //printf("---\n");

            float r1 = crealf(upperLeft) + (widthPiece * x);
            float i1 = cimagf(upperLeft) - (heightPiece * y);
            x++;
            float r2 = crealf(upperLeft) + (widthPiece * x);
            x++;
            float r3 = crealf(upperLeft) + (widthPiece * x);
            x++;
            float r4 = crealf(upperLeft) + (widthPiece * x);


            //printf("C1 r:%f, i:%f\n",r1,i1);
            //printf("C2 r:%f, i:%f\n",r2,i1);
            //printf("C3 r:%f, i:%f\n",r3,i1);
            //printf("C4 r:%f, i:%f\n",r4,i1);

            int value = testEscapeSeriesForPoint(r1, r2, r3, r4, i1, maxIterations, 0);
            int offset = (y * width + x) * 3;
            colorMapYUV(value, maxIterations, image + offset);
        }
    }

    return image;
}

