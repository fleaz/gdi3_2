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
#include <math.h>

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
colorMapYUV(__m128 index, int maxIterations, unsigned char* color)
{
    //Initialisiere Variablen
    __m128 r;
    __m128 g;
    __m128 b;

    __m128 maxIt = _mm_set_ps1(maxIterations);

    __m128 v0 = _mm_set_ps1(2.0);
    __m128 v1 = _mm_set_ps1(-1.0);
    __m128 v2 = _mm_set_ps1(0.5);
    __m128 v3 = _mm_set_ps1(0.492);
    __m128 v4 = _mm_set_ps1(0.877);

    __m128 v8 = _mm_set_ps1(255.0);
    __m128 v9 = _mm_set_ps1(0.39393);
    __m128 v10 = _mm_set_ps1(0.58081);

    //Berechnung von Y, U und V
    __m128 y = _mm_set_ps1(0.2);
    __m128 u = _mm_add_ps(v1, _mm_mul_ps(v0, _mm_div_ps(index,maxIt)));
    __m128 v = _mm_sub_ps(v2, _mm_div_ps(index,maxIt));

    //Berechnung von R, B, G
    r = _mm_mul_ps(v8, _mm_add_ps(y, _mm_div_ps(v, v4)));
    b = _mm_mul_ps(v8, _mm_add_ps(y, _mm_div_ps(u, v3)));

    g = _mm_mul_ps(v9, u);
    g = _mm_sub_ps(y, g);
    g = _mm_sub_ps(g, _mm_mul_ps(v10,v));

    //Lege Werte aus r,g,b in arrays, sowie die iterations werte
    float rgbR[4];
    _mm_store_ps(rgbR, r);
    float rgbG[4];
    _mm_store_ps(rgbG, g);
    float rgbB[4];
    _mm_store_ps(rgbB, b);

    float iterationsA[4];
    _mm_store_ps(iterationsA, index);

    //Setze Farben, für custom Wert -1 setze r=g=b=0
    for(int i = 0; i < 4; i++) {
        if(iterationsA[i] == -1) {
            color[i * 3] = 0;
            color[i * 3 + 1] = 0;
            color[i * 3 + 2] = 0;
        } else {
            color[i * 3] = (unsigned char)rgbR[i];
            color[i * 3 + 1] = (unsigned char)rgbG[i];
            color[i * 3 + 2] = (unsigned char)rgbB[i];
        }
    }
    return;
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
 *Betrag einer Komplexen Zahl, bekommt zwei sse-Register jeweils mit Real- und Imaginärteil
 *
 *Returns:
 *Betrag der vier "komplexen Zahlen" in einem sse-Register
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
__m128
testEscapeSeriesForPoint(__m128 realC, __m128 imagC, int maxIterations, complex float * last)
{

    //Initialisiere Variablen
    __m128 realZ = _mm_set_ps1(0.0);
    __m128 imagZ = _mm_set_ps1(0.0);

    __m128 tempZ = _mm_set_ps1(0.0);

    __m128 comp = _mm_set_ps1(0.0);
    __m128 radius = _mm_set_ps1(RADIUS);

    float it[] = {0.0,0.0,0.0,0.0};
    int loop = 0;

    while (loop <= maxIterations) {
        //Abfrage abs(z) <= Radius
        comp = _mm_cmple_ps (sseABS(realZ, imagZ), radius);
        int x = _mm_movemask_ps(comp);


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

        //Berechne neues Real und Imaginärteile von Z
        tempZ = _mm_sub_ps(_mm_mul_ps(realZ, realZ),_mm_mul_ps(imagZ, imagZ));
        imagZ = _mm_add_ps(_mm_mul_ps(realZ, imagZ),_mm_mul_ps(imagZ,realZ));
        realZ = tempZ;

        realZ = _mm_add_ps(realZ, realC);
        imagZ = _mm_add_ps(imagC, imagZ);

        loop++;
    }

    float rZ[4];
    _mm_store_ps (rZ, realZ);
    float iZ[4];
    _mm_store_ps (iZ, imagZ);

    //Berechne/Verwende mu für die jeweiligen Iterationswerte
    for(int i=0; i < 4; i++){
        if( it[i] < maxIterations){
            complex float complax = rZ[i] + iZ[i] * I;
            double mu = log(log((double)absComplex(complax)) / log(2.0)) / log(2.0);
            int muNeu = (int) (mu + 0.5);
            if(mu == 1) {
                it[i] = it[i] + 1 - muNeu;
            }

        }
        else{
            it[i] = -1;
        }
    }

    return _mm_set_ps(it[0], it[1], it[2], it[3]);
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


    for(int y = 0; y < height; y++) {
        for(int x = 0; x < width; x += 4) {

            int offset = (y * width + x) * 3;

            //Berechne aktuelle Stelle im Bild
            __m128 sseWidth = _mm_set_ps1(WIDTH);
            __m128 sseHeight = _mm_set_ps1(HEIGHT);

            __m128 realLow = _mm_set_ps1(crealf(lowerRight));
            __m128 imagLow = _mm_set_ps1(cimagf(lowerRight));
            __m128 realUpper = _mm_set_ps1(crealf(upperLeft));
            __m128 imagUpper = _mm_set_ps1(cimagf(upperLeft));

            __m128 xReal = _mm_set_ps(x, x+1, x+2, x+3);

            xReal = _mm_div_ps(xReal, sseWidth);
            xReal = _mm_mul_ps(xReal, _mm_sub_ps(realLow,realUpper));
            xReal = _mm_add_ps(xReal, realUpper);


            __m128 imag = _mm_set_ps1(y);

            imag = _mm_div_ps(imag, sseHeight);
            imag = _mm_mul_ps(imag, _mm_sub_ps(imagUpper, imagLow));
            imag = _mm_add_ps(imag, imagLow);


            //Rufe testEscapeSeriesForPoint auf, Rückgabewerte sind die jeweiligen Iterationswerte
            //parallel für 4 Pixel gleichzeiti
            __m128 itValues = testEscapeSeriesForPoint(xReal, imag, maxIterations, 0);

            //Es werde Farbe \o/
            colorMapYUV(itValues, maxIterations, image + offset);
        }
    }

    return image;
}

