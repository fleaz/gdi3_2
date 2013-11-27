#include "mandelbrot.h"
#include "stdio.h"
#include <stdlib.h>
#include <xmmintrin.h>
#include <math.h>


void
colorMapSSE(__m128 index, int maxIterations)
{
    __m128 r;
    __m128 g;
    __m128 b;

    __m128 maxIt = _mm_set_ps1(maxIterations);

    __m128 v0 = _mm_set_ps1(2.0);
    __m128 v1 = _mm_set_ps1(-1.0);
    __m128 v2 = _mm_set_ps1(0.5);
    __m128 v3 = _mm_set_ps1(0.492);
    __m128 v4 = _mm_set_ps1(0.877);
    __m128 v5 = _mm_set_ps1(1.704);
    __m128 v6 = _mm_set_ps1(0.509);
    __m128 v7 = _mm_set_ps1(0.194);
    __m128 v8 = _mm_set_ps1(255.0);

    __m128 less = _mm_set_ps1(-1);
    __m128 comp = _mm_cmpeq_ps(index,less);

    int x = _mm_movemask_ps(comp);

    __m128 y = _mm_set_ps1(0.2);
    __m128 u = _mm_add_ps(v1, _mm_mul_ps(v0, _mm_div_ps(index,maxIt)));
    __m128 v = _mm_sub_ps(v2, _mm_div_ps(index,maxIt));
    
    float arrY[4];
    _mm_store_ps(arrY, y);
    float arrU[4];
    _mm_store_ps(arrU, u);
    float arrV[4];
    _mm_store_ps(arrV, v);

    r = _mm_mul_ps(v8, _mm_add_ps(y, _mm_div_ps(v, v4)));
    b = _mm_mul_ps(v8, _mm_add_ps(y, _mm_div_ps(u, v3)));
    g = _mm_mul_ps(v8,_mm_sub_ps(_mm_sub_ps(_mm_mul_ps(v5,y), _mm_mul_ps(v6,r)),_mm_mul_ps(v7, b)));
    
    float rgbR[4];
    _mm_store_ps(rgbR, r);
    float rgbG[4];
    _mm_store_ps(rgbG, g);
    float rgbB[4];
    _mm_store_ps(rgbB, b);
        
    printf(" Y: %f\n", arrY[0]);
    printf(" U: %f\n", arrU[0]);
    printf(" V: %f\n", arrV[0]);

    if(x == 15) {
        for(int i = 0; i <= 9; i += 3) {
                printf("Schwarz\n");
        }
        return;
    }

    if (x % 2 == 1){
        rgbR[0] =(char) 0;
        rgbG[0] =(char) 0;
        rgbB[0] =(char) 0;
    }
    x = x >> 1;

    if (x % 2 == 1){
        rgbR[1] =(char) 0;
        rgbG[1] =(char) 0;
        rgbB[1] =(char) 0;
    }
    x = x >> 1;

    if (x % 2 == 1){
        rgbR[2] =(char) 0;
        rgbG[2] =(char) 0;
        rgbB[2] =(char) 0;
    }
    x = x >> 1;

    if (x % 2 == 1){
        rgbR[3] =(char) 0;
        rgbG[3] =(char) 0;
        rgbB[3] =(char) 0;
    }

    printf(" R: %d\n", (unsigned char) rgbR[0]);
    printf(" G: %d\n", (unsigned char) rgbG[0]);
    printf(" B: %d\n", (unsigned char) rgbB[0]);
    return;
}

void colorMapSeriel(int index, int maxIterations)
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
        g = (1.704 * y) - (0.509 * r) -( 0.194 * b);

        b = b * 255.0;
        r = r * 255.0;
        g = g * 255.0;
    
        printf(" Y: %f\n", y);
        printf(" U: %f\n", u);
        printf(" V: %f\n", v);
    }

    printf(" R: %d\n",(unsigned char) r);
    printf(" G: %d\n",(unsigned char) g);
    printf(" B: %d\n",(unsigned char) b);
}


void main(){
     int iteration = 100;
     int maxIterations = 100;

     __m128 itValues = _mm_set_ps1(iteration);
    
    printf("SSE:\n");
    colorMapSSE(itValues, maxIterations);
    printf("Seriel:\n");
    colorMapSeriel(iteration,maxIterations);
}
