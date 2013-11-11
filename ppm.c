#include <stdlib.h>
#include <stdio.h>

#include "ppm.h"

void exportPPM(const char *filename, struct PPM *image)
{
	FILE *img;
    img = fopen(filename, "w");
    if (img == NULL) {
        printf("Error saving image!\n");
        return;
    }

    fprintf(img, "P3 %d %d 255 ", image->width, image->height);;
    for(int y = image->height - 1; y >= 0; y--) {
        for(int x = 0; x < image->width * 3; x++) {
            fprintf(img, "%d\n",(int)(image->data[y * image->width * 3 + x]));
		}
	}
	fclose(img);
}
