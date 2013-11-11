#ifndef PPM_HEADER
#define PPM_HEADER

/*
 * Representation of an image.
 */
struct PPM {
	int width;
	int height;
	unsigned char* data;
};

/*
 * Saves a PPM file (3-channels, each 8-bit unsigned char) to disk into a file named path.
 * 
 * Arguments:
 *	filename - The filename including path for the image file.
 *	image - PPM structure containing width, height and actual data of an image.
 */
void
exportPPM(const char *filename, struct PPM *image);

#endif /* PPM_HEADER */
