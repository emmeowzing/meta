#ifndef IMG_H_
#define IMG_H_

#include <stdint.h>

typedef struct
{
    uint32_t width;
    uint32_t height;
    int bit_depth;
    png_bytep *datap;   // data
} rTuple;

const unsigned char signature[8] = { 0x89, 0x50, 0x4e, 0x47, 
                                     0x0d, 0x0a, 0x1a, 0x0a };

struct coord
{
    uint32_t left;
    uint32_t right;
};

typedef struct
{
    double *avV;
    double **var_covar_mat;
} varTuple;

/*
typedef int bool;
#define true 1
#define false 0
*/

#define PNG_BYTES_TO_CHECK 8
#define CHANNELS 3
#define SAMPLES 1e4

#endif // IMG_H_
