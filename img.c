/*
 * Open a PNG and write a `mahalanobisified' version.
 *
 * Probably a hundred ways this could be simplified &| improved - you're 
 * welcome to hit me up with ideas, I'd like to see how fast we can make
 * this!
 * 
 *                                    Author: Brandon Doyle (01.03.17)
 * Genius who suggested we try this in C/C++: Timothy Van Slyke
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <png.h>
#include <time.h>
#include "img.h"

void
print_var_covar(double **mat)
{
    // print the matrix nicely in the terminal
    uint32_t i, j;
    for (i = 0; i < CHANNELS; ++i)
    {
        for (j = 0; j < CHANNELS; ++j)
            printf("%5.6lf  ", mat[i][j]);

        printf("\n");
    }
}

double **
fabricate2DArray(uint32_t X, uint32_t Y)
{
    /* Buuld a 2d array of doubles */
    uint32_t i;
    double **newMat;
    newMat = malloc(X * sizeof(double *));

    for (i = 0; i < X; ++i)
    {
        newMat[i] = malloc(Y * sizeof(double));
        memset(newMat[i], 0, Y * sizeof(double));
    }

    return newMat;
}

int **
fabricate2DintArray(uint32_t X, uint32_t Y)
{
    /* Build a 2d array of ints */
    uint32_t i;
    int **newMat = malloc(X * sizeof(int *));
    
    for (i = 0; i < X; ++i)
    {
        newMat[i] = malloc(Y * sizeof(int));
        memset(newMat[i], 0, Y * sizeof(int));
    }

    return newMat;
}

void
destroy2DintArray(uint32_t X, int **arr)
{
    /* free a 2D integer array */
    uint32_t i;
    for (i = 0; i < X; ++i) free(arr[i]);
    free(arr);
}

void 
destroy2DPNGbyteArray(uint32_t X, png_byte **arr)
{
    uint32_t i;
    for (i = 0; i < X; ++i) free(arr[i]);
    free(arr);
}

void
destroy2DArray(uint32_t X, double **arr)
{
    /* free a 2D double array */
    uint32_t i;
    for (i = 0; i < X; ++i) free(arr[i]);
    free(arr);
}

void
outer_product(double **res, double *vec)
{
    /* multiply an array (vector) with itself as an outer product and store in 
       `res` (watch out for segfaults here too - make sure `len` is the correct
        value) */
    uint32_t i, j;

    // it's actually only necessary to solve for n(n+1) of the matrix entries
    // here, since $\Sigma$ is a symmetric matrix
    for (i = 0; i < CHANNELS; ++i)
        for (j = 0; j < CHANNELS; ++j)
            res[i][j] = vec[i] * vec[j];
}

void
add_mat_entrywise(double **X, double **Y, uint32_t length)
{
    /* Add two square matrices entrywise */
    uint32_t i, j;
    for (i = 0; i < length; ++i) 
        for (j = 0; j < length; ++j) 
            X[i][j] += Y[i][j];
}

int
check_PNG_signature(unsigned char *buffer)
{
    /* make sure the first PNG_BYTES_TO_CHECK bytes are the same as the byte
     * signature specified in the PNG specification.
     */
    unsigned i;
    for (i = 0; i < PNG_BYTES_TO_CHECK; ++i) 
    {
        if (buffer[i] != signature[i]) 
        {
            fprintf(stderr, "** File sig does not match PNG, received ");
            for (i = 0; i < PNG_BYTES_TO_CHECK; ++i)
                fprintf(stderr, "%.2X ", buffer[i]);
            fprintf(stderr, "\n");
            abort();
    }   }

    return 1;
}

double
determ(double **mat, uint32_t N)
{
    int i, j, k, t;
    double det = 0;
    double **m;

    switch (N)
    {
        case 0:
            fprintf(stderr, "Cannot compute det of matrix\n");
            abort();

        case 1:
            det = mat[0][0];
            break;

        case 2:
            det = mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
            break;
    
        default:
            det = 0.0;  // we'll have to sum it all up
            for (k = 0; k < N; ++k)
            {
                // reserve enough space for cofactor matrix
                m = malloc((N - 1) * sizeof(double *));
                for (i = 0; i < N - 1; ++i)
                    m[i] = malloc((N - 1) * sizeof(double));
                
                for (i = 1; i < N; ++i)
                {
                    t = 0.0;
                    for (j = 0; j < N; ++j)
                    {
                        if (j == k)
                            continue;

                        m[i - 1][t] = mat[i][j];
                        ++t;
                    }
                }
                det += pow(-1.0, 2.0 + k) * mat[0][k] * determ(m, N - 1);
                for (i = 0; i < N - 1; ++i) free(m[i]);
                free(m);
            }
    }

    return det;
}

double **
inverse3x3(double **matrix, double det)
{
    /* Hard-coded 3x3 matrix inverse, not a general solution. */
    double **inv = fabricate2DArray(CHANNELS, CHANNELS);
    double **tmp = fabricate2DArray(2, 2);
    uint32_t i, j;

    // Hideous, but it works for 3x3 matrices :P
    tmp[0][0] = matrix[1][1];  tmp[0][1] = matrix[1][2];
    tmp[1][0] = matrix[2][1];  tmp[1][1] = matrix[2][2];
    inv[0][0] = determ(tmp, 2);

    tmp[0][0] = matrix[0][2];  tmp[0][1] = matrix[0][1];
    tmp[1][0] = matrix[2][2];  tmp[1][1] = matrix[2][1];
    inv[0][1] = determ(tmp, 2);

    tmp[0][0] = matrix[0][1];  tmp[0][1] = matrix[0][2];
    tmp[1][0] = matrix[1][1];  tmp[1][1] = matrix[1][2];
    inv[0][2] = determ(tmp, 2);

    tmp[0][0] = matrix[1][2];  tmp[0][1] = matrix[1][0];
    tmp[1][0] = matrix[2][2];  tmp[1][1] = matrix[2][0];
    inv[1][0] = determ(tmp, 2);

    tmp[0][0] = matrix[0][0];  tmp[0][1] = matrix[0][2];
    tmp[1][0] = matrix[2][0];  tmp[1][1] = matrix[2][2];
    inv[1][1] = determ(tmp, 2);

    tmp[0][0] = matrix[0][2];  tmp[0][1] = matrix[0][0];
    tmp[1][0] = matrix[1][2];  tmp[1][1] = matrix[1][0];
    inv[1][2] = determ(tmp, 2);

    tmp[0][0] = matrix[1][0];  tmp[0][1] = matrix[1][1];
    tmp[1][0] = matrix[2][0];  tmp[1][1] = matrix[2][1];
    inv[2][0] = determ(tmp, 2);

    tmp[0][0] = matrix[0][1];  tmp[0][1] = matrix[0][0];
    tmp[1][0] = matrix[2][1];  tmp[1][1] = matrix[2][0];
    inv[2][1] = determ(tmp, 2);

    tmp[0][0] = matrix[0][0];  tmp[0][1] = matrix[0][1];
    tmp[1][0] = matrix[1][0];  tmp[1][1] = matrix[1][1];
    inv[2][2] = determ(tmp, 2);

    destroy2DArray(2, tmp);

    for (i = 0; i < CHANNELS; ++i) 
        for (j = 0; j < CHANNELS; ++j) 
            inv[i][j] /= det;

    return inv;
}

varTuple
var_covar(uint32_t samples, rTuple *data)
{
    /* Compute the variance covariance matrix. Works well currently for 
     * samples << pixel count.
     */
    double pixels = (double)data->width * (double)data->height;
    if ((double)samples > pixels) {
        fprintf(stderr, "** Samples must be lower than number of pixels\n");
        abort();
    }
    
    double *avV = malloc(CHANNELS * sizeof(double));
    memset(avV, 0, CHANNELS * sizeof(double));
    
    struct coord *picks = malloc(samples * sizeof(struct coord));
    uint32_t i, j;

    // build an array of random coordinates in the image; I've tested this for
    // pretty large samples and it has a distribution that's pretty level
    for (i = 0; i < samples; ++i) 
    {
        picks[i].left  = rand() % data->width;
        picks[i].right = rand() % data->height;
    }

    png_bytep row, px;
    
    // compute the average vector estimate using the array of random coords    
    for (i = 0; i < samples; ++i)
    {
        row = data->datap[picks[i].right];
        px  = &(row[picks[i].left * sizeof(int)]);
        for (j = 0; j < CHANNELS; ++j) avV[j] += (double)px[j];
    };

    for (j = 0; j < CHANNELS; ++j) avV[j] /= (double)samples;

    // forge variance-covariance matrix 
    double **outer_res = fabricate2DArray(CHANNELS, CHANNELS);
    double *diff = malloc(CHANNELS * sizeof(double));
    double **var_covar_mat = fabricate2DArray(CHANNELS, CHANNELS);
    memset(diff, 0, CHANNELS * sizeof(double));

    for (i = 0; i < samples; ++i)
    {
        // compute the difference between a pixel and the average vector
        row = data->datap[picks[i].right];
        px  = &(row[picks[i].left * sizeof(int)]);
        for (j = 0; j < CHANNELS; ++j) diff[j] = (double)px[j] - avV[j];
    
        // put the outer product in outer_res
        outer_product(outer_res, diff);

        // add this outer product to the var_covar_mat entrywise
        add_mat_entrywise(var_covar_mat, outer_res, CHANNELS);
    }

    // normalize our variance-covariance matrix
    for (i = 0; i < CHANNELS; ++i)
        for (j = 0; j < CHANNELS; ++j)
            var_covar_mat[i][j] /= (double)samples;

    // free up memory, all we return is the variance-covariance matrix
    free(picks);
    free(diff);
    destroy2DArray(CHANNELS, outer_res);

    varTuple r = { avV, var_covar_mat };

    return r;
}


//double *
//left_multiply_mat(double *left, double **mat)
//{
//    /* Perform the left matrix multiplication. */
//    uint32_t i, j;
//    double sum;
//    double *res = malloc(CHANNELS * sizeof(double));
//
//    for (i = 0; i < CHANNELS; ++i)
//    {
//        sum = 0.0;
//        for (j = 0; j < CHANNELS; ++j) sum += left[j] * mat[j][i];
//        res[i] = sum;
//    }
//
//    return res;
//}
//
//double
//right_multiply_mat(double *left, double *right)
//{
//    /* Perform the right matrix multiplication (really just a dot prod) */
//    uint32_t i;
//    double sum = 0.0;
//
//    for (i = 0; i < CHANNELS; ++i) sum += left[i] * right[i];
//
//    return sum;
//}

double
fast_mahal_metric(double *px, double **inv_covar, uint32_t N)
{
    /* Courtesy of Tim Van Slyke---more cache friendly */
    uint32_t i, j;
    double dist_sqrd = 0.0;
    
    for (i = 0; i < N; ++i)
        for (j = 0; j < N; ++j)
            dist_sqrd += px[j] * px[i] * inv_covar[i][j];

    return sqrt(dist_sqrd);
}

int **
mahal_metric(rTuple *data, varTuple *var)
{
    /* Now compute the mahalanobis distance. This is the real memory/time eater 
     * in the process.
     */

    int **dists = fabricate2DintArray(data->height, data->width);
    double *remapped_px = malloc(CHANNELS * sizeof(double));
    //double *res;
    uint32_t i, j, k;
    png_bytep row, px;
    //double sum;

    for (i = 0; i < data->height; ++i) 
    {
        row = data->datap[i];
        for (j = 0; j < data->width; ++j)
        {
            px  = &(row[j * sizeof(int)]);

            // Remap the pixel values to a double && compute difference between
            // pixel and average vector
            for (k = 0; k < CHANNELS; ++k) 
                remapped_px[k] = (double)px[k] - var->avV[k];

            // now two matrix multiplication
            //res = left_multiply_mat(remapped_px, var->var_covar_mat);
            //sum = right_multiply_mat(remapped_px, res);
            
            dists[i][j] = (int)(fast_mahal_metric(remapped_px, 
                var->var_covar_mat, CHANNELS));
        }
    }

    free(remapped_px);
    //free(res);

    return dists;
}

rTuple 
read_png_file(char *file_name) {
    /* Get PNG data - I've pieced this together by reading `example.c` from
       beginning to end */
    printf("** Reading data from %s\n", file_name);

    png_uint_32 width, height;

    uint32_t row;
    int bit_depth, color_type, interlace_type;
    unsigned char buff[PNG_BYTES_TO_CHECK];

    FILE *fp = fopen(file_name, "rb");
    if (fp == NULL) abort();

    if (fread(buff, 1, PNG_BYTES_TO_CHECK, fp) != PNG_BYTES_TO_CHECK) {
        fprintf(stderr, "** Could not read %d bytes\n", PNG_BYTES_TO_CHECK);
        abort();
    }

    check_PNG_signature(buff);
    rewind(fp);

    // create and initialize the png_struct, which will be destroyed later
    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING
        , NULL  /* Following 3 mean use stderr & longjump method */
        , NULL
        , NULL
    );
    if (!png_ptr) abort();

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) abort();

    // following I/O initialization method is required
    png_init_io(png_ptr, fp);
    png_set_sig_bytes(png_ptr, 0);  // libpng has this built in too

    // call to png_read_info() gives us all of the information from the
    // PNG file before the first IDAT (image data chunk)
    png_read_info(png_ptr, info_ptr);

    // Get header metadata now
    png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type, 
        &interlace_type, NULL, NULL);

    // Scale 16-bit images to 8-bits as accurately as possible (shouldn't be an
    // issue though, since we're working with RGB data)
#ifdef PNG_READ_SCALE_16_TO_8_SUPPORTED
    png_set_scale_16(png_ptr);
#else
    png_set_strip_16(png_ptr);
#endif

    png_set_packing(png_ptr);

    // PNGs we're working with should have a color_type RGB
    if (color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_palette_to_rgb(png_ptr);

    if (color_type == PNG_COLOR_TYPE_RGB 
        || color_type == PNG_COLOR_TYPE_GRAY 
        || color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_filler(png_ptr, 0xff, PNG_FILLER_AFTER);

    png_read_update_info(png_ptr, info_ptr);

    png_bytep *row_pointers = malloc(sizeof(png_bytep) * height);

    for (row = 0; row < height; ++row)
        row_pointers[row] = malloc(png_get_rowbytes(png_ptr, info_ptr));

    png_read_image(png_ptr, row_pointers);
    png_read_end(png_ptr, info_ptr);

    // Now clean up - the image data is in memory
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    fclose(fp);

    rTuple t = { width, height, bit_depth, row_pointers };

    return t;
}

void 
write_png_file(char *file_name, rTuple *data)
{
    /* Write data to file - this is slightly harder than reading it */
    printf("** Writing image %s\n", file_name);

    FILE *fp = fopen(file_name, "wb");
    if (fp == NULL)
        abort(); // probably would be better to add a custom error message 
                 // struct & redirect to /var/log/capture.log
    
    // Recreate the info and data pointer
    png_structp png_ptr;
    png_infop info_ptr;

    // As before we created a read struct, we now create a _write_ struct with
    // stderr and longjump methods
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) {
        fclose(fp);     // close this before aborting
        abort();
    }

    // Do the same thing for the info pointer
    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        fclose(fp);
        png_destroy_write_struct(&png_ptr,  NULL);  // get rid of former struct
        abort();
    }

    // Set error handling
    if (setjmp(png_jmpbuf(png_ptr)))
    {
        fclose(fp);
        png_destroy_write_struct(&png_ptr, &info_ptr);
        abort();
    }

    png_init_io(png_ptr, fp);

    // Now set info in the ihdr
    png_set_IHDR(png_ptr, info_ptr, data->width, data->height, data->bit_depth
        , PNG_COLOR_TYPE_GRAY          /* 1 channel of data */
        , PNG_INTERLACE_NONE           /* don't worry about it */
        , PNG_COMPRESSION_TYPE_BASE    /* default */
        , PNG_FILTER_TYPE_BASE
    );

    // Write the image data to fp
    png_set_rows(png_ptr, info_ptr, data->datap);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    png_write_end(png_ptr, NULL);

    /* Write additional comments &| descr. here (optional) */

    png_write_info(png_ptr, info_ptr); 
     
    png_destroy_write_struct(&png_ptr, &info_ptr);  // Required
    fclose(fp);
}

void
overwrite(rTuple *data, int **m)
{
    /* Overwrite the values in *data */
    uint32_t i, j;
    png_bytep row;

    // build a new framework/array to hold data
    png_bytep *row_pointers = malloc(sizeof(png_bytep) * data->height);

    // iterate over whole image and map to a 2d png_byte array
    for (i = 0; i < data->height; ++i)
    {
        row_pointers[i] = malloc(sizeof(png_byte) * data->width);
        row = row_pointers[i];
        for (j = 0; j < data->width; ++j)
        {
            row[j * sizeof(png_byte)] = m[i][j];
        }
    }

//    destroy2DPNGbyteArray(data->width, data->datap);
    free(data->datap);
    data->datap = row_pointers;
}

int 
main(int argc, char *argv[])
{
    if (argc < 2 || 4 < argc) {
        fprintf(stderr, "** Provide filename\n");
        abort();
    }

    clock_t start, end;
    char *fileName  = argv[1];
    char *writeName = argv[2];
    srand(time(NULL));
    
    start = clock();

    // get data read
    rTuple data = read_png_file(fileName);

    // Compute variance-covariance matrix
    varTuple variance_covariance = var_covar(SAMPLES, &data);

    // get the inverse of the variance-covariance
    double det = determ(variance_covariance.var_covar_mat, CHANNELS);
    if (abs(det) < 1e-4) {  // may or may not be a good bound
        fprintf(stderr, "Variance-covariance likely degenerate\n");
        abort();
    }

    double **inverse = inverse3x3(variance_covariance.var_covar_mat, det);
    destroy2DArray(CHANNELS, variance_covariance.var_covar_mat);
    variance_covariance.var_covar_mat = inverse;

    // compute mahalanobis distances now
    int **m = mahal_metric(&data, &variance_covariance);

    // write new data now
    overwrite(&data, m);
    write_png_file(writeName, &data);

    // free everything
    destroy2DPNGbyteArray(data.height, data.datap);
    free(variance_covariance.avV);
    destroy2DArray(CHANNELS, inverse);
    destroy2DintArray(1, m);

    end = clock();
    printf("** Time: %lf\n", (double)(end - start) / CLOCKS_PER_SEC);

    return 1;
}
