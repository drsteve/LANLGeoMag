#include <png.h>
#include <stdlib.h>
#include <GL/gl.h>
#include <GL/glu.h>

#define ERROR   -1
#define NOT_PNG -2


/*
 * The number of bytes to read off PNG image.  These first few bytes will be
 * used to verify that the image really is a PNG image. (We can use fewer than
 * 8 bytes, byte the validation is not as acurate. Dont use more than 8 bytes.)
 */
#define PNG_HEADER_BYTES    8 

int ReadPng( char *Filename, int *Width, int *Height, GLubyte **qImage ) {


    unsigned char   Header[8]; // 8-byte PNG Header. Used to verify image authenticity.
    png_uint_32     width, height;
    int             bit_depth, color_type, interlace_type;
    int             compression_type, filter_method;
    int             IsPng, BytesPerPixel, i, j, k, k0, k1;
    png_bytep       *row_pointers = NULL;
    GLubyte         *p, *pImage;
    FILE            *fp;
   

    if ( Filename == NULL ) {
        printf( "Filename is NULL\n" );
        return( ERROR );
    }


    /*
     *  Open the PNG file for reading. Open in binary format.  Bail out if open
     *  fails.
     */
    if ( (fp = fopen(Filename, "rb")) == NULL ) {
        printf( "Could not open file %s for reading\n", Filename );
        return( ERROR );
    }


    /*
     *  Verify that the image is truely a PNG image. Bail out if its not.
     */
    fread( Header, 1, PNG_HEADER_BYTES, fp );
    IsPng = !png_sig_cmp( Header, 0, PNG_HEADER_BYTES );
    if ( !IsPng ) {
        printf( "File (%s) does not contain a a PNG image\n", Filename );
        return( NOT_PNG );
    }




    /*
     *  Allocate png_ structs
     */
    png_structp png_ptr = png_create_read_struct( PNG_LIBPNG_VER_STRING, (png_voidp)NULL, NULL, NULL );
    if ( !png_ptr ) return (ERROR);

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
        return (ERROR);
    }


    /* 
     *  Setup the the setjmp stuff needed for returning properly when a libpng
     *  error occurs.
     */
    if ( setjmp( png_jmpbuf(png_ptr) ) ) {
        png_destroy_read_struct( &png_ptr, &info_ptr, NULL );
        if ( row_pointers ) free( row_pointers );
        return( ERROR );
    }



/* set "png_read" callback function and give source of data */
//png_set_read_fn (png_ptr, (png_voidp *)file, png_read_from_mem);

    /* read png info */
    png_init_io(png_ptr, fp);
    png_set_sig_bytes(png_ptr, PNG_HEADER_BYTES);
    png_read_info( png_ptr, info_ptr );














//    png_init_io(png_ptr, fp);
//    png_set_sig_bytes(png_ptr, PNG_HEADER_BYTES);
//    png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);



    /*
     *  Retrive info about image.
     */
    png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type, &interlace_type, &compression_type, &filter_method);


    /*
     *  Apply transformations
     */
    if (color_type == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png_ptr);
//    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) png_set_gray_1_2_4_to_8(png_ptr);
    if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS)) png_set_tRNS_to_alpha(png_ptr);
    if (bit_depth == 16) png_set_strip_16(png_ptr);
//    if (color_type & PNG_COLOR_MASK_ALPHA) png_set_strip_alpha(png_ptr);
    if (bit_depth < 8) png_set_packing(png_ptr);

    png_read_update_info( png_ptr, info_ptr );


    /*
     *  re-Retrive info about image.
     */
    png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type, &interlace_type, &compression_type, &filter_method);



    /*
     *  Determine mhow many bytes per pixel there are. (actually
     *  its how many channels/pixel there are. Is it always bytes?)
     */
    switch( color_type ) {

        case PNG_COLOR_TYPE_GRAY:
            printf("Image is a PNG_COLOR_TYPE_GRAY\n"); 
            BytesPerPixel = 1;
            break;

        case PNG_COLOR_TYPE_GRAY_ALPHA:
            printf("Image is a PNG_COLOR_TYPE_GRAY_ALPHA\n"); 
            BytesPerPixel = 2;
            break;

        case PNG_COLOR_TYPE_RGB:
            printf("Image is a PNG_COLOR_TYPE_RGB\n"); 
            BytesPerPixel = 3;
            break;

        case PNG_COLOR_TYPE_RGB_ALPHA:
            printf("Image is a PNG_COLOR_TYPE_RGB_ALPHA\n"); 
            BytesPerPixel = 4;
            break;

        default:
            printf("Unknown Color Type in PNG file\n");
            BytesPerPixel = 1; // ?
            break;

    }

                


    /*
     *  Allocate memory for OpenGL texture. We assume that we want 4 bytes/pixel here.
     */
    pImage = (GLubyte *)calloc(  width*height*4, sizeof(GLubyte) );
    printf("pImage = %p\n", pImage);


    /* 
     *  Allocate row_pointers
     */
    row_pointers = png_malloc( png_ptr, height*sizeof(png_bytep) );
    for ( i=0; i<height; i++ ) {
        row_pointers[i] = png_malloc( png_ptr, width*BytesPerPixel );
    }
    printf("width*BytesPerPixel = %d\n", (int)(width*BytesPerPixel));
    png_set_rows( png_ptr, info_ptr, row_pointers );


    /*
     *  Read the image in
     */
    png_read_image(png_ptr, row_pointers);


    /* 
     *  Finish the read and clean up
     */
    png_read_end( png_ptr, NULL );
    png_destroy_read_struct( &png_ptr, &info_ptr, NULL );
    fclose(fp);
    





    /*
     *   We want to stuff the pixel data into a GL texture image.  The texture
     *   image is basically just a single array comprised of RGBA data values.
     *   And we also need to know the width and height. Since not all PNG
     *   images are RGBA, we need to do different things for different color
     *   types...
     */
    p = &pImage[0];
    switch( color_type ) {

        case PNG_COLOR_TYPE_GRAY:
            for ( j=0; j<height; j++ ){
                for ( i=0; i<width; i++ ){
                    *p++ = row_pointers[height-1-j][i];
                    *p++ = row_pointers[height-1-j][i];
                    *p++ = row_pointers[height-1-j][i];
                    *p++ = 255;
                }
            }
            break;

        case PNG_COLOR_TYPE_GRAY_ALPHA:
            for ( j=0; j<height; j++ ){
                for ( i=0; i<width; i++ ){
                    k0 = i*BytesPerPixel;      // BytesPerPixel=2
                    k1 = i*BytesPerPixel+1;
                    *p++ = row_pointers[height-1-j][k0]; // Gray
                    *p++ = row_pointers[height-1-j][k0]; // Gray
                    *p++ = row_pointers[height-1-j][k0]; // Gray
                    *p++ = row_pointers[height-1-j][k1]; // Alpha
                }
            }
            break;

        case PNG_COLOR_TYPE_RGB:
            for ( j=0; j<height; j++ ){
                for ( i=0; i<width; i++ ){
                    k0 = i*BytesPerPixel;                                           // BytesPerPixel=3
                    for (k=k0; k < k0+BytesPerPixel; k++) *(p++) = row_pointers[height-1-j][k]; // Copy over the RGB vals
                    *(p++) = 255;                                                        // Add in a made-up A val
                }
            }
            break;

        case PNG_COLOR_TYPE_RGB_ALPHA:
            for ( j=0; j<height; j++ ){
                for ( i=0; i<width; i ++ ){
                    k0 = i*BytesPerPixel;                                           // BytesPerPixel=4
                    for (k=k0; k < k0+BytesPerPixel; k++) *p++ = row_pointers[height-1-j][k]; // Copy over the RGBA vals
                }
            }
            break;

        default:
                break;

    }
    



    *Width  = (int)width;
    *Height = (int)height;



    for ( i=0; i<height; i++ ) {
        free(row_pointers[i]);
    }
    free(row_pointers);



    *qImage = pImage;

    return( 1 );


}
