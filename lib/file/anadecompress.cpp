#include "redux/file/anadecompress.hpp"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>


int redux::file::anadecrunch32(const unsigned char* x, int32_t* array, int slice, int nx, int ny, int little_endian)
/* decompress a bit stream in x; result is n I*4 elements, put in array;
   bit slice size slice */
{
    int iq;
    int r0, r1, r2, nb;
    int j, in, i, k, ix, iy, mask;
    long long y64;
    unsigned char xq;
    union {
        int i;
        short w;
        unsigned char b[4];
    } y;
    union {
        long long l64;
        unsigned char b[8];
    } yy;

    /* begin execution */
    mask = (1 << slice) - 1;
    /* determine the # of bytes to transfer to 32 bit int for fixed portion */
    nb = (slice + 14) / 8;   /* range 1 to 5 */
    if(slice == 0)
        nb = 0;     /* but slice = 0 a special case */
    y.i = 0;
    i = 0;
    r0 = 0;
    r1 = 0;
    in = 0;
    for(iy = 0; iy < ny; iy++) {    /* start of iy (outer) loop */
        /* get the first value, 4 bytes */
        if(little_endian) {
            y.b[0] = x[i];
            y.b[1] = x[i + 1];
            y.b[2] = x[i + 2];
            y.b[3] = x[i + 3];
            iq = array[in++] = y.i;
        }
        else {
            y.b[0] = x[i + 3];
            y.b[1] = x[i + 2];
            y.b[2] = x[i + 1];
            y.b[3] = x[i];
            iq = array[in++] = y.i;
        }
        r1 = r1 + 32;
        r2 = 0;
        for(ix = 1; ix < nx; ix++) {
            /* first the fixed slice portion */
            i = r1 / 8;
            j = r1 % 8;
            if(little_endian) {
                yy.b[0] = x[i];
                if(nb > 1)
                    yy.b[1] = x[i + 1];
                if(nb > 2)
                    yy.b[2] = x[i + 2];
                if(nb > 3)
                    yy.b[3] = x[i + 3];
                if(nb > 4)
                    yy.b[4] = x[i + 4];
            }
            else {
                yy.b[7] = x[i];
                if(nb > 1)
                    yy.b[6] = x[i + 1];
                if(nb > 2)
                    yy.b[5] = x[i + 2];
                if(nb > 3)
                    yy.b[4] = x[i + 3];
                if(nb > 4)
                    yy.b[3] = x[i + 4];
            }
            /* shift and mask out the bit slice */
            r2 = (int)((yy.l64 >> j) & mask);
            /* the variable bit portion, find the first set bit */
            r1 = r1 + slice;   /* bump r1 pass the fixed part */
            i = r1 / 8;
            j = r1 % 8;
            if((xq = (x[i] >> j)) != 0) {
                /* caught it on first byte, find the bit */
                if((xq & 1) != 0)
                    r0 = 1;
                else if((xq & 2) != 0)
                    r0 = 2;
                else if((xq & 4) != 0)
                    r0 = 3;
                else if((xq & 8) != 0)
                    r0 = 4;
                else if((xq & 16) != 0)
                    r0 = 5;
                else if((xq & 32) != 0)
                    r0 = 6;
                else if((xq & 64) != 0)
                    r0 = 7;
                else if((xq & 128) != 0)
                    r0 = 8;
            }
            else {
                /* not in first byte (or part of one) checked, carry on, first count bits in
                   that first byte */
                r0 = 8 - j;
                /* check up to 4 more bytes, if not found than an error */
                for(k = i + 1; k < i + 5; k++) {
                    if((xq = x[k]) != 0) {
                        /* caught it here, find the bit and then jump from loop */
                        if((xq & 1) != 0)
                            r0 += 1;
                        else if((xq & 2) != 0)
                            r0 += 2;
                        else if((xq & 4) != 0)
                            r0 += 3;
                        else if((xq & 8) != 0)
                            r0 += 4;
                        else if((xq & 16) != 0)
                            r0 += 5;
                        else if((xq & 32) != 0)
                            r0 += 6;
                        else if((xq & 64) != 0)
                            r0 += 7;
                        else if((xq & 128) != 0)
                            r0 += 8;
                        break;
                    }
                    else {
                        r0 = r0 + 8;
                        /* add 8 bits for each all zero byte */
                        if(r0 > 32) {
                            fprintf(stderr, "DECRUNCH -- bad bit sequence, cannot continue\n");
                            fprintf(stderr, "i = %d, r1 = %d, ix= %d, iy = %d\n", i, r1, ix, iy);
                            return -1;
                        }
                    }
                }
            }
            r1 = r1 + r0;   /* update pointer */
            /* r0 even or odd determines sign of difference */
            if((r0 & 1) != 0) {
                /* positive case */
                r0 = (r0 / 2) << slice;
                iq = iq + r2;
                iq = iq + r0;
                array[in] = iq;
            }
            else {
                if(r0 == 32) {
                    /* a long one, yank out the next 33 bits and use as difference */
                    i = r1 / 8;
                    j = r1 % 8;
                    if(little_endian) {
                        yy.b[0] = x[i];
                        yy.b[1] = x[i + 1];
                        yy.b[2] = x[i + 2];
                        yy.b[3] = x[i + 3];
                        yy.b[4] = x[i + 4];
                    }
                    else {
                        yy.b[7] = x[i];
                        yy.b[6] = x[i + 1];
                        yy.b[5] = x[i + 2];
                        yy.b[4] = x[i + 3];
                        yy.b[3] = x[i + 4];
                    }
                    /* shift and mask out the 33 bit slice */
                    y64 = (yy.l64 >> j) & 0x1ffffffffLL;
                    r1 = r1 + 33;
                    /* if the top bit was set, do a sign extend, note that 64 bit arithmetic used */
                    if((y64 & 0x100000000LL) != 0)
                        y64 = y64 | 0xffffffff00000000LL;
                    y64 = y64 + (long long)array[in - 1];
                    iq = array[in] = (long)y64;
                }
                else {
                    /* minus case (normal) */
                    r0 = (-r0 / 2) << slice;
                    iq = iq + r2;
                    iq = iq + r0;
                    array[in] = iq;
                }
            }
            in = in + 1;
        }   /* end of ix loop */
        i = (r1 + 7) / 8;
        r1 = 8 * i;
    }   /* end of iy loop */
    return 1;
}   /* end of routine */

/*--------------------------------------------------------------------------*/
int redux::file::anadecrunch(const unsigned char* x, int16_t* array, int slice, int nx, int ny, int little_endian)
/* decompress a bit stream in x; result is n I*2 elements, put in array;
   bit slice size slice */
{
    short iq;
    int r0, r1, r2, r4, nb, mask;
    int j, in, i, k, ix, iy;
    unsigned char xq;
    union {
        int i;
        short w;
        unsigned char b[4];
    } y;

    /* begin execution */
    mask = (1 << slice) - 1;
    /* determine the # of bytes to transfer to 32 bit int for fixed portion */
    if(slice == 0)
        nb = 0;
    else {
        if(slice < 2)
            nb = 1;
        else {
            if(slice < 10)
                nb = 2;
            else
                nb = 3;
        }
    };
    y.i = 0;
    i = 0;
    r0 = 0;
    r1 = 0;
    in = 0;
    for(iy = 0; iy < ny; iy++) {    /* start of iy (outer) loop */
        /* get the first value */
        if(little_endian) {
            y.b[0] = x[i];
            y.b[1] = x[i + 1];
            iq = y.w;
            array[in++] = iq;
        }
        else {
            y.b[0] = x[i + 1];
            y.b[1] = x[i];
            iq = y.w;
            array[in++] = iq;
        }
        r1 = r1 + 16;
        r2 = 0;
        for(ix = 1; ix < nx; ix++) {
            /* first the fixed slice portion */
            i = r1 / 8;
            j = r1 % 8;
            if(little_endian) {
                y.b[0] = x[i];
            }
            else {
                y.b[3] = x[i];
            }
            /* test effect on timing */
            if(little_endian) {
                if(nb > 1) {
                    y.b[1] = x[i + 1];
                    if(nb > 2)
                        y.b[2] = x[i + 2];
                }
            }
            else {
                if(nb > 1) {
                    y.b[2] = x[i + 1];
                    if(nb > 2)
                        y.b[1] = x[i + 2];
                }
            }
            /* shift and mask out the bit slice */
            r2 = (y.i >> j) & mask;
            /* the variable bit portion, find the first set bit */
            r1 = r1 + slice;   /* bump r1 pass the fixed part */
            i = r1 / 8;
            j = r1 % 8;
            if((xq = x[i] >> j) != 0) {
                /* caught it on first byte, find the bit */
                if((xq & 1) != 0)
                    r0 = 1;
                else if((xq & 2) != 0)
                    r0 = 2;
                else if((xq & 4) != 0)
                    r0 = 3;
                else if((xq & 8) != 0)
                    r0 = 4;
                else if((xq & 16) != 0)
                    r0 = 5;
                else if((xq & 32) != 0)
                    r0 = 6;
                else if((xq & 64) != 0)
                    r0 = 7;
                else if((xq & 128) != 0)
                    r0 = 8;
            }
            else {
                /* not in first byte (or part of one) checked, carry on, first count bits in
                   that first byte */
                r0 = 8 - j;
                /* check up to 4 more bytes, if not found than an error */
                for(k = i + 1; k < i + 5; k++) {
                    if((xq = x[k]) != 0) {
                        /* caught it here, find the bit and then jump from loop */
                        if((xq & 1) != 0)
                            r0 += 1;
                        else if((xq & 2) != 0)
                            r0 += 2;
                        else if((xq & 4) != 0)
                            r0 += 3;
                        else if((xq & 8) != 0)
                            r0 += 4;
                        else if((xq & 16) != 0)
                            r0 += 5;
                        else if((xq & 32) != 0)
                            r0 += 6;
                        else if((xq & 64) != 0)
                            r0 += 7;
                        else if((xq & 128) != 0)
                            r0 += 8;
                        break;
                    }
                    else {
                        r0 = r0 + 8;
                        /* add 8 bits for each all zero byte */
                        if(r0 > 32) {
                            fprintf(stderr, "DECRUNCH -- bad bit sequence, cannot continue\n");
                            fprintf(stderr, "i = %d, r1 = %d, ix= %d, iy = %d\n", i, r1, ix, iy);
                            return -1;
                        }
                    }
                }
            }
            r1 = r1 + r0;   /* update pointer */
            /* r0 even or odd determines sign of difference */
            if((r0 & 1) != 0) {
                /* positive case */
                r0 = (r0 / 2) << slice;
                iq = iq + r2;
                iq = iq + r0;
                array[in] = iq;
            }
            else {
                if(r0 == 32) {
                    /* a long one, yank out the next 17 bits and use as difference */
                    i = r1 / 8;
                    j = r1 % 8;
                    if(little_endian) {
                        y.b[0] = x[i];
                        y.b[1] = x[i + 1];
                        y.b[2] = x[i + 2];
                    }
                    else {
                        y.b[3] = x[i];
                        y.b[2] = x[i + 1];
                        y.b[1] = x[i + 2];
                    }
                    /* shift and mask out the 17 bit slice */
                    r2 = (y.i >> j) & 0x1ffff;
                    r1 = r1 + 17;
                    /* if the top bit was set, do a sign extend, note that 32 bit arithmetic used */
                    if((r2 & 0x10000) != 0)
                        r2 = r2 | 0xffff0000;
                    r4 = array[in - 1];
                    r4 = r4 + r2;
                    array[in] = r4;
                    iq = r4;
                }
                else {
                    /* minus case (normal) */
                    r0 = (-r0 / 2) << slice;
                    iq = iq + r2;
                    iq = iq + r0;
                    array[in] = iq;
                }
            }
            in = in + 1;
        }   /* end of ix loop */
        i = (r1 + 7) / 8;
        r1 = 8 * i;
    }   /* end of iy loop */
    return 1;
}   /* end of routine */

/*--------------------------------------------------------------------------*/
int redux::file::anadecrunch8(const unsigned char* x, int8_t* array, int slice, int nx, int ny, int little_endian)
/* decompress a bit stream in x; result is n I*1 elements, put in array;
   bit slice size slice */
{
    uint8_t iq;
    int r0, r1, r2, r4, nb, mask;
    int j, in, i, k, ix, iy;
    uint8_t xq;
    union {
        int i;
        short w;
        unsigned char b[4];
    } y;

    /* begin execution */
    mask = (1 << slice) - 1;
    /* determine the # of bytes to transfer to 32 bit int for fixed portion */
    if(slice == 0)
        nb = 0;
    else if(slice < 2)
        nb = 1;
    else if(slice < 10)
        nb = 2;
    else
        nb = 3;
    y.i = 0;
    i = 0;
    r0 = 0;
    r1 = 0;
    in = 0;
    for(iy = 0; iy < ny; iy++) {    /* start of iy (outer) loop */
        /* get the first value */
        iq = x[i];
        array[in++] = iq;
        r1 = r1 + 8;
        r2 = 0;
        for(ix = 1; ix < nx; ix++) {
            /* first the fixed slice portion */
            i = r1 / 8;
            j = r1 % 8;
            if(little_endian) {
                y.b[0] = x[i];
                if(nb > 1)
                    y.b[1] = x[i + 1];
                if(nb > 2)
                    y.b[2] = x[i + 2];
            }
            else {
                y.b[3] = x[i];
                if(nb > 1)
                    y.b[2] = x[i + 1];
                if(nb > 2)
                    y.b[1] = x[i + 2];
            }
            /* shift and mask out the bit slice */
            r2 = (y.i >> j) & mask;
            /* the variable bit portion, find the first set bit */
            r1 = r1 + slice;   /* bump r1 pass the fixed part */
            i = r1 / 8;
            j = r1 % 8;
            if((xq = x[i] >> j) != 0) {
                /* caught it on first byte, find the bit */
                if((xq & 1) != 0)
                    r0 = 1;
                else if((xq & 2) != 0)
                    r0 = 2;
                else if((xq & 4) != 0)
                    r0 = 3;
                else if((xq & 8) != 0)
                    r0 = 4;
                else if((xq & 16) != 0)
                    r0 = 5;
                else if((xq & 32) != 0)
                    r0 = 6;
                else if((xq & 64) != 0)
                    r0 = 7;
                else if((xq & 128) != 0)
                    r0 = 8;
            }
            else {
                /* not in first byte (or part of one) checked, carry on, first count bits in
                   that first byte */
                r0 = 8 - j;
                /* check up to 4 more bytes, if not found than an error */
                for(k = i + 1; k < i + 5; k++) {
                    if((xq = x[k]) != 0) {
                        /* caught it here, find the bit and then jump from loop */
                        if((xq & 1) != 0)
                            r0 += 1;
                        else if((xq & 2) != 0)
                            r0 += 2;
                        else if((xq & 4) != 0)
                            r0 += 3;
                        else if((xq & 8) != 0)
                            r0 += 4;
                        else if((xq & 16) != 0)
                            r0 += 5;
                        else if((xq & 32) != 0)
                            r0 += 6;
                        else if((xq & 64) != 0)
                            r0 += 7;
                        else if((xq & 128) != 0)
                            r0 += 8;
                        break;
                    }
                    else {
                        r0 = r0 + 8;
                        /* add 8 bits for each all zero byte */
                        if(r0 > 32) {
                            fprintf(stderr, "DECRUNCH -- bad bit sequence, cannot continue");
                            return -1;
                        }
                    }
                }
            }
            r1 = r1 + r0;   /* update pointer */
            /* r0 even or odd determines sign of difference */
            if((r0 & 1) != 0) {
                /* positive case */
                r0 = (r0 / 2) << slice;
                iq = iq + r2;
                iq = iq + r0;
                array[in] = iq;
            }
            else {
                if(r0 == 32) {
                    /* a long one, yank out the next 9 bits and use as difference */
                    i = r1 / 8;
                    j = r1 % 8;
                    if(little_endian) {
                        y.b[0] = x[i];
                        y.b[1] = x[i + 1];
                    }
                    else {
                        y.b[3] = x[i];
                        y.b[2] = x[i + 1];
                    }
                    /* shift and mask out the 9 bit slice */
                    r2 = (y.i >> j) & 0x1ff;
                    r1 = r1 + 9;
                    /* if the top bit was set, do a sign extend, note that 32 bit arithmetic used */
                    if((r2 & 0x100) != 0)
                        r2 = r2 | 0xffffff00;
                    r4 = array[in - 1];
                    r4 = r4 + r2;
                    array[in] = r4;
                    iq = r4;
                }
                else {
                    /* minus case (normal) */
                    r0 = (-r0 / 2) << slice;
                    iq = iq + r2;
                    iq = iq + r0;
                    array[in] = iq;
                }
            }
            in = in + 1;
        }   /* end of ix loop */
        i = (r1 + 7) / 8;
        r1 = 8 * i;
    }   /* end of iy loop */
    return 1;
}   /* end of routine */

/*--------------------------------------------------------------------------*/

int redux::file::anadecrunchrun(const unsigned char* x, int16_t* array, int slice, int nx, int ny, int little_endian)
/* decompress a bit stream in x; result is n I*2 elements, put in array;
   bit slice size slice */
/* this version handles the run length encoding used in anacrunchrun */
{
    short iq;
    int r0, r1, r2, r4, nb, mask, nrun, n, nc;
    int j, in, i, k, iy;
    unsigned char xq;
    union {
        int i;
        short w;
        unsigned char b[4];
    } y;

    /* begin execution */
    mask = (1 << slice) - 1;
    /* determine the # of bytes to transfer to 32 bit int for fixed portion */
    if(slice == 0)
        nb = 0;
    else if(slice < 2)
        nb = 1;
    else if(slice < 10)
        nb = 2;
    else
        nb = 3;
    y.i = 0;
    i = 0;
    r0 = 0;
    r1 = 0;
    in = 0;
    for(iy = 0; iy < ny; iy++) {    /* start of iy (outer) loop */
        /* get the first value, assume bytes reversed */
        if(little_endian) {
            y.b[0] = x[i++];
            y.b[1] = x[i++];
            iq = y.w;
            array[in++] = iq;
        }
        else {
            y.b[1] = x[i++];
            y.b[0] = x[i++];
            iq = y.w;
            array[in++] = iq;
        }
        r1 = r1 + 16;
        r2 = 0;
        nc = nx - 1;
        while(nc > 0) {
            /* look at the next run length code */
            nrun = (int)x[i++];
            if(nrun > 127) {    /* a run of a constant difference */
                n = 255 - nrun + 2;
                nc = nc - n;
                if(little_endian) {
                    y.b[0] = x[i++];
                    y.b[1] = x[i++];
                }
                else {
                    y.b[1] = x[i++];
                    y.b[0] = x[i++];
                }
                while(n--) {
                    array[in] = array[in - 1] + y.w;
                    in++;
                }
                iq = array[in - 1];
                r1 = 8 * i;
            }
            else {      /* a literal */
                r1 = 8 * i;
                nc = nc - nrun;
                while(nrun--) {
                    /* first the fixed slice portion */
                    i = r1 / 8;
                    j = r1 % 8;
                    if(little_endian) {
                        y.b[0] = x[i];
                        /* test effect on timing */
                        if(nb > 1)
                            y.b[1] = x[i + 1];
                        if(nb > 2)
                            y.b[2] = x[i + 2];
                    }
                    else {
                        y.b[3] = x[i];
                        /* test effect on timing */
                        if(nb > 1)
                            y.b[2] = x[i + 1];
                        if(nb > 2)
                            y.b[1] = x[i + 2];
                    }
                    /* shift and mask out the bit slice */
                    r2 = (y.i >> j) & mask;
                    /* the variable bit portion, find the first set bit */
                    r1 = r1 + slice;   /* bump r1 pass the fixed part */
                    i = r1 / 8;
                    j = r1 % 8;
                    if((xq = x[i] >> j) != 0) {
                        /* caught it on first byte, find the bit */
                        if((xq & 1) != 0)
                            r0 = 1;
                        else if((xq & 2) != 0)
                            r0 = 2;
                        else if((xq & 4) != 0)
                            r0 = 3;
                        else if((xq & 8) != 0)
                            r0 = 4;
                        else if((xq & 16) != 0)
                            r0 = 5;
                        else if((xq & 32) != 0)
                            r0 = 6;
                        else if((xq & 64) != 0)
                            r0 = 7;
                        else if((xq & 128) != 0)
                            r0 = 8;
                    }
                    else {
                        /* not in first byte (or part of one) checked, carry on, first count bits in
                           that first byte */
                        r0 = 8 - j;
                        /* check up to 4 more bytes, if not found than an error */
                        for(k = i + 1; k < i + 5; k++) {
                            if((xq = x[k]) != 0) {
                                /* caught it here, find the bit and then jump from loop */
                                if((xq & 1) != 0)
                                    r0 += 1;
                                else if((xq & 2) != 0)
                                    r0 += 2;
                                else if((xq & 4) != 0)
                                    r0 += 3;
                                else if((xq & 8) != 0)
                                    r0 += 4;
                                else if((xq & 16) != 0)
                                    r0 += 5;
                                else if((xq & 32) != 0)
                                    r0 += 6;
                                else if((xq & 64) != 0)
                                    r0 += 7;
                                else if((xq & 128) != 0)
                                    r0 += 8;
                                break;
                            }
                            else {
                                r0 = r0 + 8;
                                /* add 8 bits for each all zero byte */
                                if(r0 > 32) {
                                    fprintf(stderr, "DECRUNCH -- bad bit sequence, cannot continue\n");
                                    fprintf(stderr, "i = %d, r1 = %d, iy = %d\n", i, r1, iy);
                                    return -1;
                                }
                            }
                        }
                    }
                    r1 = r1 + r0;   /* update pointer */
                    /* r0 even or odd determines sign of difference */
                    if((r0 & 1) != 0) {
                        /* positive case */
                        r0 = (r0 / 2) << slice;
                        iq = iq + r2;
                        iq = iq + r0;
                        array[in] = iq;
                    }
                    else {
                        if(r0 == 32) {
                            /* a long one, yank out the next 17 bits and use as difference */
                            i = r1 / 8;
                            j = r1 % 8;
                            if(little_endian) {
                                y.b[0] = x[i];
                                y.b[1] = x[i + 1];
                                y.b[2] = x[i + 2];
                            }
                            else {
                                y.b[3] = x[i];
                                y.b[2] = x[i + 1];
                                y.b[1] = x[i + 2];
                            }
                            /* shift and mask out the 17 bit slice */
                            r2 = (y.i >> j) & 0x1ffff;
                            r1 = r1 + 17;
                            /* if the top bit was set, do a sign extend, note that 32 bit arithmetic used */
                            if((r2 & 0x10000) != 0)
                                r2 = r2 | 0xffff0000;
                            r4 = array[in - 1];
                            r4 = r4 + r2;
                            array[in] = r4;
                            iq = r4;
                        }
                        else {
                            /* minus case (normal) */
                            r0 = (-r0 / 2) << slice;
                            iq = iq + r2;
                            iq = iq + r0;
                            array[in] = iq;
                        }
                    }
                    in = in + 1;
                }
                i = (r1 + 7) / 8;
                r1 = 8 * i;
            }
        }   /* end of ix loop */
        if(nc < 0) {
            fprintf(stderr, "bad loop in decrunchrun, nc=%d, iy=%d, in= %d\n", nc, iy, in);
            return -1;
        }

        i = (r1 + 7) / 8;
        r1 = 8 * i;
    }   /* end of iy loop */
    return 1;
}   /* end of routine */

/*--------------------------------------------------------------------------*/
int redux::file::anadecrunchrun8(const unsigned char* x, int8_t* array, int slice, int nx, int ny, int little_endian)
/* decompress a bit stream in x; result is n I*2 elements, put in array;
   bit slice size slice */
/* this version handles the run length encoding used in anacrunchrun */
{
    uint8_t iq;
    int r0, r1, r2, r4, nb, mask, nrun, n, nc;
    int j, in, i, k, iy;
    unsigned char xq;
    union {
        int i;
        short w;
        unsigned char b[4];
    } y;

    /* begin execution */
    mask = (1 << slice) - 1;
    /* determine the # of bytes to transfer to 32 bit int for fixed portion */
    if(slice == 0)
        nb = 0;
    else if(slice < 2)
        nb = 1;
    else if(slice < 10)
        nb = 2;
    else
        nb = 3;
    y.i = 0;
    i = 0;
    r0 = 0;
    r1 = 0;
    in = 0;
    for(iy = 0; iy < ny; iy++) {    /* start of iy (outer) loop */
        /* get the first value */
        iq = x[i++];
        array[in++] = iq;
        r1 = r1 + 16;
        r2 = 0;
        nc = nx - 1;
        while(nc > 0) {
            /* look at the next run length code */
            nrun = (int)x[i++];
            if(nrun > 127) {    /* a run of a constant difference */
                n = 255 - nrun + 2;
                nc = nc - n;
                if(little_endian) {
                    y.b[0] = x[i++];
                    y.b[1] = x[i++];
                }
                else {
                    y.b[1] = x[i++];
                    y.b[0] = x[i++];
                }
                while(n--) {
                    array[in] = array[in - 1] + y.w;
                    in++;
                }
                iq = array[in - 1];
                r1 = 8 * i;
            }
            else {      /* a literal */
                r1 = 8 * i;
                nc = nc - nrun;
                while(nrun--) {
                    /* first the fixed slice portion */
                    i = r1 / 8;
                    j = r1 % 8;
                    if(little_endian) {
                        y.b[0] = x[i];
                        /* test effect on timing */
                        if(nb > 1)
                            y.b[1] = x[i + 1];
                        if(nb > 2)
                            y.b[2] = x[i + 2];
                    }
                    else {
                        y.b[3] = x[i];
                        /* test effect on timing */
                        if(nb > 1)
                            y.b[2] = x[i + 1];
                        if(nb > 2)
                            y.b[1] = x[i + 2];
                    }
                    /* shift and mask out the bit slice */
                    r2 = (y.i >> j) & mask;
                    /* the variable bit portion, find the first set bit */
                    r1 = r1 + slice;   /* bump r1 pass the fixed part */
                    i = r1 / 8;
                    j = r1 % 8;
                    if((xq = x[i] >> j) != 0) {
                        /* caught it on first byte, find the bit */
                        if((xq & 1) != 0)
                            r0 = 1;
                        else if((xq & 2) != 0)
                            r0 = 2;
                        else if((xq & 4) != 0)
                            r0 = 3;
                        else if((xq & 8) != 0)
                            r0 = 4;
                        else if((xq & 16) != 0)
                            r0 = 5;
                        else if((xq & 32) != 0)
                            r0 = 6;
                        else if((xq & 64) != 0)
                            r0 = 7;
                        else if((xq & 128) != 0)
                            r0 = 8;
                    }
                    else {
                        /* not in first byte (or part of one) checked, carry on, first count bits in
                           that first byte */
                        r0 = 8 - j;
                        /* check up to 4 more bytes, if not found than an error */
                        for(k = i + 1; k < i + 5; k++) {
                            if((xq = x[k]) != 0) {
                                /* caught it here, find the bit and then jump from loop */
                                if((xq & 1) != 0)
                                    r0 += 1;
                                else if((xq & 2) != 0)
                                    r0 += 2;
                                else if((xq & 4) != 0)
                                    r0 += 3;
                                else if((xq & 8) != 0)
                                    r0 += 4;
                                else if((xq & 16) != 0)
                                    r0 += 5;
                                else if((xq & 32) != 0)
                                    r0 += 6;
                                else if((xq & 64) != 0)
                                    r0 += 7;
                                else if((xq & 128) != 0)
                                    r0 += 8;
                                break;
                            }
                            else {
                                r0 = r0 + 8;
                                /* add 8 bits for each all zero byte */
                                if(r0 > 32) {
                                    fprintf(stderr, "DECRUNCH -- bad bit sequence, cannot continue\n");
                                    fprintf(stderr, "i = %d, r1 = %d, iy = %d\n", i, r1, iy);
                                    return -1;
                                }
                            }
                        }
                    }
                    r1 = r1 + r0;   /* update pointer */
                    /* r0 even or odd determines sign of difference */
                    if((r0 & 1) != 0) {
                        /* positive case */
                        r0 = (r0 / 2) << slice;
                        iq = iq + r2;
                        iq = iq + r0;
                        array[in] = iq;
                        /* printf("r0,r2,iq = %d %d %d\n", r0,r2,iq); */
                    }
                    else {
                        if(r0 == 32) {
                            /* a long one, yank out the next 9 bits and use as difference */
                            i = r1 / 8;
                            j = r1 % 8;
                            if(little_endian) {
                                y.b[0] = x[i];
                                y.b[1] = x[i + 1];
                            }
                            else {
                                y.b[3] = x[i];
                                y.b[2] = x[i + 1];
                            }
                            /* shift and mask out the 9 bit slice */
                            r2 = (y.i >> j) & 0x1ff;
                            r1 = r1 + 9;
                            /* if the top bit was set, do a sign extend, note that 32 bit arithmetic used */
                            if((r2 & 0x100) != 0)
                                r2 = r2 | 0xffffff00;
                            r4 = array[in - 1];
                            r4 = r4 + r2;
                            array[in] = r4;
                            iq = r4;
                        }
                        else {
                            /* minus case (normal) */
                            r0 = (-r0 / 2) << slice;
                            iq = iq + r2;
                            iq = iq + r0;
                            array[in] = iq;
                        }
                    }
                    in = in + 1;
                }
                i = (r1 + 7) / 8;
                r1 = 8 * i;
            }
        }   /* end of ix loop */
        if(nc < 0) {
            fprintf(stderr, "bad loop in decrunchrun8, nc=%d, iy=%d, in= %d\n", nc, iy, in);
            return -1;
        }

        i = (r1 + 7) / 8;
        r1 = 8 * i;
    }   /* end of iy loop */
    return 1;
}   /* end of routine */

/*--------------------------------------------------------------------------*/
