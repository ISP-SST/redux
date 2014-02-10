#include "redux/file/anacompress.hpp"

#include "redux/file/fileana.hpp"
#include "redux/util/endian.hpp"

#include <cstring>
#include <cstdlib>
#include <cstdio>

using namespace redux::util;

#define BitValue(b) (1UL<<b)

uint32_t redux::file::anacrunchrun8(uint8_t* x, const uint8_t* array, int slice, int nx, int ny, uint32_t limit, int t_endian)
/* compress 8 bit array into x (a byte array) using ny blocks each of size
   nx, bit slice size slice, returns # of bytes in x */
{
    Ana::compressed_header* ch;
    const uint8_t* p;
    uint32_t nb;
    uint32_t outBytes, j, outBits;
    int r0, r3, mask, nrun, lrun, ic;
    int* dif, *d, nc, zq, yq, *dd;
    int i2, k, iy;
    union {
        int i;
        short w;
        unsigned char b[4];
    } y;

    /* begin execution */
    if(limit < 25) {
        fprintf(stderr, "limit (%d) too small in crunchrun8\n", limit);
        return -1;
    }
    limit = limit - 24; /* need 14 for header and some margin since
                   we don't check all times */
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
    /* do the compression header */
    ch = reinterpret_cast<Ana::compressed_header*>(x);
    /* important note - can't use the sizeof(struct compresshead) because it
       is 14 on some machines and rounded up to 16 on others */
    x = x + 14;
    ch->bsize = nx;
    ch->nblocks = ny;
    ch->slice_size = slice;
    ch->type = 3;
    outBytes = 0;
    outBits = 0;
    dif = (int*)malloc(nx * 4);      /* line buffer */
    for(iy = 0; iy < ny; iy++) {    /* start of iy (outer) loop */
        /* load the first value */
        x[outBytes++] = array[iy * nx];

        /* compute and store the first differences for this line */
        p = (array + nx * iy);
        nc = nx - 1;
        d = dif;
        yq = (int) * p++;
        zq = (int) * p++;
        while(nc--) {
            *d++ = zq - yq;
            yq = zq;
            zq = (int) * p++;
        }
        outBits = outBits + 8;
        p = (array + nx * iy);
        nc = nx - 1;
        d = dif;
        ic = outBytes++;   /* count position */
        outBits = outBits + 8;    /* to cover first count */
        lrun = 0;   /* literal count) */
        while(1) {
            /* look for a run */
            y.i = *d++;
            if(nc > 1) {
                while(y.i == *d) {   /* at least a run of 2 */
                    dd = d + 1;
                    nrun = 2;
                    while(nc-- > 2 && y.i == *dd) {
                        nrun++;
                        dd++;
                    }
                    /* short runs are not worth it, make the legal limit 4 */
                    if(nrun >= 4) {     /* code the run */
                        /* a previous literal ? */
                        if(lrun != 0) {
                            x[ic] = lrun;
                            outBytes = (outBits + 7) / 8;
                            lrun = 0;
                        }
                        else
                            outBytes = ic;
                        while(nrun > 128) {     /* a big one, multiple runs */
                            /* need only 2 bytes to represent run, runs can't be 17 bits */
                            if(nrun == 129) {   /* beware the dreaded 129 */
                                x[outBytes++] = 0x82;
                                nrun -= 127;
                            }
                            else {
                                x[outBytes++] = 0x81;
                                nrun -= 128;
                            }
                            if(t_endian) {
                                x[outBytes++] = y.b[3];
                                x[outBytes++] = y.b[2];
                            }
                            else {
                                x[outBytes++] = y.b[0];
                                x[outBytes++] = y.b[1];
                            }
                        }
                        if(t_endian) {   // big endian
                            x[outBytes++] = -(nrun - 1);
                            x[outBytes++] = y.b[3];
                            x[outBytes++] = y.b[2];
                        }
                        else {
                            x[outBytes++] = -(nrun - 1);
                            x[outBytes++] = y.b[0];
                            x[outBytes++] = y.b[1];
                        }
                        /* prepare for a literal and for next run check */
                        nc--;
                        if(nc <= 0)
                            goto ended_on_run;
                        lrun = 0;
                        ic = outBytes++;
                        outBits = 8 * outBytes;
                        d = dd;
                        y.i = *d++;
                    }
                    else {
                        nc = nc + nrun - 1;
                        break;
                    }
                }   /* not a run, do next literal, assume setup for literals */
            }
            else if(nc <= 0)
                break;
            nc--;
            /* increment the literal count */
            if(++lrun > 127) {   /* need a new literal run count */
                x[ic] = 127;
                /* bump to next byte boundary */
                outBytes = (outBits + 7) / 8;
                ic = outBytes++;
                outBits = 8 * outBytes;
                lrun = 1;
            }
            /* first the fixed slice portion */
            r3 = (y.i >> slice);
            outBytes = outBits >> 3;   /* byte number */
            j = outBits & 7; /* bit number */
            if(outBytes > limit)
                return -1;  /* bad news, went too far */
            /* now load nb bytes into x */
            /*low order byte of y.i is first in stream */
            if(t_endian) {   // big endian
                if(j == 0) {
                    y.i = (y.i & mask);
                    x[outBytes] = y.b[3];
                }
                else {
                    y.i = (y.i & mask) << j;
                    x[outBytes] = x[outBytes] | y.b[3];
                }
                if(nb > 1) {
                    x[outBytes + 1] = y.b[2];
                    if(nb > 2)
                        x[outBytes + 2] = y.b[1];
                }
            }
            else {
                if(j == 0) {
                    y.i = (y.i & mask);
                    x[outBytes] = y.b[0];
                }
                else {
                    y.i = (y.i & mask) << j;
                    x[outBytes] = x[outBytes] | y.b[0];
                }
                if(nb > 1) {
                    x[outBytes + 1] = y.b[1];
                    if(nb > 2)
                        x[outBytes + 2] = y.b[2];
                }
            }

            outBits = outBits + slice;    /* bump outBits pass the fixed part */
            outBytes = outBits >> 3;
            j = outBits & 7;
            /* note that r3 is the # of bits required minus 1 */
            if(r3 == 0) {
                if(j == 0) {
                    x[outBytes] = BitValue(j);
                }
                else {
                    x[outBytes] = x[outBytes] | BitValue(j);
                }
                outBits += 1;
            }
            else {
                r3 = 2 * r3;
                if(r3 < 0)
                    r3 = -r3 - 1;
                if(r3 < 31) {
                    r0 = j + r3;    /* this is the bit that needs setting offset from x[outBytes] */
                    if(r0 < 8) {
                        if(j == 0)
                            x[outBytes] = BitValue(r0);
                        else
                            x[outBytes] = x[outBytes] | BitValue(r0);
                    }
                    else {
                        if(j == 0)
                            x[outBytes] = 0;
                        j = r0 % 8;
                        if(r0 < 16)
                            x[outBytes + 1] = BitValue(j);
                        else {
                            i2 = outBytes + r0 / 8;
                            for(k = outBytes + 1; k < i2; k++)
                                x[k] = 0;
                            x[i2] = BitValue(j);
                        }
                    }
                    outBits += 1 + r3;
                }
                else {      /* big one exception, should be rare */
                    /* does not need to be efficient, used rarely */
                    if(j == 0)
                        x[outBytes] = 0;   /* gotta zero the virgins */
                    r0 = j + 31;
                    j = r0 % 8;
                    i2 = outBytes + r0 / 8;
                    for(k = outBytes + 1; k < i2; k++)
                        x[k] = 0;
                    x[i2] = BitValue(j);
                    /* recompute the difference and load 9 bits (always 2 bytes) */
                    outBits = outBits + 32;
                    outBytes = outBits / 8;
                    j = outBits % 8;
                    if(j == 0)
                        x[outBytes] = 0;   /* gotta zero the virgins */
                    y.i = ((*(d - 1)) & 0x1ff) << j;
                    if(t_endian) {   // big endian
                        x[outBytes] = x[outBytes] | y.b[3];
                        x[outBytes + 1] = y.b[2];
                    }
                    else {
                        x[outBytes] = x[outBytes] | y.b[0];
                        x[outBytes + 1] = y.b[1];
                    }
                    outBits = outBits + 9;
                }   /* end of big one exception */
            }   /* end of (r3==0) conditional */
        }   /* end of ix loop */
        /* some checks here */
        /* bump to next byte boundary */
        /* a final literal ? */
        if(lrun != 0) {
            x[ic] = lrun;
            lrun = 0;
        }
        outBytes = (outBits + 7) / 8;
    ended_on_run:
        outBits = 8 * outBytes;
    }           /* end of iy loop */
    ch->tsize = outBytes = outBytes + 14;
    /* we have to put these in a form readable by the Vax (these may be used
       by fcwrite) */
    if(t_endian) {      // big endian
        swapEndian(&(ch->tsize));
        swapEndian(&(ch->bsize));
        swapEndian(&(ch->nblocks));
    }
    free(dif);
    return outBytes;       /*return # of bytes used */
}

uint32_t redux::file::anacrunch8(uint8_t* x, const uint8_t* array, int slice, int nx, int ny, uint32_t limit, int t_endian)
/* compress 8 bit array into x (a byte array) using ny blocks each of size
   nx, bit slice size slice, returns # of bytes in x */
{
    Ana::compressed_header* ch;
    uint32_t nb, ixa, ixb;
    uint32_t outBytes, j, outBits, inIndex;
    int r0, r3, mask;
    int i2, k, iy;
    union {
        int i;
        short w;
        unsigned char b[4];
    } y;

    /* begin execution */
    if(limit < 25) {
        fprintf(stderr, "limit (%d) too small in crunch8\n", limit);
        return -1;
    }
    limit = limit - 24; /* need 14 for header and some margin since
                   we don't check all times */
    mask = (1 << slice) - 1;
    /* determine the # of bytes to transfer to 32 bit int for fixed portion */
    if(slice > 8)
        slice = 8;
    if(slice == 0)
        nb = 0;
    else if(slice < 2)
        nb = 1;
    else if(slice < 10)
        nb = 2;
    else
        nb = 3;
    y.i = 0;
    /* do the compression header */
    ch = reinterpret_cast<Ana::compressed_header*>(x);
    /* important note - can't use the sizeof(struct compresshead) because it
       is 14 on some machines and rounded up to 16 on others */
    x = x + 14;
    ch->bsize = nx;
    ch->nblocks = ny;
    ch->slice_size = slice;
    ch->type = 1;
    outBytes = 0;
    outBits = 0;
    inIndex = 0;
    for(iy = 0; iy < ny; iy++) {    /* start of iy (outer) loop */
        /* load the first value */
        x[outBytes] = array[inIndex];
        outBits = outBits + 8;
        ixa = 1 + iy * nx;
        ixb = (iy + 1) * nx;
        for(inIndex = ixa; inIndex < ixb; inIndex++) {     /* start of ix (inner) loop */
            /* first the fixed slice portion */
            y.i = (int)array[inIndex] - (int)array[inIndex - 1];
            r3 = (y.i >> slice);
            outBytes = outBits >> 3;
            j = outBits % 8;
            if(outBytes > limit)
                return -1;  /* bad news, went too far */
            /* now load nb bytes into x */
            /*low order byte of y.i is first in stream */
            if(t_endian) {   // big endian
                if(j == 0) {
                    y.i = (y.i & mask);
                    x[outBytes] = y.b[3];
                }
                else {
                    y.i = (y.i & mask) << j;
                    x[outBytes] = x[outBytes] | y.b[3];
                }
                if(nb > 1) {
                    x[outBytes + 1] = y.b[2];
                }
            }
            else {
                if(j == 0) {
                    y.i = (y.i & mask);
                    x[outBytes] = y.b[0];
                }
                else {
                    y.i = (y.i & mask) << j;
                    x[outBytes] = x[outBytes] | y.b[0];
                }
                if(nb > 1) {
                    x[outBytes + 1] = y.b[1];
                }
            }
            outBits = outBits + slice;    /* bump outBits pass the fixed part */
            outBytes = outBits >> 3;
            j = outBits % 8;
            /* note that r3 is the # of bits required minus 1 */
            if(r3 == 0) {
                if(j == 0) {
                    x[outBytes] = BitValue(j);
                }
                else {
                    x[outBytes] = x[outBytes] | BitValue(j);
                }
                outBits += 1;
            }
            else {
                r3 = 2 * r3;
                if(r3 < 0)
                    r3 = -r3 - 1;
                if(r3 < 31) {
                    r0 = j + r3;    /* this is the bit that needs setting offset from x[outBytes] */
                    if(r0 < 8) {
                        if(j == 0)
                            x[outBytes] = BitValue(r0);
                        else
                            x[outBytes] = x[outBytes] | BitValue(r0);
                    }
                    else {
                        if(j == 0)
                            x[outBytes] = 0;
                        j = r0 % 8;
                        if(r0 < 16)
                            x[outBytes + 1] = BitValue(j);
                        else {
                            i2 = outBytes + r0 / 8;
                            for(k = outBytes + 1; k < i2; k++)
                                x[k] = 0;
                            x[i2] = BitValue(j);
                        }
                    }
                    outBits += 1 + r3;
                }
                else {      /* big one exception, should be rare */
                    /* does not need to be efficient, used rarely */
                    if(j == 0)
                        x[outBytes] = 0;   /* gotta zero the virgins */
                    r0 = j + 31;
                    j = r0 % 8;
                    i2 = outBytes + r0 / 8;
                    for(k = outBytes + 1; k < i2; k++)
                        x[k] = 0;
                    x[i2] = BitValue(j);
                    /* recompute the difference and load 9 bits (always 2 bytes) */
                    outBits = outBits + 32;
                    outBytes = outBits / 8;
                    j = outBits % 8;
                    if(j == 0)
                        x[outBytes] = 0;   /* gotta zero the virgins */
                    y.i = ((array[inIndex] - array[inIndex - 1]) & 0x1ff) << j;
                    if(t_endian) {   // big endian
                        x[outBytes] = x[outBytes] | y.b[3];
                        x[outBytes + 1] = y.b[2];
                    }
                    else {
                        x[outBytes] = x[outBytes] | y.b[0];
                        x[outBytes + 1] = y.b[1];
                    }
                    outBits = outBits + 9;
                }   /* end of big one exception */
            }   /* end of (r3==0) conditional */
        }   /* end of ix loop */
        /* some checks here */
        /* bump to next byte boundary */
        outBytes = (outBits + 7) / 8;
        outBits = 8 * outBytes;
    }           /* end of iy loop */
    ch->tsize = outBytes = outBytes + 14;
    /* we have to put these in a form readable by the Vax (these may be used
       by fcwrite) */
    if(t_endian) {      // big endian
        swapEndian(&(ch->tsize));
        swapEndian(&(ch->bsize));
        swapEndian(&(ch->nblocks));
    }
    return outBytes;       /*return # of bytes used */
}               /* end of routine */

uint32_t redux::file::anacrunchrun(uint8_t* x, const int16_t* array, int slice, int nx, int ny, uint32_t limit, int t_endian)
/* compress 16 bit array into x (a byte array) using ny blocks each of size
   nx, bit slice size slice, returns # of bytes in x */
{
    Ana::compressed_header* ch;
    short* p;
    uint32_t nb;
    uint32_t outBytes, j, outBits;
    int r0, r3, mask, nrun, lrun, ic;
    int* dif, *d, nc, zq, yq, *dd;
    int i2, k, iy;
    union {
        int i;
        short w;
        unsigned char b[4];
    } y;

    /* begin execution */
    if(limit < 25) {
        fprintf(stderr, "limit (%d) too small in crunchrun\n", limit);
        return -1;
    }
    limit = limit - 24; /* need 14 for header and some margin since
                   we don't check all times */
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
    /* do the compression header */
    ch = reinterpret_cast<Ana::compressed_header*>(x);
    /* important note - can't use the sizeof(struct compresshead) because it
       is 14 on some machines and rounded up to 16 on others */
    x = x + 14;
    ch->bsize = nx;
    ch->nblocks = ny;
    ch->slice_size = slice;
    ch->type = 2;
    outBytes = 0;
    outBits = 0;
    dif = (int*)malloc(nx * 4);      /* line buffer */
    for(iy = 0; iy < ny; iy++) {    /* start of iy (outer) loop */
        /* load the first value, reverse bytes (VAX style) */
        if(t_endian) {   // big endian
            y.w = array[iy * nx];
            x[outBytes++] = y.b[1];
            x[outBytes++] = y.b[0];
        }
        else {      // little endian
            y.w = array[iy * nx];
            x[outBytes++] = y.b[0];
            x[outBytes++] = y.b[1];
        }
        /* compute and store the first differences for this line */
        p = (short int*)(array + nx * iy);
        nc = nx - 1;
        d = dif;
        yq = (int) * p++;
        zq = (int) * p++;
        while(nc--) {
            *d++ = zq - yq;
            yq = zq;
            zq = (int) * p++;
        }
        outBits = outBits + 16;
        p = (short int*)(array + nx * iy);
        nc = nx - 1;
        d = dif;
        ic = outBytes++;   /* count position */
        outBits = outBits + 8;    /* to cover first count */
        lrun = 0;   /* literal count) */
        while(1) {
            /* look for a run */
            y.i = *d++;
            if(nc > 1) {
                while(y.i == *d) {   /* at least a run of 2 */
                    dd = d + 1;
                    nrun = 2;
                    while(nc-- > 2 && y.i == *dd) {
                        nrun++;
                        dd++;
                    }
                    /* short runs are not worth it, make the legal limit 4 */
                    if(nrun >= 4) {     /* code the run */
                        /* a previous literal ? */
                        if(lrun != 0) {
                            x[ic] = lrun;
                            outBytes = (outBits + 7) / 8;
                            lrun = 0;
                        }
                        else
                            outBytes = ic;
                        while(nrun > 128) {     /* a big one, multiple runs */
                            /* need only 2 bytes to represent run, runs can't be 17 bits */
                            if(nrun == 129) {   /* beware the dreaded 129 */
                                x[outBytes++] = 0x82;
                                nrun -= 127;
                            }
                            else {
                                x[outBytes++] = 0x81;
                                nrun -= 128;
                            }
                            if(t_endian) {   // big endian
                                x[outBytes++] = y.b[3];
                                x[outBytes++] = y.b[2];
                            }
                            else {
                                x[outBytes++] = y.b[0];
                                x[outBytes++] = y.b[1];
                            }
                        }
                        if(t_endian) {   // big endian
                            x[outBytes++] = -(nrun - 1);
                            x[outBytes++] = y.b[3];
                            x[outBytes++] = y.b[2];
                        }
                        else {
                            x[outBytes++] = -(nrun - 1);
                            x[outBytes++] = y.b[0];
                            x[outBytes++] = y.b[1];
                        }
                        /* prepare for a literal and for next run check */
                        nc--;
                        if(nc <= 0)
                            goto ended_on_run;
                        lrun = 0;
                        ic = outBytes++;
                        outBits = 8 * outBytes;
                        d = dd;
                        y.i = *d++;
                    }
                    else {
                        nc = nc + nrun - 1;
                        break;
                    }
                }   /* not a run, do next literal, assume setup for literals */
            }
            else if(nc <= 0)
                break;
            nc--;
            /* increment the literal count */
            if(++lrun > 127) {   /* need a new literal run count */
                x[ic] = 127;
                /* bump to next byte boundary */
                outBytes = (outBits + 7) / 8;
                ic = outBytes++;
                outBits = 8 * outBytes;
                lrun = 1;
            }
            /* first the fixed slice portion */
            r3 = (y.i >> slice);
            outBytes = outBits >> 3;   /* byte number */
            j = outBits & 7; /* bit number */
            if(outBytes > limit)
                return -1;  /* bad news, went too far */
            /* now load nb bytes into x */
            /*low order byte of y.i is first in stream */
            if(t_endian) {   // big endian
                if(j == 0) {
                    y.i = (y.i & mask);
                    x[outBytes] = y.b[3];
                }
                else {
                    y.i = (y.i & mask) << j;
                    x[outBytes] = x[outBytes] | y.b[3];
                }
                if(nb > 1) {
                    x[outBytes + 1] = y.b[2];
                    if(nb > 2)
                        x[outBytes + 2] = y.b[1];
                }
            }
            else {
                if(j == 0) {
                    y.i = (y.i & mask);
                    x[outBytes] = y.b[0];
                }
                else {
                    y.i = (y.i & mask) << j;
                    x[outBytes] = x[outBytes] | y.b[0];
                }
                if(nb > 1) {
                    x[outBytes + 1] = y.b[1];
                    if(nb > 2)
                        x[outBytes + 2] = y.b[2];
                }
            }
            outBits = outBits + slice;    /* bump outBits pass the fixed part */
            outBytes = outBits >> 3;
            j = outBits & 7;
            /* note that r3 is the # of bits required minus 1 */
            if(r3 == 0) {
                if(j == 0) {
                    x[outBytes] = BitValue(j);
                }
                else {
                    x[outBytes] = x[outBytes] | BitValue(j);
                }
                outBits += 1;
            }
            else {
                r3 = 2 * r3;
                if(r3 < 0)
                    r3 = -r3 - 1;
                if(r3 < 31) {
                    r0 = j + r3;    /* this is the bit that needs setting offset from x[outBytes] */
                    if(r0 < 8) {
                        if(j == 0)
                            x[outBytes] = BitValue(r0);
                        else
                            x[outBytes] = x[outBytes] | BitValue(r0);
                    }
                    else {
                        if(j == 0)
                            x[outBytes] = 0;
                        j = r0 % 8;
                        if(r0 < 16)
                            x[outBytes + 1] = BitValue(j);
                        else {
                            i2 = outBytes + r0 / 8;
                            for(k = outBytes + 1; k < i2; k++)
                                x[k] = 0;
                            x[i2] = BitValue(j);
                        }
                    }
                    outBits += 1 + r3;
                }
                else {      /* big one exception, should be rare */
                    /* does not need to be efficient, used rarely */
                    if(j == 0)
                        x[outBytes] = 0;   /* gotta zero the virgins */
                    r0 = j + 31;
                    j = r0 % 8;
                    i2 = outBytes + r0 / 8;
                    for(k = outBytes + 1; k < i2; k++)
                        x[k] = 0;
                    x[i2] = BitValue(j);
                    /* recompute the difference and load 17 bits (always 3 bytes) */
                    outBits = outBits + 32;
                    outBytes = outBits / 8;
                    j = outBits % 8;
                    if(j == 0)
                        x[outBytes] = 0;   /* gotta zero the virgins */
                    y.i = ((*(d - 1)) & 0x1ffff) << j;
                    if(t_endian) {   // big endian
                        x[outBytes] = x[outBytes] | y.b[3];
                        x[outBytes + 1] = y.b[2];
                        x[outBytes + 2] = y.b[1];
                    }
                    else {
                        x[outBytes] = x[outBytes] | y.b[0];
                        x[outBytes + 1] = y.b[1];
                        x[outBytes + 2] = y.b[2];
                    }
                    outBits = outBits + 17;
                }   /* end of big one exception */
            }   /* end of (r3==0) conditional */
        }   /* end of ix loop */
        /* some checks here */
        /* bump to next byte boundary */
        /* a final literal ? */
        if(lrun != 0) {
            x[ic] = lrun;
            lrun = 0;
        }
        outBytes = (outBits + 7) / 8;
    ended_on_run:
        outBits = 8 * outBytes;
    }           /* end of iy loop */
    ch->tsize = outBytes = outBytes + 14;
    /* we have to put these in a form readable by the Vax (these may be used
       by fcwrite) */
    if(t_endian) {      // big endian
        swapEndian(&(ch->tsize));
        swapEndian(&(ch->bsize));
        swapEndian(&(ch->nblocks));
    }
    free(dif);
    return outBytes;       /*return # of bytes used */
}               /* end of routine */

uint32_t redux::file::anacrunch(uint8_t* x, const int16_t* array, int slice, int nx, int ny, uint32_t limit, int t_endian)
// compress 16 bit array into x (a byte array) using ny blocks each of size
// nx, bit slice size slice, returns # of bytes in x
{

    uint32_t outBytes, j, outBits, inIndex;
    int r0, r3, mask;
    union {
        int i;
        short w;
        unsigned char b[4];
    } y;

    if(limit < 25) {
        fprintf(stderr, "limit (%d) too small in crunch\n", limit);
        return -1;
    }
    limit -= 24;        // need 14 for header and some margin since we don't check all times
    mask = (1 << slice) - 1;
    unsigned nb;        // determine the # of bytes to transfer to 32 bit int for fixed portion

    if(slice == 0)
        nb = 0;
    else if(slice < 2)
        nb = 1;
    else if(slice < 10)
        nb = 2;
    else
        nb = 3;
    y.i = 0;
    Ana::compressed_header* ch = reinterpret_cast<Ana::compressed_header*>(x);   // the compression header

// important note: can't use the sizeof(struct compresshead) because it
//                 is 14 on some machines and rounded up to 16 on others
    x += 14;
    ch->bsize = nx;
    ch->nblocks = ny;
    ch->slice_size = slice;
    ch->type = 0;
    outBits = 0;         // outBits is the bit index in the stream...?
    inIndex = 0;         // in is the byte index in the uncompressed stream...?
    outBytes = 0;         // i is the byte index in the compressed stream...?
    int iy;

    for(iy = 0; iy < ny; ++iy) {    // start of iy (outer) loop
        y.w = array[inIndex];
        if(t_endian) {   // load the first value, reverse bytes (VAX style)
            x[outBytes] = y.b[1];
            x[outBytes + 1] = y.b[0];
        }
        else {
            x[outBytes + 1] = y.b[1];
            x[outBytes] = y.b[0];
        }
        outBits += 16;
        unsigned int iynx = iy * nx;

        for(inIndex = iynx + 1; inIndex < iynx + nx; ++inIndex) {   // start of ix (inner) loop
            y.i = array[inIndex] - array[inIndex - 1];    // first the fixed slice portion
            r3 = (y.i >> slice);
            outBytes = outBits >> 3;   // compressed data size (number of bits/8)
            j = outBits % 8;
            if(outBytes > limit)
                return -1;  // bad news: compressed data too big...
            if(j == 0) {    // now load nb bytes into x, low order byte of y.i is first in stream
                y.i = (y.i & mask);
                x[outBytes] = (uint8_t) y.i;
                if(slice > 8)
                    x[outBytes + 1] = (uint8_t)(y.i >> 8);     // since we started at bit 0, spillover to the next byte is determined as follows (and is unlikely since slice gt 8 is unusual
            }
            else {
                y.i = (y.i & mask) << j;
                x[outBytes] = x[outBytes] | (uint8_t) y.i;
                if(nb > 1) {    // spillover more likely here
                    x[outBytes + 1] = (uint8_t)(y.i >> 8);
                    if(nb > 2)
                        x[outBytes + 2] = (uint8_t)(y.i >> 16);
                }
            }
            outBits += slice;    // bump outBits pass the fixed part
            outBytes = outBits >> 3;
            j = outBits % 8;
            if(r3 == 0) {   // note that r3 is the # of bits required minus 1
                if(j == 0) {
                    x[outBytes] = 1;
                }
                else {
                    x[outBytes] = x[outBytes] | BitValue(j);
                }
                ++outBits;
            }
            else {
                r3 *= 2;
                if(r3 < 0)
                    r3 = -r3 - 1;
                if(r3 < 31) {
                    r0 = j + r3;    // this is the bit that needs setting offset from x[outBytes]
                    if(r0 < 8) {
                        if(j == 0)
                            x[outBytes] = BitValue(r0);
                        else
                            x[outBytes] = x[outBytes] | BitValue(r0);
                    }
                    else {      // note, discovered on 2/3/96 that the j==0 case not done above for the sunbow version, was OK in the umbra version, may have happened while cleaning up code?, caused extra bits to be set if x[outBytes] wasn't zero
                        if(j == 0)
                            x[outBytes] = 0;
                        j = r0 % 8;
                        if(r0 < 16) {
                            x[outBytes + 1] = BitValue(j);
                        }
                        else {
                            int i2 = outBytes + r0 / 8;
                            int k;

                            for(k = outBytes + 1; k < i2; ++k)
                                x[k] = 0;
                            x[i2] = BitValue(j);
                        }
                    }
                    outBits += r3 + 1;
                }
                else {      // big one exception, should be rare so does not need to be efficient
                    if(j == 0)
                        x[outBytes] = 0;   // gotta zero the virgins
                    r0 = j + 31;
                    j = r0 % 8;
                    int i2 = outBytes + r0 / 8;
                    int k;

                    for(k = outBytes + 1; k < i2; ++k)
                        x[k] = 0;
                    x[i2] = BitValue(j);
                    outBits += 32;   // recompute the difference and load 17 bits (always 3 bytes)
                    outBytes = outBits / 8;
                    j = outBits % 8;
                    if(j == 0)
                        x[outBytes] = 0;   // gotta zero the virgins
                    y.i = ((array[inIndex] - array[inIndex - 1]) & 0x1ffff) << j;
                    if(t_endian) {   // big endian
                        x[outBytes] = x[outBytes] | y.b[3];
                        x[outBytes + 1] = y.b[2];
                        x[outBytes + 2] = y.b[1];
                    }
                    else {
                        x[outBytes] = x[outBytes] | y.b[0];
                        x[outBytes + 1] = y.b[1];
                        x[outBytes + 2] = y.b[2];
                    }
                    outBits += 17;
                }   // end of big one exception
            }   // end of (r3==0) conditional
        }   // end of ix loop
        outBytes = (outBits + 7) / 8;   // bump to next byte boundary
        outBits = 8 * outBytes;
    }           // end of iy loop
    ch->tsize = (outBytes += 14);
    if(t_endian) {      // we have to put these in a form readable by the Vax (these may be used by fcwrite)
        swapEndian(&(ch->tsize));
        swapEndian(&(ch->bsize));
        swapEndian(&(ch->nblocks));
    }
    return outBytes;       // return # of bytes used
}

uint32_t redux::file::anacrunch32(uint8_t* x, const int32_t* array, int slice, int nx, int ny, uint32_t limit, int t_endian)
/* compress 32 bit array into x (a byte array) using ny blocks each of size
   nx, bit slice size slice, returns # of bytes in x */
{

    Ana::compressed_header* ch;
    uint32_t nb, ixa, ixb, big = 0;
    uint32_t outBytes, j, outBits, inIndex;
    int r0;
    long long r3, mask, y64;
    int i2, k, iy;
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
    if(limit < 25) {
        fprintf(stderr, "limit (%d) too small in crunch32\n", limit);
        return -1;
    }
    limit = limit - 24; /* need 14 for header and some margin since
                   we don't check all times */
    mask = (1 << slice) - 1;
    /* determine the # of bytes to transfer to 32 bit int for fixed portion */
    nb = (slice + 14) / 8;   /* range 1 to 5 */
    if(slice == 0)
        nb = 0;     /* but slice = 0 a special case */

    y.i = 0;
    /* do the compression header */
    ch = reinterpret_cast<Ana::compressed_header*>(x);
    /* important note - can't use the sizeof(struct compresshead) because it
       is 14 on some machines and rounded up to 16 on others */
    x = x + 14;
    ch->bsize = nx;
    ch->nblocks = ny;
    ch->slice_size = slice;
    ch->type = 4;
    outBytes = 0;
    outBits = 0;
    inIndex = 0;
    for(iy = 0; iy < ny; iy++) {    /* start of iy (outer) loop */
        /* load the first value, reverse bytes (VAX style) */
        if(t_endian) {   // big endian
            y.i = array[inIndex];
            x[outBytes] = y.b[3];
            x[outBytes + 1] = y.b[2];
            x[outBytes + 2] = y.b[1];
            x[outBytes + 3] = y.b[0];
        }
        else {
            y.i = array[inIndex];
            x[outBytes] = y.b[0];
            x[outBytes + 1] = y.b[1];
            x[outBytes + 2] = y.b[2];
            x[outBytes + 3] = y.b[3];
        }
        outBits = outBits + 32;
        ixa = 1 + iy * nx;
        ixb = (iy + 1) * nx;
        for(inIndex = ixa; inIndex < ixb; inIndex++) {     /* start of ix (inner) loop */
            /* first the fixed slice portion */
            y64 = (long long)array[inIndex] - (long long)array[inIndex - 1];
            r3 = (y64 >> slice);
            outBytes = outBits >> 3;
            j = outBits % 8;
            if(outBytes > limit) {
                fprintf(stderr, "went to far: limit (%d) in crunch32\n", limit);
                return -1;  /* bad news, went too far */
            }
            /* now load nb bytes into x */
            /*low order byte of y.i is first in stream */
            if(j == 0) {
                y64 = (y64 & mask);
                x[outBytes] = (uint8_t) y64;
                /* since we started at bit 0, spillover to the next byte is
                   determined as follows (and is unlikely since slice gt 8 is unusual */
                if(slice > 8)
                    x[outBytes + 1] = (uint8_t)(y64 >> 8);
                if(slice > 16)
                    x[outBytes + 2] = (uint8_t)(y64 >> 16);
                if(slice > 24)
                    x[outBytes + 3] = (uint8_t)(y64 >> 24);
            }
            else {
                y64 = (y64 & mask) << j;
                x[outBytes] = x[outBytes] | (uint8_t) y64;
                /* spillover more likely here */
                if(nb > 1)
                    x[outBytes + 1] = (uint8_t)(y64 >> 8);
                if(nb > 2)
                    x[outBytes + 2] = (uint8_t)(y64 >> 16);
                if(nb > 3)
                    x[outBytes + 3] = (uint8_t)(y64 >> 24);
                if(nb > 4)
                    x[outBytes + 4] = (uint8_t)(y64 >> 32);
            }

            outBits = outBits + slice;    /* bump outBits pass the fixed part */
            outBytes = outBits >> 3;
            j = outBits % 8;
            /* note that r3 is the # of bits required minus 1 */
            if(r3 == 0) {
                if(j == 0) {
                    x[outBytes] = 1;
                }
                else {
                    x[outBytes] = x[outBytes] | BitValue(j);
                }
                outBits++;
            }
            else {
                r3 = 2 * r3;
                if(r3 < 0)
                    r3 = -r3 - 1;
                if(r3 < 31) {
                    r0 = j + r3;    /* this is the bit that needs setting offset from x[outBytes] */
                    if(r0 < 8) {
                        if(j == 0)
                            x[outBytes] = BitValue(r0);
                        else
                            x[outBytes] = x[outBytes] | BitValue(r0);
                    }
                    else {
                        if(j == 0)
                            x[outBytes] = 0;
                        j = r0 % 8;
                        if(r0 < 16)
                            x[outBytes + 1] = BitValue(j);
                        else {
                            i2 = outBytes + r0 / 8;
                            for(k = outBytes + 1; k < i2; k++)
                                x[k] = 0;
                            x[i2] = BitValue(j);
                        }
                    }
                    outBits += 1 + r3;
                }
                else {      /* big one exception, should be rare */
                    /* does not need to be efficient, used rarely */
                    big++;
                    if(j == 0)
                        x[outBytes] = 0;   /* gotta zero the virgins */
                    r0 = j + 31;
                    j = r0 % 8;
                    i2 = outBytes + r0 / 8;
                    for(k = outBytes + 1; k < i2; k++)
                        x[k] = 0;
                    x[i2] = BitValue(j);
                    /* recompute the difference and load 33 bits (always 5 bytes) */
                    outBits = outBits + 32;
                    outBytes = outBits / 8;
                    j = outBits % 8;
                    if(j == 0)
                        x[outBytes] = 0;   /* gotta zero the virgins */
                    yy.l64 = (((long long int)array[inIndex] - (long long int)array[inIndex - 1]) & 0x1ffffffffLL) << j;
                    if(t_endian) {   // big endian
                        x[outBytes] = x[outBytes] | yy.b[7];
                        x[outBytes + 1] = yy.b[6];
                        x[outBytes + 2] = yy.b[5];
                        x[outBytes + 3] = yy.b[4];
                        x[outBytes + 4] = yy.b[3];
                    }
                    else {
                        x[outBytes] = x[outBytes] | yy.b[0];
                        x[outBytes + 1] = yy.b[1];
                        x[outBytes + 2] = yy.b[2];
                        x[outBytes + 3] = yy.b[3];
                        x[outBytes + 4] = yy.b[4];
                    }
                    outBits = outBits + 33;
                }   /* end of big one exception */
            }   /* end of (r3==0) conditional */
        }  /* end of ix loop */
        /* some checks here */
        /* bump to next byte boundary */
        outBytes = (outBits + 7) / 8;
        outBits = 8 * outBytes;
    }           /* end of iy loop */
    ch->tsize = outBytes = outBytes + 14;
    /* we have to put these in a form readable by the Vax (these may be used
       by fcwrite) */
    if(t_endian) {      // big endian
        swapEndian(&(ch->tsize));
        swapEndian(&(ch->bsize));
        swapEndian(&(ch->nblocks));
    }
    return outBytes;       /*return # of bytes used */
}               /* end of routine */
