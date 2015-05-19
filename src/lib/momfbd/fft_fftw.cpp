#include "redux/momfbd/fft.hpp"

#include <fftw3.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

//#include "types.h"
//#include "mem.h"
#include "redux/util/arrayutil.hpp"
using namespace redux::util;
using namespace redux;
//
// WARNING!!! fftw returns the FT^* and not the FT (i.e. the complex
// part has different sign).
//

struct fft_plan_str {
    fftw_plan fftf_plan, fftb_plan;
    int n1, n2;
};

struct fft_plan_str **plans = 0;

void redux::fft_init(int n1, int n2) {
    if(!plans) { // first call to fft_init
        plans = new fft_plan_str* [1];
        plans[0] = 0;
    }
    for(int p = 0; plans[p]; ++p)
        if((plans[p]->n1 == n1) && (plans[p]->n2 == n2)) return;
    int q;
    for(q = 0; plans[q]; ++q);
    fft_plan_str **tmp = new fft_plan_str* [q + 2];
    memcpy(tmp, plans, q * sizeof(fft_plan_str*));
    delete[] plans;
    plans = tmp;
    plans[q] = new fft_plan_str;
    plans[q]->n1 = n1;
    plans[q]->n2 = n2;
    complex_t **in = newArray<complex_t>(n1,n2);
    complex_t **out = newArray<complex_t>(n1,n2);
    plans[q]->fftf_plan = fftw_plan_dft_2d(n1, n2, (fftw_complex*)(*in), (fftw_complex*)(*out), FFTW_FORWARD, FFTW_MEASURE);
    plans[q]->fftb_plan = fftw_plan_dft_2d(n1, n2, (fftw_complex*)(*in), (fftw_complex*)(*out), FFTW_BACKWARD, FFTW_MEASURE);
    delArray(in);
    delArray(out);
    plans[q + 1] = 0;
}

void redux::fft_2d(complex_t **data, int nY, int nX, int isign) {
    for(int p = 0; plans[p]; ++p) {
        if((plans[p]->n1 == nY) && (plans[p]->n2 == nX)) {
            complex_t **res = newArray<complex_t>(nY,nX);
            if(isign > 0) {
                fftw_execute_dft(plans[p]->fftf_plan, (fftw_complex*)(*data), (fftw_complex*)(*res));
                for(int y = 0; y < nY; ++y)
                    for(int x = 0; x < nX; ++x) res[y][x].imag(-res[y][x].imag());
            }
            else { // inverse
                for(int y = 0; y < nY; ++y)
                    for(int x = 0; x < nX; ++x) data[y][x].imag(-data[y][x].imag());
                fftw_execute_dft(plans[p]->fftb_plan, (fftw_complex*)(*data), (fftw_complex*)(*res));
                double N = (double)(nY * nX);
                for(int y = 0; y < nY; ++y)
                    for(int x = 0; x < nX; ++x) {
                        res[y][x] /= N;
                        //res[x][y].im /= N;
                    }
            }
            memcpy(*data,*res,nY*nX*sizeof(complex_t));
            delArray(res);
            return;
        }
    }
    fprintf(stderr, "Error: dimension mismatch (%dx%d) not found in fft plans (forgot to call fft_init)?\n", nY, nX);
}

void redux::fft_done(int n1, int n2) {
    if(plans) {
        for(int p = 0; plans[p]; ++p)
            if((plans[p]->n1 == n1) && (plans[p]->n2 == n2)) {
                fftw_destroy_plan(plans[p]->fftf_plan);
                fftw_destroy_plan(plans[p]->fftb_plan);
                int q;
                for(q = 0; plans[q]; ++q);
                fft_plan_str **tmp = new fft_plan_str* [q];
                memcpy(tmp, plans, p * sizeof(fft_plan_str*));
                memcpy(tmp + p, plans + p + 1, (q - p)*sizeof(fft_plan_str*));
                delete plans[p];
                delete[] plans;
                plans = tmp;
                p = q - 2; // break the loop
            }
        if(!plans[0]) { // last call to fft_done
            delete[] plans;
            plans = 0;
        }
    }
}

void redux::fft_reorder(complex_t **f, int np) {
    int nh = np / 2;
    complex_t *buf = new complex_t [nh];
    for(int y = 0; y < nh; ++y) {
        memcpy(buf, f[y], nh * sizeof(complex_t));
        memcpy(f[y], f[y + nh] + nh, nh * sizeof(complex_t));
        memcpy(f[y + nh] + nh, buf, nh * sizeof(complex_t));
//
        memcpy(buf, f[y + nh], nh * sizeof(complex_t));
        memcpy(f[y + nh], f[y] + nh, nh * sizeof(complex_t));
        memcpy(f[y] + nh, buf, nh * sizeof(complex_t));
    }
    delete[] buf;
}


