#include "rdx.hpp"
#include "imgtools.hpp"

#include "redux/image/fouriertransform.hpp"
#include "redux/util/arraystats.hpp"

#include <iomanip>
#include <iostream>
#include <typeinfo>

#include <boost/algorithm/string.hpp>

using namespace redux::image;
using namespace redux::util;
using namespace redux;
using namespace std;


namespace {

    const char *var_names[] = { "undef", "int8", "int16", "int32", "float32",
                                "float64", "complex", "string", "struct", "double complex",
                                "int32", "int32", "uint16", "uint32", "int64", "uint64"
                              };

    const size_t var_sizes[] = { 0, sizeof (UCHAR), sizeof (IDL_INT), sizeof (IDL_LONG), sizeof (float),
                                 sizeof (double), sizeof (IDL_COMPLEX), sizeof (IDL_STRING), 0, sizeof (IDL_DCOMPLEX),
                                 sizeof (IDL_ULONG), sizeof (IDL_ULONG), sizeof (IDL_UINT), sizeof (IDL_ULONG), sizeof (IDL_LONG64), sizeof (IDL_ULONG64)
                               };

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
        IDL_INT help;
        IDL_INT indent;
        IDL_INT normalize;
        IDL_INT niter;
        IDL_INT nthreads;
        IDL_INT padding;
        IDL_INT verbose;
    } KW_RESULT;

    // NOTE:  The keywords MUST be listed in alphabetical order !!
    static IDL_KW_PAR kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { (char*) "HELP",      IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF (help) },
        { (char*) "INDENT",    IDL_TYP_INT, 1, 0,           0, (char*) IDL_KW_OFFSETOF (indent) },
        { (char*) "NORMALIZE", IDL_TYP_INT, 1, 0,           0, (char*) IDL_KW_OFFSETOF (normalize) },
        { (char*) "NITER",     IDL_TYP_INT, 1, 0,           0, (char*) IDL_KW_OFFSETOF (niter) },
        { (char*) "NTHREADS",  IDL_TYP_INT, 1, 0,           0, (char*) IDL_KW_OFFSETOF (nthreads) },
        { (char*) "PADDING",   IDL_TYP_INT, 1, 0,           0, (char*) IDL_KW_OFFSETOF (padding) },
        { (char*) "VERBOSE",   IDL_TYP_INT, 1, 0,           0, (char*) IDL_KW_OFFSETOF (verbose) },
        { NULL }
    };


    void rdx_help (void) {

        cout << "RDX HELP\n"
             "       Syntax:   struct_info,data,/KEYWORDS\n"
             "                 sz=struct_info(data,/KEYWORDS)\n"
             "       Accepted Keywords:\n"
             "              HELP              Display this help.\n"
             "              INDENT=[0-30]     Indentation for each level in the structure, default is 2." << endl;

    }

    void convolve_help (void) {

        cout << "CCONVOLVE\n"
             "       Syntax:   out=cconvolve(img,psf)\n"
             "       Accepted Keywords:\n"
             "              HELP                Display this help.\n"
             "              NORMALIZE           Normalize PDF before convolving.\n"
             "              NITER               Number of iterations in descatter loop.\n"
             "              NTHREADS            Number of threads to use for parallel FFT.\n"
             "              PADDING             Number of 0-pixels to add around the boundary for cdescatter.\n"
             "              VERBOSE={0,1,2}     Verbosity, default is 0 (only error output)." << endl;

    }

    template <typename T, typename U>
    void central_derivative (int n, T *x, U *y, double* yprime) {

        memset (yprime, 0, n * sizeof(double));

        double dx = (x[1] - x[0]);
        double der = (y[1] - y[0]) / dx;
        yprime[0] = der;

        for (int k = 1; k < n - 1; k++) {
            double dx1 = x[k + 1] - x[k];
            double der1 = (y[k + 1] - y[k]) / dx1;

            if ( (der * der1) > 0.0) {
                double lambda = (1.0 + dx1 / (dx1 + dx)) / 3.0;
                yprime[k] = (der / (lambda * der1 + (1.0 - lambda) * der)) * der1;
            }

            der = der1;
            dx = dx1;

        }

        yprime[n - 1] = der;

    }
    
}


//
// print the layout of an IDL structure and return the total data-size
//
size_t redux::dumpStruct (IDL_VPTR data, int current, int indent) {

    size_t sz = 0;

    if (data->type == IDL_TYP_STRUCT) {
        if (current < 0) {
            cout << "           TAG           TYPE                     OFFSET         SIZE           " << endl;
            sz += dumpStruct (data, indent, indent);
            cout << string (45, ' ') << "Total Size:         " << to_string (sz) << endl;
            return sz;
        }


        IDL_StructDefPtr structDef = data->value.s.sdef;
        //uint8_t *buf = data->value.s.arr->data;
        int nTags = IDL_StructNumTags (structDef);
        int count;

        for (int t = 0; t < nTags; ++t) {
            char *name = IDL_StructTagNameByIndex (structDef, t, 0, 0);
            IDL_VPTR v;
            IDL_MEMINT offset = IDL_StructTagInfoByIndex (structDef, t, 0, &v);
            string type = string (var_names[v->type]);
            count = 1;

            if (v->flags & IDL_V_ARR) {
                type.append ("(");

                for (int d = 0; d < v->value.arr->n_dim; ++d) {
                    if (d) type.append (",");

                    count *= v->value.arr->dim[d];
                    type.append (to_string ( (int) v->value.arr->dim[d]));
                }

                type.append (")");
            }

            cout.setf (std::ios::left);

            if (v->type == IDL_TYP_STRUCT) {
                cout << std::setw (25) << (string (current, ' ') + name);
                cout << std::setw (25) << type;
                //cout.setf ( std::ios::hex );
                //cout << std::setw(20) << (void*)( buf + offset );
                //cout.setf ( std::ios::dec );
                cout << std::setw (15) << to_string ( (size_t) offset) << endl;
                sz += count * dumpStruct (v, current + indent, indent);
            } else {
                sz += count * var_sizes[v->type];
                cout << std::setw (25) << (string (current, ' ') + name);
                cout << std::setw (25) << type;
                //cout.setf ( std::ios::hex );
                //cout << std::setw(20) << (void*)( buf + offset );
                //cout.setf ( std::ios::dec );
                cout << std::setw (15) << to_string ( (size_t) offset);
                cout << std::setw (15) << to_string (count * var_sizes[v->type]) << endl;   // << endl;
            }

        }

    }

    return sz;
}


void redux::printStruct (int argc, IDL_VPTR argv[], char* argk) {

    if (argc < 1) {
        rdx_help();
        return;
    }

    IDL_VPTR data = argv[0];

    KW_RESULT kw;
    kw.help = 0;
    kw.indent = 2;
    (void) IDL_KWProcessByOffset (argc, argv, argk, kw_pars, (IDL_VPTR*) 0, 255, &kw);

    if (kw.help) {
        rdx_help();
        return;
    }

    int indent = std::min (std::max ( (int) kw.indent, 0), 30);


    dumpStruct (data, -1, indent);
}


IDL_VPTR redux::structInfo (int argc, IDL_VPTR argv[], char* argk) {
    if (argc < 1) {
        rdx_help();
        return IDL_Gettmp();
    }

    IDL_VPTR data = argv[0];

    KW_RESULT kw;
    kw.help = 0;
    kw.indent = 2;
    (void) IDL_KWProcessByOffset (argc, argv, argk, kw_pars, (IDL_VPTR*) 0, 255, &kw);

    if (kw.help) {
        rdx_help();
        return IDL_Gettmp();
    }

    int indent = std::min (std::max ( (int) kw.indent, 0), 30);

    return IDL_GettmpULong (dumpStruct (data, -1, indent));
}


template <typename T, typename U = double>
void descatter (int szY, int szX, const T* img, const T* gain, const T* psf, T* out, int nthreads, int niter, int padding, int verbosity) {

    padding = std::max (std::min (std::min (szX, szY), padding), 0);

    if (verbosity > 0) {
        cout << "descatter<" << typeid (T).name() << ">    size=(" << szY << "," << szX << "), nthreads=" << nthreads
             << ", niter=" << niter << ", padding=" << padding << ", verbosity="  << verbosity << endl;
    }

    int paddedX = szX + 2 * padding;
    int paddedY = szY + 2 * padding;
    int nPixels = szY * szX;

    Array<U> paddedImg (paddedY, paddedX);          // padded array
    Array<U> subImg (paddedImg, padding, paddedY - padding - 1, padding, paddedX - padding - 1);    // centre of padded array
    paddedImg.zero();                               // zero outside of centre.

    subImg.template copyFrom<T> (psf);              // use temporarily for initializing the OTF
    FourierTransform otf (paddedImg, FT_REORDER | FT_NORMALIZE, std::max (nthreads, 1));
    ArrayStats stats;
    stats.getMinMaxMean (subImg);
    otf /= stats.sum;                               // PSF has to be normalized so that sum=1

    if (verbosity > 0) {
        cout << "                PDF-integral=" << stats.sum << ", nElements=" << paddedImg.nElements() << endl;
    }

    const T** g = makePointers (gain, szY, szX);          // 2D pointer to gain
    const T** im = makePointers (img, szY, szX);          // raw img
    U** pim = makePointers (paddedImg.ptr(), paddedY, paddedX);     // padded img
    Array<U> tmp (paddedY, paddedX);                // temporary (convolved) image
    U** cim = makePointers (tmp.ptr(), paddedY, paddedX);           // 2D pointer to convolved image

    // Initialize FT, set normalized-flag since we normalized the OTF above.
    FourierTransform ft (tmp, FT_NORMALIZE, std::max (nthreads, 1));

    for (int y = 0; y < szY; ++y) {                 // set initial value for the padded img
        for (int x = 0; x < szX; ++x) {
            pim[y + padding][x + padding] = im[y][x] / (1.0 + g[y][x] * g[y][x]);
        }
    }

    double cs = 1E10;
    double limit = nPixels * 1E-8;                // multiply with image-size and skip dividing cs with it in every iteration
    int dataSize = paddedX * paddedY * sizeof (U);
    int i (0);

    while ( (cs > limit) && (i++ < niter)) {
        memcpy (*cim, *pim, dataSize);

        for (int y = 0; y < szY; ++y) for (int x = 0; x < szX; ++x)  cim[y + padding][x + padding] *= g[y][x];

        // convolve
        //ft.set (tmp); FIXME
        ft *= otf;
        ft.directInverse (tmp);

        //
        for (int y = 0; y < szY; ++y) for (int x = 0; x < szX; ++x)  cim[y + padding][x + padding] *= g[y][x];

        // Compute ChiSq
        cs = 0.0;

        for (int y = 0; y < szY; ++y) {
            for (int x = 0; x < szX; ++x) {
                U bla = im[y][x] - cim[y + padding][x + padding];
                cs += norm (pim[y + padding][x + padding] - bla);   // norm() works for complex numbers too.
                pim[y + padding][x + padding] = bla;
            }
        }

        if (verbosity > 0) {
            cout << "cdescatter:  iteration #" << i << "  ChiSq = " << (cs / nPixels) << endl;
        }
    }

    delPointers (g);
    delPointers (im);
    delPointers (pim);
    delPointers (cim);

    subImg.template copyTo<T> (out);

}


IDL_VPTR cdescatter (int argc, IDL_VPTR* argv, char* argk) {

    if (argc < 3) {
        cout << "cdescatter: needs at least 3 arguments (image, gain and psf). " << endl;
        convolve_help();
        return IDL_GettmpInt (-1);
    }

    IDL_VPTR imageVar  = argv[0];
    IDL_VPTR gainVar  = argv[1];
    IDL_VPTR psfVar  = argv[2];

    KW_RESULT kw;
    kw.help = 0;
    kw.verbose = 0;
    kw.normalize = 0;
    kw.padding = 0;
    kw.niter = 1;
    kw.nthreads = 1;
    (void) IDL_KWProcessByOffset (argc, argv, argk, kw_pars, (IDL_VPTR*) 0, 255, &kw);

    if (kw.help) {
        convolve_help();
        return IDL_GettmpInt (-1);
    }

    IDL_ENSURE_SIMPLE (imageVar);
    IDL_ENSURE_ARRAY (imageVar);

    IDL_ENSURE_SIMPLE (gainVar);
    IDL_ENSURE_ARRAY (gainVar);

    IDL_ENSURE_SIMPLE (psfVar);
    IDL_ENSURE_ARRAY (psfVar);

    int nDims = imageVar->value.arr->n_dim;

    if (nDims != 2 || nDims != psfVar->value.arr->n_dim || nDims != gainVar->value.arr->n_dim) {
        cout << "cdescatter: img/gain/psf must all be 2-dimensional arrays." << endl;
        return IDL_GettmpInt (-1);
    }

    if (imageVar->type != psfVar->type || imageVar->type != gainVar->type) {
        cout << "cdescatter: img/gain/psf must all of the same type." << endl;
        return IDL_GettmpInt (-1);
    }

    int totalSize = 1;
    const IDL_ARRAY_DIM& dims = imageVar->value.arr->dim;

    for (int i = 0; i < nDims; ++i) {
        if (dims[i] != psfVar->value.arr->dim[i] || dims[i] != gainVar->value.arr->dim[i]) {
            cout << "cdescatter: img/gain/psf dimensions do not match: "
                 << printArray (dims, nDims, "img")
                 << printArray (gainVar->value.arr->dim, nDims, "  gain")
                 << printArray (psfVar->value.arr->dim, nDims, "  psf") << endl;
            return IDL_GettmpInt (-1);
        }

        totalSize *= dims[i];
    }

    IDL_VPTR tmp;

    if (totalSize) {
        totalSize *= var_sizes[imageVar->type];
        char* data = new char[totalSize];
        tmp = IDL_ImportArray (nDims, imageVar->value.arr->dim, imageVar->type, (UCHAR*) data, redux::util::castAndDelete<char>, NULL);
        UCHAR* imgData = imageVar->value.arr->data;
        UCHAR* gainData = gainVar->value.arr->data;
        UCHAR* psfData = psfVar->value.arr->data;
        UCHAR* tmpData = tmp->value.arr->data;

        switch (imageVar->type) {
            case IDL_TYP_FLOAT:
                descatter (dims[1], dims[0], reinterpret_cast<float*> (imgData),
                                             reinterpret_cast<float*> (gainData),
                                             reinterpret_cast<float*> (psfData),
                                             reinterpret_cast<float*> (tmpData),
                           kw.nthreads, kw.niter, kw.padding, kw.verbose);
                break;

            case IDL_TYP_DOUBLE:
                descatter (dims[1], dims[0], reinterpret_cast<double*> (imgData),
                                             reinterpret_cast<double*> (gainData),
                                             reinterpret_cast<double*> (psfData),
                                             reinterpret_cast<double*> (tmpData),
                           kw.nthreads, kw.niter, kw.padding, kw.verbose);
                break;

            case IDL_TYP_DCOMPLEX:
                descatter<complex_t, complex_t> (dims[1], dims[0], reinterpret_cast<complex_t*> (imgData),
                                                                   reinterpret_cast<complex_t*> (gainData),
                                                                   reinterpret_cast<complex_t*> (psfData),
                                                                   reinterpret_cast<complex_t*> (tmpData),
                                                 kw.nthreads, kw.niter, kw.padding, kw.verbose);
                break;

            default:
                cout << "cdescatter: not implemented for IDL-type: " << (int) imageVar->type << endl;
        }
    } else {
        cout << "cdescatter: dimension has size 0 !!" << endl;
        tmp = IDL_GettmpInt (-1);
    }

    IDL_KW_FREE;
    return tmp;

}


template <typename T>
void convolve (int szY, int szX, T* img, T* psf, T* out, int nthreads, int normalize, int verbose) {

    if (verbose > 1) {
        cout << "convolve<" << typeid (T).name() << ">    size = (" << szY << "," << szX << ")  nthreads = " << nthreads << " normalize = " << normalize << endl;
    }

    Array<T> psfA (psf, szY, szX);
    FourierTransform ft (psfA, 0, std::max (nthreads, 1));

    if (normalize) {
        ArrayStats stats;
        stats.getMinMaxMean (psfA);
        ft /= stats.sum;
    }

    Array<T> tmp = ft.convolve (Array<T> (img, szY, szX));

    memcpy (out, tmp.ptr(), szX * szY * sizeof (T));

}


IDL_VPTR cconvolve (int argc, IDL_VPTR* argv, char* argk) {

    if (argc < 2) {
        cout << "cconvolve: needs at least 2 arguments (image and psf). " << endl;
        convolve_help();
        return IDL_GettmpInt (-1);
    }

    IDL_VPTR imageVar  = argv[0];
    IDL_VPTR psfVar  = argv[1];

    KW_RESULT kw;
    kw.help = 0;
    kw.verbose = 0;
    kw.normalize = 0;
    kw.nthreads = 1;
    (void) IDL_KWProcessByOffset (argc, argv, argk, kw_pars, (IDL_VPTR*) 0, 255, &kw);

    if (kw.help) {
        convolve_help();
        return IDL_GettmpInt (-1);
    }

    IDL_ENSURE_SIMPLE (imageVar);
    IDL_ENSURE_ARRAY (imageVar);

    IDL_ENSURE_SIMPLE (psfVar);
    IDL_ENSURE_ARRAY (psfVar);

    int nDims = imageVar->value.arr->n_dim;

    if (nDims != 2 || nDims != psfVar->value.arr->n_dim) {
        cout << "cconvolve: img/psf must both be 2-dimensional arrays." << endl;
        return IDL_GettmpInt (-1);
    }

    if (imageVar->type != psfVar->type) {
        cout << "cconvolve: img/psf must both of the same type." << endl;
        return IDL_GettmpInt (-1);
    }

    int totalSize = 1;
    const IDL_ARRAY_DIM& dims = imageVar->value.arr->dim;

    for (int i = 0; i < nDims; ++i) {
        if (dims[i] != psfVar->value.arr->dim[i]) {
            cout << "cconvolve: img/psf dimensions do not match: "
                 << printArray (dims, nDims, "img")
                 << printArray (psfVar->value.arr->dim, nDims, "  psf") << endl;
            return IDL_GettmpInt (-1);
        }

        totalSize *= dims[i];
    }

    IDL_VPTR tmp;

    if (totalSize) {
        totalSize *= var_sizes[imageVar->type];
        char* data = new char[totalSize];
        tmp = IDL_ImportArray (nDims, imageVar->value.arr->dim, imageVar->type, (UCHAR*) data, redux::util::castAndDelete<char>, NULL);
        UCHAR* imgData = imageVar->value.arr->data;
        UCHAR* psfData = psfVar->value.arr->data;
        UCHAR* tmpData = tmp->value.arr->data;

        switch (imageVar->type) {
            case IDL_TYP_FLOAT:
                convolve (dims[1], dims[0], reinterpret_cast<float*> (imgData), reinterpret_cast<float*> (psfData),
                          reinterpret_cast<float*> (tmpData), kw.nthreads, kw.normalize, kw.verbose);
                break;

            case IDL_TYP_DOUBLE:
                convolve (dims[1], dims[0], reinterpret_cast<double*> (imgData), reinterpret_cast<double*> (psfData),
                          reinterpret_cast<double*> (tmpData), kw.nthreads, kw.normalize, kw.verbose);
                break;

            case IDL_TYP_DCOMPLEX:
                convolve (dims[1], dims[0], reinterpret_cast<complex_t*> (imgData), reinterpret_cast<complex_t*> (psfData),
                          reinterpret_cast<complex_t*> (tmpData), kw.nthreads, kw.normalize, kw.verbose);
                break;

            default:
                cout << "cconvolve: not implemented for IDL-type: " << (int) imageVar->type << endl;
        }
    } else {
        cout << "cconvolve: dimension has size 0 !!" << endl;
        tmp = IDL_GettmpInt (-1);
    }

    IDL_KW_FREE;
    return tmp;

}


template <typename T, typename U>
double *bezier2 (int n, T *xin, U *yin, int np, double *xp) {
    //
    // Reverse arrays?
    //
    T *x;
    U *y;
    short allocated = 0;

    if ( (xin[1] - xin[0]) > 0.0) {
        x = xin;
        y = yin;
    } else {
        x = new T [n];
        y = new U [n];
        for (int k = 0; k < n; k++) {
            x[k] = xin[n - k - 1];
            y[k] = yin[n - k - 1];
        }
        allocated = 1;
    }


    //
    // Centered derivatives
    //
    double *yp = new double [n];
    double *res = new double [np];
    central_derivative (n, x, y, yp);

    //
    // Loop intervals
    //
    for (int k = 0; k < n - 1; k++) {
        double dx = x[k + 1] - x[k];
        double cntrl = 0.5 * (y[k] + 0.5 * dx * yp[k] + y[k + 1] - dx * 0.5 * yp[k + 1]);

        for (int j = 0; j < np; j++) {
            if ( (xp[j] >= x[k]) && (xp[j] < x[k + 1])) {
                double u = (xp[j] - x[k]) / dx;
                double u1 = 1.0 - u;
                res[j] =  y[k] * u1 * u1 + y[k + 1] * u * u + 2.0 * cntrl * u * u1;
            }
        }
    }

//
    // Out-of bound values, linear extrapolation
    //
    double aa0 = (y[1] - y[0]) / (x[1] - x[0]);
    double bb0 = y[0] - aa0 * x[0];
    double aa1 = (y[n - 1] - y[n - 2]) / (x[n - 1] - y[n - 2]);
    double bb1 = y[n - 1] - aa1 * x[n - 1];

    for (int k = 0; k < np; k++) {
        //    if(xp[k] < x[0]) res[k] = y[0];
        //    if(xp[k] >= x[n-1]) res[k] = y[n-1];
        if (xp[k] < x[0]) res[k] = xp[k] * aa0 + bb0;

        if (xp[k] >= x[n - 1]) res[k] = xp[k] * aa1 + bb1;
    }

    if (allocated) {
        delete [] x;
        delete [] y;
    }

    delete [] yp;
    return res;
}

template <typename T, typename U>
double *bezier3 (int n, T *xin, U *yin, int np, double *xp) {
    //
    // Reverse arrays?
    //
    T *x;
    U *y;
    short allocated = 0;

    if ( (xin[1] - xin[0]) > 0.0) {
        x = xin;
        y = yin;
    } else {
        x = new T [n];
        y = new U [n];
        for (int k = 0; k < n; k++) {
            x[k] = xin[n - k - 1];
            y[k] = yin[n - k - 1];
        }
        allocated = 1;
    }


    //
    // Centered derivatives
    //
    double *yp = new double [n];
    double *res = new double [np];
    central_derivative (n, x, y, yp);

    //
    // Loop intervals
    //
    for (int k = 0; k < n - 1; k++) {
        double dx = x[k + 1] - x[k];
        double cntrl1 = y[k] +  dx * yp[k] / 3.0;
        double cntrl2 = y[k + 1] - dx * yp[k + 1] / 3.0;

        for (int j = 0; j < np; j++) {
            if ( (xp[j] >= x[k]) && (xp[j] < x[k + 1])) {
                double u = (xp[j] - x[k]) / dx;
                double u1 = 1.0 - u;
                res[j] =  y[k] * u1 * u1 * u1 + y[k + 1] * u * u * u + 3.0 * cntrl1 * u * u1 * u1 + 3.0 * cntrl2 * u * u * u1;
            }
        }
    }

    //
    // Out-of bound values, linear extrapolation
    //
    double aa0 = (y[1] - y[0]) / (x[1] - x[0]);
    double bb0 = y[0] - aa0 * x[0];
    double aa1 = (y[n - 1] - y[n - 2]) / (x[n - 1] - y[n - 2]);
    double bb1 = y[n - 1] - aa1 * x[n - 1];

    for (int k = 0; k < np; k++) {
        //    if(xp[k] < x[0]) res[k] = y[0];
        //    if(xp[k] >= x[n-1]) res[k] = y[n-1];
        if (xp[k] < x[0]) res[k] = xp[k] * aa0 + bb0;

        if (xp[k] >= x[n - 1]) res[k] = xp[k] * aa1 + bb1;
    }

    if (allocated) {
        delete [] x;
        delete [] y;
    }

    delete [] yp;
    return res;
}


IDL_VPTR cbezier2 (int argc, IDL_VPTR* argv, char* argk) {
    
    if (argc < 3) {
        cout << "cbezier2: needs 3 arguments. " << endl;
        return IDL_GettmpInt (-1);
    }
    
    IDL_VPTR x = argv[0];
    IDL_VPTR y = argv[1];
    IDL_VPTR xp_in = argv[2];

    KW_RESULT kw;
    kw.help = 0;
    kw.verbose = 0;
    kw.normalize = 0;
    kw.padding = 0;
    kw.niter = 1;
    kw.nthreads = 1;
    (void) IDL_KWProcessByOffset (argc, argv, argk, kw_pars, (IDL_VPTR*) 0, 255, &kw);

    if (kw.help) {
        cout << "cbezier2: no help available yet..." << endl;
        return IDL_GettmpInt (-1);
    }
    
    IDL_ENSURE_SIMPLE (x);
    IDL_ENSURE_ARRAY (x);

    IDL_ENSURE_SIMPLE (y);
    IDL_ENSURE_ARRAY (y);

    IDL_ENSURE_SIMPLE (xp_in);
    IDL_ENSURE_ARRAY (xp_in);

    int nDims = x->value.arr->n_dim;

    if (nDims != 1 || nDims != xp_in->value.arr->n_dim || nDims != y->value.arr->n_dim) {
        cout << "cbezier2: input must be 1-dimensional arrays." << endl;
        return IDL_GettmpInt (-1);
    }
    
    if (x->type < 4 || y->type < 4 || xp_in->type < 4 ||
        x->type > 5 || y->type > 5 || xp_in->type > 5) {
        cout << "cbezier2: ERROR, input data must be double or float!" << endl;
        return IDL_GettmpInt (-1);
    }

    int n =      x->value.s.arr->dim[0];
    int np = xp_in->value.s.arr->dim[0];
    double *res;
    double *xp;
    unique_ptr<double> tmpXP;
    
    if (xp_in->type == IDL_TYP_FLOAT) {
        float* ptr = reinterpret_cast<float*>(xp_in->value.arr->data);
        tmpXP.reset(new double[np]);
        xp = tmpXP.get();
        for (int k = 0; k < np; k++) xp[k] = ptr[k];
    } else {
        xp = reinterpret_cast<double*>(xp_in->value.arr->data);
    }
    
    if (x->type == IDL_TYP_FLOAT) {
        float* xPtr = reinterpret_cast<float*>(x->value.arr->data);
        if (x->type == IDL_TYP_FLOAT) {
            float* yPtr = reinterpret_cast<float*>(y->value.arr->data);
            res = bezier2(n, xPtr, yPtr, np, xp);
        } else {
            double* yPtr = reinterpret_cast<double*>(y->value.arr->data);
            res = bezier2(n, xPtr, yPtr, np, xp);
        }
    } else {
        double* xPtr = reinterpret_cast<double*>(x->value.arr->data);
        if (x->type == IDL_TYP_FLOAT) {
            float* yPtr = reinterpret_cast<float*>(y->value.arr->data);
            res = bezier2(n, xPtr, yPtr, np, xp);
        } else {
            double* yPtr = reinterpret_cast<double*>(y->value.arr->data);
            res = bezier2(n, xPtr, yPtr, np, xp);
        }
    }
    
    IDL_MEMINT dims[] = {np};
    return IDL_ImportArray (1, dims, IDL_TYP_DOUBLE, (UCHAR*)res, redux::util::castAndDelete<double>, 0);

}


IDL_VPTR cbezier3(int argc, IDL_VPTR* argv, char* argk) {
    
    if (argc < 3) {
        cout << "cbezier2: needs 3 arguments. " << endl;
        return IDL_GettmpInt (-1);
    }
    
    IDL_VPTR x = argv[0];
    IDL_VPTR y = argv[1];
    IDL_VPTR xp_in = argv[2];

    KW_RESULT kw;
    kw.help = 0;
    kw.verbose = 0;
    kw.normalize = 0;
    kw.padding = 0;
    kw.niter = 1;
    kw.nthreads = 1;
    (void) IDL_KWProcessByOffset (argc, argv, argk, kw_pars, (IDL_VPTR*) 0, 255, &kw);

    if (kw.help) {
        cout << "cbezier2: no help available yet..." << endl;
        return IDL_GettmpInt (-1);
    }
    
    IDL_ENSURE_SIMPLE (x);
    IDL_ENSURE_ARRAY (x);

    IDL_ENSURE_SIMPLE (y);
    IDL_ENSURE_ARRAY (y);

    IDL_ENSURE_SIMPLE (xp_in);
    IDL_ENSURE_ARRAY (xp_in);

    int nDims = x->value.arr->n_dim;

    if (nDims != 1 || nDims != xp_in->value.arr->n_dim || nDims != y->value.arr->n_dim) {
        cout << "cbezier2: input must be 1-dimensional arrays." << endl;
        return IDL_GettmpInt (-1);
    }
    
    if (x->type < 4 || y->type < 4 || xp_in->type < 4 ||
        x->type > 5 || y->type > 5 || xp_in->type > 5) {
        cout << "cbezier2: ERROR, input data must be double or float!" << endl;
        return IDL_GettmpInt (-1);
    }

    int n =      x->value.s.arr->dim[0];
    int np = xp_in->value.s.arr->dim[0];
    double *res;
    double *xp;
    unique_ptr<double> tmpXP;
    
    if (xp_in->type == IDL_TYP_FLOAT) {
        float* ptr = reinterpret_cast<float*>(xp_in->value.arr->data);
        tmpXP.reset(new double[np]);
        xp = tmpXP.get();
        for (int k = 0; k < np; k++) xp[k] = ptr[k];
    } else {
        xp = reinterpret_cast<double*>(xp_in->value.arr->data);
    }
    
    if (x->type == IDL_TYP_FLOAT) {
        float* xPtr = reinterpret_cast<float*>(x->value.arr->data);
        if (x->type == IDL_TYP_FLOAT) {
            float* yPtr = reinterpret_cast<float*>(y->value.arr->data);
            res = bezier3(n, xPtr, yPtr, np, xp);
        } else {
            double* yPtr = reinterpret_cast<double*>(y->value.arr->data);
            res = bezier3(n, xPtr, yPtr, np, xp);
        }
    } else {
        double* xPtr = reinterpret_cast<double*>(x->value.arr->data);
        if (x->type == IDL_TYP_FLOAT) {
            float* yPtr = reinterpret_cast<float*>(y->value.arr->data);
            res = bezier3(n, xPtr, yPtr, np, xp);
        } else {
            double* yPtr = reinterpret_cast<double*>(y->value.arr->data);
            res = bezier3(n, xPtr, yPtr, np, xp);
        }
    }
    
    IDL_MEMINT dims[] = {np};
    return IDL_ImportArray (1, dims, IDL_TYP_DOUBLE, (UCHAR*)res, redux::util::castAndDelete<double>, 0);
    
}


typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_INT sort;
    IDL_INT unique;
    IDL_INT verbose;
} STR_TO_INT_KW;


// NOTE:  The keywords MUST be listed in alphabetical order !!
static IDL_KW_PAR str_to_int_kw_pars[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "SORT",      IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*)IDL_KW_OFFSETOF2(STR_TO_INT_KW,sort) },
    { (char*) "UNIQUE",    IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*)IDL_KW_OFFSETOF2(STR_TO_INT_KW,unique) },
    { (char*) "VERBOSE",   IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*)IDL_KW_OFFSETOF2(STR_TO_INT_KW,verbose) },
    { NULL }
};

IDL_VPTR string_to_uints(int argc, IDL_VPTR* argv, char* argk) {
 
    if (argc < 1) {
        return IDL_GettmpInt(0);
    }
    
    IDL_VPTR ret;
    vector<uint32_t> uints;
    IDL_VPTR str = argv[0];
    IDL_ENSURE_SIMPLE(str);
    IDL_ENSURE_STRING(str);
    
    STR_TO_INT_KW kw;
    kw.sort = 0;
    kw.verbose = 0;
    kw.unique = 0;
    (void) IDL_KWProcessByOffset (argc, argv, argk, str_to_int_kw_pars, (IDL_VPTR*) 0, 255, &kw);
    
 
    try {
        uints = redux::util::stringToUInts<uint32_t>( IDL_VarGetString(str) );
    } catch( exception & ) { }

    if( kw.sort || kw.unique ) {
        std::sort( uints.begin(), uints.end() );
    }
    if( kw.unique ) {
        auto it = std::unique( uints.begin(), uints.end() );
        uints.resize( std::distance( uints.begin(), it ) );
    }

    if ( uints.empty() ) {
        return IDL_GettmpInt(0);
    }
    
    IDL_MEMINT dims[] = { (int)uints.size() };
    IDL_LONG* tmpData = (IDL_LONG*)IDL_MakeTempArray(IDL_TYP_LONG,1,dims,IDL_ARR_INI_NOP,&ret);
    std::copy( uints.begin(), uints.end(), tmpData );

    return ret;

}


IDL_VPTR uints_to_string(int argc, IDL_VPTR* argv, char* argk) {
    
    string ret;
    if (argc < 1) {
        return IDL_StrToSTRING( (char*)ret.c_str() );
    }
    
    IDL_VPTR uints = argv[0];
    IDL_ENSURE_SIMPLE(uints);
    IDL_ENSURE_ARRAY(uints);
    
    STR_TO_INT_KW kw;
    kw.sort = 0;
    kw.verbose = 0;
    kw.unique = 0;
    (void) IDL_KWProcessByOffset (argc, argv, argk, str_to_int_kw_pars, (IDL_VPTR*) 0, 255, &kw);
    
    vector<uint64_t> tmp;
    if( uints->type == IDL_TYP_BYTE ) {
        auto beg = uints->value.arr->data;
        std::copy( beg, beg + uints->value.arr->n_elts , back_inserter(tmp) );
    } else if (uints->type == IDL_TYP_INT ) {
        auto beg = reinterpret_cast<int16_t*>(uints->value.arr->data);
        std::copy( beg, beg + uints->value.arr->n_elts , back_inserter(tmp) );
    } else if (uints->type == IDL_TYP_LONG ) {
        auto beg = reinterpret_cast<int32_t*>(uints->value.arr->data);
        std::copy( beg, beg + uints->value.arr->n_elts , back_inserter(tmp) );
    } else  {
        cout << "uints_to_string: input array must be of type BYTE, INT or LONG." << endl;
        return IDL_GettmpInt (0);
    }
    
    if( kw.sort || kw.unique ) {
        std::sort( tmp.begin(), tmp.end() );
    }
    if( kw.unique ) {
        auto it = std::unique( tmp.begin(), tmp.end() );
        tmp.resize( std::distance( tmp.begin(), it ) );
    }

    try {
        ret = redux::util::uIntsToString( tmp );
    } catch( exception & ) { }
    
    return IDL_StrToSTRING( (char*)ret.c_str() );

}


typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_VPTR arglist;
    IDL_STRING split_chars;
} MT_KW;


// NOTE:  The keywords MUST be listed in alphabetical order !!
static IDL_KW_PAR mt_kw_pars[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "ARG_LIST",   IDL_TYP_UNDEF,  1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(MT_KW,arglist) },
    { (char*) "SPLIT_CHARS", IDL_TYP_STRING, 1, 0, 0, (char*)IDL_KW_OFFSETOF2(MT_KW,split_chars) },
    { NULL }
};

IDL_VPTR make_template(int argc, IDL_VPTR* argv, char* argk) {
    
    string ret;
    if (argc < 1) {
        return IDL_StrToSTRING( (char*)ret.c_str() );
    }
    
    IDL_VPTR strlist = argv[0];
    IDL_ENSURE_STRING( strlist )
    IDL_ENSURE_SIMPLE( strlist );
    
    if ( !(strlist->flags & IDL_V_ARR) ) {
        return strlist;
    }
    
    IDL_ARRAY* strarr = strlist->value.arr;
    IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>(strlist->value.arr->data);

    if ( strarr->n_elts == 1) {
        return IDL_StrToSTRING( strptr->s );
    }
    
    MT_KW kw;
    kw.split_chars = IDL_STATIC_STRING(".");
    
    (void) IDL_KWProcessByOffset (argc, argv, argk, mt_kw_pars, (IDL_VPTR*)0, 255, &kw);
    
    string split_chars = kw.split_chars.s;
    string tmp = strptr->s;    // first item
    
    std::vector<std::string> segments;
    boost::split(segments, tmp, boost::is_any_of(split_chars));
    
    size_t n_segments = segments.size();
    vector< set<string> > seg_list(n_segments);

    for( int i=0; i<strarr->n_elts; ++i ) {
        tmp = strptr[i].s;
        boost::split(segments, tmp, boost::is_any_of(split_chars));
        if( segments.size() != n_segments ) {
            cout << "Item did not match template: " << tmp << endl;
            continue;
        }
        for( size_t j=0; j<n_segments; ++j ) {
            seg_list[j].insert(segments[j]);
        }
    }
    
    int arg_cnt(0);
    for( size_t j=0; j<n_segments; ++j ) {
        if(seg_list[j].size() > 1 ) {
            ret += "%"+to_string(++arg_cnt);
        } else {
            ret += *(seg_list[j].begin());
        }
        if( j < n_segments-1 ) ret += split_chars[0];
    }
    
    if( kw.arglist && arg_cnt ) {
        IDL_VPTR arglist;
        IDL_MEMINT dims[] = { arg_cnt };    // x/y, im#, nMatches 
        IDL_MakeTempArray( IDL_TYP_STRING, 1, dims, IDL_ARR_INI_ZERO, &arglist );
        strptr = reinterpret_cast<IDL_STRING*>(arglist->value.arr->data);
        for( size_t j=0; j<n_segments; ++j ) {
            if(seg_list[j].size() > 1 ) {
                tmp = boost::join( seg_list[j], " " );
                IDL_StrStore(strptr++,const_cast<char*>(tmp.c_str()));
            }
        }
        IDL_VarCopy( arglist, kw.arglist );
    }
    
    return IDL_StrToSTRING( (char*)ret.c_str() );

}


void fft_reorder(int argc, IDL_VPTR* argv) {

    if( argc < 1 ) {
        cout << "fft_reorder: needs a single 2D array as input." << endl;
        return;
    }
    
    IDL_VPTR array = argv[0];
    IDL_ENSURE_SIMPLE(array);
    IDL_ENSURE_ARRAY(array);
    
    if ( array->value.arr->n_dim != 2 ) {
        cout << "fft_reorder: needs a single 2D array as input." << endl;
        return;
    }

    size_t nX = array->value.arr->dim[0];
    size_t nY = array->value.arr->dim[1];
    UCHAR* data = array->value.arr->data;
    try {
        switch(array->type) {
            case IDL_TYP_INT: FourierTransform::reorder( reinterpret_cast<int16_t*>(data), nY, nX ); break;
            case IDL_TYP_LONG: FourierTransform::reorder( reinterpret_cast<int32_t*>(data), nY, nX ); break;
            case IDL_TYP_FLOAT: FourierTransform::reorder( reinterpret_cast<float*>(data), nY, nX ); break;
            case IDL_TYP_DOUBLE: FourierTransform::reorder( reinterpret_cast<double*>(data), nY, nX ); break;
            case IDL_TYP_COMPLEX: FourierTransform::reorder( reinterpret_cast<std::complex<float>*>(data), nY, nX ); break;
            case IDL_TYP_DCOMPLEX: FourierTransform::reorder( reinterpret_cast<std::complex<double>*>(data), nY, nX ); break;
            default: FourierTransform::reorder( data, nY, nX);
        }
    } catch ( const exception& e ) {   // catch-all to avoid killing the idl session.
        cout << "fft_reorder: Error during call to FourierTransform::reorder() : " << e.what() << endl;
    }

}


extern "C" {

    int IDL_Load (void) {

        static IDL_SYSFUN_DEF2 function_addr[] = {
            { { (IDL_VPTR (*) ()) redux::structInfo}, (char*) "STRUCT_INFO", 0, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { { (IDL_VPTR (*) ()) cconvolve}, (char*) "CCONVOLVE", 2, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { { (IDL_VPTR (*) ()) cdescatter}, (char*) "CDESCATTER", 3, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { { (IDL_VPTR (*) ()) cbezier2}, (char*) "CBEZIER2", 3, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { { (IDL_VPTR (*) ()) cbezier3}, (char*) "CBEZIER3", 3, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { { (IDL_VPTR (*) ()) string_to_uints}, (char*) "STRING_TO_UINTS", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { { (IDL_VPTR (*) ()) uints_to_string}, (char*) "UINTS_TO_STRING", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { { (IDL_VPTR (*) ()) make_template}, (char*) "MAKE_TEMPLATE", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { { (IDL_VPTR (*) ()) redux::img_align}, (char*) "IMG_ALIGN", 2, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { { (IDL_VPTR (*) ()) redux::img_project}, (char*) "IMG_PROJECT", 2, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }
        };

        static IDL_SYSFUN_DEF2 procedure_addr[] = {
            { { (IDL_SYSRTN_GENERIC) redux::printStruct}, (char*) "STRUCT_INFO", 0, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { { (IDL_SYSRTN_GENERIC) fft_reorder}, (char*) "FFT_REORDER", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
        };

        return IDL_SysRtnAdd (function_addr, TRUE, IDL_CARRAY_ELTS (function_addr)) &&
               IDL_SysRtnAdd (procedure_addr, FALSE, IDL_CARRAY_ELTS (procedure_addr));

    }

}
