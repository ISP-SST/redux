#include "rdx.hpp"

#include "idlutil.hpp"
#include "imgtools.hpp"

#include "redux/util/arrayutil.hpp"
#include "redux/util/cache.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stopwatch.hpp"
#include "redux/file/fileio.hpp"
#include "redux/version.hpp"

#include <iomanip>
#include <iostream>
#include <typeinfo>

using namespace redux::file;
using namespace redux::util;
using namespace redux;
using namespace std;


namespace {

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
    double aa1 = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
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
    double aa1 = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
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
    
    KW_RESULT kw;
    kw.help = 0;
    kw.verbose = 0;
    kw.normalize = 0;
    kw.padding = 0;
    kw.niter = 1;
    kw.nthreads = 1;
    int nPlainArgs = IDL_KWProcessByOffset (argc, argv, argk, kw_pars, (IDL_VPTR*) 0, 255, &kw);

    if (nPlainArgs < 3) {
        cout << "cbezier2: needs 3 arguments. " << endl;
        return IDL_GettmpInt (-1);
    }
    
    IDL_VPTR x = argv[0];
    IDL_VPTR y = argv[1];
    IDL_VPTR xp_in = argv[2];

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
    
    KW_RESULT kw;
    kw.help = 0;
    kw.verbose = 0;
    kw.normalize = 0;
    kw.padding = 0;
    kw.niter = 1;
    kw.nthreads = 1;
    int nPlainArgs = IDL_KWProcessByOffset (argc, argv, argk, kw_pars, (IDL_VPTR*) 0, 255, &kw);

    if (nPlainArgs < 3) {
        cout << "cbezier2: needs 3 arguments. " << endl;
        return IDL_GettmpInt (-1);
    }
    
    IDL_VPTR x = argv[0];
    IDL_VPTR y = argv[1];
    IDL_VPTR xp_in = argv[2];

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


StopWatch& getTimer(void) {
    static StopWatch sw;
    return sw;
}

typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_INT help;
    IDL_VPTR seconds;
    IDL_INT cpu;
} KW_TIMER;

// NOTE:  The keywords MUST be listed in alphabetical order !!
static IDL_KW_PAR kw_timer_pars[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "CPU",         IDL_TYP_INT,   1,   IDL_KW_ZERO,   0,   (char*)IDL_KW_OFFSETOF2(KW_TIMER,cpu) },
    { (char*) "SECONDS",   IDL_TYP_UNDEF,   1,   IDL_KW_OUT|IDL_KW_ZERO,   0,   (char*)IDL_KW_OFFSETOF2(KW_TIMER,seconds) },
    { NULL }
};


string timer_info( int lvl ) {
    
    string ret = "RDX_TIC/RDX_TOC";
    if( lvl > 0 ) {
        
        ret += ((lvl > 1)?"\n":"    ");          // newline if lvl>1
        ret += "   Syntax:   rdx_tic & wait,3 & rdx_toc\n";
    
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords (for TOC):\n"
                    "      SECONDS             (output) Get number of elapsed seconds instead of just printing.\n"
                    "      CPU                 Get CPU-time (user+system) instead of wall-time. (only used together with SECONDS)\n";
        }
    
    } else ret += "\n";

    return ret;
    
}


void timer_start(int argc, IDL_VPTR* argv) {
    
    StopWatch& sw = getTimer();
    sw.start();
    
}

void timer_elapsed( int argc, IDL_VPTR* argv, char* argk ) {
    
    StopWatch& sw = getTimer();
    KW_TIMER kw;
    kw.help = 0;
    (void) IDL_KWProcessByOffset(argc, argv, argk, kw_timer_pars, (IDL_VPTR*)0, 255, &kw);
    
    if( kw.help ) {
        cout << timer_info(2) << endl;
        return;
    }

    if( kw.seconds ) {
        IDL_ALLTYPES tmp;
        tmp.f = sw.getSeconds();
        cout << "Secs: " << tmp.f << endl;;
        IDL_StoreScalar( kw.seconds, IDL_TYP_FLOAT, &tmp );
    } else {
        cout << sw.print() << endl;
    }
    
}


typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_INT help;
    IDL_INT momfbd;
} KW_SEGMENT;

// NOTE:  The keywords MUST be listed in alphabetical order !!
static IDL_KW_PAR kw_segment_pars[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "HELP",     IDL_TYP_INT,   1,   IDL_KW_ZERO,   0,   (char*)IDL_KW_OFFSETOF2(KW_SEGMENT,help) },
    { (char*) "MOMFBD",   IDL_TYP_INT,   1,   IDL_KW_ZERO,   0,   (char*)IDL_KW_OFFSETOF2(KW_SEGMENT,momfbd) },
    { NULL }
};

string segment_info( int lvl ) {
    
    string ret = "RDX_SEGMENT";
    if( lvl > 0 ) {
        
        ret += ((lvl > 1)?"\n":"        ");          // newline if lvl>1
        ret += "   Syntax:   pos = rdx_segment( first, last, size, min-overlap) \n";
        
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      MOMFBD              Interpret arguments as pixels and do segmentation identical to the MOMFBD code.\n";
        }
    
    } else ret += "\n";

    return ret;
    
}


IDL_VPTR rdx_segment( int argc, IDL_VPTR* argv, char* argk ) {
    
    KW_SEGMENT kw;
    int nPlainArgs = IDL_KWProcessByOffset(argc, argv, argk, kw_segment_pars, (IDL_VPTR*)0, 255, &kw);
    
    if( nPlainArgs < 3 ) {
        cout << segment_info(2) << endl;
        return IDL_GettmpInt(-1);
    }
    
    IDL_LONG first = IDL_LongScalar( argv[0] );
    IDL_LONG last = IDL_LongScalar( argv[1] );
    IDL_LONG sz = IDL_LongScalar( argv[2] );

    IDL_LONG minOverLap(0);
    if( nPlainArgs > 3 ) {
        minOverLap = IDL_LongScalar( argv[3] );
    }
    
    if( kw.help ) {
        cout << segment_info(2) << endl;
        return IDL_GettmpInt(-1);
    }
    
    if( kw.momfbd ) {
        first += sz/2;
        last -= sz/2;
        minOverLap = 16 + sz/4; 
    }
    
    std::vector<IDL_LONG> resV = segment<IDL_LONG>( first, last, sz, minOverLap );
    if ( resV.empty() ) {
        return IDL_GettmpInt(-1);
    }
    
    IDL_MEMINT dims[] = { (IDL_MEMINT)resV.size() };
    unique_ptr<IDL_LONG> res( new IDL_LONG[ dims[0] ] );
    memcpy( res.get(),resV.data(),dims[0]*sizeof(IDL_LONG) );
   
    return IDL_ImportArray( 1, dims, IDL_TYP_LONG, (UCHAR*)res.release(), redux::util::castAndDelete<IDL_LONG>, 0 );
    
}


string hasopencv_info( int lvl ) {
    
    string ret = "RDX_HASOPENCV";
    if( lvl > 0 ) {
        ret += "         Syntax:   has_cv = rdx_hasopencv()\n";
    } else ret += "\n";

    return ret;
    
}

IDL_VPTR rdx_hasopencv( int argc, IDL_VPTR* argv, char* argk ) {
    
#ifdef REDUX_WITH_OPENCV
    return IDL_GettmpInt(1);
#else
    return IDL_GettmpInt(0);
#endif
    
}


string filetype_info( int lvl ) {
    
    string ret = "RDX_FILETYPE";
    if( lvl > 0 ) {
        ret += "          Syntax:   fmt = rdx_filetype()\n";
    } else ret += "\n";

    return ret;
    
}

IDL_VPTR rdx_filetype( int argc, IDL_VPTR* argv, char* argk ) {
    
    if( argc < 1 ) {
        cout << filetype_info(2) << endl;
        return IDL_StrToSTRING( (char*)"fla" );
    }
    
    IDL_VPTR filenames = argv[0];
    IDL_ENSURE_SIMPLE( filenames );
    
    if( filenames->type != IDL_TYP_STRING ) {
        return IDL_StrToSTRING( (char*)"fly" );
    }

    if ( !(filenames->flags & IDL_V_ARR) ) {
        bfs::path fn( string(filenames->value.str.s) );
        Format fmt = FMT_NONE;
        try {
            if( bfs::is_regular_file(fn) ) {
                try {
                    fmt = redux::file::readFmt( fn.string() );
                } catch( const std::ios_base::failure&) {
                    // silently ignore read errors.
                    fmt = redux::file::guessFmt( fn.string() );
                }
            } else {
                fmt = redux::file::guessFmt( fn.string() );
            }
        } catch( const exception& e) {
            cout << "Failed to determine type: " << e.what() << endl;
            cout << "filename: " << fn << endl;
        }
        switch(fmt) {
            case FMT_ANA: return IDL_StrToSTRING( (char*)"ANA" );
            case FMT_FITS: return IDL_StrToSTRING( (char*)"FITS" );
            case FMT_NCDF: return IDL_StrToSTRING( (char*)"NCDF" );
            case FMT_MOMFBD: return IDL_StrToSTRING( (char*)"MOMFBD" );
            default: return IDL_StrToSTRING( (char*)"" );
        }
    } else {
        IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>(filenames->value.arr->data);
        vector<string> fileTypes;
        for( int i=0; i<filenames->value.arr->n_elts; ++i ) {
            bfs::path fn( string(strptr[i].s) );
            Format fmt = FMT_NONE;
            try {
                if( bfs::is_regular_file(fn) ) {
                    fmt = redux::file::readFmt( fn.string() );
                } else {
                    fmt = redux::file::guessFmt( fn.string() );
                }
            } catch( const exception& e) {
                cout << "Failed to determine type: " << e.what() << endl;
                cout << "filename: " << fn << endl;
            }
            switch(fmt) {
                case FMT_ANA: fileTypes.push_back("ANA"); break;
                case FMT_FITS: fileTypes.push_back("FITS"); break;
                case FMT_NCDF: fileTypes.push_back("NCDF"); break;
                case FMT_MOMFBD: fileTypes.push_back("MOMFBD"); break;
                default: fileTypes.push_back("");
            }
        }
        if( fileTypes.empty() ) return IDL_StrToSTRING( (char*)"" );
        IDL_MEMINT dims[] = { (IDL_MEMINT)fileTypes.size() };
        IDL_VPTR ret;
        IDL_STRING* retptr = (IDL_STRING*)IDL_MakeTempArray( IDL_TYP_STRING, 1, dims, IDL_ARR_INI_NOP, &ret );
        for( auto & ft: fileTypes ) {
            IDL_StrStore( retptr, (char*)ft.c_str() );
            retptr++;
        }
        return ret;
    }

    
}


extern "C" {

    int IDL_Load (void) {
        
        IDL_ExitRegister( IdlContainer::exit );
        
        IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)IdlContainer::info, (char*)"RDX", 0, 1, 0, 0 }, 0);
        IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)timer_start, (char*)"RDX_TIC", 0, 0, 0, 0 }, 0 );
        IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)timer_elapsed, (char*)"RDX_TOC", 0, 0, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 0, timer_info );
        IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)rdx_segment, (char*)"RDX_SEGMENT", 0, 4, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, segment_info );
        IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)rdx_hasopencv, (char*)"RDX_HASOPENCV", 0, 0, 0, 0 }, 1, hasopencv_info );
        IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)rdx_filetype, (char*)"RDX_FILETYPE", 0, 1, 0, 0 }, 1, filetype_info );
        
        static IDL_SYSFUN_DEF2 function_addr[] = {
            { { (IDL_VPTR (*) ()) cbezier2}, (char*) "RDX_CBEZIER2", 3, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { { (IDL_VPTR (*) ()) cbezier3}, (char*) "RDX_CBEZIER3", 3, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { { (IDL_VPTR (*) ()) redux::inpaint}, (char*) "RDX_INPAINT", 2, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { { (IDL_VPTR (*) ()) redux::img_align}, (char*) "RDX_IMG_ALIGN", 2, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { { (IDL_VPTR (*) ()) redux::img_project}, (char*) "RDX_IMG_PROJECT", 2, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { { (IDL_VPTR (*) ()) redux::img_remap}, (char*) "RDX_IMG_REMAP", 3, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
        };

        
        int ret = IdlContainer::load();
        ret &= IDL_SysRtnAdd (function_addr, TRUE, IDL_CARRAY_ELTS (function_addr));
        
        return ret;

    }

    void IDL_ResetSession(void) {
        //cout << "IDL_ResetSession" << endl;
        IdlContainer::reset();
     //  IDL_ExitUnregister(&test_exit_callback);

    }
    
}
