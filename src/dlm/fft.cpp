#include "idlutil.hpp"


#include "redux/image/descatter.hpp"
#include "redux/image/fouriertransform.hpp"
#include "redux/util/arraystats.hpp"
#include "redux/util/cache.hpp"

#include <atomic>
#include <thread>

#include <fftw3.h>

#include <boost/format.hpp>

using namespace redux::image;
using namespace redux::util;
using namespace redux;

using namespace std;

namespace {
    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
        IDL_INT help;
        IDL_INT correlate;
        IDL_INT in_place;
        IDL_INT nonormalize;
        UCHAR nthreads;
        IDL_INT verbose;
    } KW_CONV_RESULT;

    // NOTE:  The keywords MUST be listed in alphabetical order !!
    static IDL_KW_PAR kw_conv_pars[] = {
        IDL_KW_FAST_SCAN,
        { (char*) "CORRELATE",    IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(KW_CONV_RESULT,correlate) },
        { (char*) "HELP",         IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(KW_CONV_RESULT,help) },
        { (char*) "IN_PLACE",     IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(KW_CONV_RESULT,in_place) },
        { (char*) "NONORMALIZE",  IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(KW_CONV_RESULT,nonormalize) },
        { (char*) "NTHREADS",     IDL_TYP_BYTE, 1, 0,           0, (char*) IDL_KW_OFFSETOF2(KW_CONV_RESULT,nthreads) },
        { (char*) "VERBOSE",      IDL_TYP_INT, 1, 0,           0, (char*) IDL_KW_OFFSETOF2(KW_CONV_RESULT,verbose) },
        { NULL }
    };
}


void fft_memalign(int argc, IDL_VPTR* argv);


string convolve_info( int lvl ) {
    
    string ret = "RDX_CONVOLVE";
    if( lvl > 0 ) {
        
        ret += ((lvl > 1)?"\n":"       ");          // newline if lvl>1
        ret += "   Syntax:   out = rdx_convolve(img(s), psf, /KEYWORDS)\n";
    
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      CORRELATE           Compute correlation instead.\n"
                    "      IN_PLACE            Replace input with result instead of using a temporary.\n"
                    "      NONORMALIZE         Do not normalize PSF.\n"
                    "      NTHREADS            Number of threads.\n"
                    "      VERBOSE             Verbosity, default is 0 (only error output).\n";
        }
    
    } else ret += "\n";

    return ret;
    
}


IDL_VPTR convolve (int argc, IDL_VPTR* argv, char* argk) {

    KW_CONV_RESULT kw;
    kw.nthreads = thread::hardware_concurrency();
    int nPlainArgs = IDL_KWProcessByOffset (argc, argv, argk, kw_conv_pars, (IDL_VPTR*) 0, 255, &kw);
    
    if (nPlainArgs < 2) {
        printMessage( "Needs 2 arguments (image and psf).", IDL_MSG_LONGJMP );
    }

    IDL_VPTR imageVar  = argv[0];
    IDL_VPTR psfVar  = argv[1];

    IDL_ENSURE_SIMPLE(imageVar);
    IDL_ENSURE_ARRAY(imageVar);

    IDL_ENSURE_SIMPLE(psfVar);
    IDL_ENSURE_ARRAY(psfVar);

    if( kw.help ) {
        int lvl = 1;
        if( kw.verbose ) lvl++;
        printMessage( convolve_info(lvl) );
        return IDL_GettmpInt(-1);
    }

    kw.nthreads = max<UCHAR>(1, min<UCHAR>(kw.nthreads, thread::hardware_concurrency()));

    if( psfVar->value.arr->n_dim != 2 ) {
        printMessage( "PSF must be A 2D array.", IDL_MSG_LONGJMP );
    }
    IDL_LONG64 xSize = psfVar->value.arr->dim[0];
    IDL_LONG64 ySize = psfVar->value.arr->dim[1];

    if( xSize != imageVar->value.arr->dim[0] || ySize != imageVar->value.arr->dim[1] ) {
        string msg = "PSF and image(s) must be of the same size."
                   + printArray( imageVar->value.arr->dim, imageVar->value.arr->n_dim, "\n   img")
                   + printArray( psfVar->value.arr->dim, psfVar->value.arr->n_dim, "\n   psf");
        printMessage( msg, IDL_MSG_LONGJMP );
    }

    if( !xSize || !ySize ) {
        printMessage( "Input has zero-size.", IDL_MSG_LONGJMP );
    }

    size_t nImages = 1;
    IDL_LONG64 nDims = imageVar->value.arr->n_dim;
    if( nDims == 3 ) {          // multiple images
        nImages = imageVar->value.arr->dim[2];
    }
    
    size_t dataSize = xSize*ySize;
    size_t ftSize = ySize*(xSize/2+1);
    
    try {
        int fftAlign = 0; //fftw_alignment_of( (double*)ftSize );
        if( fftAlign ) fftAlign = (16 - fftAlign);

        string statusString;
        if( kw.verbose ) {
            statusString = (kw.correlate?"Correlating":"Convolving ") + to_string(nImages)
            + " images using " +to_string((int)kw.nthreads) + string(" thread") + ((kw.nthreads>1)?"s.":".");
            cout << statusString << ((kw.verbose == 1)?"\n":"") << flush;
        }
        
        FourierTransform::Plan::Ptr plan = FourierTransform::Plan::get( ySize, xSize, FourierTransform::Plan::R2C, 1 );

        size_t tmpSize = (kw.nthreads+1)*(ftSize+fftAlign);
        std::shared_ptr<double> tmpComplex( (double*)fftw_malloc(tmpSize*sizeof(fftw_complex)), fftw_free );
        fftw_complex* psfFT = reinterpret_cast<fftw_complex*>(tmpComplex.get());
        fftw_complex* ptrC = psfFT + ftSize + fftAlign;

        shared_ptr<double> psf( (double*)fftw_malloc(dataSize*sizeof(double)), fftw_free );
        double* psfData = psf.get();
        copyToRaw( psfVar, psfData );

        double norm = 1.0 / dataSize;
        if( !kw.nonormalize ) {
            double sum = 0;
            for( size_t i=0; i<dataSize; ++i ) sum += psfData[i];
            norm /= sum;
        }
        
        for( size_t i=0; i<dataSize; ++i ) psfData[i] *= norm;
        plan->forward( psfData, psfFT );
        
        if( kw.correlate ) {        // conjugate second input
            for( size_t n=0; n<ftSize; ++n ) reinterpret_cast<complex_t&>(psfFT[n]) = conj(reinterpret_cast<complex_t&>(psfFT[n]));
        }
        
        double* allocatedData = nullptr;
        double* dataPtr = nullptr;
        if( kw.in_place && (imageVar->type == IDL_TYP_DOUBLE) ) {
            fft_memalign( 1, &imageVar );   // align if needed, else no action.
            dataPtr = reinterpret_cast<double*>( imageVar->value.arr->data );
        } else if ( kw.in_place ) {
            printMessage( "Images have to be of type double when transforming in_place.", IDL_MSG_LONGJMP );
        } else {
            //dataAlign = 0;
            allocatedData = (double*)fftw_malloc(nImages*dataSize*sizeof(double));
            dataPtr = allocatedData;
            copyToRaw( imageVar, dataPtr );
        }

        atomic<size_t> imgIndex(0);
        std::vector<std::thread> threads;           // TODO: check performance, multiple threads have almost no impact !!
        for (int i=0; i<kw.nthreads; ++i) {
            threads.push_back( std::thread(
                [=,&imgIndex](){
                    size_t myIndex;
                    fftw_complex* myFT = ptrC + i*(ftSize+fftAlign);
                    //fftw_complex* myFT = fftw_alloc_complex(ftSize);
                    while( (myIndex=imgIndex.fetch_add(1)) < nImages ) {
                        double* myData = dataPtr + myIndex*dataSize;
                        plan->forward( myData, myFT );
                        for( size_t n=0; n<ftSize; ++n) reinterpret_cast<complex_t&>(myFT[n]) *= reinterpret_cast<complex_t&>(psfFT[n]);
                        plan->backward( myFT, myData );
                        FourierTransform::reorder( myData, ySize, xSize );
                        if( kw.verbose > 1 ) printProgress( statusString, (myIndex*100.0/(nImages-1)));
                    }
                    //fftw_free(myFT);
                }));
        }
        for (auto& th : threads) th.join();

        if( kw.verbose > 1 ) {
            printProgress( statusString, 100.0 );
            cout << endl;
        }
        
        if( allocatedData ) {
            return IDL_ImportArray( nDims, imageVar->value.arr->dim, IDL_TYP_DOUBLE, (UCHAR*)allocatedData, (IDL_ARRAY_FREE_CB)fftw_free, NULL );
        }
    
    } catch( const exception& e ) {
        string msg = string("Unhandled exception: ") + e.what();
        printMessage( msg, IDL_MSG_LONGJMP );
    }


    return IDL_GettmpInt(0);
    
}


void convolve_cleanup( void ) {

    Cache::clear<FourierTransform::Plan::Index, FourierTransform::Plan::Ptr >();

}


typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_INT help;
    float epsilon;
    IDL_INT in_place;
    IDL_INT normalize;
    IDL_INT niter;
    UCHAR nthreads;
    IDL_INT padding;
    IDL_INT verbose;
} KW_DESCATTER_RESULT;

// NOTE:  The keywords MUST be listed in alphabetical order !!
static IDL_KW_PAR kw_descatter_pars[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "EPSILON",   IDL_TYP_FLOAT,  1, 0,           0, (char*) IDL_KW_OFFSETOF2(KW_DESCATTER_RESULT,epsilon) },
    { (char*) "HELP",      IDL_TYP_INT,    1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(KW_DESCATTER_RESULT,help) },
    { (char*) "IN_PLACE",  IDL_TYP_INT,    1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(KW_DESCATTER_RESULT,in_place) },
    { (char*) "NITER",     IDL_TYP_INT,    1, 0,           0, (char*) IDL_KW_OFFSETOF2(KW_DESCATTER_RESULT,niter) },
    { (char*) "NORMALIZE", IDL_TYP_INT,    1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(KW_DESCATTER_RESULT,normalize) },
    { (char*) "NTHREADS",  IDL_TYP_BYTE,   1, 0,           0, (char*) IDL_KW_OFFSETOF2(KW_DESCATTER_RESULT,nthreads) },
    { (char*) "PADDING",   IDL_TYP_INT,    1, 0,           0, (char*) IDL_KW_OFFSETOF2(KW_DESCATTER_RESULT,padding) },
    { (char*) "VERBOSE",   IDL_TYP_INT,    1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(KW_DESCATTER_RESULT,verbose) },
    { NULL }
};

string rdx_descatter_info( int lvl ) {
    
    string ret = "RDX_DESCATTER";
    if( lvl > 0 ) {
        
        ret += ((lvl > 1)?"\n":"      ");          // newline if lvl>1
        ret += "   Syntax:   out = rdx_descatter(img(s), backgain, psf, /KEYWORDS)\n";
    
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      EPSILON             Exit loop when ChiSqr is smaller than epsilon. (1E-8)\n"
                    "      IN_PLACE            Replace input with result instead of using a temporary.\n"
                    "      NONORMALIZE         Do not normalize PSF.\n"
                    "      NTHREADS            Number of threads.\n"
                    "      PADDING             Padding around image to avoid periodic edge effects. (256)\n"
                    "      VERBOSE             Verbosity, default is 0 (only error output).\n";
        }
    
    } else ret += "\n";

    return ret;
    
}

IDL_VPTR rdx_descatter( int argc, IDL_VPTR* argv, char* argk ) {

    KW_DESCATTER_RESULT kw;
    kw.epsilon = 1E-8;
    kw.padding = 256;
    kw.niter = 50;
    kw.nthreads = 1; //thread::hardware_concurrency();
    int nPlainArgs = IDL_KWProcessByOffset( argc, argv, argk, kw_descatter_pars, (IDL_VPTR*) 0, 255, &kw );

    if( nPlainArgs < 3 ) {
        string msg = "Needs at least 3 arguments, image(s), backgain and PSF.";
        printMessage( msg, IDL_MSG_LONGJMP );
    }

    IDL_VPTR imageVar  = argv[0];
    IDL_VPTR gainVar  = argv[1];
    IDL_VPTR psfVar  = argv[2];

    IDL_ENSURE_SIMPLE(imageVar);
    IDL_ENSURE_ARRAY(imageVar);

    IDL_ENSURE_SIMPLE(gainVar);
    IDL_ENSURE_ARRAY(gainVar);

    IDL_ENSURE_SIMPLE(psfVar);
    IDL_ENSURE_ARRAY(psfVar);

    if( kw.help ) {
        int lvl = 1;
        if( kw.verbose ) lvl++;
        printMessage( rdx_descatter_info(lvl) );
        return IDL_GettmpInt (-1);
    }

    kw.nthreads = max<UCHAR>(1, min<UCHAR>(kw.nthreads, thread::hardware_concurrency()));
    
    if( psfVar->value.arr->n_dim != 2 || gainVar->value.arr->n_dim != 2 ) {
        printMessage( "Backgain and psf must be A 2D arrays.", IDL_MSG_LONGJMP );
    }
    IDL_LONG64 xSize = psfVar->value.arr->dim[0];
    IDL_LONG64 ySize = psfVar->value.arr->dim[1];

    if( xSize != gainVar->value.arr->dim[0] || ySize != gainVar->value.arr->dim[1] ||
        xSize != imageVar->value.arr->dim[0] || ySize != imageVar->value.arr->dim[1] ) {
        string msg = "Backgain, PSF and image(s) must be of the same size.";
        msg += printArray(imageVar->value.arr->dim, imageVar->value.arr->n_dim, "\n   img")
             + printArray(gainVar->value.arr->dim, gainVar->value.arr->n_dim, "\n  gain")
             + printArray(psfVar->value.arr->dim, psfVar->value.arr->n_dim, "\n   psf");
        printMessage( msg, IDL_MSG_LONGJMP );
    }
    
    if( !xSize || !ySize ) {
        printMessage( "Input has zero-size.", IDL_MSG_LONGJMP );
    }

    size_t nImages = 1;
    IDL_LONG64 nDims = imageVar->value.arr->n_dim;
    if( nDims == 3 ) {          // multiple images
        nImages = imageVar->value.arr->dim[2];
    }

    
    size_t paddedX = xSize + 2*kw.padding;
    size_t paddedY = ySize + 2*kw.padding;
    
    size_t posX = (paddedX-xSize)/2;
    size_t posY = (paddedY-ySize)/2;
    
    size_t dataSize = paddedX*paddedY;
    size_t ftSize = paddedY*(paddedX/2+1);

    try {
        
        shared_ptr<double> psf, gain;
        if( kw.padding ) {
            gain.reset( (double*)fftw_malloc(dataSize*sizeof(double)), fftw_free);
            memset( gain.get(), 0, dataSize*sizeof(double) );
            copyInto( gainVar, gain.get(), paddedY, paddedX, posY, posX );
            psf.reset( (double*)fftw_malloc(dataSize*sizeof(double)), fftw_free);
            memset( psf.get(), 0, dataSize*sizeof(double) );
            copyInto( psfVar, psf.get(), paddedY, paddedX, posY, posX );
        } else {
            gain = castOrCopy<double>(gainVar);
            psf = castOrCopy<double>(psfVar);
        }
        
        double* gainPtr = gain.get();
        double* psfPtr = psf.get();
            
        // normalize psf
        double sum = 0;
        for( size_t i=0; i<dataSize; ++i ) sum += psfPtr[i];
        sum = 1.0/(sum*dataSize);
        for( size_t i=0; i<dataSize; ++i ) psfPtr[i] *= sum;
        
        FourierTransform::reorder( psfPtr, paddedY, paddedX );
        
        FourierTransform::Plan::Ptr plan = FourierTransform::Plan::get( paddedY, paddedX, FourierTransform::Plan::R2C, kw.nthreads );
        std::shared_ptr<double> otf( (double*)fftw_malloc(ftSize*sizeof(fftw_complex)), fftw_free );
        fftw_complex* otfPtr = reinterpret_cast<fftw_complex*>(otf.get());
        plan->forward( psfPtr, otfPtr );
       
        string statusString;
        if( kw.verbose ) {
            statusString = "Descattering " + to_string(nImages)
            + " images using " +to_string((int)kw.nthreads) + string(" thread") + ((kw.nthreads>1)?"s.":".");
            cout << statusString << ((kw.verbose == 1)?"\n":"") << flush;
        }
        
        double* allocatedData = nullptr;
        double* dataPtr = nullptr;
        if( kw.in_place && (imageVar->type == IDL_TYP_DOUBLE) ) {
            fft_memalign( 1, &imageVar );   // align if needed, else no action.
            dataPtr = reinterpret_cast<double*>( imageVar->value.arr->data );
        } else if ( kw.in_place ) {
            printMessage( "Images have to be of type double when transforming in_place.", IDL_MSG_LONGJMP );
        } else {
            //dataAlign = 0;
            allocatedData = (double*)fftw_malloc(nImages*dataSize*sizeof(double));
            dataPtr = allocatedData;
            copyToRaw( imageVar, dataPtr );
        }

        for( size_t n=0; n<nImages; ++n) {
            double* myData = dataPtr + n*dataSize;
            redux::image::descatter( myData, ySize, xSize, gainPtr, otfPtr, paddedY, paddedX, kw.nthreads, kw.niter, kw.epsilon );
            if( kw.verbose > 1 ) printProgress( statusString, (n*100.0/(nImages-1)));
        }
        
/*        atomic<size_t> imgIndex(0);
        std::vector<std::thread> threads;           // TODO: check performance, multiple threads have almost no impact !!
        for( int i=0; i<kw.nthreads; ++i ) {
            threads.push_back( std::thread(
                [=,&imgIndex](){
                    size_t myIndex;
                    while( (myIndex=imgIndex.fetch_add(1)) < nImages ) {
                        double* myData = dataPtr + myIndex*dataSize;
                        redux::image::descatter( myData, ySize, xSize, gainPtr, otfPtr, paddedY, paddedX, kw.nthreads, kw.niter, kw.epsilon );
                        if( kw.verbose > 1 ) printProgress( statusString, (myIndex*100.0/(nImages-1)));
                    }
                }));
        }
        for (auto& th : threads) th.join();
*/
        if( kw.verbose > 1 ) {
            printProgress( statusString, 100.0 );
            cout << endl;
        }
        
        if( allocatedData ) {
            return IDL_ImportArray( nDims, imageVar->value.arr->dim, IDL_TYP_DOUBLE, (UCHAR*)allocatedData, (IDL_ARRAY_FREE_CB)fftw_free, NULL );
        }
        
    } catch( const exception& e ) {
        string msg = string("Unhandled exception: ") + e.what();
        printMessage( msg, IDL_MSG_LONGJMP );
    }
    
    return IDL_GettmpInt(0);

}


void rdx_descatter_cleanup( void ) {

    Cache::clear<FourierTransform::Plan::Index, FourierTransform::Plan::Ptr >();

}


string fft_reorder_info( int lvl ) {
    
    string ret = "RDX_REORDER";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"        ");          // newline if lvl>1
        ret += "   Syntax:   rdx_reorder,img\n";
        if( lvl > 1 ) {
            ret +=  "   Reorganize quadrants of the image, according to the fftw3 ordering\n";
        }
    } else ret += "\n";

    return ret;
    
}


void fft_reorder(int argc, IDL_VPTR* argv) {

    IDL_VPTR array = argv[0];
    IDL_ENSURE_SIMPLE(array);
    IDL_ENSURE_ARRAY(array);
    
    if ( array->value.arr->n_dim != 2 ) {
        string msg = "A single 2D array should be passed as input.";
        printMessage( msg, IDL_MSG_LONGJMP );
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
        string msg = string("Error during call to FourierTransform::reorder(): ") + e.what();
        printMessage( msg );
    }

}


string fft_memalign_info( int lvl ) {
    
    string ret = "RDX_MEMALIGN";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"       ");          // newline if lvl>1
        ret += "   Syntax:   rdx_memalign,data\n";
        if( lvl > 1 ) {
            ret +=  "   Align array such that the fftw3 library can use vectorized instructions.\n";
        }
    } else ret += "\n";

    return ret;
    
}


void fft_memalign( int argc, IDL_VPTR* argv ) {

    IDL_VPTR array = argv[0];
    IDL_ENSURE_SIMPLE(array);
    IDL_ENSURE_ARRAY(array);
    
    try {
        
        int alignIn = 0; // fftw_alignment_of( reinterpret_cast<double*>( array->value.arr->data ) );
        if( alignIn == 0 ) return;  // nothing to do

        size_t nElements = array->value.arr->n_elts;
        size_t totalSize = nElements * IDL_TypeSizeFunc( array->type );

        switch(array->type) {
            case IDL_TYP_DOUBLE: {
                double* data = reinterpret_cast<double*>( array->value.arr->data );
                double* newData = (double*)fftw_malloc(nElements*sizeof(double));
                memcpy( newData, data, totalSize );
                IDL_VPTR tmpVar = IDL_ImportArray( array->value.arr->n_dim, array->value.arr->dim, array->type, (UCHAR*)newData, (IDL_ARRAY_FREE_CB)fftw_free, NULL );
                IDL_VarCopy( tmpVar, array );
                break;
            }
            case IDL_TYP_DCOMPLEX: {
                fftw_complex* data = reinterpret_cast<fftw_complex*>( array->value.arr->data );
                fftw_complex* newData = (fftw_complex*)fftw_malloc(nElements*sizeof(fftw_complex));
                memcpy( newData, data, totalSize );
                IDL_VPTR tmpVar = IDL_ImportArray( array->value.arr->n_dim, array->value.arr->dim, array->type, (UCHAR*)newData, (IDL_ARRAY_FREE_CB)fftw_free, NULL );
                IDL_VarCopy( tmpVar, array );
                break;
            }
            default: ;          // no action
        }
        
    } catch ( const exception& e ) {   // catch-all to avoid killing the idl session.
        string msg = string("Error during memory alignment: ") + e.what();
        printMessage( msg );
    }

    
}


void fft_cleanup(void) {

}


namespace {
    static int dummy RDX_UNUSED =
    IdlContainer::registerRoutine( {{(IDL_SYSRTN_GENERIC)convolve}, (char*)"RDX_CONVOLVE", 2, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, convolve_info, convolve_cleanup ) +
    IdlContainer::registerRoutine( {{(IDL_SYSRTN_GENERIC)rdx_descatter}, (char*)"RDX_DESCATTER", 3, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, rdx_descatter_info/*, descatter_cleanup*/ ) +
    IdlContainer::registerRoutine( {{(IDL_SYSRTN_GENERIC)fft_reorder}, (char*)"RDX_REORDER", 1, 1, 0, 0 }, 0 , fft_reorder_info) +
    IdlContainer::registerRoutine( {{(IDL_SYSRTN_GENERIC)fft_memalign}, (char*)"RDX_MEMALIGN", 1, 1, 0, 0 }, 0 , fft_memalign_info);
}
