#include "idlutil.hpp"


#include "redux/constants.hpp"
#include "redux/momfbd/modes.hpp"
#include "redux/image/utils.hpp"
#include "redux/image/zernike.hpp"
#include "redux/util/cache.hpp"

using namespace redux::image;
using namespace redux::momfbd;
using namespace redux::util;
using namespace redux;

using namespace std;

namespace {

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
        IDL_VPTR area;
        IDL_VPTR center;
        IDL_INT  help;
        IDL_INT  verbose;
    } kw_pup;
    
    // NOTE:  The keywords MUST be listed in alphabetical order !!
    static IDL_KW_PAR kw_img_trans_pars[] = {
        IDL_KW_FAST_SCAN,
        { (char*) "AREA",          IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(kw_pup, area) },
        { (char*) "CENTER",        IDL_TYP_UNDEF, 1, IDL_KW_VIN|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2( kw_pup, center ) },
        { (char*) "HELP",          IDL_TYP_INT,   1, IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2( kw_pup, help ) },
        { (char*) "VERBOSE",       IDL_TYP_INT,   1, IDL_KW_ZERO, 0, (char*)IDL_KW_OFFSETOF2( kw_pup, verbose ) },
        { NULL }
    };

    string make_pupil_info( int lvl ) {
        
        string ret = "RDX_MAKE_PUPIL";
        if( lvl > 0 ) {
            ret += ((lvl > 1)?"\n":"     ");          // newline if lvl>1
            ret += "   Syntax:   pupil = rdx_make_pupil( pixels, radius [, central_obscuration_radius] )\n";
            if( lvl > 1 ) {
                ret +=  "   Accepted Keywords:\n"
                        "      AREA                (out) Return the area of the pupil.\n"
                        "      CENTER              Center of the pupil. (float, 1-2 elements, default is [pixels,pixels]/2+0.5).\n"
                        "      HELP                Display this info.\n"
                        "      VERBOSE             Verbosity, default is 0 (only error output).\n";
            }
        } else ret += "\n";

        return ret;
        
    }

}


IDL_VPTR make_pupil(int argc, IDL_VPTR* argv, char* argk) {

    kw_pup kw;
    int nPlainArgs = IDL_KWProcessByOffset( argc, argv, argk, kw_img_trans_pars, (IDL_VPTR*)0, 255, &kw );

    if( nPlainArgs < 2 ) {
        string msg = "Two arguments needed. nPixels & pupil-radius (in pixels).";
        msg += " An optional argument can be supplied to central-obscuration radius.";
        IDL_Message( IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP, msg.c_str() );
    }
    
    if( kw.help ) {
        int lvl = 1;
        if( kw.verbose ) lvl++;
        IDL_Message( IDL_M_NAMED_GENERIC, IDL_MSG_INFO, make_pupil_info(lvl).c_str() );
        return IDL_GettmpInt(0);
    }
    
    IDL_LONG nPixels = IDL_LongScalar( argv[0] );
    double radius = IDL_DoubleScalar( argv[1] );
    double innerRadius = 0.0;
    if( nPlainArgs == 3 ) {
        innerRadius = IDL_DoubleScalar( argv[2] );
    }
    
    float c = nPixels/2.0+0.5;
    PointF center( c, c );
    if( kw.center ) {
        vector<float> centv = getAsVector<float>( kw.center );
        if( centv.size() == 1 ) {    // just use as a flag which means center in output image
            centv.resize( 2, centv[0] );  // copy the value
        }
        if( centv.size() == 2 ) {
            center.x = centv[0];
            center.y = centv[1];
        }
    }
    
    IDL_VPTR tmp;
    IDL_MEMINT dims[] = { nPixels, nPixels };
    double* tmpData = (double*)IDL_MakeTempArray( IDL_TYP_DOUBLE, 2, dims, IDL_ARR_INI_NOP, &tmp );
    
    if( kw.verbose ) {
        string msg = "Generating " + printArray( dims, 2, "" ) + " pupil with radius " + to_string(radius);
        if( innerRadius > 0.0 ) msg += " and a central obscuration of radius " + to_string(innerRadius);
        if( kw.center ) msg += ", centered at " + (string)center;
        msg += ".";
        IDL_Message( IDL_M_NAMED_GENERIC, IDL_MSG_INFO, msg.c_str() );
    }
    
    auto tmp2D = reshapeArray( tmpData, nPixels, nPixels );
    double area = redux::image::makePupil( tmp2D.get(), nPixels, center, radius, innerRadius );

    if( kw.area ) {
        IDL_VPTR tmpArea = IDL_Gettmp();
        tmpArea->type = IDL_TYP_FLOAT;
        tmpArea->value.f = area;
        IDL_VarCopy( tmpArea, kw.area );
    }

    return tmp;

}


namespace {
    
    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
        IDL_INT help;
        IDL_INT firstZernike;
        IDL_INT lastZernike;
        IDL_INT normalize;
        float angle;
        float cutoff;
        IDL_VPTR pupil;
        IDL_VPTR variance;
        IDL_INT verbose;
        IDL_INT zernike;
    } MODE_KW;


    // NOTE:  The keywords MUST be listed in alphabetical order !!
    static IDL_KW_PAR mode_kw_pars[] = {
        IDL_KW_FAST_SCAN,
        { (char*) "ANGLE",     IDL_TYP_FLOAT, 1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(MODE_KW,angle) },
        { (char*) "CUTOFF",    IDL_TYP_FLOAT, 1, 0,                      0, (char*) IDL_KW_OFFSETOF2(MODE_KW,cutoff) },
        { (char*) "FIRST",     IDL_TYP_INT,   1, 0,                      0, (char*) IDL_KW_OFFSETOF2(MODE_KW,firstZernike) },
        { (char*) "HELP",      IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(MODE_KW,help) },
        { (char*) "LAST",      IDL_TYP_INT,   1, 0,                      0, (char*) IDL_KW_OFFSETOF2(MODE_KW,lastZernike) },
        { (char*) "NORMALIZE", IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(MODE_KW,normalize) },
        { (char*) "PUPIL",     IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(MODE_KW,pupil) },
        { (char*) "VARIANCE",  IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(MODE_KW,variance) },
        { (char*) "VERBOSE",   IDL_TYP_INT,   1, 0,                      0, (char*) IDL_KW_OFFSETOF2(MODE_KW,verbose) },
        { (char*) "ZERNIKE",   IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(MODE_KW,zernike) },
        { NULL }
    };


    string make_modes_info( int lvl ) {
        
        string ret = "RDX_MAKE_MODES";
        if( lvl > 0 ) {
            ret += ((lvl > 1)?"\n":"     ");          // newline if lvl>1
            ret += "   Syntax:   modes = rdx_make_modes([2,4,5,9],pixels,radius)\n";
            if( lvl > 1 ) {
                ret +=  "   Accepted Keywords:\n"
                        "      HELP                Display this info.\n"
                        "      ANGLE               Rotate the modes (degrees).\n"
                        "      CUTOFF              Smallest coefficient to consider for the Karhunen-Loeve expansion. (0.001)\n"
                        "      FIRST               First Zernike-mode to use for the Karhunen-Loeve expansion. (2)\n"
                        "      LAST                Last Zernike-mode to use for the Karhunen-Loeve expansion.(2000)\n"
                        "      NORMALIZE           Normalize the modes.\n"
                        "      PUPIL               (input/output) Use pupil, or keep the generated pupil.\n"
                        "      VARIANCE            (output) Mode variances.\n"
                        "      VERBOSE             Verbosity, default is 0 (only error output).\n"
                        "      ZERNIKE             Use Zernike modes (default is KL).\n";
            }
        } else ret += "\n";

        return ret;
        
    }

}


IDL_VPTR make_modes(int argc, IDL_VPTR* argv, char* argk) {
    
    MODE_KW kw;
    kw.firstZernike = 2;
    kw.lastZernike = 2000;
    kw.cutoff = 1E-3;
    int nPlainArgs = IDL_KWProcessByOffset (argc, argv, argk, mode_kw_pars, (IDL_VPTR*) 0, 255, &kw);
    
    if (nPlainArgs < 3) {
        string msg = "Three arguments needed: Mode-number(s), nPixels & pupil-radius (in pixels).";
        IDL_Message( IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP, msg.c_str() );
    }
    
    IDL_VPTR index = argv[0];
    IDL_ENSURE_SIMPLE(index);
    
    vector<uint16_t> modeNumbers;

    if ( (index->flags & IDL_V_ARR) ) {
        if( index->type == IDL_TYP_BYTE ) {
            auto beg = index->value.arr->data;
            std::copy( beg, beg + index->value.arr->n_elts , back_inserter(modeNumbers) );
        } else if (index->type == IDL_TYP_INT ) {
            auto beg = reinterpret_cast<int16_t*>(index->value.arr->data);
            std::copy( beg, beg + index->value.arr->n_elts , back_inserter(modeNumbers) );
        } else if (index->type == IDL_TYP_LONG ) {
            auto beg = reinterpret_cast<int32_t*>(index->value.arr->data);
            std::copy( beg, beg + index->value.arr->n_elts , back_inserter(modeNumbers) );
        } else  {
            string msg = "Mode-numbers must be of type BYTE, INT or LONG.";
            IDL_Message( IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP, msg.c_str() );
        }
    } else {
        if( index->type == IDL_TYP_STRING ) {
            modeNumbers = stringToUInts<uint16_t>( IDL_VarGetString(index) );
        } else {
            modeNumbers.push_back( IDL_LongScalar(index) );
        }
    }
    
    IDL_LONG nPixels = IDL_LongScalar(argv[1]);
    double radius = IDL_DoubleScalar(argv[2]);
    IDL_VPTR tmp;
 

    kw.angle *= dtor;       // convert degrees to radians.
    
    if( kw.zernike ) {
        kw.firstZernike = kw.lastZernike = 0;
    }
    
    int flags(Zernike::GET_RADIAL|Zernike::GET_ANGULAR);

    Pupil pupil;
    if( kw.pupil || kw.normalize ) {  // undefined has type=0
        PupilInfo pi(nPixels, radius);
        pupil = redux::util::Cache::get<PupilInfo,Pupil>(pi);
        if( pupil.empty() ) {    // this set was inserted, so it is not generated yet.
            pupil.generate( nPixels, radius );
        }
    }
    
    ModeInfo mi(kw.firstZernike, kw.lastZernike, modeNumbers, nPixels, radius, kw.angle, kw.cutoff);
    if( kw.verbose ) {
        flags |= Zernike::VERBOSE;
        string msg = "Generating modes: " + (string)mi+ ".";
        IDL_Message( IDL_M_NAMED_GENERIC, IDL_MSG_INFO, msg.c_str() );
    }

    shared_ptr<ModeSet>& modesRef = Cache::get<ModeInfo,shared_ptr<ModeSet>>(mi );
    if( ! modesRef ) modesRef.reset(new ModeSet());
    if( modesRef->empty() ) {
        if( kw.firstZernike && kw.lastZernike ) {
            modesRef->generate( nPixels, radius, kw.angle, kw.firstZernike, kw.lastZernike, modeNumbers, kw.cutoff, flags );
        } else {    // first=last=0  =>  Zernike
            kw.firstZernike = kw.lastZernike = 0;   // in case only one was 0.
            modesRef->generate( nPixels, radius, kw.angle, modeNumbers, flags );
        }
    }
    
    ModeSet modes = modesRef->clone();
    if( kw.normalize ) {
        modes.getNorms( pupil );
        modes.normalize( 1.0 );
    }
    
    vector<IDL_LONG64> dimss;
    std::copy( modes.dimensions().begin(), modes.dimensions().end(), back_inserter(dimss) );
    std::reverse( dimss.begin(), dimss.end() );
    if( *(dimss.rbegin()) == 1 ) dimss.pop_back();      // trivial dimensions should be removed in the cloning step, so this should not be needed.
    
    int nDims = dimss.size();
    IDL_MEMINT nModes = 1;
    if( nDims == 3 ) {
        nModes = dimss[2];
    }
    if( kw.variance ) {
        IDL_VPTR tmpVar;
        double* tmpData = (double*)IDL_MakeTempArray( IDL_TYP_DOUBLE, 1, &nModes, IDL_ARR_INI_NOP, &tmpVar );
        for( int i=0; i<nModes; ++i ) tmpData[i] = modes.atm_rms[i]*modes.atm_rms[i];
        IDL_VarCopy( tmpVar, kw.variance );
    }
    
    if( kw.pupil ) {
        IDL_VPTR tmpPup; // = IDL_Gettmp();
        double* tmpData =(double*)IDL_MakeTempArray( IDL_TYP_DOUBLE, 2, dimss.data(), IDL_ARR_INI_NOP, &tmpPup );
        pupil.copyTo<double>(tmpData);
        IDL_VarCopy( tmpPup, kw.pupil );
    }

    
    double* tmpData = (double*)IDL_MakeTempArray( IDL_TYP_DOUBLE, nDims, dimss.data(), IDL_ARR_INI_NOP, &tmp );
    modes.copyTo<double>(tmpData);
    
    return tmp;

}


void clear_modes(void) {
    Cache::clear<ModeInfo,ModeSet>();
    Cache::clear<PupilInfo,Pupil>();
}


string clear_modecache_info( int lvl ) {
    
    string ret = "RDX_CLEAR_MODES";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"    ");          // newline if lvl>1
        ret +=  "  Clear the cache of existing modes.\n";
        ret += "   Syntax:   rdx_clear_modes\n";
    } else ret += "\n";

    return ret;
    
}


void clear_modecache(int, IDL_VPTR, char* ) {
    clear_modes();
}


namespace {
    static int dummy RDX_UNUSED =
    IdlContainer::registerRoutine( {{(IDL_SYSRTN_GENERIC)make_pupil}, (char*)"RDX_MAKE_PUPIL", 2, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, make_pupil_info ) +
    IdlContainer::registerRoutine( {{(IDL_SYSRTN_GENERIC)make_modes}, (char*)"RDX_MAKE_MODES", 3, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, make_modes_info, clear_modes ) +
    IdlContainer::registerRoutine( {{(IDL_SYSRTN_GENERIC)clear_modecache}, (char*)"RDX_CLEAR_MODES", 0, 0, 0, 0 }, 0 , clear_modecache_info);
}

