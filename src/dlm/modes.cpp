#include "idlutil.hpp"


#include "redux/constants.hpp"
#include "redux/momfbd/modes.hpp"
#include "redux/util/cache.hpp"

using namespace redux::image;
using namespace redux::momfbd;
using namespace redux::util;
using namespace redux;

using namespace std;


string make_pupil_info( int lvl ) {
    
    string ret = "RDX_MAKE_PUPIL";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"     ");          // newline if lvl>1
        ret += "   Syntax:   pup = rdx_make_pupil(pixels,radius)\n";
    } else ret += "\n";

    return ret;
    
}


IDL_VPTR make_pupil(int argc, IDL_VPTR* argv, char* argk) {
    
    if (argc < 2) {
        cout << "rdx_make_pupil: needs 2 arguments: nPixels & pupil-radius (in pixels). " << endl;
        return IDL_GettmpInt (-1);
    }
    
    IDL_LONG nPixels = IDL_LongScalar(argv[0]);
    double radius = IDL_DoubleScalar(argv[1]);
    IDL_VPTR tmp;
 
    PupilInfo pi(nPixels, radius);
    Pupil pupil = redux::util::Cache::get<PupilInfo,Pupil>(pi);
    if( pupil.empty() ) {    // this set was inserted, so it is not generated yet.
        pupil.generate( nPixels, radius );
    }
        
    IDL_MEMINT dims[] = { (int)pupil.dimSize(1), (int)pupil.dimSize(0) };
    float* tmpData = (float*)IDL_MakeTempArray(IDL_TYP_FLOAT,2,dims,IDL_ARR_INI_NOP,&tmp);
    pupil.copyTo<float>(tmpData);
    
    return tmp;

}


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
    { (char*) "ANGLE",     IDL_TYP_FLOAT, 1, 0,                      0, (char*) IDL_KW_OFFSETOF2(MODE_KW,angle) },
    { (char*) "CUTOFF",    IDL_TYP_FLOAT, 1, 0,                      0, (char*) IDL_KW_OFFSETOF2(MODE_KW,cutoff) },
    { (char*) "FIRST",     IDL_TYP_INT,   1, 0,                      0, (char*) IDL_KW_OFFSETOF2(MODE_KW,firstZernike) },
    { (char*) "HELP",      IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(MODE_KW,help) },
    { (char*) "LAST",      IDL_TYP_INT,   1, 0,                      0, (char*) IDL_KW_OFFSETOF2(MODE_KW,lastZernike) },
    { (char*) "NORMALIZE", IDL_TYP_INT,   1, 0,                      0, (char*) IDL_KW_OFFSETOF2(MODE_KW,normalize) },
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


IDL_VPTR make_modes(int argc, IDL_VPTR* argv, char* argk) {
    
    if (argc < 3) {
        cout << "rdx_make_modes: needs 3 arguments: mode-number(s), nPixels & pupil-radius (in pixels). " << endl;
        return IDL_GettmpInt (-1);
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
            cout << "rdx_make_modes: mode-numbers must be of type BYTE, INT or LONG." << endl;
            return IDL_GettmpInt (-1);
        }
    } else {
        modeNumbers.push_back( IDL_LongScalar(index) );
    }
    
    IDL_LONG nPixels = IDL_LongScalar(argv[1]);
    double radius = IDL_DoubleScalar(argv[2]);
    IDL_VPTR tmp;
 
    MODE_KW kw;
    kw.help = 0;
    kw.verbose = 0;
    kw.normalize = 0;
    kw.firstZernike = 2;
    kw.lastZernike = 2000;
    kw.angle = 0;
    kw.pupil = kw.variance = nullptr;
    kw.cutoff = 1E-3;

    (void) IDL_KWProcessByOffset (argc, argv, argk, mode_kw_pars, (IDL_VPTR*) 0, 255, &kw);

    kw.angle *= dtor;       // convert degrees to radians.
    
    if( kw.zernike ) {
        kw.firstZernike = kw.lastZernike = 0;
    }
    Pupil pupil;
    if( kw.pupil || kw.normalize ) {  // undefined has type=0
        PupilInfo pi(nPixels, radius);
        pupil = redux::util::Cache::get<PupilInfo,Pupil>(pi);
        if( pupil.empty() ) {    // this set was inserted, so it is not generated yet.
            pupil.generate( nPixels, radius );
        }
    }

    ModeInfo mi(kw.firstZernike, kw.lastZernike, 0, nPixels, radius, kw.angle, kw.cutoff);
    ModeSet& modes = Cache::get<ModeInfo,ModeSet>(mi);

    if( modes.empty() ) {
        if( kw.firstZernike && kw.lastZernike ) {
            modes.generate( nPixels, radius, kw.angle, kw.firstZernike, kw.lastZernike, modeNumbers, kw.cutoff );
        } else {    // first=last=0  =>  Zernike
            kw.firstZernike = kw.lastZernike = 0;   // in case only one was 0.
            modes.generate( nPixels, radius, kw.angle, modeNumbers );
        }
    }

    if( kw.normalize ) {
        modes.getNorms( pupil );
        cout << "rdx_make_modes: Normalizing not properly implemented at the moment." << endl;
        // TODO: do the actual normalization of the local copy of the modes below
    }
        
    IDL_MEMINT dims[] = { (int)modes.dimSize(2), (int)modes.dimSize(1), (int)modes.dimSize(0) };
    if( kw.variance ) {
        IDL_VPTR tmpVar;
        float* tmpData = (float*)IDL_MakeTempArray(IDL_TYP_FLOAT,1,dims+2,IDL_ARR_INI_NOP,&tmpVar);
        for( int i=0; i<dims[2]; ++i ) tmpData[i] = modes.atm_rms[i]*modes.atm_rms[i];
        IDL_VarCopy( tmpVar, kw.variance );
    }
    

    if( kw.pupil ) {

        IDL_VPTR tmpPup; // = IDL_Gettmp();
        float* tmpData =(float*)IDL_MakeTempArray( IDL_TYP_FLOAT, 2, dims, IDL_ARR_INI_NOP, &tmpPup );
        pupil.copyTo<float>(tmpData);
        IDL_VarCopy( tmpPup, kw.pupil );
    }
    
    float* tmpData = (float*)IDL_MakeTempArray(IDL_TYP_FLOAT,3,dims,IDL_ARR_INI_NOP,&tmp);
    modes.copyTo<float>(tmpData);
    
    return tmp;

}


void clear_pupils(void) {
    Cache::clear<PupilInfo,Pupil>();
}


void clear_modes(void) {
    Cache::clear<ModeInfo,ModeSet>();
}


string clear_modecache_info( int lvl ) {
    
    string ret = "RDX_CLEAR_MODES";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"    ");          // newline if lvl>1
        ret += "   Syntax:   rdx_clear_modes\n";
    } else ret += "\n";

    return ret;
    
}


void clear_modecache(int argc, IDL_VPTR argv[], char* argk) {
    Cache::clear<PupilInfo,Pupil>();
    Cache::clear<ModeInfo,ModeSet>();
}


string clear_cache_info( int lvl ) {
    
    string ret = "RDX_CLEAR_CACHE";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"    ");          // newline if lvl>1
        ret += "   Syntax:   rdx_clear_cache\n";
    } else ret += "\n";

    return ret;
    
}


void clear_cache(int argc, IDL_VPTR argv[], char* argk) {
    Cache::cleanup();
}


namespace {
    static int dummy UNUSED =
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)make_pupil, (char*)"RDX_MAKE_PUPIL", 2, 2, 0, 0 }, 1, make_pupil_info, clear_pupils ) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)make_modes, (char*)"RDX_MAKE_MODES", 3, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, make_modes_info, clear_modes ) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)clear_modecache, (char*)"RDX_CLEAR_MODES", 0, 0, 0, 0 }, 0 , clear_modecache_info) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)clear_cache, (char*)"RDX_CLEAR_CACHE", 0, 0, 0, 0 }, 0 , clear_cache_info);
}

