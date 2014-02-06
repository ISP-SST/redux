#ifndef REDUX_MOMFBD_DEFINES_HPP
#define REDUX_MOMFBD_DEFINES_HPP


#define CFG_RECONSTR           0     // reconstruction mode
#define CFG_CALIBRATE          1     // calibration mode
#define CFG_FF_ONLY            2     // only flatfield and clip, then exit

#define CFG_GRADIENT_DIFF      1     // finite difference
#define CFG_GRADIENT_VOGEL     2     // Vogel

#define CFG_GETSTEP_SDSC       1     // steepest descent
#define CFG_GETSTEP_CNJG       2     // conjugate gradient
#define CFG_GETSTEP_BFGS_inv   3     // inverse BFGS?
#define CFG_GETSTEP_BFGS       4     // BFGS?
#define CFG_GETSTEP_NAMES {"unknown","steepest descent","conjugate_gradient","inverse BFGS","unknown"}

#define CFG_FPM_MEDIAN         1
#define CFG_FPM_INVDISTWEIGHT  2
#define CFG_FPM_HORINT  3
#define CFG_FPMETHOD_NAMES {"unknown","median","inverse distance weighted","horizontal interpolation", "unknown"}

#define CFG_DEFAULT           0
#define CFG_ZERNIKE           1
#define CFG_KARHUNEN_LOEVE    2
#define CFG_OTF_BASIS_FUNC    3
#define CFG_MODE_MAX          3


#define MFBD_FT_NAMES         {"NONE","ANA","FITS","MOMFBD"}
#define MFBD_FT_EXT           {"NONE","f0","fits","momfbd"}

// flags
#define MFBD_NONE                                0
#define MFBD_FT_ANA                     0x00000001
#define MFBD_FT_FITS                    0x00000002
#define MFBD_FT_MOMFBD                  0x00000004
#define MFBD_FT_MASK                    0x00000007

#define MFBD_CALIBRATE                  0x00000100
#define MFBD_DONT_MATCH_IMAGE_NUMS      0x00000200
#define MFBD_FAST_QR                    0x00000400
#define MFBD_FIT_PLANE                  0x00000800
#define MFBD_FORCE_WRITE                0x00001000
#define MFBD_FLATFIELD                  0x00002000
#define MFBD_GET_ALPHA                  0x00004000
#define MFBD_GET_COBJ                   0x00008000
#define MFBD_GET_DIVERSITY              0x00010000
#define MFBD_GET_METRIC                 0x00020000
#define MFBD_GET_MODES                  0x00040000
#define MFBD_GET_PSF                    0x00080000
#define MFBD_GET_PSF_AVG                0x00100000
#define MFBD_GET_RESIDUAL               0x00200000
#define MFBD_GLOBAL_NOISE               0x00400000
#define MFBD_NEW_CONSTRAINTS            0x00800000
#define MFBD_NO_CLIP                    0x01000000
#define MFBD_NO_CONSTRAINTS             0x02000000
#define MFBD_NO_FILTER                  0x04000000
#define MFBD_NO_RESTORE                 0x08000000
#define MFBD_SAVE_FFDATA                0x10000000
#define MFBD_SWAP                       0x20000000



#define MFBD_I08T   0
#define MFBD_I16T   1
#define MFBD_I32T   2
#define MFBD_I64T   3
#define MFBD_F32T   4
#define MFBD_F64T   5

#define MFBD_TYPE_SIZES {1,2,4,8,4,8}
#define MFBD_TYPE_NAMES {"byte","int16","int32","int64","float32","float64"}


#define DEF_BASIS             "Zernike"
#define DEF_MODES             "2-35"
#define DEF_KL_MIN_MODE       2
#define DEF_KL_MAX_MODE       2000
#define DEF_NUM_POINTS        128
#define DEF_TELESCOPE_D       0.97
#define DEF_TELESCOPE_F       47.0
#define DEF_ARCSECPERPIX      0.041
#define DEF_PIXELSIZE         9.00E-06
#define DEF_MIN_ITER          5
#define DEF_MAX_ITER          500
#define DEF_N_DONE_ITER       3
#define DEF_FTOL              1.00E-03
#define DEF_EPS               1.00E-10
#define DEF_SVD_REG           1.00E-03
#define DEF_GETSTEP           "getstep_BFGS_inv"
#define DEF_GRADIENT          "gradient_diff"
#define DEF_REG_GAMMA         1.00E-04
#define DEF_PROG_DATA_DIR     "./data"
#define DEF_STOKES_WEIGHTS    {1.0,0.0,0.0,0.0}
#define DEF_BORDER_CLIP       100
#define DEF_DATE_OBS          "N/A"
#define DEF_MAX_LOCAL_SHIFT   5
#define DEF_MODE_START        5
#define DEF_MODE_STEP         5
#define DEF_NF                1.0
#define DEF_FILE_TYPE         CFG_FT_ANA
#define DEF_DATA_TYPE         "SHORT"
#define DEF_ANGLE             0.0
#define DEF_WEIGHT            1.0
#define DEF_FPMETHOD          "invdistweight"

#include <map>
#include <string>
namespace redux {
    namespace momfbd {
        enum FileType { FT_NONE = MFBD_NONE, FT_ANA = MFBD_FT_ANA, FT_FITS = MFBD_FT_FITS, FT_MOMFBD = MFBD_FT_MOMFBD };
        const std::map<FileType, std::string> FileTypeNames = {{FT_NONE, ""},
            {FT_ANA, "ANA"},
            {FT_FITS, "FITS"},
            {FT_MOMFBD, "MOMFBD"}
        };
        const std::map<FileType, std::string> FileTypeExtensions = {{FT_NONE, ""},
            {FT_ANA, "f0"},
            {FT_FITS, "fits"},
            {FT_MOMFBD, "momfbd"}
        };
    }
}


#endif  // REDUX_MOMFBD_DEFINES_HPP
