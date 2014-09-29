#include "momfbd.hpp"

#include "redux/file/filemomfbd.hpp"
#include "redux/math/functions.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/endian.hpp"
#include "redux/util/stringutil.hpp"

#include <iostream>
#include <map>

using namespace redux::file;
using namespace redux::math;
using namespace redux::util;
using namespace redux;
using namespace std;

//a=momfbd_read('/home/devel/bigtest.momfbd')
//a=momfbd_read('/scratch/Data.local/camXIX_im22Apr2008.1804856..1805202.momfbd')
//a=momfbd_read('/scratch/tomas/indata_set1.1..2.momfbd',verbose=2,/all)
//a=momfbd_read('/scratch/pit/2010.10.05/BLUE/Results/Ca-WB.00014.momfbd',verbose=2,/all)
//a=momfbd_read('/scratch/roald/sims/momfbd/r0_5_new/image_t.-100.lc0.80..116.momfbd',verbose=2,/all)
//a=momfbd_read('/scratch/roald/2009/2009.05.23/6302_diskcenter_newmomfbd/camXVIII_im23May2009.fp3.28875.35732.lc0.1887067..1887831.momfbd',verbose=2,/all)
//for q=0,30 do for o=5,10 do for i=0,6 do for j=0,6 do tvscl,a.patch(i+o,j+o).img,128*i,128*j
//for i=0,36 do tvscale,a.mode(*,*,i)
//for i=0,181 do tvscale,a.patch(1,1).psf(*,*,i)
//o=9 & for q=0,181 do for i=0,6 do for j=0,6 do tvscl,a.patch(i+o,j+o).obj(*,*,q),128*i,128*j
//o=9 & for q=0,181 do for i=0,6 do for j=0,6 do tvscl,a.patch(i+o,j+o).psf(*,*,q),128*i,128*j
//img=momfbd_mozaic( a.patch.img, a.patch(*,0).xl, a.patch(*,0).xh, a.patch(0,*).yl, a.patch(0,*).yh, /clip )
//o=10 & for f=0,181 do for i=0,6 do for j=0,6 do tvscl,a.patch(i+o,j+o).obj(*,*,f),128*i,128*j
//setenv IDL_DLM_PATH '/home/hillberg/lib/dlm-x86_64:<IDL_DEFAULT>'
//setenv IDL_DLM_PATH '/home/hillberg/lib/dlm-old:<IDL_DEFAULT>'
//http://idlastro.gsfc.nasa.gov/idl_html_help/Structure_Variables.html

#define LOAD_IMG        1
#define LOAD_PSF        2
#define LOAD_OBJ        4
#define LOAD_RES        8
#define LOAD_ALPHA     16
#define LOAD_DIV       32
#define LOAD_PATCH     63       // convenience for flagging any "per-patch" data.
#define LOAD_MODES     64
#define LOAD_NAMES    128
#define LOAD_ALL      255

namespace {



    const char *var_names[] = { "undef", "int8", "int16", "int32", "float32",
                                "float64", "complex", "string", "struct", "double complex",
                                "int32", "int32", "uint16", "uint32", "int64", "uint64"
                              };

    const size_t var_sizes[] = { 0, sizeof ( UCHAR ), sizeof ( IDL_INT ), sizeof ( IDL_LONG ), sizeof ( float ),
                                 sizeof ( double ), sizeof ( IDL_COMPLEX ), sizeof ( IDL_STRING ), 0, sizeof ( IDL_DCOMPLEX ),
                                 sizeof ( IDL_ULONG ), sizeof ( IDL_ULONG ), sizeof ( IDL_UINT ), sizeof ( IDL_ULONG ), sizeof ( IDL_LONG64 ), sizeof ( IDL_ULONG64 )
                               };


//
// print the layout of an IDL structure and return the total data-size
//
    size_t dumpStruct ( IDL_VPTR data, int current, int indent ) {

        size_t sz = 0;
        if ( data->type == IDL_TYP_STRUCT ) {

            if ( current < 0 ) {
                cout << "              TAG        TYPE                     ADDRESS             OFFSET         SIZE           " << endl;
                sz += dumpStruct ( data, indent, indent );
                cout << string ( 50, ' ' ) << "Total Size:                        " << to_string ( sz ) << endl;
                return sz;
            }


            IDL_StructDefPtr structDef = data->value.s.sdef;
            uint8_t *buf = data->value.s.arr->data;
            int nTags = IDL_StructNumTags ( structDef );
            int count;

            for ( int t = 0; t < nTags; ++t ) {

                char *name = IDL_StructTagNameByIndex ( structDef, t, 0, 0 );
                IDL_VPTR v;
                IDL_MEMINT offset = IDL_StructTagInfoByIndex ( structDef, t, 0, &v );
                string type = string ( var_names[v->type] );
                count = 1;
                if ( v->flags & IDL_V_ARR ) {
                    type.append ( "(" );
                    for ( int d = 0; d < v->value.arr->n_dim; ++d ) {
                        if ( d ) type.append ( "," );
                        count *= v->value.arr->dim[d];
                        type.append ( to_string ( ( int ) v->value.arr->dim[d] ) );
                    }
                    type.append ( ")" );
                }

                if ( v->type == IDL_TYP_STRUCT ) {
                    cout << string ( current, ' ' ) <<alignLeft ( name, 25 - current ) << alignLeft ( type, 25 );
                    cout << alignLeft ( hexString ( buf + offset ), 20 ) << alignLeft ( to_string ( ( size_t ) offset ), 15 ) << endl;
                    sz += count * dumpStruct ( v, current + indent, indent );
                } else {
                    sz += count * var_sizes[v->type];
                    cout << string ( current, ' ' ) << alignLeft ( name, 25 - current ) << alignLeft ( type, 25 );
                    cout << alignLeft ( hexString ( buf + offset ), 20 ) << alignLeft ( to_string ( ( size_t ) offset ), 15 );
                    cout << alignLeft ( to_string ( count * var_sizes[v->type] ), 15 ) << endl;
                }

            }

        }

        return sz;
    }

    
    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
        IDL_INT all;
        IDL_INT alpha;
        IDL_INT check;
        IDL_INT clip;
        IDL_INT div;
        IDL_INT help;
        IDL_INT img;
        IDL_INT margin;
        IDL_INT modes;
        IDL_INT names;
        IDL_INT obj;
        IDL_INT psf;
        IDL_INT res;
        IDL_INT transpose;
        IDL_INT verbose;
    } KW_RESULT;

    static IDL_KW_PAR kw_pars[] = {   IDL_KW_FAST_SCAN,       // NOTE:  The keywords MUST be listed in alphabetical order !!
        { ( char* ) "ALL",           IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( all ) },
        { ( char* ) "ALPHA",         IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( alpha ) },
        { ( char* ) "CHECK",         IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( check ) },
        { ( char* ) "CLIP",          IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( clip ) },
        { ( char* ) "DIV",           IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( div ) },
        { ( char* ) "HELP",          IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( help ) },
        { ( char* ) "IMG",           IDL_TYP_INT, 1, 0,                           0, ( char* ) IDL_KW_OFFSETOF ( img ) },
        { ( char* ) "MARGIN",        IDL_TYP_INT, 1, 0,                           0, ( char* ) IDL_KW_OFFSETOF ( margin ) },
        { ( char* ) "MODES",         IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( modes ) },
        { ( char* ) "NAMES",         IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( names ) },
        { ( char* ) "OBJ",           IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( obj ) },
        { ( char* ) "PSF",           IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( psf ) },
        { ( char* ) "RES",           IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( res ) },
        { ( char* ) "TRANSPOSE",     IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( transpose ) },
        { ( char* ) "VERBOSE",       IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( verbose ) },
        { NULL }
    };

    const uint8_t alignTo = 4;  // Aligning to 4-bit boundary (IDL structures are aligned to 4-byte boundaries).


    // Define a structure to contain the datablock + a vector with filenames
    struct MomfdContainer {
        int32_t nNames;
        IDL_STRING* ptr;
    };


    struct Overlaps {
        Overlaps ( int16_t n = 0, int16_t e = 0, int16_t s = 0, int16_t w = 0 ) : north ( n ), east ( e ), south ( s ), west ( w ) {}
        int16_t xPixels, yPixels, north, east, south, west;
        bool operator< ( const Overlaps& rhs ) const {
            return north < rhs.north || east < rhs.east || south < rhs.south || west < rhs.west;
        }
        bool operator== ( const Overlaps& rhs ) {
            return north == rhs.north && east == rhs.east && south == rhs.south && west == rhs.west;
        }
    };


    void cleanupFloat ( UCHAR* arg ) {
        delete[] ( ( float* ) arg );
    }


    void cleanupMomfbd ( UCHAR *arg ) {
        // The actual datablock allocated for the data begins at (arg - sizeof(struct MomfdContainer))
        // The preceeding structure just serves to point to the filename-array.
        uint8_t *data = ( uint8_t* ) ( arg - sizeof ( struct MomfdContainer ) );
        struct MomfdContainer* info = ( struct MomfdContainer* ) data;

        // free the name-array, if present.
        if ( info->nNames ) {
            for ( int i = 0; i < info->nNames; ++i ) {
                delete[] info->ptr[ i ].s;
            }
        }

        // delete version/date/time strings
        IDL_STRING *str = reinterpret_cast<IDL_STRING*> ( arg );
        delete[] str[0].s;
        delete[] str[1].s;
        delete[] str[2].s;

        // free the data-block
        delete[] data;

    }


    void momfbd_help ( void ) {

        cout << "MOMFBD DLM usage:\n";
        cout << "    Syntax:          struct = momfbd_read('/path/to/file/data.momfbd',/KEYWORDS)\n";
        cout << "                     img = momfbg_mozaic( a.patch.img, a.patch(*,0).xl, a.patch(*,0).xh, a.patch(0,*).yl, a.patch(0,*).yh, /clip )\n";
        cout << "    Accepted Keywords:\n";
        cout << "        /HELP                    Print this info.\n";
        cout << "        /CHECK                   Just parse the file and show the structure, don't load data.\n";
        cout << "        /IMG                     Load Image data. (set by default)\n";
        cout << "        /PSF                     Load PSF. (if present in file)\n";
        cout << "        /OBJ                     Load convolved object data. (if present in file)\n";
        cout << "        /RES                     Load residuals. (if present in file)\n";
        cout << "        /ALPHA                   Load Alpha. (if present in file)\n";
        cout << "        /DIV                     Load Phase diversity data. (if present in file)\n";
        cout << "        /MODES                   Load Modes. (if present in file)\n";
        cout << "        /NAMES                   Load Filenames used in reconstruction. (if present in file)\n";
        cout << "        /ALL                     Load all data from file.\n";
        cout << "        VERBOSE={0,1,2}          Verbosity, default is 0 (no output)." << endl;

    }


    size_t appendTag ( vector<IDL_STRUCT_TAG_DEF>& tags, const char* name, IDL_MEMINT* dims, void* type, UCHAR flags = 0 ) {
        IDL_STRUCT_TAG_DEF tmp;
        tmp.name = const_cast<char*> ( name );
        tmp.dims = dims;
        tmp.type = type;
        tmp.flags = flags;
        tags.push_back ( tmp );
        return tags.size();
    }


    void createPatchTags ( vector<IDL_STRUCT_TAG_DEF>& tags, uint8_t loadMask, const FileMomfbd* const info ) {

        appendTag ( tags, "XL", 0, ( void* ) IDL_TYP_INT );
        appendTag ( tags, "XH", 0, ( void* ) IDL_TYP_INT );
        appendTag ( tags, "YL", 0, ( void* ) IDL_TYP_INT );
        appendTag ( tags, "YH", 0, ( void* ) IDL_TYP_INT );

        IDL_MEMINT* tmpDims = new IDL_MEMINT[2];
        tmpDims[0] = 1;
        tmpDims[1] = info->nChannels;
        appendTag ( tags, "DX", tmpDims, ( void* ) IDL_TYP_INT );

        tmpDims = new IDL_MEMINT[2];
        tmpDims[0] = 1;
        tmpDims[1] = info->nChannels;
        appendTag ( tags, "DY", tmpDims, ( void* ) IDL_TYP_INT );

        tmpDims = new IDL_MEMINT[2];
        tmpDims[0] = 1;
        tmpDims[1] = info->nChannels;
        appendTag ( tags, "NIMG", tmpDims, ( void* ) IDL_TYP_INT );

        const FileMomfbd::PatchInfo& tmpPatch = info->patches ( 0 );

        if ( tmpPatch.imgPos && ( loadMask & LOAD_IMG ) ) {
            tmpDims = new IDL_MEMINT[3];
            tmpDims[0] = 2;
            tmpDims[1] = tmpPatch.nPixelsX;
            tmpDims[2] = tmpPatch.nPixelsY;
            appendTag ( tags, "IMG", tmpDims, ( void* ) IDL_TYP_FLOAT );
        }
        if ( tmpPatch.npsf && ( loadMask & LOAD_PSF ) ) {
            tmpDims = new IDL_MEMINT[4];
            tmpDims[0] = 3;
            tmpDims[1] = tmpPatch.nPixelsX;
            tmpDims[2] = tmpPatch.nPixelsY;
            tmpDims[3] = tmpPatch.npsf;
            appendTag ( tags, "PSF", tmpDims, ( void* ) IDL_TYP_FLOAT );
        }
        if ( tmpPatch.nobj && ( loadMask & LOAD_OBJ ) ) {
            tmpDims = new IDL_MEMINT[4];
            tmpDims[0] = 3;
            tmpDims[1] = tmpPatch.nPixelsX;
            tmpDims[2] = tmpPatch.nPixelsY;
            tmpDims[3] = tmpPatch.nobj;
            appendTag ( tags, "OBJ", tmpDims, ( void* ) IDL_TYP_FLOAT );
        }
        if ( tmpPatch.nres && ( loadMask & LOAD_RES ) ) {
            tmpDims = new IDL_MEMINT[4];
            tmpDims[0] = 3;
            tmpDims[1] = tmpPatch.nPixelsX;
            tmpDims[2] = tmpPatch.nPixelsY;
            tmpDims[3] = tmpPatch.nres;
            appendTag ( tags, "RES", tmpDims, ( void* ) IDL_TYP_FLOAT );
        }
        if ( tmpPatch.nalpha && ( loadMask & LOAD_ALPHA ) ) {
            tmpDims = new IDL_MEMINT[3];
            tmpDims[0] = 2;
            tmpDims[1] = tmpPatch.nm;
            tmpDims[2] = tmpPatch.nalpha;
            appendTag ( tags, "ALPHA", tmpDims, ( void* ) IDL_TYP_FLOAT );
        }
        if ( tmpPatch.ndiv && ( loadMask & LOAD_DIV ) ) {
            tmpDims = new IDL_MEMINT[4];
            tmpDims[0] = 3;
            tmpDims[1] = tmpPatch.nphx;
            tmpDims[2] = tmpPatch.nphy;
            tmpDims[3] = tmpPatch.ndiv;
            appendTag ( tags, "DIV", tmpDims, ( void* ) IDL_TYP_FLOAT );
        }

        appendTag ( tags, 0, 0, 0 );           // tag-list has to be terminated with a "zero" entry

    }


    void createTags ( vector<IDL_STRUCT_TAG_DEF>& tags, uint8_t loadMask, const FileMomfbd* const info ) {

        appendTag ( tags, "VERSION", 0, ( void* ) IDL_TYP_STRING );
        appendTag ( tags, "TIME", 0, ( void* ) IDL_TYP_STRING );
        appendTag ( tags, "DATE", 0, ( void* ) IDL_TYP_STRING );

        IDL_MEMINT* tmpDims;
        if ( ( loadMask & LOAD_MODES ) && info->nModes ) {
            if ( info->nPH ) {
                tmpDims = new IDL_MEMINT[3];
                tmpDims[0] = 2;
                tmpDims[1] = info->nPH;
                tmpDims[2] = info->nPH;
                appendTag ( tags, "PUPIL", tmpDims, ( void* ) IDL_TYP_FLOAT );

                tmpDims = new IDL_MEMINT[4];
                tmpDims[0] = 3;
                tmpDims[1] = info->nPH;
                tmpDims[2] = info->nPH;
                tmpDims[3] = info->nModes;
                appendTag ( tags, "MODE", tmpDims, ( void* ) IDL_TYP_FLOAT );
            }
        }

        if ( loadMask && info->nChannels ) {
            tmpDims = new IDL_MEMINT[4];
            tmpDims[0] = 3;
            tmpDims[1] = info->nChannels;
            tmpDims[2] = 2;
            tmpDims[3] = 2;
            appendTag ( tags, "CLIP", tmpDims, ( void* ) IDL_TYP_INT );
        }

        if ( loadMask & LOAD_PATCH && info->nPatchesX > 0 && info->nPatchesY > 0 ) {
            vector<IDL_STRUCT_TAG_DEF> patchTags;
            createPatchTags ( patchTags, loadMask, info );
            // Append the "patch substructure" to the IDL structure
            IDL_StructDefPtr ds = IDL_MakeStruct ( 0, patchTags.data() );
            tmpDims = new IDL_MEMINT[3];
            tmpDims[0] = 2;
            tmpDims[1] = info->nPatchesX;
            tmpDims[2] = info->nPatchesY;
            appendTag ( tags, "PATCH", tmpDims, ( void* ) ds );
            for ( auto & it : patchTags ) { // delete the "dims" array for the tags that has them
                if ( it.dims ) {
                    delete[] it.dims;
                    it.dims = 0;
                }
            }

        }
        // Append file-names to the IDL structure
        if ( ( loadMask & LOAD_NAMES ) && info->nFileNames ) {
            tmpDims = new IDL_MEMINT[2];
            tmpDims[0] = 1;
            tmpDims[1] = info->nFileNames;
            appendTag ( tags, "NAME", tmpDims, ( void* ) IDL_TYP_STRING );
        }

        appendTag ( tags, 0, 0, 0 );           // tag-list has to be terminated with a "zero" entry

    }


    size_t getPatchSize ( const FileMomfbd* const info, uint8_t loadMask ) {

        size_t patchSize  = 4 * sizeof ( IDL_INT );                                 // xl,yl,xh,yl
        patchSize += 3 * info->nChannels * sizeof ( IDL_INT );                      // dx,dy,nim
        while ( patchSize % alignTo ) patchSize++;                                  // pad if not on boundary

        const FileMomfbd::PatchInfo& tmpPatch = info->patches ( 0 );

        size_t nFloats = 0;
        if ( ( loadMask & LOAD_IMG ) && tmpPatch.imgPos ) {
            nFloats += tmpPatch.nPixelsX * tmpPatch.nPixelsY;
        }
        if ( ( loadMask & LOAD_PSF ) && tmpPatch.npsf ) {
            nFloats += tmpPatch.npsf * tmpPatch.nPixelsX * tmpPatch.nPixelsY;
        }
        if ( ( loadMask & LOAD_OBJ ) && tmpPatch.nobj ) {
            nFloats += tmpPatch.nobj * tmpPatch.nPixelsX * tmpPatch.nPixelsY;
        }
        if ( ( loadMask & LOAD_RES ) && tmpPatch.nres ) {
            nFloats += tmpPatch.nPixelsX * tmpPatch.nPixelsY * tmpPatch.nres;
        }
        if ( ( loadMask & LOAD_ALPHA ) && tmpPatch.nalpha &&  tmpPatch.nm ) {
            nFloats += tmpPatch.nm * tmpPatch.nalpha * tmpPatch.nPixelsY;
        }
        if ( tmpPatch.ndiv && ( loadMask & LOAD_DIV ) ) {
            nFloats += tmpPatch.nphx * tmpPatch.nphy * tmpPatch.ndiv;
        }
        patchSize += nFloats * sizeof ( float );

        return patchSize;

    }

    void loadData ( ifstream& file, char* data, uint8_t loadMask, const FileMomfbd* const info, int verbosity ) {

        // zero the names-ptr
        MomfdContainer* container = reinterpret_cast<MomfdContainer*> ( data );
        container->nNames = 0;
        container->ptr = nullptr;

        char* ptr = data + sizeof ( MomfdContainer );
        float* fPtr;
        IDL_INT* tmpInt;

        // File version, date, time
        IDL_STRING *str = reinterpret_cast<IDL_STRING*> ( ptr );

        // Version string
        str[0].slen = info->versionString.length();
        str[0].s = new char[str->slen + 1];
        str[0].s[str[0].slen] = '\0';
        info->versionString.copy ( str[0].s, str[0].slen );
        str[0].stype = 0; //IDL_V_DYNAMIC;    // flag as dynamic (container will be deleted when destructed)

        // Time string
        str[1].slen = info->timeString.length();
        str[1].s = new char[str[1].slen + 1];
        str[1].s[str[1].slen] = '\0';
        info->timeString.copy ( str[1].s, str[1].slen );
        str[1].stype = 0; //IDL_V_DYNAMIC;    // flag as dynamic (container will be deleted when destructed)

        // Date string
        str[2].slen = info->dateString.length();
        str[2].s = new char[str[2].slen + 1];
        str[2].s[str[2].slen] = '\0';
        info->dateString.copy ( str[2].s, str[2].slen );
        str[2].stype = 0; //IDL_V_DYNAMIC;    // flag as dynamic (container will be deleted when destructed)
        
        ptr += 3*sizeof(IDL_STRING);

        file.clear();

        // pupil & modes
        if ( ( loadMask & LOAD_MODES ) && info->nModes ) {

            // pupil
            if ( info->phOffset ) {
                file.seekg ( info->phOffset );
                fPtr = reinterpret_cast<float*> ( ptr );
                ptr += readOrThrow ( file, fPtr, info->nPH * info->nPH, "MomfbdData:pupil" );
                if ( info->swapNeeded ) {
                    swapEndian ( fPtr, info->nPH * info->nPH );
                }
            }

            // modes
            if ( info->modesOffset ) {
                // Load data
                file.seekg ( info->modesOffset );
                fPtr = reinterpret_cast<float*> ( ptr );
                ptr += readOrThrow ( file, fPtr, info->nModes * info->nPH * info->nPH, "MomfbdData:modes" );
                if ( info->swapNeeded ) {
                    swapEndian ( fPtr, info->nPH * info->nPH );
                }
            }

        }

        if ( loadMask && info->nChannels ) {
            //NOTE: clipStartY & clipEndX swapped compared to file-order
            for ( int i = 0; i < info->nChannels; ++i ) {
                ( ( IDL_INT* ) ptr ) [ i + 0 * info->nChannels ] = info->clipStartX.get() [ i ];
                ( ( IDL_INT* ) ptr ) [ i + 1 * info->nChannels ] = info->clipStartY.get() [ i ];
                ( ( IDL_INT* ) ptr ) [ i + 2 * info->nChannels ] = info->clipEndX.get() [ i ];
                ( ( IDL_INT* ) ptr ) [ i + 3 * info->nChannels ] = info->clipEndY.get() [ i ];
            }
            ptr += 4 * info->nChannels * sizeof ( IDL_INT );
        }

        if ( loadMask & LOAD_PATCH && info->nPatchesX > 0 && info->nPatchesY > 0 ) {
            const FileMomfbd::PatchInfo* tmpPatch;
            size_t tmpSize;

            if ( verbosity > 1 ) {
                cout << "Total patches: " << info->nPatchesX << " x " << info->nPatchesY << endl;
            }

            for ( int y = 0; y < info->nPatchesY; ++y ) {
                for ( int x = 0; x < info->nPatchesX; ++x ) {

                    tmpPatch = info->patches.ptr ( x, y );
                    size_t nxny;
                    if ( tmpPatch ) {

                        if ( verbosity > 1 ) {
                            cout << "Loading patch (" << x << "," << y << ")   \r" << flush;
                        }
                        // region
                        tmpInt = ( IDL_INT* ) ptr;
                        for ( int i = 0; i < 4; ++i ) {
                            tmpInt[i] = tmpPatch->region[i];
                        }
                        ptr += 4 * sizeof ( IDL_INT );

                        tmpInt = ( IDL_INT* ) ptr;
                        for ( int i = 0; i < info->nChannels; ++i ) {
                            tmpInt[ i + 0 * info->nChannels ] = tmpPatch->dx.get() [ i ];
                            tmpInt[ i + 1 * info->nChannels ] = tmpPatch->dy.get() [ i ];
                            tmpInt[ i + 2 * info->nChannels ] = tmpPatch->nim.get() [ i ];
                        }
                        ptr += 3 * info->nChannels * sizeof ( IDL_INT );

                        while ( ( size_t ) ptr % alignTo ) ptr++;

                        nxny = tmpPatch->nPixelsX * tmpPatch->nPixelsX;
                        if ( ( loadMask & LOAD_IMG ) && tmpPatch->imgPos ) {
                            file.seekg ( tmpPatch->imgPos );
                            fPtr = reinterpret_cast<float*> ( ptr );
                            ptr += readOrThrow ( file, fPtr, nxny, "MomfbdPatch:img" );
                            if ( info->swapNeeded ) {
                                swapEndian ( fPtr, nxny );
                            }
                            // file data-order: [x,y]
                            // IDL has fast-index to the left, so a transpose is needed to keep index order [x,y]
                            transpose ( fPtr, tmpPatch->nPixelsY, tmpPatch->nPixelsX );
                        }

                        if ( ( loadMask & LOAD_PSF ) && tmpPatch->npsf ) {
                            file.seekg ( tmpPatch->psfPos );
                            fPtr = reinterpret_cast<float*> ( ptr );
                            ptr += readOrThrow ( file, fPtr, tmpPatch->npsf * nxny, "MomfbdPatch:psf" );
                            if ( info->swapNeeded ) {
                                swapEndian ( fPtr, tmpPatch->npsf * nxny );
                            }
                            // file data-order: [x,y]
                            // IDL has fast-index to the left, so a transpose is needed to keep index order [x,y]
                            for ( int i = 0; i < tmpPatch->npsf; ++i )
                                transpose ( ( fPtr + i * nxny ), tmpPatch->nPixelsY, tmpPatch->nPixelsX );
                        }

                        if ( ( loadMask & LOAD_OBJ ) && tmpPatch->nobj ) {
                            file.seekg ( tmpPatch->objPos );
                            fPtr = reinterpret_cast<float*> ( ptr );
                            ptr += readOrThrow ( file, fPtr, tmpPatch->nobj * nxny, "MomfbdPatch:obj" );
                            if ( info->swapNeeded ) {
                                swapEndian ( fPtr, tmpPatch->nobj * nxny );
                            }
                            // file data-order: [x,y]
                            // IDL has fast-index to the left, so a transpose is needed to keep index order [x,y]
                            for ( int i = 0; i < tmpPatch->nobj; ++i )
                                transpose ( ( fPtr + i * nxny ), tmpPatch->nPixelsY, tmpPatch->nPixelsX );
                        }

                        if ( ( loadMask & LOAD_RES ) && tmpPatch->nres ) {
                            file.seekg ( tmpPatch->resPos );
                            fPtr = reinterpret_cast<float*> ( ptr );
                            ptr += readOrThrow ( file, fPtr, tmpPatch->nres * nxny, "MomfbdPatch:res" );
                            if ( info->swapNeeded ) {
                                swapEndian ( fPtr, tmpPatch->nres * nxny );
                            }
                            // file data-order: [x,y]
                            // IDL has fast-index to the left, so a transpose is needed to keep index order [x,y]
                            for ( int i = 0; i < tmpPatch->nres; ++i )
                                transpose ( ( fPtr + i * nxny ), tmpPatch->nPixelsY, tmpPatch->nPixelsX );
                        }

                        if ( ( loadMask & LOAD_ALPHA ) && tmpPatch->nalpha ) {
                            tmpSize = tmpPatch->nalpha * tmpPatch->nm;
                            file.seekg ( tmpPatch->alphaPos );
                            fPtr = reinterpret_cast<float*> ( ptr );
                            ptr += readOrThrow ( file, fPtr, tmpSize, "MomfbdPatch:alpha" );
                            if ( info->swapNeeded ) {
                                swapEndian ( fPtr, tmpSize );
                            }
                            // No transpose is needed to get index order [nmodes,nalpha]
                        }

                        if ( ( loadMask & LOAD_DIV ) && tmpPatch->ndiv ) {
                            tmpSize = tmpPatch->ndiv * tmpPatch->nphy * tmpPatch->nphx;
                            file.seekg ( tmpPatch->diversityPos );
                            fPtr = reinterpret_cast<float*> ( ptr );
                            ptr += readOrThrow ( file, fPtr, tmpSize, "MomfbdPatch:res" );
                            if ( info->swapNeeded ) {
                                swapEndian ( fPtr, tmpSize );
                            }
                            // file data-order: [x,y]
                            // IDL has fast-index to the left, so a transpose is needed to keep index order [x,y]
                            tmpSize = tmpPatch->nphy * tmpPatch->nphx;
                            for ( int i = 0; i < tmpPatch->ndiv; ++i )
                                transpose ( ( fPtr + i * tmpSize ), tmpPatch->nphy, tmpPatch->nphx );

                        }

                    }   //  if( tmpPatch )

                }   // x-loop

            }   // y-loop

        }   // if(LOAD_PATCH)


        if ( ( loadMask & LOAD_NAMES ) && info->filenameOffset ) {
            file.seekg ( info->filenameOffset );

            int32_t nameLength;
            vector<char> tmpStr ( 1024, 0 );
            int i = 0;
            str = reinterpret_cast<IDL_STRING*> ( ptr );
            container->ptr = str;

            while ( i < info->nFileNames ) {
                readOrThrow ( file, &nameLength, 1, "MomfbdData:namelength #" + to_string ( i ) );
                if ( info->swapNeeded ) swapEndian ( nameLength );
                tmpStr.reserve ( nameLength );
                readOrThrow ( file, & ( tmpStr[0] ), nameLength,  "MomfbdData:name #" + to_string ( i ) );
                //namesPtr->push_back( string( tmpStr.begin(), tmpStr.begin() + nameLength ) );
                str[i].slen = nameLength;
                str[i].s = new char[nameLength + 1]; //const_cast<char*>(namesPtr->rbegin()->c_str());
                memcpy ( str[i].s, tmpStr.data(), nameLength );
                str[i].s[nameLength] = 0;
                str[i].stype = 0;

                i++;
            }
            container->nNames = i;

        }


    }

    Array<float> getWeights ( const Overlaps& overlaps, int margin, int ny, int nx ) {

        static std::map<Overlaps, Array<float> > map;
        typedef std::map<Overlaps, Array<float>>::iterator mit;

        mit it = map.find ( overlaps );
        if ( it != map.end() ) {
            return it->second;
        }

        Array<float> tmp ( ny, ny );
        memset ( tmp.ptr(), 0, nx * ny * sizeof ( float ) );

        for ( int x = margin; x < nx - margin; ++x ) {
            for ( int y = margin; y < ny - margin; ++y ) {
                tmp ( y, x ) = 1;
            }
        }

        std::shared_ptr<double> blend;
        int sz = 0;

        if ( overlaps.north > 0 ) {
            blend.reset ( new double[overlaps.north], [] ( double * p ) {
                delete[] p;
            } );
            sz = overlaps.north;
            hann ( blend.get(), overlaps.north );
            for ( int x = margin; x < nx - margin; ++x ) {
                for ( int y = 0; y < overlaps.north; ++y ) {
                    tmp ( ny - margin - y - 1, x ) *= blend.get() [y];
                }
            }
        }

        if ( overlaps.south > 0 ) {
            if ( overlaps.south != sz ) {
                blend.reset ( new double[overlaps.south], [] ( double * p ) {
                    delete[] p;
                } );
                sz = overlaps.south;
                hann ( blend.get(), overlaps.south );
            }
            for ( int x = margin; x < nx - margin; ++x ) {
                for ( int y = 0; y < overlaps.south; ++y ) {
                    tmp ( margin + y, x ) *= blend.get() [y];
                }
            }
        }

        if ( overlaps.east > 0 ) {
            if ( overlaps.east != overlaps.south ) {
                blend.reset ( new double[overlaps.east], [] ( double * p ) {
                    delete[] p;
                } );
                sz = overlaps.east;
                hann ( blend.get(), overlaps.east );
            }
            for ( int y = margin; y < ny - margin; ++y ) {
                for ( int x = 0; x < overlaps.east; ++x ) {
                    tmp ( y, nx - margin - x - 1 ) *= blend.get() [x];
                }
            }
        }

        if ( overlaps.west > 0 ) {
            if ( overlaps.west != sz ) {
                blend.reset ( new double[overlaps.west], [] ( double * p ) {
                    delete[] p;
                } );
                hann ( blend.get(), overlaps.west );
            }
            for ( int y = margin; y < ny - margin; ++y ) {
                for ( int x = 0; x < overlaps.west; ++x ) {
                    tmp ( y, margin + x ) *= blend.get() [x];
                }
            }
        }

        map.insert ( std::make_pair ( overlaps, tmp ) );
        return tmp;

    }


}   // anon namespace

//
// Routine to read .momfbd files.
// ex. a=momfbd_read('/scratch/roald/sims/momfbd/r0_5_new/image_t.-100.lc0.80..116.momfbd',verbose=2,/all)
//
IDL_VPTR redux::momfbd_read ( int argc, IDL_VPTR* argv, char* argk ) {

    if ( argc < 1 ) {
        cout << "MOMFBD_MOZAIC: takes 1 argument." << endl;
        momfbd_help();
        return IDL_Gettmp();
    }

    IDL_VPTR fileName = argv[0];

    KW_RESULT kw;
    kw.verbose = 0;
    kw.help = 0;
    kw.img = 1;
    ( void ) IDL_KWProcessByOffset ( argc, argv, argk, kw_pars, ( IDL_VPTR* ) 0, 255, &kw );

    if ( kw.help ) {
        momfbd_help();
        return IDL_Gettmp();
    }

    IDL_ENSURE_SIMPLE ( fileName );
    IDL_ENSURE_STRING ( fileName );

    char *name = IDL_VarGetString ( fileName );
    int verbosity = std::min ( std::max ( ( int ) kw.verbose, 0 ), 8 );
    int checkData = kw.check;
    uint8_t loadMask = 0;
    if ( kw.img )   loadMask |= LOAD_IMG;
    if ( kw.psf )   loadMask |= LOAD_PSF;
    if ( kw.obj )   loadMask |= LOAD_OBJ;
    if ( kw.res )   loadMask |= LOAD_RES;
    if ( kw.alpha ) loadMask |= LOAD_ALPHA;
    if ( kw.div )   loadMask |= LOAD_DIV;
    if ( kw.modes ) loadMask |= LOAD_MODES;
    if ( kw.names ) loadMask |= LOAD_NAMES;
    if ( kw.all )   loadMask |= LOAD_ALL;

    IDL_KW_FREE;

    if ( checkData && !verbosity ) verbosity = 1;            // When checking we want some output...

    if ( verbosity > 0 ) {
        cout << "Loading file: \"" << name << "\"" << endl;
    }

    ifstream file;
    std::shared_ptr<FileMomfbd> info ( new FileMomfbd() );

    try {
        file.open ( name );
        info->read ( file );
    } catch ( exception& e ) {
        cout << "Failed to read info from Momfbd file: " << name << endl << "Reason: " << e.what() << endl;
        return IDL_Gettmp();
    }

    vector<IDL_STRUCT_TAG_DEF> tags;
    size_t patchSize;
    if ( checkData ) {
        createTags ( tags, LOAD_ALL, info.get() );
        patchSize = 0;      // no need to allocate a big chunk of memory just to print the structure.
    } else {
        createTags ( tags, loadMask, info.get() );
        patchSize = getPatchSize ( info.get(), loadMask );
    }
    IDL_StructDefPtr myStruct = IDL_MakeStruct ( 0, tags.data() );              // Generate the IDL structure defined above
    // Clean up the "dims" array for the tags that has them
    for ( auto & it : tags ) {
        if ( it.dims ) {
            delete[] it.dims;
            it.dims = 0;
        }
    }

    // Calculate size of data to load.
    size_t totalSize = 3 * sizeof ( IDL_STRING );                               // VERSION - TIME - DATE
    totalSize += 4 * info->nChannels * sizeof ( IDL_INT );                      // clip-values for each channel
    if ( ( loadMask & LOAD_MODES ) && info->nModes ) {
        totalSize += info->nPH * info->nPH * sizeof ( float );                  // Pupil Data
        totalSize += info->nModes * info->nPH * info->nPH * sizeof ( float );   // Mode Data
    }
    totalSize += info->nPatchesX * info->nPatchesY * patchSize;
    while ( totalSize % alignTo ) totalSize++;                                  // pad if not on boundary

    if ( info->nFileNames && ( loadMask & LOAD_NAMES ) ) {                          // if filenames are stored and loaded
        totalSize += info->nFileNames * sizeof ( IDL_STRING );
    }

    // Allocate the datablock needed.
    std::unique_ptr<char> data ( new char [ totalSize + sizeof ( MomfdContainer )] ); //, [](char *p) { delete[] p; });
    loadData ( file, data.get(), loadMask, info.get(), verbosity );

    // Import the datablock to the structure specified above.
    IDL_MEMINT dims[] = {1};
    IDL_VPTR v = IDL_ImportArray ( 1, dims, IDL_TYP_STRUCT, ( UCHAR* ) data.get() + sizeof ( MomfdContainer ), cleanupMomfbd, myStruct );

    
    // Dump structure layout if requested
    if ( verbosity > 1 ) {
        dumpStruct ( v, -1, 2 );
    }

    if ( checkData ) {
        return IDL_Gettmp();        // "data" will be freed on return
    }

    // release the datablock from the RAII container to prevent de-allocation on return.
    data.release();
    return v;


}


void img_clip ( Array<float>& img ) {

    int imgSizeY = img.dimSize(0);
    int imgSizeX = img.dimSize(1);

    vector<double> colSums(imgSizeX,0);
    vector<double> rowSums(imgSizeY,0);
    
    size_t firstX(0), lastX(imgSizeX-1), firstY(0), lastY(imgSizeY-1);

    for( int x=0; x < imgSizeX; ++x ) {
        for( int y=0; y < imgSizeY; ++y ) {
            colSums[x] += img(y,x);
            rowSums[y] += img(y,x);
        }
    }
    
    while(firstX < lastX && !colSums[firstX]) ++firstX;
    while(lastX && !colSums[lastX]) --lastX;
    while(firstY < lastY && !rowSums[firstY]) ++firstY;
    while(lastY && !rowSums[lastY]) --lastY;
    
    if(firstX < lastX && firstY < lastY) {
        img.setLimits(firstY,lastY,firstX,lastX);
        Array<float> tmp(img,firstY,lastY,firstX,lastX);
        img = tmp;
    }

}

IDL_VPTR redux::momfbd_mozaic ( int argc, IDL_VPTR *argv, char *argk ) {

    if ( argc < 5 ) {
        cout << "MOMFBD_MOZAIC: takes 5 arguments." << endl;
        momfbd_help();
        return IDL_Gettmp();
    }


    IDL_VPTR img_in = argv[0];
    int patchSizeX  = img_in->value.s.arr->dim[0];
    int patchSizeY  = img_in->value.s.arr->dim[1];
    int nPatchesX = img_in->value.s.arr->dim[2];
    int nPatchesY = img_in->value.s.arr->dim[3];

    KW_RESULT kw;
    kw.margin = std::max ( patchSizeX, patchSizeY ) / 8; // number of pixels to cut from the edges of each path,
                                                         // default is 12.5% (to conform with the old version)
    kw.verbose = 0;
    ( void ) IDL_KWProcessByOffset ( argc, argv, argk, kw_pars, ( IDL_VPTR* ) 0, 1, &kw );
    int do_clip = kw.clip;
    int margin = kw.margin;
    int transpose = kw.transpose;
    int verbosity = std::min ( std::max ( (int)kw.verbose, 0 ), 8 );
    IDL_KW_FREE;

    if ( verbosity > 0 ) {
        cout << "Mozaic:  nPatches = (" << nPatchesX << "," << nPatchesY << ")" << endl;
        cout << "        patchSize = (" << patchSizeX << "," << patchSizeY << ")" << endl;
    }

    
    if ( margin > patchSizeX >> 1 || margin > patchSizeY >> 1 ) {
        cout << "Margin is too big, nothing will be left." << endl;
        return IDL_Gettmp();
    }

    int16_t* patchesFirstX = ( int16_t* ) ( argv[1]->value.s.arr->data );
    int16_t* patchesLastX = ( int16_t* ) ( argv[2]->value.s.arr->data );
    int16_t* patchesFirstY = ( int16_t* ) ( argv[3]->value.s.arr->data );
    int16_t* patchesLastY = ( int16_t* ) ( argv[4]->value.s.arr->data );
    
    if ( verbosity > 0 ) {
        cout << "      " << printArray(patchesFirstX,nPatchesX,"firstPixelX") << endl;
        cout << "       " << printArray(patchesLastX,nPatchesX,"lastPixelX") << endl;
        cout << "      " << printArray(patchesFirstY,nPatchesY,"firstPixelY") << endl;
        cout << "       " << printArray(patchesLastY,nPatchesY,"lastPixelY") << endl;
    }

    Array<struct Overlaps> overlaps ( nPatchesY, nPatchesX ); // How many pixels of overlap between the patches, in order N,E,S,W
    // calculate the overlaps.
    for ( int y = 0; y < nPatchesY; ++y ) {
        for ( int x = 0; x < nPatchesX; ++x ) {
            overlaps ( y, x ) = Overlaps (
                                    y < ( nPatchesY - 1 ) ? patchesLastY[y] - patchesFirstY[y + 1] + 1 - 2 * margin : 0,    // north
                                    x < ( nPatchesX - 1 ) ? patchesLastX[x] - patchesFirstX[x + 1] + 1 - 2 * margin : 0,    // east
                                    y ? patchesLastY[y - 1] - patchesFirstY[y] + 1 - 2 * margin : 0,                        // south
                                    x ? patchesLastX[x - 1] - patchesFirstX[x] + 1 - 2 * margin : 0                         // west
                                );
        }
    }


    Array<float> img ( ( float* ) ( img_in->value.s.arr->data ), nPatchesY, nPatchesX, patchSizeY, patchSizeX );


    int imgFirstX = patchesFirstX[0];
    int imgLastX = patchesLastX[0];
    for ( int x = 0; x < nPatchesX; ++x ) {
        if ( patchesFirstX[x] < imgFirstX ) imgFirstX = patchesFirstX[x];
        if ( patchesLastX[x] > imgLastX ) imgLastX = patchesLastX[x];
    }
    for ( int x = 0; x < nPatchesX; ++x ) {
        patchesFirstX[x] -= imgFirstX;
        patchesLastX[x] -= imgFirstX;
    }
    imgLastX -= imgFirstX;
    imgFirstX = 0;
    int imgFirstY = patchesFirstY[0];
    int imgLastY = patchesLastY[0];
    for ( int y = 0; y < nPatchesY; ++y ) {
        if ( patchesFirstY[y] < imgFirstY ) imgFirstY = patchesFirstY[y];
        if ( patchesLastY[y] > imgLastY ) imgLastY = patchesLastY[y];
    }
    for ( int y = 0; y < nPatchesY; ++y ) {
        patchesFirstY[y] -= imgFirstY;
        patchesLastY[y] -= imgFirstY;
    }
    imgLastY -= imgFirstY;
    imgFirstY = 0;
    //
    int imgSizeX = imgLastX - imgFirstX + 1;
    int imgSizeY = imgLastY - imgFirstY + 1;

    Array<float> pic ( imgSizeY, imgSizeX );
    memset ( pic.ptr ( 0 ), 0, imgSizeX * imgSizeY * sizeof ( float ) );

    for ( int yp = 0; yp < nPatchesY; ++yp ) {
        for ( int xp = 0; xp < nPatchesX; ++xp ) {
            Array<float> patch ( img, yp, yp, xp, xp, 0, patchSizeY-1, 0, patchSizeX-1 );
            Array<float> weights = getWeights ( overlaps ( yp, xp ), margin, patchSizeY, patchSizeX );
            for ( int y = patchesFirstY[yp] + margin; y <= patchesLastY[yp] - margin; ++y ) {
                for ( int x = patchesFirstX[xp] + margin; x <= patchesLastX[xp] - margin; ++x ) {
                    pic ( y, x ) += patch ( 0, 0, y - patchesFirstY[yp], x - patchesFirstX[xp] ) * weights ( y - patchesFirstY[yp], x - patchesFirstX[xp] );
                }
            }
        }
    }

    if ( transpose ) {
        if ( verbosity > 0 ) {
            cout << "       Transposing image." << endl;
        }
        redux::util::transpose( pic.ptr(0), pic.dimSize(0), pic.dimSize(1) );
        pic.permuteDimensions({0,1});
       // std::swap(dims[0],dims[1]);
    }

    if ( do_clip ) {
        if ( verbosity > 0 ) {
            cout << "       Clipping image." << endl;
        }
        img_clip ( pic );
    }

    IDL_MEMINT dims[] = { (int)pic.dimSize(1), (int)pic.dimSize(0) };
    IDL_VPTR v = IDL_ImportArray ( 2, dims, IDL_TYP_FLOAT, ( UCHAR* ) ( pic.cloneData() ), cleanupFloat, 0 );

    return v;

}




extern "C" {

    int IDL_Load ( void ) {

        static IDL_SYSFUN_DEF2 function_addr[] = {
            { { ( IDL_VPTR ( * ) () ) redux::momfbd_read}, ( char* ) "MOMFBD_READ", 0, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
            { { ( IDL_VPTR ( * ) () ) redux::momfbd_mozaic}, ( char* ) "MOMFBD_MOZAIC", 0, 5, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
        };

        // static IDL_SYSFUN_DEF2 procedure_addr[] = {
        //    { { ( IDL_SYSRTN_GENERIC ) redux::fzread}, ( char* ) "FZREAD", 0, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
        //    { { ( IDL_SYSRTN_GENERIC ) redux::fzwrite}, ( char* ) "FZWRITE", 0, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
        //    { { ( IDL_SYSRTN_GENERIC ) redux::fcwrite}, ( char* ) "FCWRITE", 0, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
        //   };

        /* Register our routine. The routines must be specified exactly the same as in testmodule.dlm. */
        return IDL_SysRtnAdd ( function_addr, TRUE, IDL_CARRAY_ELTS ( function_addr ) );
        // &&  IDL_SysRtnAdd ( procedure_addr, FALSE, IDL_CARRAY_ELTS ( procedure_addr ) );

    }

}

