#include "momfbd.hpp"

#include "redux/file/filemomfbd.hpp"
#include "redux/math/functions.hpp"
#include "redux/util/array.hpp"
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
                cout << "           TAG           TYPE                     OFFSET         SIZE           " << endl;
                sz += dumpStruct ( data, indent, indent );
                cout << string ( 45, ' ' ) << "Total Size:         " << to_string ( sz ) << endl;
                return sz;
            }


            IDL_StructDefPtr structDef = data->value.s.sdef;
            //uint8_t *buf = data->value.s.arr->data;
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

            cout.setf ( std::ios::left );
            if ( v->type == IDL_TYP_STRUCT ) {
                cout << std::setw(25) << (string(current, ' ') + name);
                cout << std::setw(25) << type;
                //cout.setf ( std::ios::hex );
                //cout << std::setw(20) << (void*)( buf + offset );
                //cout.setf ( std::ios::dec );
                cout << std::setw(15) << to_string ( (size_t)offset ) << endl;
                sz += count * dumpStruct ( v, current + indent, indent );
            } else {
                sz += count * var_sizes[v->type];
                cout << std::setw(25) << (string(current,' ')+name);
                cout << std::setw(25) << type;
                //cout.setf ( std::ios::hex );
                //cout << std::setw(20) << (void*)( buf + offset );
                //cout.setf ( std::ios::dec );
                cout << std::setw(15) << to_string ( (size_t)offset );
                cout << std::setw(15) << to_string ( count * var_sizes[v->type] ) << endl;
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
        { ( char* ) "IMG",           IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( img ) },
        { ( char* ) "MARGIN",        IDL_TYP_INT, 1, IDL_KW_ZERO,                 0, ( char* ) IDL_KW_OFFSETOF ( margin ) },
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
        cout << "    Syntax example:    a = momfbd_read('/path/to/file/data.momfbd',/img,verbose=2)\n";
        cout << "                     img = momfbd_mozaic( a.patch.img, a.patch(*,0).xl, a.patch(*,0).xh, a.patch(0,*).yl, a.patch(0,*).yh, /clip )\n";
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

        appendTag ( tags, "YL", 0, ( void* ) IDL_TYP_INT );                 // fast index first => in the IDL struct X & Y are interchanged compared to the c++ names
        appendTag ( tags, "YH", 0, ( void* ) IDL_TYP_INT );
        appendTag ( tags, "XL", 0, ( void* ) IDL_TYP_INT );
        appendTag ( tags, "XH", 0, ( void* ) IDL_TYP_INT );
        if ( info->version >= 20110714.0) {
            appendTag ( tags, "OFFY", 0, ( void* ) IDL_TYP_INT );
            appendTag ( tags, "OFFX", 0, ( void* ) IDL_TYP_INT );
        }
        
        IDL_MEMINT* tmpDims = new IDL_MEMINT[2];
        tmpDims[0] = 1;
        tmpDims[1] = info->nChannels;
        appendTag ( tags, "DY", tmpDims, ( void* ) IDL_TYP_INT );           // fast index first => in the IDL struct X & Y are interchanged compared to the c++ names   

        tmpDims = new IDL_MEMINT[2];
        tmpDims[0] = 1;
        tmpDims[1] = info->nChannels;
        appendTag ( tags, "DX", tmpDims, ( void* ) IDL_TYP_INT );

        tmpDims = new IDL_MEMINT[2];
        tmpDims[0] = 1;
        tmpDims[1] = info->nChannels;
        appendTag ( tags, "NIMG", tmpDims, ( void* ) IDL_TYP_INT );

        const FileMomfbd::PatchInfo& tmpPatch = info->patches ( 0 );

        if ( tmpPatch.imgPos && ( loadMask & MOMFBD_IMG ) ) {
            tmpDims = new IDL_MEMINT[3];
            tmpDims[0] = 2;
            tmpDims[1] = tmpPatch.nPixelsY;     // fast index first
            tmpDims[2] = tmpPatch.nPixelsX;
            appendTag ( tags, "IMG", tmpDims, ( void* ) IDL_TYP_FLOAT );
        }
        if ( tmpPatch.npsf && ( loadMask & MOMFBD_PSF ) ) {
            tmpDims = new IDL_MEMINT[4];
            tmpDims[0] = 3;
            tmpDims[1] = tmpPatch.nPixelsY;     // fast index first
            tmpDims[2] = tmpPatch.nPixelsX;
            tmpDims[3] = tmpPatch.npsf;
            appendTag ( tags, "PSF", tmpDims, ( void* ) IDL_TYP_FLOAT );
        }
        if ( tmpPatch.nobj && ( loadMask & MOMFBD_OBJ ) ) {
            tmpDims = new IDL_MEMINT[4];
            tmpDims[0] = 3;
            tmpDims[1] = tmpPatch.nPixelsY;     // fast index first
            tmpDims[2] = tmpPatch.nPixelsX;
            tmpDims[3] = tmpPatch.nobj;
            appendTag ( tags, "OBJ", tmpDims, ( void* ) IDL_TYP_FLOAT );
        }
        if ( tmpPatch.nres && ( loadMask & MOMFBD_RES ) ) {
            tmpDims = new IDL_MEMINT[4];
            tmpDims[0] = 3;
            tmpDims[1] = tmpPatch.nPixelsY;     // fast index first
            tmpDims[2] = tmpPatch.nPixelsX;
            tmpDims[3] = tmpPatch.nres;
            appendTag ( tags, "RES", tmpDims, ( void* ) IDL_TYP_FLOAT );
        }
        if ( tmpPatch.nalpha && ( loadMask & MOMFBD_ALPHA ) ) {
            tmpDims = new IDL_MEMINT[3];
            tmpDims[0] = 2;
            tmpDims[1] = tmpPatch.nm;
            tmpDims[2] = tmpPatch.nalpha;
            appendTag ( tags, "ALPHA", tmpDims, ( void* ) IDL_TYP_FLOAT );
        }
        if ( tmpPatch.ndiv && ( loadMask & MOMFBD_DIV ) ) {
            tmpDims = new IDL_MEMINT[4];
            tmpDims[0] = 3;
            tmpDims[1] = tmpPatch.nphy;     // fast index first
            tmpDims[2] = tmpPatch.nphx;
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

        if ( loadMask && info->nChannels ) {
            tmpDims = new IDL_MEMINT[4];
            tmpDims[0] = 3;
            tmpDims[1] = info->nChannels;
            tmpDims[2] = 2;
            tmpDims[3] = 2;
            appendTag ( tags, "CLIP", tmpDims, ( void* ) IDL_TYP_INT );
        }
        
        if ( loadMask & MOMFBD_MODES ) {
            if ( info->version >= 20110714.0 ) {
        cout << "BLÖBLÖ" << endl;
                appendTag ( tags, "PIX2CF", 0, ( void* ) IDL_TYP_FLOAT );
                appendTag ( tags, "CF2PIX", 0, ( void* ) IDL_TYP_FLOAT );
            }
            if ( info->nPH ) {
                tmpDims = new IDL_MEMINT[3];
                tmpDims[0] = 2;
                tmpDims[1] = info->nPH;
                tmpDims[2] = info->nPH;
                appendTag ( tags, "PUPIL", tmpDims, ( void* ) IDL_TYP_FLOAT );

                if ( info->nModes ) {
                    tmpDims = new IDL_MEMINT[4];
                    tmpDims[0] = 3;
                    tmpDims[1] = info->nPH;
                    tmpDims[2] = info->nPH;
                    tmpDims[3] = info->nModes;
                    appendTag ( tags, "MODE", tmpDims, ( void* ) IDL_TYP_FLOAT );
                }
            }
        }

        if ( loadMask & MOMFBD_PATCH && info->nPatchesX > 0 && info->nPatchesY > 0 ) {
            vector<IDL_STRUCT_TAG_DEF> patchTags;
            createPatchTags ( patchTags, loadMask, info );
            // Append the "patch substructure" to the IDL structure
            IDL_StructDefPtr ds = IDL_MakeStruct ( 0, patchTags.data() );
            tmpDims = new IDL_MEMINT[3];
            tmpDims[0] = 2;
            tmpDims[1] = info->nPatchesY;       // fast index first
            tmpDims[2] = info->nPatchesX;
            appendTag ( tags, "PATCH", tmpDims, ( void* ) ds );
            for ( auto & it : patchTags ) { // delete the "dims" array for the tags that has them
                if ( it.dims ) {
                    delete[] it.dims;
                    it.dims = 0;
                }
            }

        }
        // Append file-names to the IDL structure
        if ( ( loadMask & MOMFBD_NAMES ) && info->nFileNames ) {
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
        if ( ( loadMask & MOMFBD_IMG ) && tmpPatch.imgPos ) {
            nFloats += tmpPatch.nPixelsX * tmpPatch.nPixelsY;
        }
        if ( ( loadMask & MOMFBD_PSF ) && tmpPatch.npsf ) {
            nFloats += tmpPatch.npsf * tmpPatch.nPixelsX * tmpPatch.nPixelsY;
        }
        if ( ( loadMask & MOMFBD_OBJ ) && tmpPatch.nobj ) {
            nFloats += tmpPatch.nobj * tmpPatch.nPixelsX * tmpPatch.nPixelsY;
        }
        if ( ( loadMask & MOMFBD_RES ) && tmpPatch.nres ) {
            nFloats += tmpPatch.nPixelsX * tmpPatch.nPixelsY * tmpPatch.nres;
        }
        if ( ( loadMask & MOMFBD_ALPHA ) && tmpPatch.nalpha &&  tmpPatch.nm ) {
            nFloats += tmpPatch.nm * tmpPatch.nalpha * tmpPatch.nPixelsY;
        }
        if ( tmpPatch.ndiv && ( loadMask & MOMFBD_DIV ) ) {
            nFloats += tmpPatch.nphx * tmpPatch.nphy * tmpPatch.ndiv;
        }
        patchSize += nFloats * sizeof ( float );

        return patchSize;

    }

    
    void loadData ( ifstream& file, char* data, uint8_t loadMask, FileMomfbd* info, int verbosity ) {

        // zero the list of filenames
        MomfdContainer* container = reinterpret_cast<MomfdContainer*> ( data );
        container->nNames = 0;
        container->ptr = nullptr;

        char* ptr = data + sizeof ( MomfdContainer );

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
        
        
        if ( loadMask && info->nChannels ) {
            //NOTE: clipStartY & clipEndX swapped compared to file-order
            for ( int i = 0; i < info->nChannels; ++i ) {
                ( ( IDL_INT* ) ptr ) [ i + 0 * info->nChannels ] = info->clipStartY.get() [ i ];
                ( ( IDL_INT* ) ptr ) [ i + 1 * info->nChannels ] = info->clipStartX.get() [ i ];
                ( ( IDL_INT* ) ptr ) [ i + 2 * info->nChannels ] = info->clipEndY.get() [ i ];
                ( ( IDL_INT* ) ptr ) [ i + 3 * info->nChannels ] = info->clipEndX.get() [ i ];
            }
            ptr += 4 * info->nChannels * sizeof ( IDL_INT );
        }

        if ( ( loadMask & MOMFBD_MODES ) && info->version >= 20110714.0 ) {
        cout << "BLÄBLÄ" << endl;
            *((float*)ptr) = info->pix2cf;
            ptr += sizeof(float);
            *((float*)ptr) = info->cf2pix;
            ptr += sizeof(float);
        }
        
        // Load the data
        ptr = info->load( file, ptr, loadMask, verbosity, 4 );
        
        // add list of filenames, if requested
        if( loadMask & MOMFBD_NAMES ) {
            int i = 0;
            str = reinterpret_cast<IDL_STRING*> ( ptr );
            container->ptr = str;
            for(auto &fn: info->fileNames) {
                size_t nameLength = fn.length();
                str[i].slen = nameLength;
                str[i].s = new char[nameLength + 1];
                strncpy ( str[i].s, fn.c_str(), nameLength );
                str[i].s[nameLength] = 0;
                str[i++].stype = 0;
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
        cout << "MOMFBD_READ: takes 1 argument." << endl;
        momfbd_help();
        return IDL_GettmpInt(0);    // return a dummy
    }

    IDL_VPTR fileName = argv[0];

    KW_RESULT kw;
    kw.verbose = 0;
    kw.help = 0;
    ( void ) IDL_KWProcessByOffset ( argc, argv, argk, kw_pars, ( IDL_VPTR* ) 0, 255, &kw );

    if ( kw.help ) {
        momfbd_help();
        return IDL_GettmpInt(0);    // return a dummy
    }

    IDL_ENSURE_SIMPLE ( fileName );
    IDL_ENSURE_STRING ( fileName );

    char *name = IDL_VarGetString ( fileName );
    int verbosity = std::min ( std::max ( ( int ) kw.verbose, 0 ), 8 );
    int checkData = kw.check;
    
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
        return IDL_GettmpInt(-1);    // return a dummy
    }
        
    
    uint8_t loadMask = 0;
    
    if ( kw.img && (info->dataMask & MOMFBD_IMG) )     loadMask |= MOMFBD_IMG;
    if ( kw.psf && (info->dataMask & MOMFBD_PSF) )     loadMask |= MOMFBD_PSF;
    if ( kw.obj && (info->dataMask & MOMFBD_OBJ) )     loadMask |= MOMFBD_OBJ;
    if ( kw.res && (info->dataMask & MOMFBD_RES) )     loadMask |= MOMFBD_RES;
    if ( kw.alpha && (info->dataMask & MOMFBD_ALPHA) ) loadMask |= MOMFBD_ALPHA;
    if ( kw.div && (info->dataMask & MOMFBD_DIV) )     loadMask |= MOMFBD_DIV;
    if ( kw.modes && (info->dataMask & MOMFBD_MODES) ) loadMask |= MOMFBD_MODES;
    if ( kw.names && (info->dataMask & MOMFBD_NAMES) ) loadMask |= MOMFBD_NAMES;
    if ( kw.all )   loadMask = info->dataMask;

    IDL_KW_FREE;
    
    if( !loadMask && !checkData ) {
        if ( info->dataMask & MOMFBD_IMG ) loadMask = info->dataMask;
    }

    IDL_MEMINT dims[] = {1};
    IDL_VPTR v;
    
    if ( checkData ) {
        vector<IDL_STRUCT_TAG_DEF> allTags;
        createTags ( allTags, info->dataMask, info.get() );
        IDL_StructDefPtr allStruct = IDL_MakeStruct ( 0, allTags.data() );
        // Clean up the "dims" array for the tags that has them
        for ( auto & it : allTags ) {
            if ( it.dims ) {
                delete[] it.dims;
                it.dims = 0;
            }
        }
        v = IDL_ImportArray ( 1, dims, IDL_TYP_STRUCT, 0, 0, allStruct );
        dumpStruct ( v, -1, 2 );
        return IDL_GettmpInt(-1);    // return a dummy
    }

    vector<IDL_STRUCT_TAG_DEF> tags;
    size_t patchSize = getPatchSize ( info.get(), loadMask );
    createTags ( tags, loadMask, info.get() );
    
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
    if ( ( loadMask & MOMFBD_MODES ) && info->nModes ) {
        if ( info->version >= 20110714.0 ) {
        cout << "BLUBLU" << endl;
            totalSize += 2*sizeof ( float );                                    // pix2cf, cf2pix
        }
        totalSize += info->nPH * info->nPH * sizeof ( float );                  // Pupil Data
        totalSize += info->nModes * info->nPH * info->nPH * sizeof ( float );   // Mode Data
    }
    totalSize += info->nPatchesX * info->nPatchesY * patchSize;
    while ( totalSize % alignTo ) totalSize++;                                  // pad if not on boundary

    if ( info->nFileNames && ( loadMask & MOMFBD_NAMES ) ) {                          // if filenames are stored and loaded
        totalSize += info->nFileNames * sizeof ( IDL_STRING );
    }

    // Allocate the datablock needed.
    std::unique_ptr<char> data ( new char [ totalSize + sizeof ( MomfdContainer )] ); //, [](char *p) { delete[] p; });
    
    v = IDL_ImportArray ( 1, dims, IDL_TYP_STRUCT, ( UCHAR* ) data.get() + sizeof ( MomfdContainer ), cleanupMomfbd, myStruct );

    // Dump structure layout if requested
    if ( verbosity > 1 ) {
        dumpStruct ( v, -1, 2 );
    }

    loadData ( file, data.get(), loadMask, info.get(), verbosity );

    // release the datablock from the RAII container to prevent de-allocation on return.
    data.release();

    return v;


}


void redux::momfbd_write( int argc, IDL_VPTR* argv, char* argk ) {

    if( argc < 2 ) {
        cout << "MOMFBD_WRITE: needs (at least) 2 arguments." << endl;
        momfbd_help();
        return;
    }
    
    IDL_VPTR momfbdStruct  = argv[0];
    IDL_VPTR fileName  = argv[1];

    KW_RESULT kw;
    kw.help = 0;
    kw.verbose = 0;
    ( void ) IDL_KWProcessByOffset( argc, argv, argk, kw_pars, ( IDL_VPTR* )0, 255, &kw );
    
    if( kw.help ) {
        momfbd_help();
        return;
    }

    IDL_ENSURE_STRUCTURE( momfbdStruct );

    IDL_ENSURE_SIMPLE( fileName );
    IDL_ENSURE_STRING( fileName );

    char *name = IDL_VarGetString ( fileName );
    int verbosity = std::min ( std::max ( ( int ) kw.verbose, 0 ), 8 );
    uint8_t writeMask = 0;
    bool defaultMask(false);
    if ( kw.img )   writeMask |= MOMFBD_IMG;
    if ( kw.psf )   writeMask |= MOMFBD_PSF;
    if ( kw.obj )   writeMask |= MOMFBD_OBJ;
    if ( kw.res )   writeMask |= MOMFBD_RES;
    if ( kw.alpha ) writeMask |= MOMFBD_ALPHA;
    if ( kw.div )   writeMask |= MOMFBD_DIV;
    if ( kw.modes ) writeMask |= MOMFBD_MODES;
    if ( kw.names ) writeMask |= MOMFBD_NAMES;
    if ( kw.all )   writeMask |= MOMFBD_ALL;

    IDL_KW_FREE;
    
    if( !writeMask ) {  // write all by default
        writeMask |= MOMFBD_ALL;
        defaultMask = true;
    }

    if ( verbosity > 0 ) {
        cout << "Writing file: \"" << name << "\"  writeMask=" << bitString(writeMask) << endl;
    }
    
    FileMomfbd* infoPtr = new FileMomfbd();
    std::shared_ptr<FileMomfbd> info ( infoPtr );
    
    IDL_StructDefPtr structDef = momfbdStruct->value.s.sdef;
    
    IDL_VPTR tagPtr;
    IDL_INT* intPtr;
    IDL_STRING* stringPtr;
    IDL_MEMINT tagOffset;
    string tag,type;
    int nTags = IDL_StructNumTags ( structDef );
    for ( int t = 0; t < nTags; ++t ) {

            tag = IDL_StructTagNameByIndex ( structDef, t, 0, 0 );
            tagOffset = IDL_StructTagInfoByIndex ( structDef, t, 0, &tagPtr );

            //cout << "TAG #" << t << "   " << tag << "    offset = " << tagOffset << "    data = " << (int64_t)momfbdStruct->value.arr->data << endl;

            if(tagOffset < 0) continue;
            
            if(tagPtr->type == IDL_TYP_STRING) {
                stringPtr = (IDL_STRING*)(momfbdStruct->value.arr->data + tagOffset);
                if(contains("VERSION",tag,true)){
                    infoPtr->versionString = stringPtr->s;
                } else if(contains("DATE",tag,true)){
                    infoPtr->dateString = stringPtr->s;
                } else if(contains("TIME",tag,true)){
                    infoPtr->timeString = stringPtr->s;
                } else if((writeMask & MOMFBD_NAMES) && contains("NAME",tag,true)){
                    if ( tagPtr->flags & IDL_V_ARR && tagPtr->value.arr->n_dim==1) {
                        for ( int n = 0; n < infoPtr->nFileNames; ++n ) {
                            infoPtr->fileNames.push_back(string(stringPtr[n].s));
                        }
                        infoPtr->nFileNames = infoPtr->fileNames.size();
                    }
                } 
            } else if(tagPtr->type == IDL_TYP_INT) {
                intPtr = (IDL_INT*)(momfbdStruct->value.arr->data + tagOffset);
                if((writeMask) && contains("CLIP",tag,true)) {    // dimensions: (nChannels,2,2)
                    if ( tagPtr->flags & IDL_V_ARR && tagPtr->value.arr->n_dim==3) {
                        int32_t nChannels = infoPtr->nChannels = tagPtr->value.arr->dim[0];
                        infoPtr->clipStartX.reset ( new int16_t[ nChannels ], [] ( int16_t * p ) { delete[] p; } );
                        infoPtr->clipEndX.reset ( new int16_t[ nChannels ], [] ( int16_t * p ) { delete[] p; } );
                        infoPtr->clipStartY.reset ( new int16_t[ nChannels ], [] ( int16_t * p ) { delete[] p; } );
                        infoPtr->clipEndY.reset ( new int16_t[ nChannels ], [] ( int16_t * p ) { delete[] p; } );

                        // fast index first => in the IDL struct X & Y are interchanged compared to the c++ names
                        for ( int i = 0; i < nChannels; ++i ) {
                            infoPtr->clipStartY.get() [ i ] = intPtr[ i + 0 * nChannels ];
                            infoPtr->clipStartX.get() [ i ] = intPtr[ i + 1 * nChannels ];
                            infoPtr->clipEndY.get() [ i ] = intPtr[ i + 2 * nChannels ];
                            infoPtr->clipEndX.get() [ i ] = intPtr[ i + 3 * nChannels ];
                        }
                    }
                }
            } else if(tagPtr->type == IDL_TYP_FLOAT) {  // pupil, mode
                if(writeMask & MOMFBD_MODES) {
                    if(contains("PUPIL",tag,true)) {    // dimensions: (nPH,nPH)
                        if ( tagPtr->flags & IDL_V_ARR && tagPtr->value.arr->n_dim==2) {
                            cout << "PUPIL tag exists -> writing" << endl;
                            infoPtr->nPH = tagPtr->value.arr->dim[0];
                            infoPtr->phOffset = tagOffset;
                        }
                    } else if(contains("MODE",tag,true)) {    // dimensions: (nPH,nPH, nModes)
                        if ( tagPtr->flags & IDL_V_ARR && tagPtr->value.arr->n_dim==3) {
                            cout << "MODE tag exists -> writing" << endl;
                            infoPtr->nModes = tagPtr->value.arr->dim[2];
                            infoPtr->modesOffset = tagOffset;
                        }
                    }
                }
            } else if(tagPtr->type == IDL_TYP_STRUCT) { // patches
                if((writeMask & MOMFBD_PATCH) && contains("PATCH",tag,true)){
                    if ( tagPtr->flags & IDL_V_ARR && tagPtr->value.arr->n_dim==2) {
                        infoPtr->nPatchesY = tagPtr->value.arr->dim[0];     // fast index first => in the IDL struct X & Y are interchanged compared to the c++ names
                        infoPtr->nPatchesX = tagPtr->value.arr->dim[1];
                        int64_t offset = 0;
                        infoPtr->patches.resize(infoPtr->nPatchesX,infoPtr->nPatchesY);
                        for ( int x = 0; x < infoPtr->nPatchesX; ++x ) {
                            for ( int y = 0; y < infoPtr->nPatchesY; ++y ) {
                                //cout << "patch:  offset = " << offset << endl;
                                FileMomfbd::PatchInfo* patch = infoPtr->patches.ptr ( x, y );
                                patch->offset = tagOffset + offset;
                                IDL_StructDefPtr subStruct = tagPtr->value.s.sdef;
                                int64_t subOffset = IDL_StructTagInfoByName(subStruct, (char*)"YL", 0, 0);     // fast index first => in the IDL struct X & Y are interchanged compared to the c++ names
                                if( subOffset >= 0 ) {
                                    patch->region[0] = *((IDL_INT*)(momfbdStruct->value.arr->data + patch->offset + subOffset));
                                } else cout << "patch (" << x << "," << y << ")  has no YL tag." << endl;
                                subOffset = IDL_StructTagInfoByName(subStruct, (char*)"YH", 0, 0);
                                if( subOffset >= 0 ) {
                                    patch->region[1] = *((IDL_INT*)(momfbdStruct->value.arr->data + patch->offset + subOffset));
                                } else cout << "patch (" << x << "," << y << ")  has no YH tag." << endl;
                                subOffset = IDL_StructTagInfoByName(subStruct, (char*)"XL", 0, 0);             // fast index first => in the IDL struct X & Y are interchanged compared to the c++ names
                                if( subOffset >= 0 ) {
                                    patch->region[2] = *((IDL_INT*)(momfbdStruct->value.arr->data + patch->offset + subOffset));
                                } else cout << "patch (" << x << "," << y << ")  has no XL tag." << endl;
                                subOffset = IDL_StructTagInfoByName(subStruct, (char*)"XH", 0, 0);
                                if( subOffset >= 0 ) {
                                    patch->region[3] = *((IDL_INT*)(momfbdStruct->value.arr->data + patch->offset + subOffset));
                                } else cout << "patch (" << x << "," << y << ")  has no XH tag." << endl;
                                
                                patch->nPixelsX = patch->region[1] - patch->region[0] + 1;
                                patch->nPixelsY = patch->region[3] - patch->region[2] + 1;

                                IDL_VPTR tmpPtr;
                                subOffset = IDL_StructTagInfoByName(subStruct, (char*)"DY", IDL_MSG_INFO|IDL_MSG_ATTR_NOPRINT, &tmpPtr);   // fast index first => in the IDL struct X & Y are interchanged compared to the c++ names
                                if( subOffset >= 0 ) {
                                    int32_t nChannels = tmpPtr->value.arr->dim[0];
                                    intPtr = (IDL_INT*)(momfbdStruct->value.arr->data + patch->offset + subOffset);
                                    if( nChannels > 0 ) {
                                        patch->dx.reset ( new int32_t[ nChannels ], [] ( int32_t * p ) { delete[] p; } );
                                        for ( int i = 0; i < nChannels; ++i ) {
                                            patch->dx.get() [ i ] = intPtr[ i ];
                                        }
                                        if( !infoPtr->nChannels ) infoPtr->nChannels = nChannels;
                                        else {
                                            if( nChannels != infoPtr->nChannels ) cout << "patch (" << x << "," << y << ")  nC mismatch !!!" << endl;
                                        }
                                    }
                                }
                                subOffset = IDL_StructTagInfoByName(subStruct, (char*)"DX", IDL_MSG_INFO|IDL_MSG_ATTR_NOPRINT, &tmpPtr);
                                if( subOffset >= 0 ) {
                                    int32_t nChannels = tmpPtr->value.arr->dim[0];
                                    intPtr = (IDL_INT*)(momfbdStruct->value.arr->data + patch->offset + subOffset);
                                    if( nChannels > 0 ) {
                                        patch->dy.reset ( new int32_t[ nChannels ], [] ( int32_t * p ) { delete[] p; } );
                                        for ( int i = 0; i < nChannels; ++i ) {
                                            patch->dy.get() [ i ] = intPtr[ i ];
                                        }
                                        if( !infoPtr->nChannels ) infoPtr->nChannels = nChannels;
                                        else {
                                            if( nChannels != infoPtr->nChannels ) cout << "patch (" << x << "," << y << ")  nC mismatch !!!" << endl;
                                        }
                                    }
                                }
                                
                                subOffset = IDL_StructTagInfoByName(subStruct, (char*)"NIMG", IDL_MSG_INFO|IDL_MSG_ATTR_NOPRINT, &tmpPtr);
                                if( subOffset >= 0 ) {
                                    int32_t nChannels = tmpPtr->value.arr->dim[0];
                                    intPtr = (IDL_INT*)(momfbdStruct->value.arr->data + patch->offset + subOffset);
                                    if( nChannels > 0 ) {
                                        patch->nim.reset ( new int32_t[ nChannels ], [] ( int32_t * p ) { delete[] p; } );
                                        for ( int i = 0; i < nChannels; ++i ) {
                                            patch->nim.get() [ i ] = intPtr[ i ];
                                        }
                                        if( !infoPtr->nChannels ) infoPtr->nChannels = nChannels;
                                        else {
                                            if( nChannels != infoPtr->nChannels ) cout << "patch (" << x << "," << y << ")  nC mismatch !!!" << endl;
                                        }
                                    }
                                }
                                
                                if(writeMask & MOMFBD_IMG) {
                                    subOffset = IDL_StructTagInfoByName(subStruct, (char*)"IMG", IDL_MSG_INFO|IDL_MSG_ATTR_NOPRINT, &tmpPtr);
                                    if( subOffset >= 0 ) {
                                        patch->imgPos = patch->offset + subOffset;
                                    } else if ( defaultMask ) {
                                        writeMask &= ~MOMFBD_IMG;   // no data, so don't try to write it later
                                    } else {
                                        cout << "Attempting to write IMG data without it being present in the structure." << endl;
                                    }
                                }

                                if(writeMask & MOMFBD_PSF) {
                                    subOffset = IDL_StructTagInfoByName(subStruct, (char*)"PSF", IDL_MSG_INFO|IDL_MSG_ATTR_NOPRINT, &tmpPtr);
                                    if( subOffset >= 0 ) {    // dimensions: (xPixels,yPixels,nPSF(=nFiles) )
                                        patch->psfPos = patch->offset + subOffset;
                                        patch->npsf = tmpPtr->value.arr->dim[2];
                                    } else if ( defaultMask ) {
                                        writeMask &= ~MOMFBD_PSF;   // no data, so don't try to write it later
                                    } else {
                                        cout << "Attempting to write PSF data without it being present in the structure." << endl;
                                    }
                                }

                                if(writeMask & MOMFBD_OBJ) {
                                    subOffset = IDL_StructTagInfoByName(subStruct, (char*)"OBJ", IDL_MSG_INFO|IDL_MSG_ATTR_NOPRINT, &tmpPtr);
                                    if( subOffset >= 0 ) {    // dimensions: (xPixels,yPixels,nOBJ(=nFiles) )
                                        patch->objPos = patch->offset + subOffset;
                                        patch->nobj = tmpPtr->value.arr->dim[2];
                                    } else if ( defaultMask ) {
                                        writeMask &= ~MOMFBD_OBJ;   // no data, so don't try to write it later
                                    } else {
                                        cout << "Attempting to write OBJ data without it being present in the structure." << endl;
                                    }
                                }

                                if(writeMask & MOMFBD_RES) {
                                    subOffset = IDL_StructTagInfoByName(subStruct, (char*)"RES", IDL_MSG_INFO|IDL_MSG_ATTR_NOPRINT, &tmpPtr);
                                    if( subOffset >= 0 ) {    // dimensions: (xPixels,yPixels,nRES(=nFiles) )
                                        cout << "RES tag exists -> writing" << endl;
                                        patch->resPos = patch->offset + subOffset;
                                        patch->nres = tmpPtr->value.arr->dim[2];
                                    } else if ( defaultMask ) {
                                        writeMask &= ~MOMFBD_RES;   // no data, so don't try to write it later
                                    } else {
                                        cout << "Attempting to write RES data without it being present in the structure." << endl;
                                    }
                                }

                                if(writeMask & MOMFBD_ALPHA) {
                                    subOffset = IDL_StructTagInfoByName(subStruct, (char*)"ALPHA", IDL_MSG_INFO|IDL_MSG_ATTR_NOPRINT, &tmpPtr);
                                    if( subOffset >= 0 ) {    // dimensions: (nModes,nALPHA(=nFiles) )
                                        patch->alphaPos = patch->offset + subOffset;
                                        patch->nm = tmpPtr->value.arr->dim[0];
                                        patch->nalpha = tmpPtr->value.arr->dim[1];
                                    } else if ( defaultMask ) {
                                        writeMask &= ~MOMFBD_ALPHA;   // no data, so don't try to write it later
                                    } else {
                                        cout << "Attempting to write ALPHA data without it being present in the structure." << endl;
                                    }
                                }

                                if(writeMask & MOMFBD_DIV) {
                                    subOffset = IDL_StructTagInfoByName(subStruct, (char*)"DIV", IDL_MSG_INFO|IDL_MSG_ATTR_NOPRINT, &tmpPtr);
                                    if( subOffset >= 0 ) {    // dimensions: (nModes,nDIV(=nFiles) )
                                        patch->diversityPos = patch->offset + subOffset;
                                        patch->nphx = tmpPtr->value.arr->dim[0];
                                        patch->nphy = tmpPtr->value.arr->dim[1];
                                        patch->ndiv = tmpPtr->value.arr->dim[2];
                                    } else if ( defaultMask ) {
                                        writeMask &= ~MOMFBD_DIV;   // no data, so don't try to write it later
                                    } else {
                                        cout << "Attempting to write DIV data without it being present in the structure." << endl;
                                    }
                                }
                                patch-> nChannels = infoPtr->nChannels;
                                offset += tagPtr->value.arr->elt_len;
                                
                            }   // patches x-loop
                        }       // patches y-loop

                    }           // if nDims == 2
                }               // PATCH
                
            }                   // struct

        }       // loop over tags
    
    infoPtr->write(name, (const char*)momfbdStruct->value.arr->data, writeMask, verbosity);
    
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
        cout << "             clip = " << (do_clip?"YES":"NO") << endl;
        cout << "           margin = " << margin << endl;
        cout << "        transpose = " << (transpose?"YES":"NO") << endl;
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

        static IDL_SYSFUN_DEF2 procedure_addr[] = {
            { {( IDL_SYSRTN_GENERIC )redux::momfbd_write}, ( char* )"MOMFBD_WRITE", 0, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
        };

        /* Register our routine. The routines must be specified exactly the same as in testmodule.dlm. */
        return IDL_SysRtnAdd ( function_addr, TRUE, IDL_CARRAY_ELTS ( function_addr )) &&
               IDL_SysRtnAdd ( procedure_addr, FALSE, IDL_CARRAY_ELTS ( procedure_addr ));

    }

}

