#include "redux/python/common.hpp"

#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "redux/file/filemomfbd.hpp"

#include <fstream>

namespace redux {

    namespace python {

        namespace py = pybind11;
        using namespace redux::file;
        using namespace std;


        static py::object read_floats( ifstream& in, int64_t pos, size_t count, vector<ssize_t> shape ) {
            if( !count || pos <= 0 ) return py::none();
            py::array_t<float> arr( shape );
            in.seekg( pos );
            in.read( reinterpret_cast<char*>( arr.mutable_data() ), count*sizeof( float ) );
            if( !in ) throw runtime_error( "read_momfbd: failed to read " + to_string( count ) + " floats at offset " + to_string( pos ) );
            return arr;
        }


        static py::dict read_momfbd( const string& path ) {

            require_file( path );

            FileMomfbd f( path );      // constructor parses the header + per-patch info
            if( f.swapNeeded ) {
                // Reading a file on different endian than creation is not well tested/implemented -> throw
                throw invalid_argument( "read_momfbd: file requires byte-swapping, not supported yet" );
            }

            py::dict result;
            result["version"] = f.version;
            result["nPH"] = f.nPH;
            result["nModes"] = f.nModes;
            result["nPatchesX"] = f.nPatchesX;
            result["nPatchesY"] = f.nPatchesY;
            result["region"] = py::make_tuple( f.region[0], f.region[1], f.region[2], f.region[3] );

            py::list fileNames;
            for( auto& fn: f.fileNames ) fileNames.append( fn );
            result["fileNames"] = fileNames;

            ifstream in( path, ios::binary );

            result["pupil"] = py::none();
            result["modes"] = py::none();
            if( ( f.dataMask & MOMFBD_MODES ) && f.nPH > 0 ) {
                result["pupil"] = read_floats( in, f.phOffset, (size_t)f.nPH*f.nPH, { f.nPH, f.nPH } );
                if( f.nModes > 0 && f.modesOffset >= 0 ) {
                    size_t nPix = (size_t)f.nPH*f.nPH;
                    py::array_t<float> modes( vector<ssize_t>{ f.nModes, f.nPH, f.nPH } );
                    float* dst = modes.mutable_data();
                    in.seekg( f.modesOffset );
                    for( int32_t i=0; i<f.nModes; ++i ) {
                        in.seekg( 1, ios::cur );     // legacy format: 1-byte skip before each mode block
                        in.read( reinterpret_cast<char*>( dst+i*nPix ), nPix*sizeof( float ) );
                    }
                    if( !in ) throw runtime_error( "read_momfbd: failed to read modes block" );
                    result["modes"] = modes;
                }
            }

            py::list patches;
            for( int py_=0; py_<f.nPatchesY; ++py_ ) {
                for( int px_=0; px_<f.nPatchesX; ++px_ ) {
                    auto& patch = f.patches( py_, px_ );
                    size_t nxny = (size_t)patch.nPixelsX*patch.nPixelsY;
                    py::dict pd;
                    pd["py"] = py_;
                    pd["px"] = px_;
                    pd["region"] = py::make_tuple( patch.region[0], patch.region[1], patch.region[2], patch.region[3] );
                    pd["offx"] = patch.offx;
                    pd["offy"] = patch.offy;
                    pd["nChannels"] = patch.nChannels;
                    pd["img"] = read_floats( in, patch.imgPos, nxny, { patch.nPixelsY, patch.nPixelsX } );
                    pd["alpha"] = read_floats( in, patch.alphaPos, (size_t)patch.nalpha*patch.nm, { patch.nalpha, patch.nm } );
                    pd["obj"] = read_floats( in, patch.objPos, (size_t)patch.nobj*nxny, { patch.nobj, patch.nPixelsY, patch.nPixelsX } );
                    pd["psf"] = read_floats( in, patch.psfPos, (size_t)patch.npsf*nxny, { patch.npsf, patch.nPixelsY, patch.nPixelsX } );
                    pd["res"] = read_floats( in, patch.resPos, (size_t)patch.nres*nxny, { patch.nres, patch.nPixelsY, patch.nPixelsX } );
                    patches.append( pd );
                }
            }
            result["patches"] = patches;

            return result;
        }

        void init_momfbd( py::module_& m ) {
            m.def( "read_momfbd", &read_momfbd,
                   "Read a .momfbd result file: global pupil/modes + per-patch img/alpha/obj/psf/res.",
                   py::arg( "path" ) );
        }

    }  // namespace python

}  // namespace redux
