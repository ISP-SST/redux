#include "redux/python/common.hpp"

#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "redux/file/fileio.hpp"
#include "redux/util/stringutil.hpp"

#include <boost/date_time/posix_time/posix_time.hpp>

namespace redux {

    namespace python {

        namespace py = pybind11;
        namespace bpx = boost::posix_time;
        using namespace redux::file;
        using namespace std;


        // Maps redux's internal type-id (FileMeta::getIDLType(), shared with the IDL DLM) to a NumPy dtype.
        static py::dtype get_numpy_type( size_t rdx_type ) {
            switch( rdx_type ) {
                case 1:  return py::dtype::of<uint8_t>();
                case 2:  return py::dtype::of<int16_t>();
                case 3:  return py::dtype::of<int32_t>();
                case 4:  return py::dtype::of<float>();
                case 5:  return py::dtype::of<double>();
                case 12: return py::dtype::of<uint16_t>();
                case 13: return py::dtype::of<uint32_t>();
                default: return py::dtype::of<float>();
            }
        }


        static py::dict readdata( py::object input, bool header = false, bool date_beg = false, bool raw = false, bool all = false ) {

            py::dict result;
            vector<string> existingFiles;

            if( py::isinstance<py::str>( input ) ) {
                string f = input.cast<string>();
                if( boost::filesystem::exists( f ) ) existingFiles.push_back( f );
            } else if( py::isinstance<py::list>( input ) ) {
                for( auto item: input.cast<py::list>() ) {
                    string f = item.cast<string>();
                    if( boost::filesystem::exists( f ) ) existingFiles.push_back( f );
                }
            }

            if( existingFiles.empty() ) {
                result["data"] = py::none();
                return result;
            }

            FileMeta::Ptr myMeta = getMeta( existingFiles[0] );
            if( !myMeta ) {
                result["data"] = py::none();
                return result;
            }

            if( header ) {
                myMeta->getAverageTime();
                myMeta->getEndTime();
                myMeta->getStartTime();
                vector<string> hdrTexts = myMeta->getText( raw );
                py::list cards;
                int nTexts = hdrTexts.size();
                if( all == 0 ) nTexts = 1;
                if( nTexts ) {
                    size_t charCount( 0 );
                    for( int i=0; i<nTexts; ++i ) charCount += hdrTexts[i].size();
                    if( charCount%80 ) {
                        // skip returning header cards for headers which aren't a multiple of 80 bytes long.
                        result["header"] = py::list();
                    } else {
                        for( int i=0; i<nTexts; ++i ) {
                            string hdr = hdrTexts[i];
                            while( !hdr.empty() ) {
                                string tmp = hdr.substr( 0, 80 );
                                cards.append( tmp );
                                hdr.erase( 0, 80 );
                            }
                        }
                        result["header"] = cards;
                    }
                } else {
                    result["header"] = cards;
                }
            }

            if( date_beg ) {
                py::list dates;
                vector<bpx::ptime> date_beg = myMeta->getStartTimes();
                if( date_beg.empty() ) date_beg.push_back( myMeta->getStartTime() );
                for( const auto& t: date_beg ) {
                    dates.append( bpx::to_iso_extended_string( t ) );
                }
                result["date_beg"] = dates;
            }

            size_t nDims = myMeta->nDims();
            if( !nDims ) {
                result["data"] = py::none();
                return result;
            }
            vector<ssize_t> dims;
            for( size_t i=0; i<nDims; ++i ) dims.push_back( myMeta->dimSize( i ) );

            py::dtype dt = get_numpy_type( myMeta->getIDLType() );
            py::array data = py::array( dt, dims );
            readFile( existingFiles[0], (char*)data.mutable_data(), myMeta );
            result["data"] = data;

            return result;
        }

        void init_io( py::module_& m ) {
            m.def( "readdata", &readdata,
                   "Read a file (FITS, ANA, momfbd). Direct port of readdata() from filetools.cpp.",
                   py::arg( "input" ),
                   py::arg( "header" ) = false,
                   py::arg( "date_beg" ) = false,
                   py::arg( "raw" ) = false,
                   py::arg( "all" ) = false );
        }

    }  // namespace python

}  // namespace redux
