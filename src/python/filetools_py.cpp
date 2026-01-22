#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "redux/file/fileio.hpp"
#include "redux/util/stringutil.hpp"

#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>


namespace py = pybind11;
using namespace redux::file;
using namespace redux::util;
using namespace std;
namespace bfs = boost::filesystem;
namespace bpx = boost::posix_time;

// Helper to map Redux internal types to NumPy descriptors
py::dtype get_numpy_type(size_t rdx_type) {
    switch (rdx_type) {
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

py::dict py_readdata(py::object input, bool header = false, bool date_beg = false, bool raw = false, bool all=false) {

    py::dict result;
    vector<string> existingFiles;

    // --- Input Parsing (Matching original logic for input files) ---
    if (py::isinstance<py::str>(input)) {
        string f = input.cast<string>();
        if (bfs::exists(f)) existingFiles.push_back(f);
    } else if (py::isinstance<py::list>(input)) {
        for (auto item : input.cast<py::list>()) {
            string f = item.cast<string>();
            if (bfs::exists(f)) existingFiles.push_back(f);
        }
    }

    if (existingFiles.empty()) {
        result["data"] = py::none();
        return result;
    }

    try {
        // --- Call getMeta (Line 278) ---
        FileMeta::Ptr myMeta = getMeta(existingFiles[0]);
        if (!myMeta) {
            result["data"] = py::none();
            return result;
        }

        // --- Exact Time calls (Lines 280-282) ---
        // --- Header Processing (Lines 285-307) ---
        // Using the literal call you specified: getText(raw)
        if( header ) {
            myMeta->getAverageTime();
            myMeta->getEndTime();
            myMeta->getStartTime();
            vector<string> hdrTexts = myMeta->getText( raw );
            py::list cards;
            int nTexts = hdrTexts.size();
            if( all == 0 ) nTexts = 1;
            if( nTexts ) {
                size_t charCount(0);
                for( int i=0; i<nTexts; ++i ) charCount += hdrTexts[i].size();
                if( charCount%80 ) {
                    throw logic_error("Header text-size is not a multiple of 80.");
                }
                for( int i=0; i<nTexts; ++i ) {
                    string hdr = hdrTexts[i];
                    while( !hdr.empty() ) {
                        string tmp = hdr.substr( 0, 80 );
                        cards.append(tmp);
                        hdr.erase( 0, 80 );
                    }
                }
            }
            result["header"] = cards;
        }

        // --- Date Beg Processing (Lines 350-370) ---
        if (date_beg) {
            py::list dates;
            vector<bpx::ptime> date_beg = myMeta->getStartTimes();
            if( date_beg.empty() ) date_beg.push_back( myMeta->getStartTime() );
            for (const auto& t : date_beg) {
                dates.append(bpx::to_iso_extended_string(t));
            }
            result["date_beg"] = dates;
        }
        size_t nDims = myMeta->nDims();
        if( !nDims ) {
            cerr << "rdx_readdata: file contains no data." << endl;
            return result;
        }
        vector<ssize_t> dims;

        for( size_t i=0; i<nDims; ++i ) {
            dims.push_back( myMeta->dimSize(i) );
        }
        //std::reverse( dims.begin(), dims.end() );

        py::dtype dt = get_numpy_type( myMeta->getIDLType() );
        py::array data = py::array(dt, dims);
        readFile( existingFiles[0], (char*)data.mutable_data(), myMeta );
        result["data"] = data;

    } catch (const std::exception& e) {
        throw std::runtime_error("C++ Exception in pyredux: " + std::string(e.what()));
    }

    return result;
}

PYBIND11_MODULE(pyredux, m) {
    m.def("readdata", &py_readdata,
          "Direct port of readdata from filetools.cpp",
          py::arg("input"),
          py::arg("header") = false,
          py::arg("date_beg") = false,
          py::arg("raw") = false,
          py::arg("all") = false);
}
