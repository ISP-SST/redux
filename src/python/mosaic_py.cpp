#include "redux/python/common.hpp"

#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "redux/image/utils.hpp"
#include "redux/util/arrayutil.hpp"


namespace redux {

    namespace python {

        namespace py = pybind11;
        using namespace std;


        static py::array_t<double> mozaic( py::array_t<double, py::array::c_style | py::array::forcecast> patches,
                                            vector<int32_t> posX, vector<int32_t> posY,
                                            int32_t blend = -1, int32_t margin = -1, bool transpose = false ) {

            auto buf = patches.request();
            if( buf.ndim != 3 ) throw invalid_argument( "mozaic: patches must be a 3D stack (nPatches, pRows, pCols)" );
            size_t nPatches = buf.shape[0], pRows = buf.shape[1], pCols = buf.shape[2];
            if( posX.size() != nPatches || posY.size() != nPatches ) {
                throw invalid_argument( "mozaic: posX/posY length must match the number of patches" );
            }

            int32_t minX = *min_element( posX.begin(), posX.end() );
            int32_t maxX = *max_element( posX.begin(), posX.end() ) - minX;
            for( auto& v: posX ) v -= minX;
            int32_t minY = *min_element( posY.begin(), posY.end() );
            int32_t maxY = *max_element( posY.begin(), posY.end() ) - minY;
            for( auto& v: posY ) v -= minY;

            if( margin < 0 ) margin = std::max( pCols, pRows )/8;
            if( blend < 0 ) blend = ( std::min( pCols, pRows ) - 2*margin )/3;
            if( margin > (int32_t)(pCols/2) || margin > (int32_t)(pRows/2) ) {
                throw invalid_argument( "mozaic: margin is too large, nothing would be left" );
            }

            size_t imgCols = maxX+pCols, imgRows = maxY+pRows;

            auto patchPtrs = redux::util::reshapeArray( (const double*)buf.ptr, nPatches, pRows, pCols );

            redux::util::Array<double> outArr( imgRows, imgCols );
            auto outView = outArr.reshape( imgRows, imgCols );

            redux::image::mozaic( outView.get(), imgRows, imgCols, patchPtrs.get(), nPatches, pRows, pCols,
                                   posY.data(), posX.data(), blend, margin, transpose );

            py::array_t<double> out( vector<ssize_t>{ (ssize_t)imgRows, (ssize_t)imgCols } );
            std::copy( outArr.get(), outArr.get()+imgRows*imgCols, out.mutable_data() );
            return out;
        }

        void init_mosaic( py::module_& m ) {
            m.def( "mozaic", &mozaic,
                   "Stitch a stack of patches (nPatches, pRows, pCols) into a single image, with apodized blending. "
                   "posX/posY give each patch's relative position. "
                   "blend/margin default to patchSize/8 and remainder/3, same as reduxd uses.",
                   py::arg( "patches" ), py::arg( "posX" ), py::arg( "posY" ),
                   py::arg( "blend" ) = -1, py::arg( "margin" ) = -1, py::arg( "transpose" ) = false );
        }

    }  // namespace python

}  // namespace redux
