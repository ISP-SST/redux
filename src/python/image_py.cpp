#include "redux/python/common.hpp"

#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <pybind11/complex.h>

#include "redux/image/fouriertransform.hpp"
#include "redux/image/pupil.hpp"
#include "redux/image/utils.hpp"
#include "redux/image/zernike.hpp"
#include "redux/momfbd/modes.hpp"


namespace redux {

    namespace python {

        namespace py = pybind11;
        using namespace redux::image;
        using namespace std;


        static py::tuple calculatePupilSize( double wavelength, uint32_t nPixels, double telescopeDiameter, double arcSecsPerPixel ) {
            double freqCutoff, pupilRadius;
            uint16_t nPupilPixels;
            Pupil::calculatePupilSize( freqCutoff, pupilRadius, nPupilPixels, wavelength, nPixels, telescopeDiameter, arcSecsPerPixel );
            return py::make_tuple( freqCutoff, pupilRadius, nPupilPixels );
        }


        static py::array_t<double> generatePupil( uint16_t pixels, double pupilRadius ) {
            Pupil p( pixels, pupilRadius );
            py::array_t<double> arr( vector<ssize_t>{ pixels, pixels } );
            std::copy( p.get(), p.get()+p.nElements(), arr.mutable_data() );
            return arr;
        }


        static py::tuple NollToNM( unsigned int j ) {
            uint16_t n; int16_t m;
            Zernike::NollToNM( j, n, m );
            return py::make_tuple( n, m );
        }


        static double calcCovariance( unsigned int i, unsigned int j ) {
            return Zernike::calcCovariance( i, j );
        }


        static py::array_t<double> getZernike( unsigned int nPixels, float radius, unsigned int nollIndex, float angle = 0.0f, bool normalize = true ) {
            py::array_t<double> arr( vector<ssize_t>{ (ssize_t)nPixels, (ssize_t)nPixels } );
            int flags = Zernike::GET_ZERNIKE | ( normalize ? Zernike::NORMALIZE : 0 );
            Zernike::getZernike( arr.mutable_data(), nPixels, radius, angle, nollIndex, flags );
            return arr;
        }


        static py::array_t<double> apodize( py::array_t<double, py::array::c_style | py::array::forcecast> input, size_t rowBlend, size_t colBlend, size_t rowMargin = 0, size_t colMargin = 0 ) {
            auto buf = input.request();
            if( buf.ndim != 2 ) throw invalid_argument( "apodize: input must be 2D" );
            size_t nRows = buf.shape[0], nCols = buf.shape[1];
            redux::util::Array<double> arr( nRows, nCols );
            std::copy( (double*)buf.ptr, (double*)buf.ptr+nRows*nCols, arr.get() );
            redux::image::apodizeInPlace( arr, rowBlend, colBlend, rowMargin, colMargin );
            py::array_t<double> out( vector<ssize_t>{ (ssize_t)nRows, (ssize_t)nCols } );
            std::copy( arr.get(), arr.get()+nRows*nCols, out.mutable_data() );
            return out;
        }


        static py::array_t<double> generateModeSet( uint16_t pixels, double radius, double angle, vector<uint16_t> modeNumbers, bool kl = true, uint16_t firstMode = 0, uint16_t lastMode = 0, double cutoff = 0.0, bool normalize = true ) {
            if( modeNumbers.empty() ) throw invalid_argument( "generate_modeset: modeNumbers must not be empty" );
            redux::momfbd::ModeList modes;
            for( auto mn: modeNumbers ) {
                modes.push_back( redux::momfbd::ModeID( mn, kl ? redux::momfbd::KARHUNEN_LOEVE : redux::momfbd::ZERNIKE ) );
            }
            redux::momfbd::ModeSet ms;
            int flags = normalize ? Zernike::NORMALIZE : 0;
            if( kl ) {
                uint16_t lm = lastMode ? lastMode : *max_element( modeNumbers.begin(), modeNumbers.end() );
                ms.generate( pixels, radius, angle, firstMode, lm, modes, cutoff, flags );
            } else {
                ms.generate( pixels, radius, angle, modes, flags );
            }
            size_t n = modeNumbers.size();
            py::array_t<double> out( vector<ssize_t>{ (ssize_t)n, (ssize_t)pixels, (ssize_t)pixels } );
            std::copy( ms.get(), ms.get()+ms.nElements(), out.mutable_data() );
            return out;
        }


        static py::array_t<double> psfFromPhase( uint16_t pupilPixels, double pupilRadius, py::array_t<double, py::array::c_style | py::array::forcecast> phase, bool normalizeSum = true ) {
            auto buf = phase.request();
            if( buf.ndim != 2 || (size_t)buf.shape[0] != pupilPixels || (size_t)buf.shape[1] != pupilPixels ) {
                throw invalid_argument( "psf_from_phase: phase must be pupilPixels x pupilPixels" );
            }
            Pupil pupil( pupilPixels, pupilRadius );
            size_t otfSize = 2*pupilPixels;
            vector<complex_t> PF( otfSize*otfSize, complex_t( 0, 0 ) );
            const double* pupPtr = pupil.get();
            const double* phiPtr = static_cast<double*>( buf.ptr );
            for( const auto& idx: pupil.pupilInOTF ) {
                size_t pi = idx.first, oi = idx.second;
                PF[oi] = std::polar( pupPtr[pi], phiPtr[pi] );
            }
            FourierTransform ft( otfSize, otfSize, FULLCOMPLEX );
            vector<complex_t> ftOut( otfSize*otfSize );
            ft.ft( PF.data(), ftOut.data(), FULLCOMPLEX );
            vector<double> psf( otfSize*otfSize );
            for( size_t i=0; i<psf.size(); ++i ) psf[i] = std::norm( ftOut[i] );
            FourierTransform::reorder( psf.data(), otfSize, otfSize );
            if( normalizeSum ) {
                double sum = 0;
                for( double v: psf ) sum += v;
                for( double& v: psf ) v /= sum;
            }
            py::array_t<double> out( vector<ssize_t>{ (ssize_t)otfSize, (ssize_t)otfSize } );
            std::copy( psf.begin(), psf.end(), out.mutable_data() );
            return out;
        }


        static py::array_t<double> convolve( py::array_t<double, py::array::c_style | py::array::forcecast> image, py::array_t<double, py::array::c_style | py::array::forcecast> kernel ) {
            auto ibuf = image.request();
            auto kbuf = kernel.request();
            if( ibuf.ndim != 2 || kbuf.ndim != 2 ) throw invalid_argument( "convolve: inputs must be 2D" );
            if( ibuf.shape[0] != kbuf.shape[0] || ibuf.shape[1] != kbuf.shape[1] ) {
                throw invalid_argument( "convolve: image and kernel must be the same size" );
            }
            size_t ny = ibuf.shape[0], nx = ibuf.shape[1];
            redux::util::Array<double> kArr( ny, nx );
            std::copy( (double*)kbuf.ptr, (double*)kbuf.ptr+ny*nx, kArr.get() );
            FourierTransform ft( kArr, REORDER_FT|NORMALIZE_FT );
            redux::util::Array<double> inArr( ny, nx ), outArr( ny, nx );
            std::copy( (double*)ibuf.ptr, (double*)ibuf.ptr+ny*nx, inArr.get() );
            ft.convolve( inArr.get(), outArr.get() );
            py::array_t<double> out( vector<ssize_t>{ (ssize_t)ny, (ssize_t)nx } );
            std::copy( outArr.get(), outArr.get()+ny*nx, out.mutable_data() );
            return out;
        }


        static py::array_t<complex_t> calcImgFt( uint16_t pupilPixels,
                py::array_t<double, py::array::c_style|py::array::forcecast> image ) {
            size_t otfSize = 2*pupilPixels, otfSize2 = otfSize*otfSize;
            auto buf = image.request();
            if( buf.ndim != 2 || (size_t)buf.shape[0] != otfSize || (size_t)buf.shape[1] != otfSize ) {
                throw invalid_argument( "calc_img_ft: image must be (2*pupilPixels) x (2*pupilPixels)" );
            }
            // same FourierTransform config + single-arg ft() call as SubImage::imgFT
            // uses in initialize().
            FourierTransform imgFT( otfSize, otfSize, REORDER_FT|FULLCOMPLEX );
            imgFT.ft( static_cast<double*>( buf.ptr ) );
            py::array_t<complex_t> out( vector<ssize_t>{ (ssize_t)otfSize, (ssize_t)otfSize } );
            std::copy( imgFT.get(), imgFT.get()+otfSize2, out.mutable_data() );
            return out;
        }


        // Inverse of calc_img_ft. Follows Object::restorePatch's avgObjFT.ift().
        static py::array_t<double> calcImgIft( uint16_t pupilPixels,
                py::array_t<complex_t, py::array::c_style|py::array::forcecast> imgft ) {
            size_t otfSize = 2*pupilPixels, otfSize2 = otfSize*otfSize;
            auto buf = imgft.request();
            if( buf.ndim != 2 || (size_t)buf.shape[0] != otfSize || (size_t)buf.shape[1] != otfSize ) {
                throw invalid_argument( "calc_img_ift: imgft must be (2*pupilPixels) x (2*pupilPixels)" );
            }
            FourierTransform ft( otfSize, otfSize, FULLCOMPLEX|REORDER_FT );
            std::copy( static_cast<complex_t*>(buf.ptr), static_cast<complex_t*>(buf.ptr)+otfSize2, ft.get() );
            py::array_t<double> out( vector<ssize_t>{ (ssize_t)otfSize, (ssize_t)otfSize } );
            ft.ift( out.mutable_data() );
            return out;
        }


        static py::dict karhunenLoeveExpansion( uint16_t firstMode, uint16_t lastMode ) {
            const auto& kle = Zernike::karhunenLoeveExpansion( firstMode, lastMode );
            py::dict result;
            for( const auto& entry: kle ) {
                py::list weights;
                for( const auto& w: entry.second->zernikeWeights ) {
                    weights.append( py::make_tuple( w.first, w.second ) );
                }
                result[py::int_( entry.first )] = py::make_tuple( entry.second->covariance, weights );
            }
            return result;
        }


        void init_image( py::module_& m ) {
            m.def( "calculate_pupil_size", &calculatePupilSize,
                   "(freqCutoff, pupilRadiusInPixels, nPupilPixels) for the given optical setup.",
                   py::arg( "wavelength" ), py::arg( "nPixels" ), py::arg( "telescopeDiameter" ), py::arg( "arcSecsPerPixel" ) );

            m.def( "generatePupil", &generatePupil,
                   "A circular pupil amplitude array (pixels x pixels), via redux's Pupil class.",
                   py::arg( "pixels" ), py::arg( "pupilRadius" ) );

            m.def( "zernike_noll_to_nm", &NollToNM,
                   "Noll index -> (radial order n, azimuthal order m).", py::arg( "nollIndex" ) );

            m.def( "zernike_covariance", &calcCovariance,
                   "Noll cross-covariance of two Zernike modes (Kolmogorov statistics, D/r0=1).",
                   py::arg( "nollIndex1" ), py::arg( "nollIndex2" ) );

            m.def( "zernike_mode", &getZernike,
                   "A single normalized Zernike mode image (pixels x pixels).",
                   py::arg( "nPixels" ), py::arg( "radius" ), py::arg( "nollIndex" ),
                   py::arg( "angle" ) = 0.0f, py::arg( "normalize" ) = true );

            m.def( "kl_expansion", &karhunenLoeveExpansion,
                   "Karhunen-Loeve basis over [firstMode, lastMode]: {noll_index: (covariance, [(noll_index, weight), ...])}.",
                   py::arg( "firstMode" ), py::arg( "lastMode" ) );

            m.def( "apodize", &apodize,
                   "Raised-cosine edge taper (same as redux::image::apodizeInPlace applies to real channel images).",
                   py::arg( "input" ), py::arg( "rowBlend" ), py::arg( "colBlend" ),
                   py::arg( "rowMargin" ) = 0, py::arg( "colMargin" ) = 0 );

            m.def( "generate_modeset", &generateModeSet,
                   "A complete mode basis (nModes x pixels x pixels), via redux's own ModeSet::generate().",
                   py::arg( "pixels" ), py::arg( "radius" ), py::arg( "angle" ), py::arg( "modeNumbers" ),
                   py::arg( "kl" ) = true, py::arg( "firstMode" ) = 0, py::arg( "lastMode" ) = 0,
                   py::arg( "cutoff" ) = 0.0, py::arg( "normalize" ) = true );

            m.def( "psf_from_phase", &psfFromPhase,
                   "PSF (sum-normalized by default) from a phase map over the pupil, using redux's exact "
                   "embedding/FFT convention (ported from the validated gen_synthetic_psf.cpp).",
                   py::arg( "pupilPixels" ), py::arg( "pupilRadius" ), py::arg( "phase" ), py::arg( "normalizeSum" ) = true );

            m.def( "convolve", &convolve,
                   "FFT-based convolution of two same-size arrays, via redux's own FourierTransform::convolve().",
                   py::arg( "image" ), py::arg( "kernel" ) );

            m.def( "calc_img_ft", &calcImgFt,
                   "imgFT from a raw (2*pupilPixels)x(2*pupilPixels) image, via the exact "
                   "FourierTransform config + ft() call SubImage::initialize() uses.",
                   py::arg( "pupilPixels" ), py::arg( "image" ) );

            m.def( "calc_img_ift", &calcImgIft,
                   "Inverse of calc_img_ft -- exact mirror of Object::restorePatch's avgObjFT.ift() call.",
                   py::arg( "pupilPixels" ), py::arg( "imgft" ) );
        }

    }  // namespace python

}  // namespace redux
