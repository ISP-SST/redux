#ifndef REDUX_IMAGE_FOURIERTRANSFORM_HPP
#define REDUX_IMAGE_FOURIERTRANSFORM_HPP

#include "redux/types.hpp"
#include "redux/util/array.hpp"
#include "redux/util/stringutil.hpp"

#include <fftw3.h>

namespace redux {

    namespace image {
        
        
        // TODO: optimize 
        // TODO: support >2 dimensions.

        enum FT_FLAGS { FT_REORDER=1,           //!< Re-order transform to have the 0 frequency at the center (FFTW expects it at (0,0)
                        FT_NORMALIZE,           //!< Normalize transform
                        FT_FULLCOMPLEX=4        //!< Force full complex format (default is to use "r2c half-complex" format for real input data)
                      };

        class FourierTransform : public redux::util::Array<complex_t> {

        public:
            
            struct Plan {
                typedef std::shared_ptr<const Plan> Ptr;
                enum TYPE { R2C=1, C2C };
                struct Index {
                    Index(const std::vector<size_t>& dims, TYPE t, uint8_t nt);
                    bool operator<( const Index& rhs ) const;
                    TYPE tp;
                    uint8_t nThreads;
                    std::vector<size_t> sizes;
                } id;
                fftw_plan forward_plan, backward_plan;
                Plan( const Index& );
                ~Plan();
                void init( void );
                bool operator<( const Plan& rhs ) const { return (id < rhs.id); };
            };

            static Plan::Ptr getPlan(const std::vector<size_t>& dims, Plan::TYPE tp, uint8_t nThreads=1);
            
            FourierTransform();
            FourierTransform(size_t ySize, size_t xSize, int flags=0, uint8_t nThreads=1);
            FourierTransform( const FourierTransform&);
            template <typename T>
            FourierTransform( const redux::util::Array<T>&, int flags=0, uint8_t nThreads=1);
            template <typename T>
            FourierTransform( const T*, size_t ySize, size_t xSize, int flags=0, uint8_t nThreads=1);
            
            template <typename T>
            void reset( const T*, size_t ySize, size_t xSize, int flags=0, uint8_t nThreads=1);
            template <typename T>
            void reset( const redux::util::Array<T>&in, int flags=0, uint8_t nThreads=1) { reset(in.get(),in.dimSize(0),in.dimSize(1),flags,nThreads); }
            
            template <typename T>
            void directInverse( T* );
            template <typename T>
            void directInverse( redux::util::Array<T>& out ) { directInverse(out.get()); }
            template <typename T>
            void inv( redux::util::Array<T>&, int flags=0 ) const;
            
            template <typename T>
            redux::util::Array<T> correlate( const redux::util::Array<T>& ) const;
            
            void autocorrelate( void );
            void autocorrelate( double scale );
            template <typename T>
            static void autocorrelate( T*, size_t, size_t );
            template <typename T>
            static void autocorrelate( redux::util::Array<T>& in ) { autocorrelate(in.get(), in.dimSize(0), in.dimSize(1)); }
            template <typename T>
            static void autocorrelate( const redux::util::Array<T>&, redux::util::Array<T>& );

            
            redux::util::Array<double> power( void ) const;
            double noise(int mask=0, double cutoff=0) const;
            
            template <typename T>
            void convolveInPlace( redux::util::Array<T>& in, int flags=0  ) const;

            template <typename T>
            redux::util::Array<T> convolve( const redux::util::Array<T>& in ) const {
                int flags = FT_NORMALIZE;
                auto dims = in.dimensions();
                if (halfComplex) {            // for half-complex, the last dimension has half size (n/2+1)
                    dims.back() >>= 1;
                    dims.back() += 1;
                } else flags |= FT_FULLCOMPLEX;
                
                if (dims != dimensions()) {
                    return Array<T>();
                }
                redux::util::Array<double> tmp(in.dimensions());
                if (centered) {
                    flags |= FT_REORDER;
                }
                FourierTransform inFT (in, flags);
                inFT *= *this;
                inFT.directInverse(tmp);
                return std::move(tmp.copy<T>());
            }
            
            static void normalize( FourierTransform& );
            void normalize(void) { normalize(*this); };
            void setThreads(uint8_t nT) { nThreads = nT; };
            double nInputElements(void) const { return inputSize; };

            void reorder(void);
            template <typename T>
            static void reorderInto( const T* inSmall, size_t inSizeY, size_t inSizeX, T* outBig, size_t outSizeY, size_t outSizeX);
            template <typename T>
            static void reorder( T*, size_t, size_t );
            template <typename T>
            static void reorder( redux::util::Array<T>& in ) { reorder(in.get(), in.dimSize(0), in.dimSize(1)); }

            FourierTransform reordered(void) const;
            
            void resize(size_t ySize, size_t xSize, int flags=0, uint8_t nT=1);
            
            operator const redux::util::Array<complex_t>&() const { return reinterpret_cast<const redux::util::Array<complex_t>&>(*this); }
            const FourierTransform& operator*=( const FourierTransform& rhs );
            const FourierTransform& operator=( const FourierTransform& rhs );

            template <typename T>
            const FourierTransform& operator=( const T& val ) { Array<complex_t>::operator=(val); return *this; }
            template <typename T>
            const FourierTransform& operator*=( const T& val ) { Array<complex_t>::operator*=(val); return *this; }
            

            void ft( double* );         //!< NOTE use of these requires FT and input to be properly sized already !!
            void ft( complex_t* );
            void ift( double* );        //!> NOTE fftw:c2r is destructive, so only use these direct functions in one-shot contexts.
            void ift( complex_t* );
            
            
            void init(void);
            void init(const double*, size_t, size_t);
            void init(const complex_t*, size_t, size_t);

        private:
            Plan::Ptr plan;

            bool centered;                     // true if the FT is centered, output from FFTW is NOT centered
            bool halfComplex;
            bool normalized;
            uint8_t nThreads;
            double inputSize;
            
        };
       
  //      template <>
  //      redux::util::Array<double> FourierTransform::convolve( const redux::util::Array<double>& ) const;

    }   // image

}   // redux


#endif  // REDUX_IMAGE_FOURIERTRANSFORM_HPP
