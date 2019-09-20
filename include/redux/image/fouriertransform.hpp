#ifndef REDUX_IMAGE_FOURIERTRANSFORM_HPP
#define REDUX_IMAGE_FOURIERTRANSFORM_HPP

#include "redux/types.hpp"
#include "redux/util/array.hpp"
#include "redux/util/stringutil.hpp"

#include <mutex>

namespace redux {

    namespace image {
        
        
        // TODO: optimize 
        // TODO: support >2 dimensions.

        enum FT_FLAGS { FT_REORDER=1,           //!< Re-order transform to have the 0 frequency at the center (FFTW expects it at (0,0)
                        FT_NORMALIZE,           //!< Normalize transform
                        FT_FULLCOMPLEX=4        //!< Force full complex format (default is to use "r2c half-complex" format for real input data)
                      };

        class FourierTransform : public redux::util::Array<complex_t> {
            
            struct PlansContainer {
                PlansContainer() { fftw_init_threads(); };
                ~PlansContainer(){ fftw_cleanup_threads(); };
                std::mutex mtx;
            };

        public:
            struct Plan
#ifdef RDX_TRACE_MEM
            : public redux::util::TraceObject<Plan>
#endif
            {
                typedef std::shared_ptr<const Plan> Ptr;
                enum TYPE { R2C=1, C2C };
                struct Index {
                    Index(const std::vector<size_t>& dims, TYPE t, uint8_t nt);
                    Index( size_t sizeY, size_t sizeX, TYPE t, uint8_t nt);
                    bool operator<( const Index& rhs ) const;
                    TYPE tp;
                    uint8_t nThreads;
                    std::vector<size_t> sizes;
                } id;
                fftw_plan forward_plan, backward_plan;
                explicit Plan( const Index& );
                ~Plan();
                static Plan::Ptr get(const std::vector<size_t>& dims, Plan::TYPE tp, uint8_t nThreads=1);
                static Plan::Ptr get(size_t sizeY, size_t sizeX, Plan::TYPE tp, uint8_t nThreads=1);
                static void clear(void);
                static PlansContainer pc;
                void init( void );
                void forward( double* in, fftw_complex* out ) const ;
                void backward( fftw_complex* in, double* out ) const;

                bool operator<( const Plan& rhs ) const { return (id < rhs.id); };
            };

            
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
            
            void conj( void );
            void norm( void );

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
            size_t nInputElements(void) const { return inputSize; };

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
            FourierTransform& operator*=( const FourierTransform& rhs );
            FourierTransform& operator=( const FourierTransform& rhs );

            template <typename T>
            const FourierTransform& operator=( const T& val ) { Array<complex_t>::operator=(val); return *this; }
            template <typename T>
            const FourierTransform& operator*=( const T& val ) { Array<complex_t>::operator*=(val); return *this; }
            
            void conjugate(void);
            void ft( double* in, complex_t* out );         //!< NOTE use of these requires FT and input to be properly sized already !!
            void ft( double* in );         //!< NOTE use of these requires FT and input to be properly sized already !!
            void ft( complex_t* in );
            void ift( complex_t* in );
            void ift( complex_t* in, double* out );
            void getIFT( double* out );        //!> NOTE fftw:c2r is destructive, so only use these direct functions in one-shot contexts.
            void getIFT( complex_t* out ) const;
            void getFT( complex_t* out ) const;
            
            Plan::Ptr getPlan(void) { return plan; };
            
            void init(void);
            template <typename T> void init(const T*, size_t, size_t);

        private:
            Plan::Ptr plan;

            bool centered;                     // true if the FT is centered, output from FFTW is NOT centered
            bool halfComplex;
            bool normalized;
            uint8_t nThreads;
            size_t inputSize;
            
        };
       
        template <>
        inline void FourierTransform::init( const std::complex<double>* in, size_t ySize, size_t xSize ) {
            if( (ySize != dimSize(0)) || (xSize != dimSize(1))) {
                redux::util::Array<complex_t>::resize(ySize,xSize);
            }

            plan = Plan::get({ySize, xSize}, Plan::C2C, nThreads);
            complex_t* dataPtr = const_cast<complex_t*>(in);    // fftw takes non-const, even if input is not modified.
            fftw_execute_dft (plan->forward_plan, reinterpret_cast<fftw_complex*> (dataPtr), reinterpret_cast<fftw_complex*> (ptr()));
        }
        
        template <>
        inline void FourierTransform::init( const std::complex<float>* in, size_t ySize, size_t xSize ) {
            if( (ySize != dimSize(0)) || (xSize != dimSize(1))) {
                redux::util::Array<complex_t>::resize(ySize,xSize);
            }

            plan = Plan::get({ySize, xSize}, Plan::C2C, nThreads);
            size_t nElements = ySize*xSize;
            complex_t* tmpData = new complex_t[ySize*xSize];
            std::copy(in, in+nElements, tmpData);
            fftw_execute_dft (plan->forward_plan, reinterpret_cast<fftw_complex*> (tmpData), reinterpret_cast<fftw_complex*> (ptr()));
            delete[] tmpData;
        }

    }   // image

}   // redux


#endif  // REDUX_IMAGE_FOURIERTRANSFORM_HPP
