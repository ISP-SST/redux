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

        enum FT_FLAGS { FT_REORDER=1, FT_NORMALIZE };

        class FourierTransform : public redux::util::Array<complex_t> {

        public:
            
            struct Plan : public std::enable_shared_from_this<const Plan>  {
                typedef std::shared_ptr<const Plan> Ptr;
                enum TYPE { R2C=1, C2C } tp;
                fftw_plan forward_plan, backward_plan;
                std::vector<size_t> sizes;
                Plan( const std::vector<size_t>& dims, TYPE tp=R2C );
                //template <typename ...S> Plan( S ...s ) : Plan({s...}) {}
                ~Plan();
                void init( void );
                Ptr shared(void) const { return shared_from_this(); };
                bool operator<( const Plan& rhs ) const { if(tp == rhs.tp) return (sizes < rhs.sizes); return tp < rhs.tp; };
            };

            static Plan::Ptr getPlan(const std::vector<size_t>& dims, Plan::TYPE tp);
            //template <typename ...S> static const Plan& getPlan( S ...s ) { return getPlan({s...}); }
            
            FourierTransform() {}

            template <typename T>
            FourierTransform( const redux::util::Array<T>&, int flags=0);
            
            template <typename T>
            void inv( redux::util::Array<T>& );
            
            template <typename T>
            redux::util::Array<T> correlate( const redux::util::Array<T>& ) const;
            
            void autocorrelate( void );
            redux::util::Array<double> power( void ) const;
            
            void convolveInPlace( redux::util::Array<double>& in, int flags=0 ) const;
            void convolveInPlace( redux::util::Array<complex_t>& in, int flags=0 ) const;
            void convolveInPlaceHC( redux::util::Array<complex_t>& in, int flags=0 ) const;
            redux::util::Array<double> convolve( const redux::util::Array<double>& in, int flags=0 ) const;
            redux::util::Array<complex_t> convolve( const redux::util::Array<complex_t>& in, int flags=0 ) const;
            
            //template <typename T> void convolveInPlace( redux::util::Array<T>&, int flags=0 ) const;
            //template <typename T> redux::util::Array<T> convolve( const redux::util::Array<T>&, int flags=0 ) const;
            
            static void normalize( FourierTransform& );
            void normalize() { normalize(*this); };

            template <typename T>
            static void reorder( redux::util::Array<T>& );

            void reorder(void) { reorder(*this); centered = !centered; };

        private:
            void init( redux::util::Array<double>&);
            void init( redux::util::Array<complex_t>&);

            Plan::Ptr plan;

            bool centered;
            bool halfComplex;
            bool normalized;
            bool reOrdered;
            double nInputElements;
            
        };
       
  //      template <>
  //      redux::util::Array<double> FourierTransform::convolve( const redux::util::Array<double>& ) const;

    }   // image

}   // redux


#endif  // REDUX_IMAGE_FOURIERTRANSFORM_HPP
