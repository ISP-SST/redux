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
            
            struct Plan {
                fftw_plan forward_plan, backward_plan;
                std::vector<size_t> sizes;
                Plan( const std::vector<size_t>& dims );
                template <typename ...S> Plan( S ...s ) : Plan({s...}) {}
                ~Plan();
                void init( void );
                bool operator<( const Plan& rhs ) const { return (sizes < rhs.sizes); };
            };

            static const Plan& getPlan(const std::vector<size_t>& dims);
            template <typename ...S> static const Plan& getPlan( S ...s ) { return getPlan({s...}); }
            
            FourierTransform() {}

            template <typename T>
            FourierTransform( const redux::util::Array<T>&, int flags=0);

            template <typename T>
            redux::util::Array<T> convolve( const redux::util::Array<T>& ) const;
            
            template <typename T>
            static void normalize( redux::util::Array<T>& );

            template <typename T>
            static void reorder( redux::util::Array<T>& );


        };
       
        template <>
        redux::util::Array<double> FourierTransform::convolve( const redux::util::Array<double>& ) const;

    }   // image

}   // redux


#endif  // REDUX_IMAGE_FOURIERTRANSFORM_HPP
