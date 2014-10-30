#ifndef REDUX_MATH_DIFFERENTIATE_HPP
#define REDUX_MATH_DIFFERENTIATE_HPP

#include <functional>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/function.hpp>



namespace redux {

    namespace math {
        
        const double DEFAULT_DIFF_STEP_SIZE = 1e-3;
        
        typedef boost::numeric::ublas::matrix<double> Matrix;

        class NumericDifferentiator {
            typedef boost::numeric::ublas::vector<double> Vector;
            //typedef boost::numeric::ublas::matrix<double> Matrix;
        public:
            explicit NumericDifferentiator ( boost::function<Vector ( const Vector & ) > );
            NumericDifferentiator ( boost::function<Vector ( const Vector & ) >, double stepsize );
            NumericDifferentiator ( boost::function<Vector ( const Vector & ) >, const Vector &stepsizes );

            Matrix operator() ( const Vector & ) const;
        private:
            Vector m_Stepsizes;
            boost::function<Vector ( const Vector & ) > f;
        };
        
        
        template <class IN, class OUT=IN>
        class Differentiate {
            
        public:
            explicit Differentiate ( std::function<OUT(const IN&)> f ) : m_Stepsizes( 1, DEFAULT_DIFF_STEP_SIZE ), f( f ) {};
            Differentiate ( std::function<OUT(const IN&)> f, double stepsize) : m_Stepsizes( 1, stepsize ), f( f ) {};
            Differentiate ( std::function<OUT(const IN&)> f, const std::vector<double>& stepsizes ) : m_Stepsizes( stepsizes ), f( f ) {
                if(m_Stepsizes.size() < 1) m_Stepsizes.resize( 1, 1E-3 );
            };
            void setSteps( const std::vector<double>& stepsizes) { m_Stepsizes = stepsizes; };
            void setSteps( const double& step, size_t n=1 ) { m_Stepsizes.resize(std::max(n,1UL),step); };
            Matrix operator() ( const IN& ) const;
        private:
            std::vector<double> m_Stepsizes;
            std::function<OUT(const IN&)> f;
        };

    }

}

#endif  // REDUX_MATH_DIFFERENTIATE_HPP
