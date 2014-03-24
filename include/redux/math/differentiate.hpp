#ifndef REDUX_MATH_DIFFERENTIATE_HPP
#define REDUX_MATH_DIFFERENTIATE_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/function.hpp>

namespace redux {

    namespace math {

        class NumericDifferentiator {
            typedef boost::numeric::ublas::vector<double> Vector;
            typedef boost::numeric::ublas::matrix<double> Matrix;
        public:
            explicit NumericDifferentiator ( boost::function<Vector ( const Vector & ) > );
            NumericDifferentiator ( boost::function<Vector ( const Vector & ) >, double stepsize );
            NumericDifferentiator ( boost::function<Vector ( const Vector & ) >, const Vector &stepsizes );

            Matrix operator() ( const Vector & ) const;
        private:
            Vector m_Stepsizes;
            boost::function<Vector ( const Vector & ) > f;
        };

    }

}

#endif  // REDUX_MATH_DIFFERENTIATE_HPP
