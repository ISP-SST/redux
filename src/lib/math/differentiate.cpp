#include "redux/math/differentiate.hpp"


using namespace redux::math;

namespace {
    static const double DEFAULT_STEP_SIZE = 1e-3;
}

NumericDifferentiator::NumericDifferentiator( boost::function<Vector( const Vector & ) > f )
    : m_Stepsizes( 1 ), f( f ) {
    m_Stepsizes( 0 ) = DEFAULT_STEP_SIZE;
}

NumericDifferentiator::NumericDifferentiator( boost::function<Vector( const Vector & ) > f, double stepsize )
    : m_Stepsizes( 1 ), f( f ) {
    m_Stepsizes( 0 ) = stepsize;
}

NumericDifferentiator::NumericDifferentiator( boost::function<Vector( const Vector & ) > f, const Vector &stepsizes )
    : m_Stepsizes( stepsizes ), f( f ) {
}

NumericDifferentiator::Matrix NumericDifferentiator::operator()( const Vector &x ) const {
    Vector fx = f( x );
    size_t n = x.size();
    size_t m = fx.size();
    Matrix jacobian( fx.size(), x.size() );
    Vector dx1( n );
    Vector dx2( n );

    for( size_t i = 0; i < n; ++i ) {
        dx1.clear();
        dx2.clear();
        double h = ( m_Stepsizes.size() == 1 ) ? m_Stepsizes( 0 ) : m_Stepsizes( i );
        dx1( i ) -= h;
        dx2( i ) += h;
        Vector fx1 = f( x + dx1 );
        Vector fx2 = f( x + dx2 );

        for( size_t j = 0; j < m; ++j ) {
            jacobian( j, i ) = ( fx2( j ) - fx1( j ) ) / ( 2 * h );
        }
    }

    return jacobian;
}

