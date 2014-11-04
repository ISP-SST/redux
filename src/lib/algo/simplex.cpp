#include "redux/algo/simplex.hpp"


using namespace redux::algo;
using namespace Eigen;

#include <iostream>
using namespace std;

namespace {

    const double EPSILON = 1.0E-10;

}

Simplex::Simplex( MatrixXd& a, VectorXd& b, VectorXd& c ) {

    m_nConstraints = b.rows();
    m_nVariables = c.rows();

    m_Tableu.resize( m_nConstraints + 1, m_nConstraints + m_nVariables + 1 );

    for( int i = 0; i < m_nConstraints; i++ )  {
        for( int j = 0; j < m_nVariables; j++ ) {
            m_Tableu( i, j ) = a( i, j );
        }
    }
    for( int i = 0; i < m_nConstraints; i++ ) m_Tableu( i, m_nVariables + i ) = 1.0;
    for( int j = 0; j < m_nVariables; j++ ) m_Tableu( m_nConstraints, j )   = c( j );
    for( int i = 0; i < m_nConstraints; i++ ) m_Tableu( i, m_nConstraints + m_nVariables ) = b( i );


    std::cout << "Tableu:" << endl << m_Tableu << endl;
    solve();
    std::cout << "Tableu2:" << endl << m_Tableu << endl;

}


// pivot on entry (p, q) using Gauss-Jordan elimination
void Simplex::pivot( int p, int q ) {

    // everything but row p and column q
    for( int i = 0; i <= m_nConstraints; i++ )
        for( int j = 0; j <= m_nConstraints + m_nVariables; j++ )
            if( i != p && j != q ) m_Tableu( i, j ) -= m_Tableu( p, j ) * m_Tableu( i, q ) / m_Tableu( p, q );

    // zero out column q
    for( int i = 0; i <= m_nConstraints; i++ )
        if( i != p ) m_Tableu( i, q ) = 0.0;

    // scale row p
    for( int j = 0; j <= m_nConstraints + m_nVariables; j++ )
        if( j != q ) m_Tableu( p, j ) /= m_Tableu( p, q );
    m_Tableu( p, q ) = 1.0;
}


// run simplex algorithm starting from initial BFS
void Simplex::solve( void ) {
    while( true ) {

        // find entering column q
        int q = bland();
        if( q == -1 ) break; // optimal

        // find leaving row p
        int p = minRatioRule( q );
        if( p == -1 ) {
            cout << "Linear program is unbounded" << endl;
            return;
            //throw new ArithmeticException("Linear program is unbounded");
        }
        // pivot
        pivot( p, q );

    }
}


// return optimal objective value
double Simplex::value( void ) {
    return -m_Tableu( m_nConstraints, m_nConstraints + m_nVariables );
}


// lowest index of a non-basic column with a positive cost
int Simplex::bland( void ) {
    for( int j = 0; j < m_nConstraints + m_nVariables; j++ )
        if( m_Tableu( m_nConstraints, j ) > 0 ) return j;
    return -1;  // optimal
}


// index of a non-basic column with most positive cost
int Simplex::dantzig( void ) {
    int q = 0;
    for( int j = 1; j < m_nConstraints + m_nVariables; j++ )
        if( m_Tableu( m_nConstraints, j ) > m_Tableu( m_nConstraints, q ) ) q = j;

    if( m_Tableu( m_nConstraints, q ) <= 0 ) return -1; // optimal
    else return q;
}


// find row p using min ratio rule (-1 if no such row)
int Simplex::minRatioRule( int q ) {
    int p = -1;
    for( int i = 0; i < m_nConstraints; i++ ) {
        if( m_Tableu( i, q ) <= 0 ) continue;
        else if( p == -1 ) p = i;
        else if( ( m_Tableu( i, m_nConstraints + m_nVariables ) / m_Tableu( i, q ) ) < ( m_Tableu( p, m_nConstraints + m_nVariables ) / m_Tableu( p, q ) ) ) p = i;
    }
    return p;
}
