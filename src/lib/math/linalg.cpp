#include "redux/math/linalg.hpp"

#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>

#include <numeric>
#include <iostream>

using namespace redux::math;
using namespace redux::util;
using namespace std;

void redux::math::normalize_columns( double* A, size_t rows, size_t cols ) {
    if( rows == 0 || cols == 0 ) return;
    unique_ptr<double[]> tmp( new double[cols] );
    double* col_norms = tmp.get();
    memset( col_norms, 0, cols*sizeof(double) );
    for( size_t i(0); i<rows; ++i ) {
        double* rat = A + i*cols;
        std::transform( col_norms, col_norms+cols, rat, col_norms, [&]( double& a, double& b ){ return a + b*b; } );
    }
    std::transform( col_norms, col_norms+cols, col_norms, [&]( double& a ){ return 1.0/sqrt(a); } );
    for( size_t i(0); i<rows; ++i ) {
        double* rat = A + i*cols;
        std::transform( rat, rat+cols, col_norms, rat, [&]( double& a, double& b ){ return a*b; } );
    }
}


void redux::math::normalize_columns( redux::util::Array<double>& A ) {
    const vector<size_t>& adims = A.dimensions();
    if( adims.size() != 2 ) {
        cout << "redux::math::normalize_columns:: size of A is not ok: " << printArray(adims,"adims") << endl;
        return;
    }
    normalize_columns( A.get(), adims[0], adims[1] );
}


void redux::math::normalize_rows( double* A, size_t rows, size_t cols ) {
    if( rows == 0 || cols == 0 ) return;
    for( size_t i(0); i<rows; ++i ) {
        double* rat = A + i*cols;
        double norm(0);
        norm = std::accumulate( rat, rat+cols, norm, [&]( double& a, double& b ){ return a+b*b; } );
        if( norm > 0.0 ) norm = 1.0/sqrt(norm);
        else norm = 0;
        std::transform( rat, rat+cols, rat, [&]( double& a ){ return a*norm; } );
    }
    
}


void redux::math::normalize_rows( redux::util::Array<double>& A ) {
    const vector<size_t>& adims = A.dimensions();
    if( adims.size() != 2 ) {
        cout << "redux::math::normalize_rows:: size of A is not ok: " << printArray(adims,"adims") << endl;
        return;
    }
    normalize_rows( A.get(), adims[0], adims[1] );
}
        

void redux::math::gauss_jordan( double* A, size_t rows, size_t cols ) {
    
    cout << "gauss_jordan: nRows = " << rows << "   nCols = " << cols << endl;
    bool transposed(false);

    if( rows > cols ) {
        transpose( A, rows, cols );
        std::swap( rows, cols );
        transposed = true;
    }

    cout << "gauss_jordan:  " << __LINE__ << "  nRows = " << rows << "    nCols = " << cols << endl;

    unique_ptr<double[]> tmp( new double[cols] );
    double* tmpRow = tmp.get();
       
    // forward elimination with pivoting
    size_t rank(0);
    for( size_t j(0); j<rows; ++j ) {          
        // find maximum of column j
        double col_max = fabs( A[j*cols+j] );
        size_t row_with_max = j;
        for( size_t i(j+1); i<rows; ++i ) {
            double aij = fabs( A[i*cols+j] );
            if( aij > col_max ) {
                col_max = aij;
                row_with_max = i;
            }
        }
        if( row_with_max != j ) { // swap rows
            memcpy( tmpRow, A+j*cols, cols*sizeof(double) );
            memcpy( A+j*cols, A+row_with_max*cols, cols*sizeof(double) );
            memcpy( A+row_with_max*cols, tmpRow, cols*sizeof(double) );
        }
        double* ratj = A + j*cols + j;
        double ajj = A[j*cols+j];
        if( fabs(ajj) > 0.0 ) {
            rank++;
            size_t count = cols-j;
            for( size_t i(j+1); i<rows; ++i ) {
                double* rati = A + i*cols + j;
                double aij = A[i*cols+j] / ajj;
                if( fabs(aij) > 0.0 ) {
                    std::transform( rati, rati+count, ratj, rati, [&]( double& a, double& b ){ return a-aij*b; } );
                }
            }
        }
    }

cout << " L = " << __LINE__ << "   rank = " << rank << endl;
    if( rank < 2 ) return;
    
    // reverse elimination
    for( size_t j(rank-1); j>=1; --j ) {          
        double* ratj = A + j*cols + j;
        double ajj = A[j*cols+j];
        if( fabs(ajj) > 0.0 ) {
            size_t count = cols-j;
            for( size_t i(0); i<j; ++i ) {
                double* rati = A + i*cols + j;
                double aij = A[i*cols+j] / ajj;
                if( fabs(aij) > 0.0 ) {
                    std::transform( rati, rati+count, ratj, rati, [&]( double& a, double& b ){ return a-aij*b; } );
                }
            }
        }
    }

    
    if( transposed ) transpose( A, rows, cols );
        
}


void redux::math::gauss_jordan( Array<double>& A ) {
    const vector<size_t>& adims = A.dimensions();
    if( adims.size() != 2 ) {
        cout << "redux::math::gauss_jordan:: size of A is not ok: " << printArray(adims,"adims") << endl;
        return;
    }
    gauss_jordan( A.get(), adims[0], adims[1] );
}

void redux::math::lu_decomp( double* A, int rows, int cols ) {
    
    size_t N = std::max( rows, cols );
    bool do_transpose(false);
    if( rows > cols ) do_transpose = true;
    
    if( do_transpose ) transpose( A, rows, cols );
    
    gsl_matrix* a = gsl_matrix_alloc( N, N );
    gsl_permutation* p = gsl_permutation_alloc(N);

    memcpy( a->data, A, rows*cols*sizeof(double) );
    
    int signum;
    gsl_linalg_LU_decomp( a, p, &signum );

    memcpy( A, a->data, rows*cols*sizeof(double) );

    if( do_transpose ) transpose( A, cols, rows );
    
    gsl_matrix_free(a);
    gsl_permutation_free(p);
    
}

   
void redux::math::lu_decomp( Array<double>& A ) {
    const vector<size_t>& adims = A.dimensions();
    if( adims.size() != 2 ) {
        cout << "redux::math::lu_decomp:: size of A is not ok: " << printArray(adims,"adims") << endl;
        return;
    }
    lu_decomp( A.get(), adims[0], adims[1]);
}
   
   
void redux::math::qr_decomp(const double* A, int rows, int cols, double* Q, double* R) {
    
    gsl_matrix* a = gsl_matrix_alloc(rows, cols);
    gsl_vector* tau = gsl_vector_alloc(std::min(rows,cols));
    
    memcpy(a->data,A,rows*cols*sizeof(double));

    gsl_linalg_QR_decomp(a, tau);
    
    gsl_matrix q,r;
    q.data = Q;
    r.data = R;
    q.size1 = r.size1 = rows;
    q.tda = q.size2 = rows;
    r.tda = r.size2 = cols;

    gsl_linalg_QR_unpack (a, tau, &q, &r);

    
    gsl_vector_free(tau);
    gsl_matrix_free(a);
    
}
   
   
void redux::math::qr_decomp(const Array<double>& A, Array<double>& Q, Array<double>& R) {
    
    const vector<size_t>& adims =  A.dimensions();
    
    if( adims.size() != 2 /*|| adims[0] >= adims[1]*/ ) {
        cout << "redux::math::qr_decomp():: size of A is not ok: " << printArray(adims,"adims") << endl;
        return;
    }
    Q.resize(adims[0],adims[0]);
    R.resize(adims[0],adims[1]);
    qr_decomp(A.get(),adims[0],adims[1],Q.get(),R.get());
}


void redux::math::qr_decomp_pivot(const double* A, int rows, int cols, double* Q, double* R, gsl_permutation* p) {
    
    gsl_matrix* a = gsl_matrix_alloc(rows, cols);
    gsl_vector* tau = gsl_vector_alloc(std::min(rows,cols));
    gsl_vector* tmp = gsl_vector_alloc(std::max(rows,cols));
    
    int signum;
    
    memcpy(a->data,A,rows*cols*sizeof(double));

    gsl_linalg_QRPT_decomp(a, tau, p, &signum, tmp);
    
    gsl_matrix q,r;
    q.data = Q;
    r.data = R;
    q.size1 = r.size1 = rows;
    q.tda = q.size2 = rows;
    r.tda = r.size2 = cols;

    gsl_linalg_QR_unpack (a, tau, &q, &r);

    
    gsl_vector_free(tau);
    gsl_matrix_free(a);
    
}


void redux::math::qr_decomp_pivot(const double* A, int rows, int cols, double* Q, double* R) {
    gsl_permutation* p = gsl_permutation_alloc(std::max(rows,cols));
    //gsl_permutation_init(p);
    qr_decomp_pivot(A, rows, cols, Q, R, p);
    gsl_permutation_free(p);
}


void redux::math::qr_decomp_pivot(const Array<double>& A, Array<double>& Q, Array<double>& R, gsl_permutation* p) {
    const vector<size_t>& adims =  A.dimensions();
    
    if( adims.size() != 2 || adims[0] >= adims[1] ) {
        cout << "redux::math::qr_decomp_pivot():: size of A is not ok: " << printArray(adims,"adims") << endl;
        return;
    }
    Q.resize(adims[0],adims[0]);
    R.resize(adims[0],adims[1]);
    qr_decomp_pivot(A.get(),adims[0],adims[1],Q.get(),R.get(),p);
    
}
   
   
void redux::math::cod_decomp( const double* A, int rows, int cols, double* Q, double* R, double* Z ) {

    gsl_matrix* a = gsl_matrix_alloc(rows, cols);
    gsl_vector* tauQ = gsl_vector_alloc(std::min(rows,cols));
    gsl_vector* tauZ = gsl_vector_alloc(std::min(rows,cols));
    gsl_vector* work = gsl_vector_alloc( cols );
    gsl_permutation* p = gsl_permutation_alloc( cols );
    size_t rank;

    memcpy(a->data,A,rows*cols*sizeof(double));

    gsl_linalg_COD_decomp( a, tauQ, tauZ, p, &rank, work );

    gsl_matrix q,r,z;
    q.data = Q;
    r.data = R;
    z.data = Z;
    q.size1 = r.size1 = rows;
    q.tda = q.size2 = rows;
    r.tda = r.size2 = cols;
    z.tda = z.size1 = z.size2 = cols;

    gsl_linalg_COD_unpack( a, tauQ, tauZ, rank, &q, &r, &z );

    gsl_permutation_free(p);
    gsl_vector_free(work);
    gsl_vector_free(tauZ);
    gsl_vector_free(tauQ);
    gsl_matrix_free(a);
    
}
   
   
void redux::math::cod_decomp( const Array<double>& A, Array<double>& Q, Array<double>& R, Array<double>& Z ) {
    
    const vector<size_t>& adims =  A.dimensions();
    
    if( adims.size() != 2 /*|| adims[0] >= adims[1]*/ ) {
        cout << "redux::math::cod_decomp():: size of A is not ok: " << printArray(adims,"adims") << endl;
        return;
    }
    Q.resize( adims[0], adims[0] );
    R.resize( adims[0], adims[1] );
    Z.resize( adims[1], adims[1] );
    cod_decomp( A.get(), adims[0], adims[1], Q.get(), R.get(), Z.get() );
}


void redux::math::svd(double* A_U, int rows, int cols, double* S, double* V) {
 
    if( rows < cols ) {
        transpose(A_U,rows,cols);
        svd(A_U,cols,rows,S,V);
        return;
    }
    
    gsl_matrix a_u_v[2];
    gsl_vector s;
    
    memset(a_u_v,0,2*sizeof(gsl_matrix));
    memset(&s,0,sizeof(gsl_vector));
    
    a_u_v[0].data = A_U;
    a_u_v[1].data = V;
    s.data = S;

    a_u_v[0].size1 = rows;
    a_u_v[0].size2 = a_u_v[0].tda = a_u_v[1].size1 = a_u_v[1].size2 = a_u_v[1].tda = s.size = cols;
    s.stride = 1;
   
    gsl_vector* work = gsl_vector_alloc(cols);
    
    gsl_linalg_SV_decomp( a_u_v, a_u_v+1, &s, work );
    
    gsl_vector_free(work);
    
}
