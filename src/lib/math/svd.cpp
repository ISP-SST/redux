#include "redux/math/svd.hpp"

#include "redux/util/datautil.hpp"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <iostream>
using namespace redux::math;
using namespace redux::util;
using namespace std;
void redux::math::svd(double* A_U, int rows, int cols, double* S, double* V) {
 
    if( rows < cols ) {
        //cout << "." << flush;
        transpose(A_U,rows,cols);
        //cout << "-" << flush;
        svd(A_U,cols,rows,S,V);
        //cout << "|" << flush;
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
   
/*        
            gsl_matrix AA;
    gsl_matrix AAA;
    gsl_matrix UU;
    gsl_matrix VV;
    gsl_vector* SS;
    gsl_vector* work;

    SS = gsl_vector_alloc( nn );
    work = gsl_vector_alloc( nn );

    AA.block = AAA.block = UU.block = VV.block = 0;
    AA.owner = AAA.owner = UU.owner = VV.owner = 0;

    AA.size1 = AAA.size1 = UU.size1 = UU.size2 = UU.tda = mm;
    AA.size2 = AAA.size2 = VV.size1 = VV.size2 = nn;
    AA.tda = AAA.tda = VV.tda = nn;
    AA.data = *dataP2;
    AAA.data = *sigmaP2;
    UU.data = *uP2;
    VV.data = *vP2;
    {
        boost::timer::auto_cpu_timer timer;
        for( size_t n = 0; n < 1000000; ++n ) {
            memcpy( *sigmaP2, *dataP2, 30 * sizeof( double ) );
            gsl_linalg_SV_decomp( &AAA, &VV, SS, work );
            //svdcmp( sigmaP, mm, nn, uP, vP );
        }
        cout << "gsl_linalg_SV_decomp:   " << endl;
    }
  */      