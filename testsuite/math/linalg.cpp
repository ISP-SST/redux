#include "redux/math/linalg.hpp"

#include "redux/util/array.hpp"
#include "redux/util/stringutil.hpp"

#include <boost/test/unit_test.hpp>
//#include <boost/timer/timer.hpp>

using namespace redux::math;
using namespace redux::util;

using namespace std;

namespace testsuite {

    namespace math {

        void qrTest(void) {
            
            const int ROWS = 512;
            const int COLS = 512;

            Array<double> A(ROWS, COLS);
            Array<double> Q(ROWS, ROWS);
            Array<double> R(ROWS, COLS);
            double** AA = newArray<double>(ROWS, COLS);
            //double** Q = newArray<double>(ROWS, ROWS);
            //double** R = newArray<double>(ROWS, COLS);
            // set Lotkin matrix
            // a_rc = 1 (r = 0) or 1/(r+c-1) (r > 0) (i.e. Hilbert matrix with first row set to 1)
            for(int c = 0; c < COLS; ++c) A(0,c) = 1;
            for(int r = 1; r < ROWS; ++r) {
                for(int c = 0; c < COLS; ++c) {
                    A(r,c) = 1.0 / (double)(r + c + 1);
                }
            }

        #if 0
            // set Frank matrix
            // a_ij = DIM - min(i,j) + 1
            for(int r = 1; r < ROWS; ++r) {
                for(int c = 0; c < COLS; ++c) {
                    A(r, c) = ROWS - std::min(r, c) + 1;
                }
            }
        #endif
            memcpy(*AA,A.get(),ROWS*COLS*sizeof(double));
            if(ROWS < 20 && COLS < 20) {
                printf("Matrix(%d,%d)\n", ROWS, COLS);
                for(int r = 0; r < ROWS; ++r) {
                    printf("%3d: ", r);
                    for(int c = 0; c < COLS; ++c) {
                        cout << alignRight(to_string(A(r,c)), 8) << "  ";
                    }
                    printf("\n");
                }
                printf("\n");
            }
            /*double** q;
            {
                //boost::timer::auto_cpu_timer timer;
                q = legacy::qr(AA, ROWS, COLS, 0);
                printf("Q(%d,%d)\n", ROWS, ROWS);
                if(ROWS < 20 && COLS < 20) {
                    for(int r = 0; r < ROWS; ++r) {
                        printf("%3d: ", r);
                        for(int c = 0; c < ROWS; ++c) {
                            cout << alignRight(to_string(q[r][c]), 8) << "  ";
                        }
                        printf("\n");
                    }
                    printf("\n");
                }

            }*/

            //memcpy(*AA,*A,ROWS*COLS*sizeof(double));
            {
                //boost::timer::auto_cpu_timer timer;
                qr_decomp(A.get(), ROWS, COLS, Q.get(), R.get());
                if(ROWS < 20 && COLS < 20) {
                    printf("Q2(%d,%d)\n", ROWS, ROWS);
                    for(int r = 0; r < ROWS; ++r) {
                        printf("%3d: ", r);
                        for(int c = 0; c < ROWS; ++c) {
                            cout << alignRight(to_string(Q(r,c)), 8) << "  ";
                        }
                        printf("\n");
                    }
                    printf("\n");
                }
            }

            //cout << printArray(AA,0,ROWS-1,0,COLS-1,"R1") << endl;
            //cout << printArray(q,0,ROWS-1,0,ROWS-1,"Q1") << endl;
            //cout << printArray(A,0,ROWS-1,0,COLS-1,"A") << endl;
            //cout << printArray(Q,0,ROWS-1,0,ROWS-1,"Q") << endl;
            //cout << printArray(R,0,ROWS-1,0,COLS-1,"R") << endl;


            //delArray(A);
            //delArray(R);
            //delArray(Q);
            //delArray(q);
            delArray(AA);

        }
        
        using namespace boost::unit_test;
        void add_linalg_tests( test_suite* ts ) {
            
            ts->add( BOOST_TEST_CASE_NAME( &qrTest, "QR decomposition" ) );

        }

    }

}

