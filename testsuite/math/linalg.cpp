#include "redux/math/linalg.hpp"

#include "redux/util/array.hpp"
#include "redux/util/stringutil.hpp"

#include "catch.hpp"

using namespace redux::math;
using namespace redux::util;

using namespace std;

namespace testsuite {

    namespace math {

        TEST_CASE( "QR decomposition of a Lotkin matrix", "[math]" ) {

            const int ROWS = 512;
            const int COLS = 512;

            Array<double> A(ROWS, COLS);
            Array<double> Q(ROWS, ROWS);
            Array<double> R(ROWS, COLS);

            for(int c = 0; c < COLS; ++c) A(0,c) = 1;
            for(int r = 1; r < ROWS; ++r) {
                for(int c = 0; c < COLS; ++c) {
                    A(r,c) = 1.0 / (double)(r + c + 1);
                }
            }

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

            CHECK_NOTHROW( qr_decomp(A.get(), ROWS, COLS, Q.get(), R.get()) );

            if(ROWS < 20 && COLS < 20) {
                printf("Q(%d,%d)\n", ROWS, ROWS);
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

    }

}
