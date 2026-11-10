#include "redux/math/interval.hpp"

#include "catch.hpp"

using namespace redux::math;

namespace testsuite {

    namespace math {

        TEMPLATE_TEST_CASE( "Interval class", "[math]", uint8_t, int, size_t, double ) {

                using T = TestType;
                Interval<T> interval(3, 9);
                CHECK( interval.start() == 3 );
                CHECK( interval.end() == 9 );

                CHECK( interval.contains(5) );
                CHECK_FALSE( interval.contains(1) );
                CHECK_FALSE( interval.contains(20) );

                CHECK( interval.surrounds(Interval<T>(4, 8)) );
                CHECK_FALSE( interval.inside(Interval<T>(4, 8)) );
                CHECK( interval.inside(Interval<T>(1, 10)) );
                CHECK_FALSE( interval.surrounds(Interval<T>(1, 10)) );

                CHECK( interval.intersects(Interval<T>(0, 3)) );
                CHECK( interval.intersects(Interval<T>(9, 15)) );
                CHECK( interval.intersects(Interval<T>(1, 15)) );
                CHECK_FALSE( interval.intersects(Interval<T>(10, 15)) );
                CHECK_FALSE( interval.intersects(Interval<T>(0, 2)) );

                CHECK( (interval == Interval<T>(3, 9)) );
                CHECK_FALSE( (interval == Interval<T>(2, 9)) );

                CHECK_FALSE( (interval != Interval<T>(3, 9)) );
                CHECK( (interval != Interval<T>(2, 9)) );

                CHECK( (interval > Interval<T>(0, 2)) );
                CHECK_FALSE( (interval > Interval<T>(10, 20)) );
                CHECK_FALSE( (interval < Interval<T>(0, 2)) );
                CHECK( (interval < Interval<T>(10, 20)) );

                CHECK_THROWS_AS( interval.join(Interval<T>(10, 20)), std::logic_error );
                CHECK_THROWS_AS( interval.joined(Interval<T>(10, 20)), std::logic_error );
                CHECK_THROWS_AS( interval.intersection(Interval<T>(10, 20)), std::logic_error );

                CHECK( (interval.from(5) == Interval<T>(5, 9)) );
                CHECK( (interval.until(5) == Interval<T>(3, 5)) );

                CHECK_THROWS_AS( interval.from(10), std::logic_error );
                CHECK_THROWS_AS( interval.until(10), std::logic_error );

        }

    }

}
