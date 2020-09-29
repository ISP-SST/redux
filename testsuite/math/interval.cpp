#include "redux/math/interval.hpp"

#include <boost/test/unit_test.hpp>

using namespace redux::math;

namespace testsuite {

    namespace math {
        
        template <typename T>
        void typedIntervalTest(void) {

            Interval<T> interval(3, 9);
            BOOST_CHECK_EQUAL(interval.start(), 3);
            BOOST_CHECK_EQUAL(interval.end(), 9);

            BOOST_CHECK(interval.contains(5));
            BOOST_CHECK(!interval.contains(1));
            BOOST_CHECK(!interval.contains(20));

            BOOST_CHECK(interval.surrounds(Interval<T>(4, 8)));
            BOOST_CHECK(!interval.inside(Interval<T>(4, 8)));
            BOOST_CHECK(interval.inside(Interval<T>(1, 10)));
            BOOST_CHECK(!interval.surrounds(Interval<T>(1, 10)));

            BOOST_CHECK(interval.intersects(Interval<T>(0, 3)));
            BOOST_CHECK(interval.intersects(Interval<T>(9, 15)));
            BOOST_CHECK(interval.intersects(Interval<T>(1, 15)));
            BOOST_CHECK(!interval.intersects(Interval<T>(10, 15)));
            BOOST_CHECK(!interval.intersects(Interval<T>(0, 2)));

            BOOST_CHECK((interval == Interval<T>(3, 9)));
            BOOST_CHECK(!(interval == Interval<T>(2, 9)));

            BOOST_CHECK(!(interval != Interval<T>(3, 9)));
            BOOST_CHECK((interval != Interval<T>(2, 9)));

            BOOST_CHECK((interval > Interval<T>(0, 2)));
            BOOST_CHECK(!(interval > Interval<T>(10, 20)));
            BOOST_CHECK(!(interval < Interval<T>(0, 2)));
            BOOST_CHECK((interval < Interval<T>(10, 20)));

            BOOST_CHECK_THROW(interval.join(Interval<T>(10, 20)), std::logic_error);
            BOOST_CHECK_THROW(interval.joined(Interval<T>(10, 20)), std::logic_error);
            BOOST_CHECK_THROW(interval.intersection(Interval<T>(10, 20)), std::logic_error);

            BOOST_CHECK((interval.from(5) == Interval<T>(5, 9)));
            BOOST_CHECK((interval.until(5) == Interval<T>(3, 5)));

            BOOST_CHECK_THROW(interval.from(10), std::logic_error);
            BOOST_CHECK_THROW(interval.until(10), std::logic_error);

        }

        void intervalTest(void) {
            
            // try for signed/unsigned/floats
            typedIntervalTest<uint8_t>();
            typedIntervalTest<int>();
            typedIntervalTest<size_t>();
            typedIntervalTest<double>();

        }

        
        using namespace boost::unit_test;
        void add_interval_tests( test_suite* ts ) {
            
            ts->add( BOOST_TEST_CASE_NAME( &intervalTest, "Interval class" ) );

        }

    }

}

