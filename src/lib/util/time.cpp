#include "redux/util/time.hpp"

//#include "redux/util/math.hpp"

#include <cstdlib>


using namespace redux::util;
using namespace std;


timespec redux::util::operator+ ( const timespec& a, const timespec& b ) {
    struct timespec c;
    c.tv_sec = a.tv_sec + b.tv_sec;
    c.tv_nsec = a.tv_nsec + b.tv_nsec;

    if ( c.tv_nsec > 999999999 ) {
        c.tv_sec++;
        c.tv_nsec -= 1000000000;
    }

    return c;
}


timespec& redux::util::operator+= ( timespec& a, const timespec& b ) {
    if ( &a == &b )
        return a;

    a.tv_sec  += b.tv_sec;
    a.tv_nsec += b.tv_nsec;

    if ( a.tv_nsec > 999999999 ) {
        a.tv_sec++;
        a.tv_nsec -= 1000000000;
    }

    return a;
}


timespec redux::util::operator- ( const timespec& a, const timespec& b ) {

    struct timespec c;
    c.tv_sec = abs( a.tv_sec - b.tv_sec );
    c.tv_nsec = abs( a.tv_nsec - b.tv_nsec );
    return c;
}


timespec redux::util::operator+ ( const timespec& a, const long& b ) {
    struct timespec c;
    c.tv_sec = abs( a.tv_sec + b );
    c.tv_nsec = a.tv_nsec;
    return c;
}


timespec& redux::util::operator+= ( timespec& a, const long& b ) {
    a.tv_sec = abs( a.tv_sec + b );
    return a;
}


timespec redux::util::operator- ( const timespec& a, const long& b ) {
    struct timespec c;
    c.tv_sec = abs( a.tv_sec - b );
    c.tv_nsec = a.tv_nsec;
    return c;
}


/*timespec operator+ ( const timespec& a, const double& b ) {
    struct timespec c;
    c.tv_sec = a.tv_sec;
    c.tv_nsec = ABS(a.tv_sec+(int)(b));
    return c;
}*/


bool redux::util::operator< ( const timespec& a, const timespec& b ) {
    if ( a.tv_sec < b.tv_sec )
        return true;
    else
        if ( a.tv_sec > b.tv_sec )
            return false;
        else
            return ( a.tv_nsec < b.tv_nsec );
}


bool redux::util::operator> ( const timespec& a, const timespec& b ) {
    if ( a.tv_sec > b.tv_sec )
        return true;
    else
        if ( a.tv_sec < b.tv_sec )
            return false;
        else
            return ( a.tv_nsec > b.tv_nsec );
}


bool redux::util::operator<= ( const timespec& a, const timespec& b ) {
    if ( a.tv_sec < b.tv_sec )
        return true;
    else
        if ( a.tv_sec > b.tv_sec )
            return false;
        else
            return ( a.tv_nsec <= b.tv_nsec );
}


bool redux::util::operator>= ( const timespec& a, const timespec& b ) {
    if ( a.tv_sec > b.tv_sec )
        return true;
    else
        if ( a.tv_sec < b.tv_sec )
            return false;
        else
            return ( a.tv_nsec >= b.tv_nsec );
}
