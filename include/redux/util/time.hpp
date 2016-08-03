#ifndef REDUX_UTIL_TIME_HPP
#define REDUX_UTIL_TIME_HPP

#include <sys/time.h>

namespace redux {
    
    namespace util {


        /*!  @ingroup util
        *  @{
        */

        timespec operator+ ( const timespec& a, const timespec& b );
        timespec& operator+= ( timespec& a, const timespec& b );
        timespec operator- ( const timespec& a, const timespec& b );
        timespec operator+ ( const timespec& a, const long& b );
        timespec& operator+= ( timespec& a, const long& b );
        timespec operator- ( const timespec& a, const long& b );
    //timespec operator+ ( const timespec& a, const double& b );
        bool operator< ( const timespec& a, const timespec& b );
        bool operator> ( const timespec& a, const timespec& b );
        bool operator<= ( const timespec& a, const timespec& b );
        bool operator>= ( const timespec& a, const timespec& b );



        /*! @} */


    } // end namespace util

} // end namespace redux


#endif // REDUX_UTIL_TIME_HPP

