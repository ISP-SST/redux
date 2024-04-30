#ifndef REDUX_VERSION_HPP
#define REDUX_VERSION_HPP

#include "redux/types.hpp"

#include <string>

namespace redux {
    
    /*! @defgroup redux REDUX
     *  @{
     */

    /*!  @file      version.hpp
     *   @brief     Functions for retrieving revision numbers/strings in some organized form.
     *   @name      Version/Revision information
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2013
     */

    /*!
     *   @fn         extern int getVersionNumber (void)
     *   @brief      Returns an integer representing the version
     *   @details    The different version numbers are packed into a 64-bit integer that can be used
     *               to compare versions. E.g. v. 1.2.3-9 would be returned as the integer 4328718345
     *               (= (1<<32)+(2<<24)+(3<<16)+9)
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2013
     */
    extern uint64_t getVersionNumber( void );


    /*!
     *   @fn         extern std::string getVersionString(void)
     *   @brief      Returns a human readable version string, e.g. "1.2.3"
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2013
     */
    extern std::string getVersionString( void );
    extern std::string getVersionString( uint64_t );


    /*!
     *   @fn         extern std::string getLongVersionString ( bool )
     *   @brief      Returns a string with version info and commit message.
     *   @details    E.g. If the build was made 7 commits after the last version-tag, the
     *               result would be something like: "1.2.3 (commit: 7 - Added implementation of getLongVersionString)"
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2013
     */
    extern std::string getLongVersionString( bool includeMessage=true );

    /*! @} */

}

#endif // REDUX_VERSION_HPP
