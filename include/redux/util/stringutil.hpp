#ifndef REDUX_UTIL_STRINGUTIL_HPP
#define REDUX_UTIL_STRINGUTIL_HPP

#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>

namespace redux {

    namespace util {


        /*!  @ingroup util
         *  @{
         */

        
        /*!  @file      stringutil.hpp
         *   @brief     Collection of functions for string-manipulation
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2013
         */


        /*! @fn bool onlyDigits( const std::string &s );
         *  @fn bool onlyAlpha( const std::string &s );
         *  @brief Functions for checking string-content
         */
        bool onlyDigits( const std::string &s );
        bool onlyAlpha( const std::string &s );
        bool onlyAlnum( const std::string &s );
        bool onlyHex( const std::string &s );
        bool isInteger( const std::string &s );
        bool isHex( const std::string &s );


        /*! @fn std::string alignCenter( const std::string& s, size_t n=20, unsigned char c=' ' )
         *  @brief Append char 'c' to both sides of s, and form a block of width n
         *  @param s input string
         *  @param n block width
         *  @param c fill-character
         */
        std::string alignCenter( const std::string& s, size_t n = 20, unsigned char c = ' ' );

        /*! @fn std::string alignLeft( const std::string& s, size_t n=20, unsigned char c=' ' )
         *  @brief Append char 'c' to the right of s, and form a block of width n
         *  @param s input string
         *  @param n block width
         *  @param c fill-character
         */
        std::string alignLeft( const std::string& s, size_t n = 20, unsigned char c = ' ' );

        /*! @fn std::string alignRight( const std::string& s, size_t n=20, unsigned char c=' ' )
        *  @brief Append char 'c' to the left of s, and form a block of width n
        *  @param s input string
        *  @param n block width
        *  @param c fill-character
        */
        std::string alignRight( const std::string& s, size_t n = 20, unsigned char c = ' ' );

        std::string getUname( __uid_t id=0 );
        std::string cleanPath( std::string path, std::string base="" );

        /*! @fn std::string hexString( const T& v, bool prefix=true )
         *  @brief Converts integer types to a hexadecimal std::string
         */
        template <typename T>
        std::string hexString( const T& v, bool prefix = true ) {

            std::ostringstream oss;
            if( prefix ) {
                oss << "0x";
            }
            oss << ( std::hex ) << ( std::noshowbase ) << ( uint64_t )v;
            return oss.str();

        }


        /*! @fn template <typename T> std::string bitString( T var )
         *  @brief Print the bit-content of a variable as an std::string, e.g bitString((int)4) returns "00000000 00000000 00000000 00000100"
         */
        template <typename T> std::string bitString( T var ) {

            std::string tmp;
            uint8_t* begin = ( uint8_t* )( &var );
            uint8_t* end = begin + ( sizeof( T ) - 1 );

            while( end >= begin ) {
                for( uint8_t j( 128 ); j; j >>= 1 ) {
                    tmp += ( *end & j ) ? '1' : '0';
                }
                if( end > begin ) {
                    tmp += ' ';
                }
                end--;
            }

            return tmp;
        }


        /*! @fn std::string bitString( T* var, size_t n=1 )
         *  @brief Print the bit-content of an array of data as an std::string \n
         *  e.g uint8_t a[] = {0,1,2,3} => bitString(aa,4) returns "00000011 00000010 00000001 00000000"
         *  @param var pointer to input
         *  @param n number of elements to iterate over
         */
        template <typename T> std::string bitString( T* var, size_t n ) {

            std::string tmp;
            uint8_t* begin = reinterpret_cast<uint8_t*>( var );
            uint8_t* end = begin + n * sizeof( T ) - 1;

            while( end >= begin ) {
                for( uint8_t j( 128 ); j; j >>= 1 ) {
                    tmp += ( *end & j ) ? '1' : '0';
                }
                if( end > begin ) {
                    tmp += ' ';
                }
                end--;
            }

            return tmp;
        }


        /*! @fn std::string printBits( T* var, size_t n=1 )
         *  @brief Generate a std::string showing the memory layout of var
         *  @param var pointer to input
         *  @param n number of elements to be printed
         */
        template <typename T> inline std::string printBits( T* var, size_t n = 1 ) {

            std::string tmp( "     MSB <<  7  6  5  4  3  2  1  0  << LSB\n" );
            unsigned char* ptr = reinterpret_cast<unsigned char*>( var );
            unsigned char* end = ptr + n * sizeof( T ) - 1;

            while( end >= ptr )  {
                bool alpha = ( *end > 31 && *end < 126 );
                tmp += alignRight( "Byte ", 7 ) + alignLeft( std::to_string( ( int )( end - ptr ) ), 5 );

                for( int j = 128; j; j >>= 1 ) {
                    tmp += alignCenter( ( *end & j ) ? "1" : "0", 3 );
                }

                tmp += alignRight( std::to_string( ( int ) * end ), 4 ) + std::string( " (" );
                tmp += hexString( ( ( int ) * end ) >> 4, false )
                       + hexString( ( int )( *end ) & 0xF, false );

                if( alpha ) {
                    tmp += std::string( "," ) + std::string( ( char* )end, 1 );
                }

                tmp += std::string( ")\n" );
                end--;
            }

            tmp += alignCenter( std::string( 24, '-' ), 48 ) + std::string( "\n" );

            return tmp;
        }

        
        /*! @fn std::string printBits( const T& var, size_t n=1 )
         *  @brief Generate a std::string showing the memory layout of var
         *  @param var input
         *  @param n number of elements to be printed
         */
        template <typename T> inline std::string printBits( const T& var, size_t n = 1 ) {

            std::string tmp( "     MSB <<  7  6  5  4  3  2  1  0  << LSB\n" );
            unsigned char* ptr = ( unsigned char* )( &var );
            unsigned char* end = ptr + n * sizeof( T ) - 1;

            while( end >= ptr )  {
                bool alpha = ( *end > 31 && *end < 126 );
                tmp += alignRight( "Byte ", 7 ) + alignLeft( std::to_string( ( int )( end - ptr ) ), 5 );

                for( int j = 128; j; j >>= 1 ) {
                    tmp += alignCenter( ( *end & j ) ? "1" : "0", 3 );
                }

                tmp += alignRight( std::to_string( ( int ) * end ), 4 ) + std::string( " (" );
                tmp += hexString( ( ( int ) * end ) >> 4, false )
                       + hexString( ( int )( *end ) & 0xF, false );

                if( alpha ) {
                    tmp += std::string( "," ) + std::string( ( char* )end, 1 );
                }

                tmp += std::string( ")\n" );

                end--;
            }

            tmp += alignCenter( std::string( 24, '-' ), 48 ) + std::string( "\n" );

            return tmp;
        }

        
        template <typename T>
        std::string printBlock( T* data, size_t n, int d = 3, const std::string& delimiters = "[]", const char separator = ',' ) {

            std::ostringstream oss;
            size_t mid = delimiters.length() >> 1;
            if( mid ) {
                oss << delimiters.substr( 0, mid );
            }

            oss << std::setprecision( d );
            bool notFirst( true );
            for( size_t i = 0; i < n; ++i ) {
                if( notFirst ) {
                    oss << separator;
                }
                oss << data[i];
                notFirst = true;
            }

            oss << delimiters[1];
            if( mid ) {
                oss << delimiters.substr( mid );
            }

            return oss.str();
        }


        //@{
        /*! @brief Print the content of a vector as data = [[ 1, 2], [3, 4], ..., [n-1, n] ]
         *  @details Simple way to output a vector in a form that can be copy/pasted into IDL etc.
         *  @param data reference to input or pointer to data
         *  @param n number of elements to print (for pointer version)
         *  @param s name to be printed
         *  @param d width/precision parameter to be passed on to toString()
         */
        template <typename T>
        std::string printArray( T* data, size_t n, const char* s = "vector", int d = 3 ) {

            std::ostringstream oss;
            oss.precision(d);
            oss << s << "=[";
            bool separator( false );
            for( size_t i = 0; i < n; ++i ) {
                if( separator ) {
                    oss << ", ";
                }
                oss << data[i];
                separator = true;
            }

            oss << "]";

            return oss.str();
        }

        
        template <typename T>
        inline std::string printArray( const T& data, const char* s = "vector", int d = 3 ) {

            std::ostringstream oss;
            oss << std::setprecision( d ) << s << "=[";
            bool separator( false );
            for( auto & it : data ) {
                if( separator ) {
                    oss << ", ";
                }
                oss << it;
                separator = true;
            }

            oss << "]";

            return oss.str();
        }
        //@}



        /*! @} */


    }

}

#endif // REDUX_UTIL_STRINGUTIL_HPP
