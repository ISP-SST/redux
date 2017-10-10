#ifndef REDUX_UTIL_STRINGUTIL_HPP
#define REDUX_UTIL_STRINGUTIL_HPP

#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <iostream>
#include <typeinfo>
#include <iomanip>

#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

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
        bool contains(const std::string & haystack, const std::string & needle, bool ignoreCase=false);
        std::string replace_n(std::string input, const std::string& loc, const std::string& replace, size_t n=1);
        bool nocaseLess(const std::string& lhs, const std::string& rhs);
        
        bool isRelative( const std::string& );
        inline bool isRelative( const bfs::path &p ) { return isRelative( p.string() ); }
        
        std::vector<std::set<std::string>> make_template( const std::vector<std::string>& list, std::string& out, std::string split_chars="." );

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

        std::string getUname( __uid_t id = 0 );
        std::string cleanPath( std::string path, std::string base = "" );
        
        void printProgress(const std::string& text, float progress);
        
        template <typename T=uint32_t>
        std::vector<T> stringToUInts(const std::string& str);
        template <typename T>
        std::string uIntsToString(const std::vector<T>& ints);

        /*! @brief Converts integer types to a hexadecimal std::string
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


        /*! @name bitString
         *  @brief Return the bit-content of a variable as an std::string.
         */
        //@{
        /*!
         *  @param var pointer to input
         *  @code 
         *  int64_t foo = (1L<<60);
         *  cout << bitString(foo) << endl;
         *  @endcode
         *  will output
         *  @code
         *  00010000 00000000 00000000 00000000 00000000 00000000 00000000 00000000
         *  @endcode
         */
        template <typename T> std::string bitString( const T& var ) {

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


        /*! 
         *  @param var pointer to input
         *  @param n number of elements to iterate over
         *  @code 
         *  char foo[] = {0,0,0,0,0,0,0,16};
         *  cout << bitString(foo,8) << endl;
         *  @endcode
         *  will output
         *  @code
         *  00010000 00000000 00000000 00000000 00000000 00000000 00000000 00000000
         *  @endcode
         */
        template <typename T> std::string bitString( const T* var, size_t n ) {

            std::string tmp;
            const uint8_t* begin = reinterpret_cast<const uint8_t*>( var );
            const uint8_t* end = begin + n * sizeof( T ) - 1;

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
        //@}

        /*! @name
         *  @brief Generate a std::string showing the memory layout of var
         *  @param var pointer to input
         *  @param n number of elements to be printed
         *  @code 
         *  int64_t foo = (1L<<60);
         *  cout << printBits(foo) << endl;
         *  cout << endl << bitString(foo) << endl;
         *  @endcode
         *  will output
         *  @code
         *     MSB <<  7  6  5  4  3  2  1  0  << LSB
         *  Byte 7     0  0  0  1  0  0  0  0   16 (10)
         *  Byte 6     0  0  0  0  0  0  0  0    0 (00)
         *  Byte 5     0  0  0  0  0  0  0  0    0 (00)
         *  Byte 4     0  0  0  0  0  0  0  0    0 (00)
         *  Byte 3     0  0  0  0  0  0  0  0    0 (00)
         *  Byte 2     0  0  0  0  0  0  0  0    0 (00)
         *  Byte 1     0  0  0  0  0  0  0  0    0 (00)
         *  Byte 0     0  0  0  0  0  0  0  0    0 (00)
         *            ------------------------            
         *  @endcode
         */
        template <typename T> inline std::string printBits( T* var, size_t n = 1 ) {

            std::ostringstream oss;
            oss << "     MSB <<  7  6  5  4  3  2  1  0  << LSB\n";
            const unsigned char* ptr = reinterpret_cast<const unsigned char*>( var );
            const unsigned char* end = ptr + n * sizeof( T ) - 1;
            
            std::ios_base::fmtflags ff = oss.flags();
            
            while( end >= ptr )  {
                bool alpha = ( *end > 31 && *end < 126 );
                oss << "  Byte ";
                oss.setf ( std::ios::left );
                oss << std::setw(5) << std::to_string( ( int )( end - ptr ) );
                oss.flags( ff );

                for( int j = 128; j; j >>= 1 ) {
                    oss << (( *end & j ) ? " 1 " : " 0 ");
                }

                oss << " " << std::setw(5) << std::to_string( ( int ) * end );
                oss << std::string( " (" );
                oss << hexString( ( ( int ) * end ) >> 4, false ) << hexString( ( int )( *end ) & 0xF, false );

                if( alpha ) {
                    oss << ",'" << std::string( ( char* )end, 1 ) << "'";
                }

                oss << ")\n";
                end--;
            }

            oss <<  "            ------------------------            \n";

            return oss.str();

        }


        /*! @fn std::string printBits( const T& var, size_t n=1 )
         *  @brief Generate a std::string showing the memory layout of var
         *  @param var input
         *  @param n number of elements to be printed
         */
        template <typename T> inline std::string printBits( const T& var, size_t n = 1 ) {

            std::ostringstream oss;
            oss << "     MSB <<  7  6  5  4  3  2  1  0  << LSB\n";

            unsigned char* ptr = ( unsigned char* )( &var );
            unsigned char* end = ptr + n * sizeof( T ) - 1;
            
            std::ios_base::fmtflags ff = oss.flags();
            
            while( end >= ptr )  {
                bool alpha = ( *end > 31 && *end < 126 );
                oss << "  Byte ";
                oss.setf ( std::ios::left );
                oss << std::setw(5) << std::to_string( ( int )( end - ptr ) );
                oss.flags ( ff );

                for( int j = 128; j; j >>= 1 ) {
                    oss << (( *end & j ) ? " 1 " : " 0 ");
                }

                oss << " " << std::setw(5) << std::to_string( ( int ) * end ) << std::string( " (" );
                oss << hexString( ( ( int ) * end ) >> 4, false ) << hexString( ( int )( *end ) & 0xF, false );

                if( alpha ) {
                    oss << ",'" << std::string( ( char* )end, 1 ) << "'";
                }

                oss << ")\n";

                end--;
            }

            oss <<  "            ------------------------            \n";

            return oss.str();

        }


        /*! @name printArray
         *  @details Simple way to output arrays in a form that can be copy/pasted into IDL/Matlab etc.
         */
        //@{
        /*! @brief Print the content of a matrix as name = [[ 1, 2], [3, 4], ..., [n-1, n]]
         *  @param data Pointer to data (second order)
         *  @param firstY First row index.
         *  @param lastY  Last row index.
         *  @param firstX First column index.
         *  @param lastX  Last column index.
         *  @param name Name to be printed
         *  @param delimiters List with two delimiters (column/row separators). The square brackets are always added around each row.
         *  @param d width/precision parameter to be passed on to stringstream.precision()
         *  @code
         *  cout << printArray(myData,2,4,2,4,"foo", {",",";\n      "}) << endl;
         *  @endcode
         *  will output something like
         *  @code
         *  foo = [[-0.433, 0.852, -0.294];
         *        [-0.477, -0.493, -0.727];
         *        [-0.764, -0.175, 0.62]]
         *  @endcode
         */
        template <typename T>
        std::string printArray( T** data, size_t firstY, size_t lastY, size_t firstX, size_t lastX, const std::string& name = "matrix", const std::vector<std::string>& delimiters = {",",","}, int d = 3 ) {
            std::ostringstream oss;
            oss.precision( d );
            if( name.empty() ) oss << "[";
            else oss << name << "=[";
            bool rowSeparator( false );
            for( size_t i = firstY; i <= lastY; ++i ) {
                bool elementSeparator( false );
                    if( rowSeparator ) {
                        oss << delimiters[1] << " ";
                    } else {
                        rowSeparator = true;
                    }
                    oss << "[";
                for( size_t j = firstX; j <= lastX; ++j ) {
                    if( elementSeparator ) {
                        oss << delimiters[0] << " ";
                    } else elementSeparator = true;
                    oss << data[i][j];
                }
                oss << "]";
            }
            oss << "]";
            return oss.str();
        }

        /*! @brief Print the content of a matrix as name = [[ 1, 2], [3, 4], ..., [n-1, n]]
         *  @param data Pointer to data (second order)
         *  @param m First (slow/row) dimension size.
         *  @param n Second (fast/column) dimension size.
         *  @param name Name to be printed
         *  @param delimiters List with two delimiters (column/row separators). The square brackets are always added around each row.
         *  @param d width/precision parameter to be passed on to stringstream.precision()
         *  @code
         *  cout << printArray(myData,3,3,"foo", {",",";\n      "}) << endl;
         *  @endcode
         *  will output something like
         *  @code
         *  foo = [[-0.433, 0.852, -0.294];
         *        [-0.477, -0.493, -0.727];
         *        [-0.764, -0.175, 0.62]]
         *  @endcode
         */
        template <typename T>
        std::string printArray( T** data, size_t m, size_t n, const std::string& name = "matrix", const std::vector<std::string>& delimiters = {",",","}, int d = 3 ) {
            std::ostringstream oss;
            oss.precision( d );
            if( name.empty() ) oss << "[";
            else oss << name << "=[";
            bool rowSeparator( false );
            for( size_t i = 0; i < m; ++i ) {
                bool elementSeparator( false );
                    if( rowSeparator ) {
                        oss << delimiters[1] << " ";
                    } else {
                        rowSeparator = true;
                    }
                    oss << "[";
                for( size_t j = 0; j < n; ++j ) {
                    if( elementSeparator ) {
                        oss << delimiters[0] << " ";
                    } else elementSeparator = true;
                    oss << data[i][j];
                }
                oss << "]";
            }
            oss << "]";
            return oss.str();
        }


        /*! @brief Print the content of an array as name = [ 1, 2, 3, ..., n ]
         *  @param data Pointer to data (second order)
         *  @param n Size of array.
         *  @param name Name to be printed
         *  @param d width/precision parameter to be passed on to stringstream.precision()
         *  @code
         *  cout << printArray(myData,9,"foo", 5) << endl;
         *  @endcode
         *  will output something like
         *  @code
         *  foo = [-0.43326, 0.85281, -0.29473, -0.47733, -0.49376, -0.72729, -0.76411, -0.17500, 0.62737]
         *  @endcode
         */
        template <typename T>
        std::string printArray( T* data, size_t n, const std::string& name = "vector", int d = 3 ) {

            std::ostringstream oss;
            oss.precision( d );
            if( name.empty() ) oss << "[";
            else oss << name << "=[";
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


        /*! @brief Print the content of an STL container (or anything that has begin()/end() iterator initializers as name = [ 1, 2, 3, ..., n ]
         *  @param data Reference to container
         *  @param name Name to be printed
         *  @param d width/precision parameter to be passed on to stringstream.precision()
         *  @code
         *  cout << printArray(myData,"foo") << endl;
         *  @endcode
         *  will output something like
         *  @code
         *  foo = [-0.433, 0.881, -0.293, -0.433, -0.493, -0.729, -0.764, -0.175, 0.627]
         *  @endcode
         */
        template <typename T, typename U>
        inline std::string printArray( const std::map<T,U>& data, const std::string& name = "vector", int d = 3 ) {

            std::ostringstream oss;
            oss << std::setprecision( d );
            if( !name.empty() ) oss << name << "=[";
            bool separator( false );
            for( auto & element : data ) {
                if( separator ) {
                    oss << ",";
                }
                //if( std::is_integral<T>::value ) {
                //    oss << +element.first;   // trick to promote char to int before output to avoid interpretation as character
                //} else {
                    oss << element.first;                    
                //}
                oss << "->";
                if( std::is_integral<U>::value ) {
                    oss << +element.second;   // trick to promote char to int before output to avoid interpretation as character
                } else {
                    oss << element.second;                    
                }
                separator = true;
            }

            if( !name.empty() ) oss << "]";

            return oss.str();
        }
        template <typename T>
        inline std::string printArray( const T& data, const std::string& name = "vector", int d = 3 ) {

            std::ostringstream oss;
            oss << std::setprecision( d );
            if( !name.empty() ) oss << name << "=[";
            bool separator( false );
            for( auto & element : data ) {
                if( separator ) {
                    oss << ",";
                }
                //if( std::is_integral<decltype(it)>::value ) {
                //    oss << +it;   // FIXME for some reason it=string passes as integer...weird
                //} else {
                    oss << element;
                //}
                separator = true;
            }

            if( !name.empty() ) oss << "]";

            return oss.str();
        }
        //@}


        enum StringColor { BLACK = 30, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE,
                           BBLACK= 90, BRED, BGREEN, BYELLOW, BBLUE, BMAGENTA, BCYAN, BWHITE
        };
        /* ! @fn std::string colorString( const std::string& t, color col = WHITE )
        *  @brief Function to wrap console color-codes around a std::string.
        *  @param t input variable
        *  @param col color
        *  @param n width/precision parameter to be passed on to toString()
        */
        std::string colorString( const std::string& in, StringColor col=WHITE );

        /*! @fn string tvToString ( const timeval& a, bool millis )
        *  @brief Convert timeval to a string of the form: "HH:MM:SS.mmm"
        *  @param a timeval structure.
        *  @param millis Should milliseconds be included ?
        *  @returns "HH:MM:SS.mmm"
        */
        std::string tvToString( const timeval&, bool millis = false );
                
                
        /*! @fn string tsToString ( const timespec& a, bool millis )
        *  @brief Convert timeval to a string of the form: "HH:MM:SS.mmm"
        *  @param a timespec structure.
        *  @param millis Should milliseconds be included ?
        *  @returns "HH:MM:SS.mmm"
        */
        std::string tsToString( const timespec&, bool millis = false );


        /*! @} */


    }

}

#endif // REDUX_UTIL_STRINGUTIL_HPP
