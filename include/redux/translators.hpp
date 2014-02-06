#ifndef REDUX_TRANSLATORS_HPP
#define REDUX_TRANSLATORS_HPP

#include "redux/momfbd/defines.hpp"
#include "redux/util/stringutil.hpp"

#include <algorithm>
#include <limits>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>


// property_tree custom translator for bool
template <typename Ch, typename Traits, typename Alloc>
class BoolTranslator {
public:
    typedef std::basic_string<Ch, Traits, Alloc> internal_type;
    typedef bool external_type;

    // Converts a string to bool
    boost::optional<external_type> get_value( const internal_type& str ) const {
        if( !str.empty() ) {
            using boost::algorithm::iequals;

            if( iequals( str, "true" ) || iequals( str, "yes" ) || str == "1" )
                return boost::optional<external_type>( true );
            else
                return boost::optional<external_type>( false );
        }
        else
            return boost::optional<external_type>( true );    // empty value => flag, so should default to true.
    }

    // Converts a bool to string
    boost::optional<internal_type> put_value( const external_type& b ) {
        return boost::optional<internal_type>( b ? "true" : "" );
    }
};

// property_tree custom translator for FileType
template <typename Ch, typename Traits, typename Alloc>
class FileTypeTranslator {
public:
    typedef std::basic_string<Ch, Traits, Alloc> internal_type;
    typedef redux::momfbd::FileType external_type;

    // Converts a string to FileType
    external_type get_value( const internal_type& str ) const {
        using namespace redux::momfbd;
        if( !str.empty() ) {
            using boost::algorithm::iequals;
            if( iequals( str, "ANA" ) )
                return FT_ANA ;
            else if( iequals( str, "FITS" ) )
                return FT_FITS;
            else if( iequals( str, "MOMFBD" ) )
                return FT_MOMFBD;
            else
                return FT_NONE;
        }
        else
            return FT_NONE;
    }

    // Converts a FileType to string
    internal_type put_value( const external_type& b ) {
        using namespace redux::momfbd;
        switch( b ) {
            case FT_ANA: return "ANA";
            case FT_FITS: return "FITS";
            case FT_MOMFBD: return "MOMFBD";
            default: return "";
        }
    }
};

// property_tree custom translator for vector
template <typename Ch, typename Traits, typename Alloc, typename T>
class VectorTranslator {
    typedef std::basic_string<Ch, Traits, Alloc> internal_type;
    typedef typename std::vector<T> external_type;
    void parseSegment( external_type& out, const internal_type& elem ) {
        size_t n = std::count( elem.begin(), elem.end(), '-' );
        if( n == 0 ) {
            out.push_back( boost::lexical_cast<T>( elem ) );
            return;
        }
        else if( n == 1 ) {
            n = elem.find_first_of( '-' );
            T first = boost::lexical_cast<T>( elem.substr( 0, n ) );
            T last = boost::lexical_cast<T>( elem.substr( n + 1 ) );
            while( first <= last ) out.push_back( first++ );
        }
    }
public:

    // Converts a string to a vector
    boost::optional<external_type> get_value( const internal_type& str ) {
        external_type result;
        if( !str.empty() ) {
            std::vector<std::string> tok;
            boost::split(tok, str, boost::is_any_of(",") );
            if( !std::numeric_limits<T>::is_signed ) {  // only expand hyphen ('-') for unsigned integers
                for( auto & it : tok ) {
                    parseSegment( result, it );
                }
            }
            else {
                std::transform( tok.begin(), tok.end(), back_inserter( result ), boost::lexical_cast<T, internal_type> );
            }
        }
        return boost::optional<external_type>( result );
    }

    // Converts a vector to a comma-separated string
    boost::optional<internal_type> put_value( const external_type& vec ) {
        std::ostringstream oss;
        if( !vec.empty() ) {
            std::copy( vec.begin(), vec.end() - 1, std::ostream_iterator<T>( oss, "," ) );
            oss << vec.back();
        }
        return oss.str();
    }
};

// specialization to handle a vector<FileType>
template <typename Ch, typename Traits, typename Alloc>
class VectorTranslator<Ch, Traits, Alloc, redux::momfbd::FileType> {
    typedef std::basic_string<Ch, Traits, Alloc> internal_type;
    typedef typename std::vector<redux::momfbd::FileType> external_type;

    FileTypeTranslator<Ch, Traits, Alloc> elementTranslator;
    void parseSegment( external_type& out, const internal_type& elem ) {
        redux::momfbd::FileType tmp = elementTranslator.get_value( elem );
        if( tmp )
            out.push_back( elementTranslator.get_value( elem ) );
    }
public:

    // Converts a string to a vector
    boost::optional<external_type> get_value( const internal_type& str ) {
        external_type result;
        if( !str.empty() ) {
            std::vector<std::string> tok;
            boost::split(tok, str, boost::is_any_of(",") );
            for( auto & it : tok ) {
                parseSegment( result, it );
            }
        }
        return boost::optional<external_type>( result );
    }

    // Converts a vector to a comma-separated string
    boost::optional<internal_type> put_value( const external_type& vec ) {
        std::ostringstream oss;
        for( auto it = vec.begin(); it < vec.end() - 1; it++ ) {
            oss << elementTranslator.put_value( vec ) << ",";
        }
//         if( !vec.empty() ) {
//             std::copy( vec.begin(), vec.end() - 1, std::ostream_iterator<T>( oss, "," ) );
//             oss << vec.back();
//         }
        return oss.str();
    }
};

namespace boost {
    namespace property_tree {

        template <typename Ch, typename Traits, typename Alloc>
        struct translator_between<std::basic_string<Ch, Traits, Alloc>, bool> {
            typedef BoolTranslator<Ch, Traits, Alloc> type;
        };

        template <typename Ch, typename Traits, typename Alloc>
        struct translator_between<std::basic_string<Ch, Traits, Alloc>, redux::momfbd::FileType> {
            typedef FileTypeTranslator<Ch, Traits, Alloc> type;
        };

        template<typename Ch, typename Traits, typename Alloc, typename T>
        struct translator_between<std::basic_string< Ch, Traits, Alloc >, std::vector<T>> {
            typedef VectorTranslator<Ch, Traits, Alloc, T> type;
        };

    } // namespace property_tree
} // namespace boost

#endif  // REDUX_TRANSLATORS_HPP
