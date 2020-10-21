#ifndef REDUX_TRANSLATORS_HPP
#define REDUX_TRANSLATORS_HPP

#include "redux/momfbd/config.hpp"
#include "redux/util/stringutil.hpp"

#include <algorithm>
#include <limits>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/range/algorithm_ext/erase.hpp>

namespace redux {

    struct BoolTranslator {

        typedef std::string internal_type;

        // Converts a string to bool (i.e. when calling get<bool>(key) on a ptree)
        boost::optional<bool> get_value(const internal_type& str) const {
            if(!str.empty()) {
                using boost::algorithm::iequals;

                if(iequals(str, "true") || iequals(str, "yes") || str == "y" || str == "1") {
                    return boost::optional<bool>(true);
                }
                else {
                    return boost::optional<bool>(false);
                }
            }
            else {
                return boost::optional<bool>(true);      // empty value => flag, so should default to true.
            }
        }

        // Converts a bool to string (i.e. when calling put<bool>(key) on a ptree)
        internal_type put_value(const bool& b) {
            if(b) return "1";
            return "0";
        }
    };

// property_tree custom translator for FileType
    struct FileTypeTranslator {

        typedef std::string internal_type;
        typedef redux::momfbd::FileType external_type;

        // Converts a string to FileType
        external_type get_value(const internal_type& str) const {
            using namespace redux::momfbd;
            if(!str.empty()) {
                using boost::algorithm::iequals;
                if(iequals(str, "ANA"))
                    return FT_ANA ;
                else if(iequals(str, "FITS"))
                    return FT_FITS;
                else if(iequals(str, "MOMFBD"))
                    return FT_MOMFBD;
                else
                    return FT_NONE;
            }
            else
                return FT_NONE;
        }

        // Converts a FileType to string
        internal_type put_value(const external_type& b) {
            using namespace redux::momfbd;
            switch(b) {
                case FT_ANA: return "ANA";
                case FT_FITS: return "FITS";
                case FT_MOMFBD: return "MOMFBD";
                default: return "";
            }
        }
    };

    
    // custom translator for ModeBase
    struct ModeBaseTranslator {

        typedef std::string internal_type;
        typedef redux::momfbd::ModeBase external_type;

        // Converts a string to FileType
        external_type get_value( const internal_type& str ) const {
            using namespace redux::momfbd;
            if( !str.empty() ) {
                using boost::algorithm::iequals;
                if( iequals(str, "Zernike") ) {
                    return ZERNIKE;
                } else if( iequals(str, "Karhunen-Loeve") ) {
                    return KARHUNEN_LOEVE;
                }
            }
            return MB_NONE;
        }

        // Converts a FileType to string
        internal_type put_value( const external_type& b ) {
            using namespace redux::momfbd;
            switch(b) {
                case ZERNIKE: return "Zernike";
                case KARHUNEN_LOEVE: return "Karhunen-Loeve";
                default: return "";
            }
        }
    };

    
    struct DiversityValueTranslator {

        typedef std::string internal_type;
        typedef redux::momfbd::DiversityValue external_type;

        // Converts a string to DiversityValue
        external_type get_value(const internal_type& str) const {
            using namespace redux::momfbd;
            external_type ret;
            return ret;
        }

        // Converts a FileType to string
        internal_type put_value(const external_type& b) {
            using namespace redux::momfbd;
            internal_type ret;
            return ret;
        }
    };

    struct ModeIDTranslator {

        typedef std::string internal_type;
        typedef redux::momfbd::ModeID external_type;

        // Converts a string to ModeID
        external_type get_value(const internal_type& str) const {
            using namespace redux::momfbd;
            external_type ret;
            return ret;
        }

        // Converts a FileType to string
        internal_type put_value(const external_type& b) {
            using namespace redux::momfbd;
            internal_type ret;
            return ret;
        }
    };



    // property_tree custom translator for vectors
    template <typename T>
    class VectorTranslator {

        typedef std::string internal_type;
        typedef typename std::vector<T> external_type;

        void parseSegment(external_type& out, const internal_type& elem) {
            size_t n1 = std::count(elem.begin(), elem.end(), '-');
            size_t n2 = std::count(elem.begin(), elem.end(), ':');
            if( (n1+n2) == 0 ) {
                out.push_back(boost::lexical_cast<T>(elem));
                return;
            } else if( (n1+n2) == 1 ) {
                size_t n = elem.find_first_of(":-");
                T first = boost::lexical_cast<T>(elem.substr(0, n));
                T last = boost::lexical_cast<T>(elem.substr(n + 1));
                bool rev(false);
                if( first > last ) {
                    std::swap(first,last);
                    rev = true;
                }
                while(first <= last) out.push_back(first++);
                if( rev ) std::reverse( out.begin(), out.end() );
            } else if( n2 == 2 ) {
                size_t n = elem.find_first_of(":");
                size_t nn = elem.find_first_of(":", n + 1);
                T first = boost::lexical_cast<T>(elem.substr(0, n));
                int step = boost::lexical_cast<T>(elem.substr(n + 1, nn - n - 1));
                T last = boost::lexical_cast<T>(elem.substr(nn + 1));
                bool rev(false);
                if( first > last ) {
                    std::swap(first,last);
                    if( step < 0 ) step = abs(step);
                    rev = true;
                }
                while(first <= last) {
                    out.push_back(first);
                    first += step;
                }
                if( rev ) std::reverse( out.begin(), out.end() );
            }
        }


    public:

        // Converts a string to a vector
        external_type get_value(const internal_type& str) {
            external_type result;
            if(!str.empty()) {
                std::vector<std::string> tokens;
                boost::split(tokens, str, boost::is_any_of(", "));
                if(std::numeric_limits<T>::is_integer && !std::numeric_limits<T>::is_signed) {    // only expand '-' & ':' for unsigned integers
                    for(auto & token : tokens) {
                        if( token.empty() ) continue;
                        external_type res;
                        parseSegment( res, token );
                        result.insert( result.end(), res.begin(), res.end() );
                    }
                }
                else {
                    std::transform(tokens.begin(), tokens.end(), back_inserter(result), boost::lexical_cast<T, internal_type>);
                }
            }
            return result;
        }

        // Converts a vector to a comma-separated string
        internal_type put_value(const external_type& vec) {
            std::ostringstream oss;
            if(!vec.empty()) {
                if( std::numeric_limits<T>::is_integer && !std::numeric_limits<T>::is_signed ) {    // only expand '-' & ':' for unsigned integers
                    auto first = vec.begin();
                    while(first != vec.end()) {
                        
                        auto last = first + 1;
                        if(last == vec.end()) { // last element in vector
                            oss << *first;
                            break;
                        }

                        int step = *last - *first;
                        T next = *last + step;
                        size_t count = 2;
                        while(++last != vec.end() && *last == next) {
                            next += step;
                            count++;
                        }
                        if(count > 2) {
                            oss << *first;
                            if( step == 1 ) oss << "-";
                            else  oss << ":";
                            if( step != 1 ) oss << step << ":";
                            oss << *(last-1);
                        } else {
                            oss << *first;
                            last--;
                        }
                        if(last != vec.end()) oss << ",";
                        first = last;

                    }
                }
                else {
                    std::copy(vec.begin(), vec.end() - 1, std::ostream_iterator<T>(oss, ","));
                    oss << vec.back();
                }
            }
            return oss.str();
        }
    };


    // property_tree custom translator for vectors
    template <>
    class VectorTranslator<std::string> {

        typedef std::string internal_type;
        typedef typename std::vector<std::string> external_type;

    public:

        // Converts a string to a vector
        external_type get_value(const internal_type& str) {
            external_type result;
            if(!str.empty()) {
                //std::vector<std::string> tok;
                boost::split(result, str, boost::is_any_of(","));
                //for( auto &it: tok ) result.push_back(it);
                //std::transform(tok.begin(), tok.end(), back_inserter(result));
            }
            return result;
        }

        // Converts a vector to a comma-separated string
        internal_type put_value(const external_type& vec) {
            std::ostringstream oss;
            if(!vec.empty()) {
                std::copy(vec.begin(), vec.end() - 1, std::ostream_iterator<std::string>(oss, ","));
                oss << vec.back();
            }
            return oss.str();
        }
    };


    // specialization to handle a vector<FileType>
    template <>
    class VectorTranslator<redux::momfbd::FileType> {

        typedef std::string internal_type;
        typedef typename std::vector<redux::momfbd::FileType> external_type;

        FileTypeTranslator elementTranslator;
        void parseSegment(external_type& out, const internal_type& elem) {
            redux::momfbd::FileType tmp = elementTranslator.get_value(elem);
            if(tmp)
                out.push_back(elementTranslator.get_value(elem));
        }
    public:

        // Converts a string to a vector
        boost::optional<external_type> get_value(const internal_type& str) {
            external_type result;
            if(!str.empty()) {
                std::vector<std::string> tokens;
                boost::split(tokens, str, boost::is_any_of(","));
                for(auto & token : tokens) {
                    parseSegment(result, token);
                }
            }
            return boost::optional<external_type>(result);
        }

        // Converts a vector to a comma-separated string
        boost::optional<internal_type> put_value(const external_type& vec) {
            std::ostringstream oss;
            for(auto it = vec.begin(); it != vec.end() - 1; it++) {
                oss << elementTranslator.put_value(*it) << ",";
            }
            //         if( !vec.empty() ) {
            //             std::copy( vec.begin(), vec.end() - 1, std::ostream_iterator<T>( oss, "," ) );
            //             oss << vec.back();
            //         }
            return oss.str();
        }
    };

    // specialization to handle a vector<DiversityValue>
    template <>
    class VectorTranslator<redux::momfbd::DiversityValue> {

        typedef std::string internal_type;
        typedef typename std::vector<redux::momfbd::DiversityValue> external_type;

        DiversityValueTranslator elementTranslator;
        void parseSegment( external_type& out, const internal_type& in ) {
            using namespace std;
            using namespace redux::momfbd;
            string tmpIn = in;
            double tmpD = 1.0;
            bool physical(false);
            if( tmpIn.find( "mm" ) != string::npos ) {
                physical = true;
                tmpD = 1.00E-03;
            } else if( tmpIn.find( "cm" ) != string::npos ) {
                physical = true;
                tmpD = 1.00E-02;
            }
            // we extracted the suffix/units above, so now we can delete the letters and extract the numbers.
            boost::remove_erase_if( tmpIn, boost::is_any_of("cm \""));
            out.push_back( DiversityValue ( tmpD*boost::lexical_cast<double>( tmpIn ), physical ) );
        }
    public:

        // Converts a string to a vector
        boost::optional<external_type> get_value(const internal_type& str) {
            external_type result;
            if(!str.empty()) {
                std::vector<std::string> tokens;
                boost::split(tokens, str, boost::is_any_of(","));
                for(auto & token : tokens) {
                    parseSegment(result, token);
                }
            }
            return boost::optional<external_type>(result);
        }

        // Converts a vector to a comma-separated string
        boost::optional<internal_type> put_value( const external_type& divs ) {
            std::ostringstream oss;
            for( const auto& d: divs ) {
                if( ! oss.str().empty() ) {
                    oss << ",";
                }
                double c = d.coefficient;
                if( d.physical ) c *= 1.0E3;       // mm by default
                oss << c;
                if( d.physical ) {
                    oss << " mm";
                }
            }
            return oss.str();
        }
    };

    
    // Translator to handle a vector<ModeID>
    class ModeListTranslator {

        typedef std::string internal_type;
        typedef typename redux::momfbd::ModeList external_type;

        void parseSegment( external_type& out, internal_type in ) {
            using namespace std;
            using namespace redux::momfbd;
            external_type tmpOut;
            string tmpIn = in;
            size_t n1 = count( tmpIn.begin(), tmpIn.end(), '-' );
            size_t n2 = count( tmpIn.begin(), tmpIn.end(), ':' );
            size_t nz = count_if( tmpIn.begin(), tmpIn.end(), [](const char& a) { return (a == 'z' || a == 'Z'); });
            size_t nk = count_if( tmpIn.begin(), tmpIn.end(), [](const char& a) { return (a == 'k' || a == 'K'); });
            if( nz && nk ) {
                throw logic_error("Conflicting mode-types in specified mode segment: \"" + in + "\"");
            }
            ModeBase tp(MB_NONE);
            if( nz ) tp = ZERNIKE;
            if( nk ) tp = KARHUNEN_LOEVE;
            boost::remove_erase_if( tmpIn, boost::is_any_of("ZzKk"));
            if( (n1+n2) == 0 ) {
                tmpOut.push_back( ModeID( boost::lexical_cast<uint16_t>( tmpIn ), tp ) );
            } else if( (n1+n2) == 1 ) {
                size_t n = tmpIn.find_first_of(":-");
                uint16_t first = boost::lexical_cast<uint16_t>( tmpIn.substr(0, n));
                uint16_t last = boost::lexical_cast<uint16_t>( tmpIn.substr(n + 1));
                bool rev(false);
                if( first > last ) {
                    std::swap( first, last );
                    rev = true;
                }
                while( first <= last ) tmpOut.push_back( ModeID( boost::lexical_cast<uint16_t>(first++), tp ) );
                if( rev ) std::reverse( tmpOut.begin(), tmpOut.end() );
            } else if( n2 == 2 ) {
                size_t n = tmpIn.find_first_of(":");
                size_t nn = tmpIn.find_first_of(":", n + 1);
                uint16_t first = boost::lexical_cast<uint16_t>( tmpIn.substr(0, n));
                int step = boost::lexical_cast<uint16_t>( tmpIn.substr(n + 1, nn - n - 1));
                uint16_t last = boost::lexical_cast<uint16_t>( tmpIn.substr(nn + 1));
                bool rev(false);
                if( first > last ) {
                    std::swap( first,last );
                    if( step < 0 ) step = abs(step);
                    rev = true;
                }
                while( first <= last ) {
                    tmpOut.push_back(ModeID( boost::lexical_cast<uint16_t>(first), tp ));
                    first += step;
                }
                if( rev ) std::reverse( tmpOut.begin(), tmpOut.end() );
            }
            out.insert( out.end(), tmpOut.begin(), tmpOut.end() );
        }
    public:

        // Converts a string to a vector
        boost::optional<external_type> get_value( const internal_type& str ) {
            external_type modeList;
            if(!str.empty()) {
                std::vector<std::string> tokens;
                boost::split(tokens, str, boost::is_any_of(","));
                for( const auto& token : tokens) {
                    parseSegment( modeList, token );
                }
            }
            return boost::optional<external_type>( modeList );
        }

        // Converts a vector to a comma-separated string
        boost::optional<internal_type> put_value( const external_type& modeList ) {
            using namespace std;
            using namespace redux::momfbd;
            string ret;
            if( !modeList.empty() ) {
                ostringstream oss;
                auto first = modeList.begin();
                while( first != modeList.end() ) {
                    oss.str("");
                    oss.clear();
                    ModeBase type = first->type;
                    auto last = first + 1;
                    if( last == modeList.end() ) { // last element in vector
                        oss << first->mode;
                        ++first;
                    } else {
                        if( last->type != type ) {
                            oss << first->mode;
                        } else {
                            int step = static_cast<int>(last->mode) - static_cast<int>(first->mode);
                            uint16_t next = static_cast<uint16_t>(last->mode + step);
                            size_t count(2);
                            while( (++last != modeList.end()) && (last->mode == next) && (last->type == type) ) {
                                next += step;
                                count++;
                            }
                            if( count > 2 ) {
                                oss << first->mode;
                                if( step == 1 ) oss << "-";
                                else  oss << ":";
                                if( step != 1 ) oss << step << ":";
                                oss << (last-1)->mode;
                            } else {
                                oss << first->mode;
                                --last;
                            }
                        }
                        if( last != modeList.end() ) oss << ",";
                        first = last;
                    }
                    if( (modeList.defaultType == MB_NONE) || (type != modeList.defaultType) ) {
                        if( type == KARHUNEN_LOEVE ) ret += "K";
                        else ret += "Z";
                    }
                    ret += oss.str();
                }
            }
            if( ret.back() == ',' ) ret.pop_back();
            return ret;
        }
    };

}


namespace boost {

    namespace property_tree {

        template <> struct translator_between<std::string, bool> { typedef redux::BoolTranslator type; };
        template <> struct translator_between<std::string, redux::momfbd::FileType> { typedef redux::FileTypeTranslator type; };
        
        template <> struct translator_between<std::string, redux::momfbd::ModeBase> { typedef redux::ModeBaseTranslator type; };
        template <> struct translator_between<std::string, redux::momfbd::ModeList> { typedef redux::ModeListTranslator type; };
        
        template<typename T>
        struct translator_between<std::string, std::vector<T>> {
            typedef redux::VectorTranslator<T> type;
        };


    } // namespace property_tree
} // namespace boost

#endif  // REDUX_TRANSLATORS_HPP
