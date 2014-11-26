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

namespace redux {

    struct BoolTranslator {

        typedef std::string internal_type;

        // Converts a string to bool (i.e. when calling get<bool>(key) on a ptree)
        boost::optional<bool> get_value(const internal_type& str) const {
            if(!str.empty()) {
                using boost::algorithm::iequals;

                if(iequals(str, "true") || iequals(str, "yes") || str == "1") {
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


    // property_tree custom translator for vectors
    template <typename T>
    class VectorTranslator {

        typedef std::string internal_type;
        typedef typename std::vector<T> external_type;

        void parseSegment(external_type& out, const internal_type& elem) {
            size_t n = std::count(elem.begin(), elem.end(), '-');
            n += std::count(elem.begin(), elem.end(), ':');
            if(n == 0) {
                out.push_back(boost::lexical_cast<T>(elem));
                return;
            }
            else if(n == 2) {
                n = elem.find_first_of(":-");
                size_t nn = elem.find_first_of(":-", n + 1);
                T first = boost::lexical_cast<T>(elem.substr(0, n));
                T step = boost::lexical_cast<T>(elem.substr(n + 1, nn - n - 1));
                T last = boost::lexical_cast<T>(elem.substr(nn + 1));
                while(first <= last) {
                    out.push_back(first);
                    first += step;
                }
            }
            else if(n == 1) {
                n = elem.find_first_of(":-");
                T first = boost::lexical_cast<T>(elem.substr(0, n));
                T last = boost::lexical_cast<T>(elem.substr(n + 1));
                while(first <= last) out.push_back(first++);
            }
        }


    public:

        // Converts a string to a vector
        external_type get_value(const internal_type& str) {
            external_type result;
            if(!str.empty()) {
                std::vector<std::string> tok;
                boost::split(tok, str, boost::is_any_of(","));
                if(std::numeric_limits<T>::is_integer && !std::numeric_limits<T>::is_signed) {    // only expand '-' & ':' for unsigned integers
                    for(auto & it : tok) {
                        parseSegment(result, it);
                    }
                }
                else {
                    std::transform(tok.begin(), tok.end(), back_inserter(result), boost::lexical_cast<T, internal_type>);
                }
            }
            return external_type(result);
        }

        // Converts a vector to a comma-separated string
        internal_type put_value(const external_type& vec) {
            std::ostringstream oss;
            if(!vec.empty()) {
                if( std::numeric_limits<T>::is_integer && !std::numeric_limits<T>::is_signed ) {    // only expand '-' & ':' for unsigned integers
                    auto first = vec.begin();
                    while(first < vec.end()) {
                        
                        auto last = first + 1;
                        if(last == vec.end()) { // last element in vector
                            oss << *first;
                            break;
                        }

                        int step = *last - *first;
                        T next = *last + step;
                        size_t count = 2;
                        while(++last < vec.end() && *last == next) {
                            next += step;
                            count++;
                        }
                        if(count > 2) {
                            oss << *first << ":";
                            if( step != 1 ) oss << step << ":";
                            oss << *(last-1);
                        } else {
                            oss << *first;
                            last--;
                        }
                        if(last < vec.end()) oss << ",";
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
                std::vector<std::string> tok;
                boost::split(tok, str, boost::is_any_of(","));
                for(auto & it : tok) {
                    parseSegment(result, it);
                }
            }
            return boost::optional<external_type>(result);
        }

        // Converts a vector to a comma-separated string
        boost::optional<internal_type> put_value(const external_type& vec) {
            std::ostringstream oss;
            for(auto it = vec.begin(); it < vec.end() - 1; it++) {
                oss << elementTranslator.put_value(*it) << ",";
            }
            //         if( !vec.empty() ) {
            //             std::copy( vec.begin(), vec.end() - 1, std::ostream_iterator<T>( oss, "," ) );
            //             oss << vec.back();
            //         }
            return oss.str();
        }
    };

}


namespace boost {

    namespace property_tree {

        template <> struct translator_between<std::string, bool> { typedef redux::BoolTranslator type; };
        template <> struct translator_between<std::string, redux::momfbd::FileType> { typedef redux::FileTypeTranslator type; };

        template<typename T>
        struct translator_between<std::string, std::vector<T>> {
            typedef redux::VectorTranslator<T> type;
        };


    } // namespace property_tree
} // namespace boost

#endif  // REDUX_TRANSLATORS_HPP
