#ifndef REDUX_FILE_FILEMETA_HPP
#define REDUX_FILE_FILEMETA_HPP

#include "redux/util/stringutil.hpp"

#include <exception>
#include <functional>
#include <iostream>
#include <map>

#include <boost/property_tree/ptree.hpp>

namespace bpt = boost::property_tree;



namespace redux {

    namespace file {

        /*! @ingroup file FileIO
         *  @{
         */

        struct FileMeta {
            struct cmp : public std::binary_function<std::string, std::string, bool> {
                bool operator()(const std::string& lhs, const std::string& rhs) const {
                    return redux::util::nocaseLess(lhs, rhs);
                }
            };

            FileMeta ( void ) {}
            virtual ~FileMeta ( void ) = default;
            
            virtual std::string getText( void ) { return ""; }
            std::string& at(const std::string& key) { return list[key]; };
            //const std::string& at (const std::string& key) const { return list[key]; };
            void clearList(void) { list.clear(); }
            bool emplace(const std::string& key, const std::string& val) { auto ret = list.emplace(key,val); return ret.second; }
            size_t erase(const std::string& key) { return list.erase(key); }
            
            std::string& operator[](const std::string& key) { return list[key]; };
            
            std::map<std::string,std::string, cmp> list;
            
        };



        template <typename T, typename U>
        inline size_t readOrThrow ( T& strm, U* out, size_t nElements, const std::string& msg = "" ) {
            size_t nBytes = nElements*sizeof ( U );
            strm.read ( reinterpret_cast<char*> ( out ), nBytes );
            if ( !strm.good() ) {
                throw std::ios_base::failure ( "Read failed: " + msg );
            }

            return nBytes;
        }



        template <typename T, typename U>
        inline size_t writeOrThrow ( T& strm, const U* data, size_t nElements, const std::string& msg = "" ) {
            size_t nBytes = nElements*sizeof ( U );
            strm.write ( reinterpret_cast<const char*> ( data ), nBytes );
            if ( !strm.good() ) {
                throw std::ios_base::failure ( "Write failed: " + msg );
            }

            return nBytes;
        }




        /*! @} */

    }

}

#endif // REDUX_FILE_FILEMETA_HPP
