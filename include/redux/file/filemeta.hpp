#ifndef REDUX_FILE_FILEMETA_HPP
#define REDUX_FILE_FILEMETA_HPP

#include "redux/util/stringutil.hpp"

#include <exception>
#include <functional>
#include <iostream>
#include <map>

#include <boost/property_tree/ptree.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace bpt = boost::property_tree;
namespace bpx = boost::posix_time;



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
            virtual size_t getNumberOfFrames(void) { return 1; }
            virtual bpx::ptime getStartTime(void) { return bpx::ptime(); };
            virtual bpx::ptime getEndTime(void) { return bpx::ptime(); };
            virtual bpx::ptime getAverageTime(void) { return bpx::ptime(); };
            virtual bpx::time_duration getExposureTime(void) { return bpx::time_duration(); };
            
            virtual size_t dataSize(void) { return 0; };
            virtual size_t dimSize(size_t) { return 0; };
            virtual uint8_t elementSize(void) { return 0; };
            virtual uint8_t nDims(void) { return 0; }
            virtual size_t nElements(void) { return 0; };
            virtual int getIDLType(void) { return 0; };
            
            virtual double getMinMaxMean( const char* data, double* Min=nullptr, double* Max=nullptr ) = 0;
            
            /*std::string& at(const std::string& key) { return list[key]; };
            //const std::string& at (const std::string& key) const { return list[key]; };
            void clearList(void) { list.clear(); }
            bool emplace(const std::string& key, const std::string& val) { auto ret = list.emplace(key,val); return ret.second; }
            size_t erase(const std::string& key) { return list.erase(key); }
            
            std::string& operator[](const std::string& key) { return list[key]; };
            
            std::map<std::string, std::string, cmp> list;*/
            
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
