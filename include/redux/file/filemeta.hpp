#ifndef REDUX_FILE_FILEMETA_HPP
#define REDUX_FILE_FILEMETA_HPP

#include <exception>
#include <functional>
#include <iostream>
#include <map>

#include <boost/date_time/posix_time/posix_time.hpp>

namespace bpx = boost::posix_time;



namespace redux {

    namespace file {

        /*! @ingroup file FileIO
         *  @{
         */

        struct FileMeta {

            FileMeta ( void ) {}
            virtual ~FileMeta ( void ) = default;
            
            virtual std::vector<std::string> getText( bool raw=false ) { return std::vector<std::string>(1,""); }
            virtual size_t getNumberOfFrames(void) { return 1; }
            virtual bpx::ptime getStartTime(void) { return bpx::ptime(); };
            virtual bpx::ptime getEndTime(void) { return bpx::ptime(); };
            virtual bpx::ptime getAverageTime(void) { return bpx::ptime(); };
            virtual bpx::time_duration getExposureTime(void) { return bpx::time_duration(); };
            
            virtual std::vector<bpx::ptime> getStartTimes(void) { return std::vector<bpx::ptime>(); };
            
            virtual size_t dataSize(void) { return 0; };
            virtual size_t dimSize(size_t) { return 0; };
            virtual uint8_t elementSize(void) { return 0; };
            virtual uint8_t nDims(void) { return 0; }
            virtual size_t nElements(void) { return 0; };
            virtual int getIDLType(void) { return 0; };
            virtual int getFormat(void) { return 0; };
            
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
