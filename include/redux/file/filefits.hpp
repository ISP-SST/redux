#ifndef REDUX_FILE_FILEFITS_HPP
#define REDUX_FILE_FILEFITS_HPP

#ifdef REDUX_WITH_FITS

#include "redux/file/filemeta.hpp"
#include "redux/image/image.hpp"
#include "redux/util/array.hpp"
#include "redux/util/arrayutil.hpp"

#include <fitsio.h>

namespace redux {

    namespace file {

        /*! @ingroup file
         *  @{
         */

        /*! Container for reading/writing FITS files.
        *
        */
        struct Fits : public redux::file::FileMeta {

            enum Magic { MAGIC_FITS = 0x504d4953 }; // = "PMIS"
            enum TypeIndex { FITS_NOTYPE = 0,
                             FITS_BYTE,          // uint8
                             FITS_WORD,          // int16
                             FITS_INT,           // int32
                             FITS_FLOAT,
                             FITS_DOUBLE,
                             FITS_COMPLEX,
                             FITS_STRING,
                             FITS_DCOMPLEX=9,
                             FITS_UWORD=12,
                             FITS_UINT,
                             FITS_LONG,
                             FITS_ULONG };
            static const uint8_t typeSizes[];   // = { 0, 1, 2, 4, 4, 8, 8, 0, 0, 16 };
            
            Fits( void );
            Fits( const std::string& );
            ~Fits();

            void close( void );
            void read( const std::string& );

            void write( std::ofstream& );

            std::string getText( void ) {
                std::string ret;
                for( auto& k: primaryHDU.keywords ) ret += (k + "\n");
                return ret;
            }
            
            size_t getNumberOfFrames(void);
            bpx::ptime getStartTime(void);
            bpx::ptime getEndTime(void);
            bpx::ptime getAverageTime(void);
            bpx::time_duration getExposureTime(void);
            
            size_t dataSize(void);
            size_t dimSize(size_t);
            uint8_t elementSize(void);
            uint8_t nDims(void) { return primaryHDU.nDims; }
            size_t nElements(void);
            int getIDLType(void);
            
            double getMinMaxMean( const char* data, double* Min=nullptr, double* Max=nullptr );

            struct hdu {
                int bitpix;
                int nDims;
                int dataType;
                size_t elementSize;
                size_t nElements;
                std::vector<int> dims;
                std::vector<std::string> keywords;
            } primaryHDU;
            
            fitsfile* fitsPtr_;
            int status_;

            /*! @name Read
             *  @brief Load a FITS file into a data block
             */
            //@{
            static void read( std::shared_ptr<redux::file::Fits>& hdr, char* data );
            template <typename T>
            static void read( const std::string& filename, redux::util::Array<T>& data, std::shared_ptr<redux::file::Fits>& hdr );
            template <typename T>
            static void read( const std::string& filename, redux::image::Image<T>& data, bool metaOnly=false );
            //@}
            
            /*! @name Write
             *  @brief Write data block into an FITS file.
             */
            //@{
            static void write( const std::string& filename, const char* data, std::shared_ptr<redux::file::Fits> hdr, bool compress = false, int slice=5 );
            template <typename T>
            static void write( const std::string& filename, const redux::util::Array<T>& data, std::shared_ptr<redux::file::Fits> hdr=0, int sliceSize=0 );
            template <typename T>
            static void write( const std::string& filename, const redux::image::Image<T>& image, int sliceSize=0 );
            template <typename T>
            static void write( const std::string& filename, const T* data, size_t n=1 );
            template <typename T>
            static void write( const std::string& filename, const std::vector<T>& v ) { write(filename,v.data(),v.size()); }
           //@}
            

        };

        /*! @} */

    } // end namespace file

} // end namespace redux

#endif  // REDUX_WITH_FITS

#endif // REDUX_FILE_FILEFITS_HPP
