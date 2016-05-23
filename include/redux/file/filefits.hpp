#ifndef REDUX_FILE_FILEFITS_HPP
#define REDUX_FILE_FILEFITS_HPP

#include "redux/file/filemeta.hpp"
#include "redux/image/image.hpp"
#include "redux/util/array.hpp"
#include "redux/util/arrayutil.hpp"


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
            enum TypeIndex { FITS_BYTE = 0,      // uint8
                             FITS_WORD,          // int16
                             FITS_LONG,          // int32
                             FITS_FLOAT,
                             FITS_DOUBLE,
                             FITS_LONGLONG,      // int64_t unsupported
                             FITS_COMPLEX=8,     // experimental support (=8 to match IDL's DCOMPLEX)
                             FITS_UNDEF=255 };
            static const uint8_t typeSizes[];   // = { 1, 2, 4, 4, 8, 8, 0, 0, 16 };

            typedef std::shared_ptr<Fits> Ptr;
            
            Fits( void );
            Fits( const std::string& );

            void read( std::ifstream& );
            void read( const std::string& );

            void write( std::ofstream& );

            std::string getText( void ) {
                return m_Header.txt + m_ExtendedHeader;
            }
            
            size_t getNumberOfFrames(void);
            bpx::ptime getStartTime(void);
            bpx::ptime getEndTime(void);
            bpx::ptime getAverageTime(void);
            bpx::time_duration getExposureTime(void);
            
            size_t dataSize(void);
            size_t dimSize(size_t);
            uint8_t elementSize(void);
            uint8_t nDims(void) { return m_Header.ndim; }
            size_t nElements(void);
            int getIDLType(void);
            
            double getMinMaxMean( const char* data, double* Min=nullptr, double* Max=nullptr );

            struct raw_header {                    // first block for ana files
                uint32_t synch_pattern;
                uint8_t subf;
                uint8_t source;
                uint8_t nhb;
                uint8_t datyp;
                uint8_t ndim;
                uint8_t free1;
                uint8_t cbytes[4];
                uint8_t free[178];
                uint32_t dim[16];
                char txt[256];
            } m_Header;

            struct compressed_header {
                uint32_t tsize, nblocks, bsize;
                uint8_t slice_size, type;
            } m_CompressedHeader;

            std::string m_ExtendedHeader;
            size_t hdrSize;

            static void readCompressed( std::ifstream& file, char* data, size_t nElements, const Fits* hdr );
            static void readUncompressed( std::ifstream& file, char* data, size_t nElements, const Fits* hdr );
            static int compressData( std::shared_ptr<uint8_t>& out, const char* data, int nElements, const std::shared_ptr<Fits>& hdr, int slice );



            /*! @name Read
             *  @brief Load an FITS file into a data block
             */
            //@{
            static void read( const std::string& filename, char* data, std::shared_ptr<redux::file::Fits>& hdr );
            template <typename T>
            static void read( const std::string& filename, redux::util::Array<T>& data, std::shared_ptr<redux::file::Fits>& hdr );
            template <typename T>
            static void read( const std::string& filename, redux::image::Image<T>& data );
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


#endif // REDUX_FILE_FILEFITS_HPP
