#ifndef REDUX_FILE_FILEANA_HPP
#define REDUX_FILE_FILEANA_HPP

#include "redux/file/fileio.hpp"
#include "redux/file/filemeta.hpp"
#include "redux/image/image.hpp"
#include "redux/util/array.hpp"

#include <fstream>


namespace redux {

    namespace file {

        /*! @ingroup file
         *  @{
         */

        /*! Container for reading/writing ANA files.
        *
        */
        struct Ana : public redux::file::FileMeta {

            enum Magic { MAGIC_ANA = 0x5555aaaa, MAGIC_ANAR = 0xaaaa5555 };
            enum TypeIndex { ANA_BYTE = 0,      // uint8
                             ANA_WORD,          // int16
                             ANA_LONG,          // int32
                             ANA_FLOAT,
                             ANA_DOUBLE,
                             ANA_LONGLONG,      // int64_t unsupported
                             ANA_COMPLEX=8,     // experimental support (=8 to match IDL's DCOMPLEX)
                             ANA_UNDEF=255 };
            static const uint8_t typeSizes[];   // = { 1, 2, 4, 4, 8, 8, 0, 0, 16 };

            typedef std::shared_ptr<Ana> Ptr;
            
            Ana( void );
            explicit Ana( const std::string& );
            ~Ana();

            void read( std::ifstream& );
            void read( const std::string& );

            void write( std::ofstream& );

            std::vector<std::string> getText( bool raw=false ) override {
                std::string hdr = m_Header.txt + m_ExtendedHeader;
                if( raw ) return std::vector<std::string>( 1, hdr );
                return makeFitsHeader(hdr);
            }
            
            std::vector<std::string> makeFitsHeader( const std::string& );
            
            size_t getNumberOfFrames(void) override;
            bpx::ptime getStartTime(void) override;
            bpx::ptime getEndTime(void) override;
            bpx::ptime getAverageTime(void) override;
            bpx::time_duration getExposureTime(void) override;
            std::vector<bpx::ptime> getStartTimes(void) override;
            
            size_t dataSize(void) override;
            size_t dimSize(size_t) override;
            uint8_t elementSize(void) override;
            uint8_t nDims(void) override { return m_Header.ndim; }
            size_t nElements(void) override;
            int getIDLType(void) override;
            int getFormat(void) override { return FMT_ANA; };
            
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

            static void readCompressed( std::ifstream& file, char* data, size_t nElements, const Ana* hdr );
            static void readUncompressed( std::ifstream& file, char* data, size_t nElements, const Ana* hdr );
            static int compressData( std::shared_ptr<uint8_t>& out, const char* data, int nElements, const std::shared_ptr<Ana>& hdr, int slice );



            /*! @name Read
             *  @brief Load an ANA file into a data block
             */
            //@{
            static void read( const std::string& filename, char* data, std::shared_ptr<redux::file::Ana>& hdr );
            template <typename T>
            static void read( const std::string& filename, redux::util::Array<T>& data, std::shared_ptr<redux::file::Ana>& hdr );
            template <typename T>
            static void read( const std::string& filename, redux::image::Image<T>& data, bool metaOnly=false );
            //@}
            
            /*! @name Write
             *  @brief Write data block into an ANA file.
             */
            //@{
            static void write( const std::string& filename, const char* data, std::shared_ptr<redux::file::Ana> hdr, bool compress = false, int slice=5 );
            template <typename T>
            static void write( const std::string& filename, const redux::util::Array<T>& data, std::shared_ptr<redux::file::Ana> hdr=0, int sliceSize=0 );
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


#endif // REDUX_FILE_FILEANA_HPP
