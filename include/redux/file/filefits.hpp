#ifndef REDUX_FILE_FILEFITS_HPP
#define REDUX_FILE_FILEFITS_HPP

#ifdef RDX_WITH_FITS

#include "redux/file/fileio.hpp"
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
            explicit Fits( const std::string& );
            ~Fits();

            void close( void );
            void read( const std::string& );

            void write( std::ofstream& );

            std::vector<std::string> getText( bool ) override;
            
            template <typename T>
            static std::string makeValue( const T& t, bool quote=true );
            static std::string makeValue( const char* s, bool quote=true );
            template <typename T>
            static T getValue( const std::string& v );
            
            
            template <typename T>
            static std::string makeCard( const std::string& key, const T& value=T(), const std::string& comment="", bool allowLong=false );
            static std::string makeCard( const std::string& key, const char* value, const std::string& comment="", bool allowLong=false );
            static std::string makeCard( const std::string& key, const std::string& comment="", bool allowLong=false );   // version for cards without value, i.e. "END", or comments.
            static std::string getCard( const std::vector<std::string>& hdr, std::string key );
            static void addCard( std::vector<std::string>& hdr, const std::string& card );
            static void removeCards( std::vector<std::string>& hdr, std::string key );
            static void insertCard( std::vector<std::string>& hdr, const std::string& card, size_t location=std::string::npos );
            static void insertCardAfter( std::vector<std::string>& hdr, const std::string& card, std::string after );
            static void insertCardBefore( std::vector<std::string>& hdr, const std::string& card, std::string before );
            static bool updateCard( std::vector<std::string>& hdr, size_t location, const std::string& card );
            static bool updateCard( std::vector<std::string>& hdr, std::string key, const std::string& card );
            static bool updateCard( std::vector<std::string>& hdr, const std::string& card );
            static bool emplaceCard( std::vector<std::string>& hdr, const std::string& key, const std::string& card );
            static bool emplaceCard( std::vector<std::string>& hdr, const std::string& card );
            
            static void splitCard( const std::string& card, std::string& key, std::string& value, std::string& comment );
            
            static std::vector<std::string> splitLongString( const std::string& card );
            static std::string mergeLongString( std::vector<std::string> cards );
            static void splitLongStrings( std::vector<std::string>& cards );
            static void mergeLongStrings( std::vector<std::string>& cards );

            
            template <typename T>
            static T getValue( const std::vector<std::string>& hdr, std::string key );
            template <typename T>
            std::vector<T> getTableArray( std::string key );
            
            size_t getNumberOfFrames(void) override;
            bpx::ptime getStartTime(void) override;
            bpx::ptime getEndTime(void) override;
            bpx::ptime getAverageTime(void) override;
            bpx::time_duration getExposureTime(void) override;
            std::vector<bpx::ptime> getStartTimes(void) override;
            std::vector<size_t> getFrameNumbers(void) override;
           
            size_t dataSize(void) override;
            size_t dimSize(size_t) override;
            uint8_t elementSize(void) override;
            uint8_t nDims(void) override { return primaryHDU.nDims; }
            size_t nElements(void) override;
            int getIDLType(void) override;
            
            double getMinMaxMean( const char* data, double* Min=nullptr, double* Max=nullptr );
            int getFormat(void) override { return FMT_FITS; };

            struct hdu {
                hdu( void ) : bitpix(0), nDims(0), dataType(0), elementSize(0), nElements(0) {}
                virtual ~hdu() {}
                int bitpix;
                int nDims;
                int dataType;               // data type as defined in cfitsio
                size_t elementSize;         // element size (in bytes) = abs(bitpix/8)
                size_t nElements;
                std::vector<int> dims;
                std::vector<std::string> cards;
                virtual void dummy(void)=0;
            };
            
            struct image_hdu : public hdu {
                image_hdu() : dHDU(0) {}
                void dummy(void) override {};
                int dHDU;       // index to hdu containing data, e.g. compressed tile image.
            };
            
            struct ascii_hdu : public hdu {
                void dummy(void) override {};
                struct table_info_t {
                    int columnStart;            // = TBCOL, offset where this column starts
                    std::string columnName;     // = TTYPEn, name of this data-column
                    std::string columnFormat;   // = TFORM, Fortran ISO 2004 format string
                    std::string columnUnit;     // = TUNIT, physical unit of the data
                };
                uint16_t nColumns;              // = TFIELDS, number of columns in this table
                std::string name;               // = EXTNAME
                std::vector<table_info_t> table_info;
                redux::util::Array<char> data;
            };
            
            struct binary_hdu : public hdu {
                void dummy(void) override {};
                std::vector<std::string> data;
            };
            
            struct image_hdu primaryHDU;
            std::vector<std::shared_ptr<struct hdu>> extHDUs;
            
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

#endif  // RDX_WITH_FITS

#endif // REDUX_FILE_FILEFITS_HPP
