#ifndef REDUX_FILE_FILEMOMFBD_HPP
#define REDUX_FILE_FILEMOMFBD_HPP

#include "redux/file/fileinfo.hpp"
#include "redux/image/image.hpp"
#include "redux/util/array.hpp"
#include "redux/util/arrayutil.hpp"

#include <fstream>

#define MOMFBD_IMG        1
#define MOMFBD_PSF        2
#define MOMFBD_OBJ        4
#define MOMFBD_RES        8
#define MOMFBD_ALPHA     16
#define MOMFBD_DIV       32
#define MOMFBD_PATCH     63       // convenience for flagging any "per-patch" data.
#define MOMFBD_MODES     64
#define MOMFBD_NAMES    128
#define MOMFBD_ALL      255


namespace redux {

    namespace file {

        struct FileMomfbd : public redux::file::FileInfo {

            typedef std::shared_ptr<FileMomfbd> Ptr;
                       
            struct PatchInfo {
                
                PatchInfo() : /*region({0,0,0,0}),*/ nChannels(0), offx(0), offy(0), nPixelsX(0), nPixelsY(0), npsf(0), nobj(0), nres(0), nalpha(0), nm(0),
                     ndiv(0), nphx(0), nphy(0), offset(0), imgPos(0), psfPos(0), objPos(0), resPos(0), alphaPos(0), diversityPos(0) {};

                int32_t region[4];         //!< FirstX, LastX, FirstY, LastY
                int32_t nChannels, offx, offy;

                std::shared_ptr<int32_t> nim, dx, dy;

                int32_t nPixelsX, nPixelsY;
                int32_t npsf, nobj, nres, nalpha, nm;
                int32_t ndiv, nphx, nphy;

                int64_t offset, imgPos, psfPos, objPos, resPos, alphaPos, diversityPos;

                uint8_t parse( std::ifstream& file, const bool& swapNeeded, const float& version );
                void load ( std::ifstream& file, char*& ptr, const bool& swapNeeded, const float& version, uint8_t loadMask=MOMFBD_ALL, int verbosity=0, uint8_t alignTo=4 ) const;
                void write( std::ofstream& file, const char* data, uint8_t writeMask=MOMFBD_ALL );

            };

            FileMomfbd( void );
            FileMomfbd( const std::string& );
            
            void clear(void);

            void read( std::ifstream& );
            void read( const std::string& );

            void write( std::ofstream&, const char* data, uint8_t writeMask=MOMFBD_ALL, int verbosity=0 );
            void write( const std::string&, const char* data, uint8_t writeMask=MOMFBD_ALL, int verbosity=0 );
            
            char* load ( std::ifstream& file, char* data, uint8_t loadMask, int verbosity=0, uint8_t alignTo=4 );
            
            std::string getText( void ) {
                return "";
            }

            float version, pix2cf, cf2pix;
            std::string dateString, timeString, versionString;
            std::vector<std::string> fileNames;

            std::shared_ptr<int16_t> clipStartX, clipEndX, clipStartY, clipEndY;

            int32_t startX, endX, startY, endY;
            int32_t nChannels, nFileNames;
            int32_t nPH, nModes;
            int32_t nPatchesX, nPatchesY, nPoints;

            int64_t phOffset, modesOffset, filenameOffset;

            size_t  patchDataSize;
            size_t  headerSize;
            
            uint8_t dataMask;
            bool swapNeeded;                    //!< File & current system have different endianess
            
            redux::util::Array<PatchInfo> patches;

        };

        std::shared_ptr<FileMomfbd> readMomfbdInfo( const std::string& );

    }
}

#endif // REDUX_FILE_FILEMOMFBD_HPP
