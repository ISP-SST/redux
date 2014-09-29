#ifndef REDUX_FILE_FILEMOMFBD_HPP
#define REDUX_FILE_FILEMOMFBD_HPP

#include "redux/file/fileinfo.hpp"
#include "redux/image/image.hpp"
#include "redux/util/array.hpp"
#include "redux/util/arrayutil.hpp"

#include <fstream>


namespace redux {

    namespace file {

        struct FileMomfbd : public redux::file::FileInfo {

            typedef std::shared_ptr<FileMomfbd> Ptr;
                       
            struct PatchInfo {

                int32_t nChannels;
                int32_t region[4];         //!< FirstX, LastX, FirstY, LastY

                std::shared_ptr<int32_t> nim, dx, dy;

                int32_t nPixelsX, nPixelsY;
                int32_t npsf, nobj, nres, nalpha, nm;
                int32_t ndiv, nphx, nphy;

                std::streampos offset, imgPos, psfPos, objPos, resPos, alphaPos, diversityPos;


            };

            FileMomfbd( void );
            FileMomfbd( const std::string& );

            void read( std::ifstream& );
            void read( const std::string& );

            void write( std::ofstream& );
            
            std::string getText( void ) {
                return "";
            }

            float version;
            std::string dateString, timeString, versionString;

            std::shared_ptr<int16_t> clipStartX, clipEndX, clipStartY, clipEndY;

            int32_t startX, endX, startY, endY;
            int32_t nChannels, nFileNames;
            int32_t nPH, nModes;
            int32_t nPatchesX, nPatchesY, nPoints;

            std::streampos phOffset;
            std::streampos modesOffset;
            std::streampos filenameOffset;

            size_t  patchDataSize;
            size_t  headerSize;
            
            bool swapNeeded;                    //!< File & current system have different endianess
            
            redux::util::Array<PatchInfo> patches;

        };

        std::shared_ptr<FileMomfbd> readMomfbdInfo( const std::string& );

    }
}

#endif // REDUX_FILE_FILEMOMFBD_HPP
