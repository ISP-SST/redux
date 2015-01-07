#ifndef REDUX_MOMFBD_CONFIG_HPP
#define REDUX_MOMFBD_CONFIG_HPP

#include <boost/asio.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/thread.hpp>

namespace po = boost::program_options;
namespace bpt = boost::property_tree;


namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */

        enum ModeBase { ZERNIKE=1, KARHUNEN_LOEVE };
        
        enum FileType { FT_NONE=0, FT_ANA, FT_FITS, FT_MOMFBD=4, FT_MASK=7 };
        extern const std::map<FileType, std::string> FileTypeNames;
        extern const std::map<FileType, std::string> FileTypeExtensions;
        
        enum DataType  { DT_I08T=0, DT_I16T, DT_I32T, DT_I64T, DT_F32T, DT_F64T };
        enum SaveFlags { SF_SAVE_ALPHA=1, SF_SAVE_COBJ, SF_SAVE_DIVERSITY=4, SF_SAVE_METRIC=8,
                         SF_SAVE_MODES=16, SF_SAVE_PSF=32, SF_SAVE_PSF_AVG=64, SF_SAVE_RESIDUAL=128,
                         SF_SAVE_NAMES=256, SF_SAVE_FFDATA=512 };
        enum RunFlags  { RF_CALIBRATE=1, RF_DONT_MATCH_IMAGE_NUMS, RF_FAST_QR=4, RF_FIT_PLANE=8,
                         RF_FLATFIELD=16, RF_GLOBAL_NOISE=32, RF_NEW_CONSTRAINTS=64, RF_NO_CLIP=128,
                         RF_NO_CONSTRAINTS=256, RF_NO_FILTER=512, RF_FORCE_WRITE=1024, RF_SWAP=2048 };

        enum FillpixMethod { FPM_MEDIAN=1, FPM_INVDISTWEIGHT, FPM_HORINT };
        extern const std::map<std::string, int> fillpixMap;
        enum GradientMethod { GM_DIFF=1, GM_VOGEL };
        extern const std::map<std::string, int> gradientMap;
        enum GetstepMethod { GSM_SDSC=1, GSM_CNJG, GSM_BFGS, GSM_BFGS_inv };
        extern const std::map<std::string, int> getstepMap;
        


        
        /*!
         * Settings with channel scope
         */
        struct ChannelCfg {

            ChannelCfg();
            virtual ~ChannelCfg();

            void parseProperties( bpt::ptree&, const ChannelCfg& defaults=ChannelCfg() );
            void getProperties( bpt::ptree&, const ChannelCfg& defaults=ChannelCfg() ) const;

            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            
            bool check(void);
            
            bool operator==(const ChannelCfg&) const;
            
            
            /********* Hardware **********/
            float arcSecsPerPixel;                  //!< Image scale   (default: 0, has to be specified)
            float pixelSize;                        //!< Physical size of pixels (default: 10\mu)
            float rotationAngle;                    //!< Tilt of this camera relative to "anchor channel" (default: 0)
            /*****************************/
            
            /**** Numerical settings *****/
            float noiseFudge;                       //!< Noise weight (default: 1)
            float weight;                           //!< Weight for this channel in calculations (default: 1)
            /*****************************/

            /********  Diversity  ********/
            // TODO: reorganize
            std::vector<double> diversity;
            std::vector<uint32_t> diversityOrders;
            std::vector<uint32_t> diversityTypes;
            /*****************************/
            
            /******* Data settings *******/
            std::vector<int16_t> alignClipl;         //!< Crop images to this region {firstX,lastX,firstY,lastY}, (default: none, has to be specified)
            uint16_t borderClip;                    //!< Disregard this many pixels at the edge when calculating statistics (default: 10)
            uint16_t maxLocalShift;                 //!< How much are the patches allowed to be shifted (default: 5 pixels)
            uint8_t incompletel;                    //!< Some files might not exist, just skip those.
            /*****************************/

            /************ Input **********/
            std::string imageDataDirl;              //!< Where the data is located
            std::string imageTemplate;
            std::string darkTemplate;
            std::string gainFile;
            std::string responseFile;
            std::string backgainFile;
            std::string psfFile;
            std::string mmFile;                     //!< Modulation matrix
            uint8_t mmRowl;                          //!< Number of rows in modulation matrix
            uint8_t mmWidthl;                        //!< Number of cols in modulation matrix
            std::string xOffsetFile,yOffsetFile;    //!< Alignment offsets (from pinhole-calibration)
            uint32_t imageNumberOffset;             //!< Add this offset to each image number in this channel
            std::vector<uint32_t> imageNumbers;     //!< Use these numbers together with the template to generate file-list
            std::vector<uint32_t> wfIndex;          //!< Identify wavefront, used to group/constrain simultaneous images if image-numbers can't be used.
            std::vector<uint32_t> darkNumbersl;     //!< Use these numbers together with the template to generate file-list
            std::vector<float> stokesWeightsl;
            /*****************************/
           
            
            
        };

        
        /*!
         * Settings with object scope
         */
        struct ObjectCfg : ChannelCfg {     // inherit channel so settings can also be specified in a wider scope (e.g. imageNumbers at any level)

            ObjectCfg();
            virtual ~ObjectCfg();

            void parseProperties( bpt::ptree&, const ObjectCfg& defaults=ObjectCfg() );
            void getProperties( bpt::ptree&, const ObjectCfg& defaults=ObjectCfg() ) const;

            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            
            bool check(void);
            
            const ObjectCfg& operator=(const ChannelCfg&);
            bool operator==(const ObjectCfg&) const;
            
            
            /******* Data settings *******/
            uint16_t saveMask;              //!< (default: 0)
            uint16_t nPatchesX, nPatchesY;
            uint16_t patchSize;             //!< (default: 128)
            uint16_t pupilSize;             //!< (default: 64)
            std::vector<uint16_t> subImagePosX, subImagePosY;
            std::string outputFileName;
            std::string pupilFile;
            float wavelength;               //!< (default: 0, has to be specified)
            /*****************************/
            
        };

        
        /*!
         * Settings with global scope
         */
        struct GlobalCfg : ObjectCfg {     // inherit object (and channel)

            GlobalCfg();
            ~GlobalCfg();

            void parseProperties( bpt::ptree& );
            void getProperties( bpt::ptree& ) const;

            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            
            bool check(void);
            
            const GlobalCfg& operator=(const ObjectCfg&);
            const GlobalCfg& operator=(const ChannelCfg&);
            bool operator==(const GlobalCfg&) const;
            
            uint16_t runFlags;

            /*********** Modes ***********/
            uint8_t modeBasis;          //!< Which basis to use for the fitting
            uint16_t klMinMode;         //!< First Zernike-mode to be considered in Karhunen-Loève expansion
            uint16_t klMaxMode;         //!< Last Zernike-mode to be considered in Karhunen-Loève expansion
            float klCutoff;             //!< If the expansion-coefficient is smaller than this, it will be ignored.
            uint16_t nInitialModes;     //!< How many modes to use in the first iteration
            uint16_t nModeIncrement;    //!< How many modes to add in each iteration
            std::vector<uint16_t> modeNumbers;  //!< Which modes to use
            /*****************************/
            
            /********* Hardware **********/
            float telescopeD;           //!< Telescope diameter
            float telescopeF;           //!< Telescope focal length  (or Ratio?)
            /*****************************/

            /**** Numerical settings *****/
            uint16_t minIterations;
            uint16_t maxIterations;
            uint8_t fillpixMethod;
            uint8_t gradientMethod;
            uint8_t getstepMethod;
            float badPixelThreshold;
            float FTOL;
            float EPS;
            float reg_gamma;       // not used atm.? def: 1E-4
            /*****************************/

            /******* Data settings *******/
            uint8_t outputFileType;
            uint8_t outputDataType;
            uint32_t sequenceNumber;
            std::string observationTime;
            std::string observationDate;
            std::string tmpDataDir;
            std::vector<std::string> outputFiles;
            /*****************************/

            
        };

        /*! @} */
        
    }

}

#endif  // REDUX_MOMFBD_CONFIG_HPP
