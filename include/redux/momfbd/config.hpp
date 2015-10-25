#ifndef REDUX_MOMFBD_CONFIG_HPP
#define REDUX_MOMFBD_CONFIG_HPP

#include <map>

#include <boost/property_tree/ptree.hpp>
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
            
            operator std::string() const;           //!< cast to string (for easy printing of configuration)

            virtual void parseProperties( bpt::ptree&, const ChannelCfg& defaults=ChannelCfg() );
            virtual void getProperties( bpt::ptree&, const ChannelCfg& defaults=ChannelCfg() ) const;

            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            
            bool operator==(const ChannelCfg&) const;
            
            
            /********* Hardware **********/
            double rotationAngle;                    //!< Tilt of this camera relative to "anchor channel" (default: 0)
            /*****************************/
            
            /**** Numerical settings *****/
            double noiseFudge;                       //!< Noise weight (default: 1)
            double weight;                           //!< Weight for this channel in calculations (default: 1)
            /*****************************/

            /********  Diversity  ********/
            // TODO: reorganize
            std::vector<double> diversity;
            std::vector<uint32_t> diversityModes;  //!< mode numbers
            std::vector<uint32_t> diversityTypes;   //!< mode types (Zernike/KL)
            /*****************************/
            
            /******* Data settings *******/
            std::vector<int16_t> alignClip;         //!< Crop images to this region {firstX,lastX,firstY,lastY}, (default: none, has to be specified)
            uint16_t borderClip;                    //!< Disregard this many pixels at the edge when calculating statistics (default: 10)
            uint8_t incomplete;                     //!< Some files might not exist, just skip those.
            std::vector<uint16_t> subImagePosX, subImagePosY;
            /*****************************/

            /************ Input **********/
            std::string imageDataDir;               //!< Where the data is located
            std::string imageTemplate;
            std::string darkTemplate;
            std::string gainFile;
            std::string responseFile;
            std::string backgainFile;
            std::string psfFile;
            std::string mmFile;                     //!< Modulation matrix
            uint8_t mmRow;                          //!< Number of rows in modulation matrix
            uint8_t mmWidth;                        //!< Number of cols in modulation matrix
            std::string xOffsetFile,yOffsetFile;    //!< Alignment offsets (from pinhole-calibration)
            uint32_t imageNumberOffset;             //!< Add this offset to each image number in this channel
            std::vector<uint32_t> imageNumbers;     //!< Use these numbers together with the template to generate file-list
            std::vector<uint32_t> wfIndex;          //!< Identify wavefront, used to group/constrain simultaneous images if image-numbers can't be used.
            std::vector<uint32_t> darkNumbers;      //!< Use these numbers together with the template to generate file-list
            std::vector<float> stokesWeights;
            /*****************************/
           
            
            
        };
        

        inline std::ostream& operator<<( std::ostream &strm, const ChannelCfg &obj ) {
            strm << static_cast<std::string>(obj);
            return strm;
        }

        
        /*!
         * Settings with object scope
         */
        struct ObjectCfg : ChannelCfg {     // inherit channel so settings can also be specified in a wider scope (e.g. imageNumbers at any level)

            ObjectCfg();
            virtual ~ObjectCfg();
            
            virtual void parseProperties( bpt::ptree&, const ChannelCfg& defaults=ObjectCfg() );
            virtual void getProperties( bpt::ptree&, const ChannelCfg& defaults=ObjectCfg() ) const;

            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            
            const ObjectCfg& operator=(const ChannelCfg&);
            bool operator==(const ObjectCfg&) const;
            
            /********* Hardware **********/
            double telescopeF;                       //!< Telescope focal length  (or Ratio?)
            double arcSecsPerPixel;                  //!< Image scale   (default: 0, has to be specified)
            double pixelSize;                        //!< Physical size of pixels (default: 10\mu)
            double alphaToPixels,pixelsToAlpha;      //!< Conversion factors for the tilt-modes. Derived from arcSecsPerPixel & telescopeD.
            double alphaToDefocus,defocusToAlpha;    //!< Conversion factors for the focus-mode. Derived from telescopeD. (defocus in meters)
            /*****************************/
            
            /******* Data settings *******/
            uint16_t maxLocalShift;                 //!< How much are the patches allowed to be shifted (default: 5 pixels)
            uint16_t minimumOverlap;                //!< Desired width of blending zone in pixels (default: 16 pixels)
            uint16_t patchSize;                     //!< (default: 128)
            uint16_t pupilPixels;                   //!< (default: 64)
            uint16_t saveMask;                      //!< (default: 0)
            std::string outputFileName;
            std::string modeFile;
            std::string pupilFile;
            double wavelength;                      //!< (default: 0, has to be specified)
            /*****************************/
            
        };

        
        /*!
         * Settings with global scope
         */
        struct GlobalCfg : ObjectCfg {     // inherit object (and channel)

            GlobalCfg();
            ~GlobalCfg();

            virtual void parseProperties( bpt::ptree&, const ChannelCfg& defaults=GlobalCfg() );
            virtual void getProperties( bpt::ptree&, const ChannelCfg& defaults=GlobalCfg() ) const;

            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            
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
            double telescopeD;           //!< Telescope diameter
            /*****************************/

            /**** Numerical settings *****/
            uint16_t minIterations;
            uint16_t maxIterations;
            uint16_t targetIterations;  //!< Exit loop after this many successful (i.e. improving) iterations
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
            std::string outputDir;                  //!< Where the output goes (defaults to current directory of jsub)
            std::vector<std::string> outputFiles;
            /*****************************/

            
        };

        /*! @} */
        
    }

}

#endif  // REDUX_MOMFBD_CONFIG_HPP
