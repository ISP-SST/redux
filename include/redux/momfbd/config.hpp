#ifndef REDUX_MOMFBD_CONFIG_HPP
#define REDUX_MOMFBD_CONFIG_HPP

#include "redux/logging/logger.hpp"
#include "redux/momfbd/modes.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/util/trace.hpp"

#include <map>
#include <ostream>
#include <string>
#include <vector>

#include <boost/property_tree/ptree.hpp>
namespace bpt = boost::property_tree;


namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */

        enum FileType { FT_NONE=0, FT_ANA, FT_FITS, FT_MOMFBD=4, FT_MASK=7 };
        extern const std::map<FileType, std::string> FileTypeNames;
        extern const std::map<FileType, std::string> FileTypeExtensions;
        
        enum DataType  { DT_I08T=0, DT_I16T, DT_I32T, DT_I64T, DT_F32T, DT_F64T };
        enum SaveFlags { SF_SAVE_ALPHA=1, SF_SAVE_COBJ, SF_SAVE_DIVERSITY=4, SF_SAVE_METRIC=8,
                         SF_SAVE_MODES=16, SF_SAVE_PSF=32, SF_SAVE_PSF_AVG=64, SF_SAVE_RESIDUAL=128,
                         SF_SAVE_NAMES=256, SF_SAVE_FFDATA=512 };
        enum RunFlags  { RF_CALIBRATE=1, RF_DONT_MATCH_IMAGE_NUMS, RF_FAST_QR=4, RF_FIT_PLANE=8,
                         RF_FLATFIELD=16, RF_GLOBAL_NOISE=32, RF_NEW_CONSTRAINTS=64, RF_NO_CLIP=128,
                         RF_NO_CONSTRAINTS=256, RF_NO_FILTER=512, RF_FORCE_WRITE=1024, RF_NOSWAP=2048,
                         RF_OLD_NS=4096, RF_SORT_MODES=8192 };
        enum NormType { NORM_NONE=0, NORM_OBJ_MAX_MEAN, NORM_OBJ_MAX_MEDIAN, NORM_OBJ_MEDIAN_MEDIAN };
        
        struct cicomp {  // case-insensitive comparator for the maps below.
            bool operator() ( const std::string& a,const std::string& b ) const { return redux::util::nocaseLess(a,b); }
        };

        enum FillpixMethod { FPM_MEDIAN=1, FPM_INVDISTWEIGHT, FPM_HORINT };
        extern const std::map<std::string, int, cicomp> fillpixMap;
        enum GradientMethod { GM_DIFF=1, GM_VOGEL };
        extern const std::map<std::string, int, cicomp> gradientMap;
        enum GetstepMethod { GSM_SDSC=1, GSM_CNJG, GSM_BFGS, GSM_BFGS_inv };
        extern const std::map<std::string, int, cicomp> getstepMap;
        extern const std::map<std::string, int, cicomp> normMap;
        
        
        /*!
         * Settings with channel scope
         */
        struct ChannelCfg
#ifdef RDX_TRACE_PROC
            : public redux::util::TraceObject<ChannelCfg>
#endif
        {

            ChannelCfg();
            virtual ~ChannelCfg();
            
            operator std::string() const;           //!< cast to string (for easy printing of configuration)

            virtual void parseProperties( bpt::ptree&, redux::logging::Logger&, const ChannelCfg& defaults=ChannelCfg() );
            virtual void getProperties( bpt::ptree&, const ChannelCfg& defaults=ChannelCfg(), bool showAll=false ) const;
            
            virtual uint64_t size(void) const;
            virtual uint64_t pack(char*) const;
            virtual uint64_t unpack(const char*, bool);
            
            bool operator==(const ChannelCfg&) const;
            
            
            /********* Hardware **********/
            double rotationAngle;                   //!< Rotation of this camera relative to "anchor channel" (default: 0)
            /*****************************/
            
            /**** Numerical settings *****/
            double noiseFudge;                      //!< Noise weight (default: 1)
            double weight;                          //!< Weight for this channel in calculations (default: 1)
            /*****************************/

            /********  Diversity  ********/
            // TODO: reorganize
            ModeBase diversityBasis;                //!< Which basis to use as default for unmarked modes.
            ModeList diversityModes;                //!< List of diversity mode types/numbers
            std::vector<DiversityValue> diversityValues;    //!< List of weights/values for the diversity modes. 
            std::string diversityModeFile;          //!< File containining modes to be used for the diversity.
            bool divModeFileNormalize;              //!< flag to disable normalization of modes supplied in file (default is enabled/true)
            bool noRestore;                         //!< Exclude this channel in the final reconstruction/deconvolution step (i.e. only use it during the fitting)
            /*****************************/
            
            /******* Data settings *******/
            std::vector<float> alignMap;            //!< Coefficients for the projective transformation (from reference channel) into this channel (default: none)
            std::vector<int16_t> alignClip;         //!< Crop images to this region {firstX,lastX,firstY,lastY}, (default: none, has to be specified)
            uint16_t borderClip;                    //!< Disregard this many pixels at the edge when calculating statistics (default: 100)
            uint8_t incomplete;                     //!< Some files might not exist, just skip those.
            std::vector<uint16_t> subImagePosXY, subImagePosX, subImagePosY;    //!< Patch coordinates
            std::vector<uint16_t> discard;          //!< Skip this many frames in the beginning/end of input data files.
            /*****************************/

            /************ Input **********/
            std::string imageDataDir;               //!< Where the data is located
            std::string imageTemplate;              //!< Filename template for the images, IMAGE_NUM(S) will be inserted into this template to generate filenames.
            std::string darkTemplate;               //!< Filename template for the darks, DARK_NUM will be inserted into this template to generate filenames.
            std::string gainFile;                   //!< Gain file to use.      (will be applied *after* backscatter correction)
            std::string responseFile;               //!< File with CCD response (will be applied *before* backscatter correction)
            std::string backgainFile;               //!< File with the back-gain to use for the backscatter correction.
            std::string psfFile;                    //!< File with the PSF to use for the backscatter correction.
            std::string mmFile;                     //!< Modulation matrix  (for including Stokes demodulation in the fitting. **N.B. Not fully implemented/tested!**)
            uint8_t mmRow;                          //!< Number of rows in modulation matrix
            uint8_t mmWidth;                        //!< Number of cols in modulation matrix
            std::string xOffsetFile,yOffsetFile;    //!< Alignment offsets (from pinhole-calibration)
            uint32_t imageNumberOffset;             //!< Add this offset to each image number in this channel
            std::vector<uint32_t> fileNumbers;      //!< Use these numbers together with the template to generate file-list
            std::vector<uint32_t> waveFrontList;    //!< Identify wavefront, used to group/constrain simultaneous images if image-numbers can't be used.
            std::vector<uint32_t> darkNumbers;      //!< Use these numbers together with the template to generate file-list
            std::vector<float> stokesWeights;       //!< Weights to use for the different Stokes components if .
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
            
            virtual void parseProperties( bpt::ptree&, redux::logging::Logger&, const ChannelCfg& defaults=ObjectCfg() ) override;
            virtual void getProperties( bpt::ptree&, const ChannelCfg& defaults=ObjectCfg(), bool showAll=false ) const override;

            uint64_t size(void) const override;
            uint64_t pack(char*) const override;
            uint64_t unpack(const char*, bool) override;
            
            const ObjectCfg& operator=(const ChannelCfg&);
            bool operator==(const ObjectCfg&) const;
            
            /********* Hardware **********/
            double telescopeF;                       //!< Telescope focal length
            double arcSecsPerPixel;                  //!< Image scale   (default: 0, has to be specified)
            double pixelSize;                        //!< Physical size of pixels (default: 10\f$\mu\f$)
            /*****************************/
            
            /******* Data settings *******/
            uint16_t maxLocalShift;                 //!< How much are the patches allowed to be shifted (default: 5 pixels)
            uint16_t minimumOverlap;                //!< Desired width of blending zone in pixels (default: 16 pixels)
            uint16_t patchSize;                     //!< (default: 128)
            uint16_t pupilPixels;                   //!< (default: 64)
            uint16_t saveMask;                      //!< (default: 0)
            std::string outputFileName;
            std::string initFile;
            std::string modeFile;
            std::string pupilFile;
            double wavelength;                      //!< (default: 0, has to be specified)
            bool modeFileNormalize;                 //!< flag to disable normalization of modes supplied in file (default is enabled/true)
            bool traceObject;                       //!< specifies that this object should be used as a reference for spatial distortions.
            /*****************************/
            
        };

        
        /*!
         * Settings with global scope
         */
        struct GlobalCfg : ObjectCfg {     // inherit object (and channel)

            GlobalCfg();
            ~GlobalCfg();

            virtual void parseProperties( bpt::ptree&, redux::logging::Logger&, const ChannelCfg& defaults=GlobalCfg() ) override;
            virtual void getProperties( bpt::ptree&, const ChannelCfg& defaults=GlobalCfg(), bool showAll=false ) const override;

            uint64_t size(void) const override;
            uint64_t pack(char*) const override;
            uint64_t unpack(const char*, bool) override;
            
            const GlobalCfg& operator=(const ObjectCfg&);
            const GlobalCfg& operator=(const ChannelCfg&);
            bool operator==(const GlobalCfg&) const;
            
            uint16_t runFlags;

            /*********** Modes ***********/
            ModeBase modeBasis;         //!< Which basis to use for the fitting
            uint16_t klMinMode;         //!< First Zernike-mode to be considered in Karhunen-Loève expansion
            uint16_t klMaxMode;         //!< Last Zernike-mode to be considered in Karhunen-Loève expansion
            float klCutoff;             //!< If the expansion-coefficient is smaller than this, it will be ignored.
            uint16_t nInitialModes;     //!< How many modes to use in the first iteration
            uint16_t nModeIncrement;    //!< How many modes to add in each iteration
            uint16_t nModes;
            ModeList modeList;        //!< Which modes to use
            /*****************************/
            
            /********* Hardware **********/
            double telescopeD;           //!< Telescope diameter
            double telescopeCO;          //!< Telescope central obscuration diameter
            /*****************************/

            /**** Numerical settings *****/
            uint16_t minIterations;
            uint16_t maxIterations;
            uint16_t targetIterations;  //!< Exit loop after this many successful (i.e. improving) iterations
            uint8_t fillpixMethod;
            uint8_t gradientMethod;
            uint8_t getstepMethod;
            uint8_t normType;
            int16_t apodizationSize;
            float badPixelThreshold;
            float FTOL;
            float EPS;
            float reg_alpha;
            float graddiff_step;            //!< step-length when calculating numerical derivative
            bool trace;                     //!< specifies that this object should be used as a reference for spatial distortions.
            /*****************************/

            /******* Data settings *******/
            uint8_t outputFileType;
            uint8_t outputDataType;
            uint32_t sequenceNumber;        //!< Not used
            std::string observationTime;    //!< Not used
            std::string observationDate;    //!< Date-tag to be placed in the header of the ooutput files.
            std::string tmpDataDir;
            std::vector<std::string> outputFiles;
            std::vector<std::string> initFiles;
            /*****************************/

            
        };

        /*! @} */
        
    }

}

#endif  // REDUX_MOMFBD_CONFIG_HPP
