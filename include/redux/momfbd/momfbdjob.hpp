#ifndef REDUX_MOMFBD_MOMFBDJOB_HPP
#define REDUX_MOMFBD_MOMFBDJOB_HPP

#include "redux/momfbd/objectcfg.hpp"
#include "redux/momfbd/patch.hpp"
#include "redux/job.hpp"
#include "redux/util/array.hpp"

#include <map>

namespace redux {

    namespace momfbd {

        /*! @defgroup momfbd MOMFBD
         *  @{
         */

        extern const std::map<std::string, int> getstepMap;
        extern const std::map<std::string, int> gradientMap;
        extern const std::map<std::string, int> fillpixMap;

        class ChannelCfg;
        /*! @brief Class containing the configuration settings for a MOMFBD job.
         *
         */
        class MomfbdJob : public Job {

        public:
            static size_t jobType;
            static int getFromMap( std::string str, const std::map<std::string, int>& m );
            static void maybeOverride( bool value, uint32_t& set, uint32_t flag );

            MomfbdJob( void );
            ~MomfbdJob( void );

            uint64_t unpackParts(const char* ptr, std::vector<Part::Ptr>&, bool);
            
            void parseProperties( po::variables_map& vm, bpt::ptree& tree );
            bpt::ptree getPropertyTree( bpt::ptree* root = nullptr );

            size_t size( void ) const;
            uint64_t pack( char* ) const;
            uint64_t unpack( const char*, bool );

            size_t getParts( WorkInProgress&, uint8_t );
            void ungetParts( WorkInProgress& );
            void returnParts( WorkInProgress& );

            void init(void);
            void cleanup(void);
            bool run( WorkInProgress&, boost::asio::io_service&, boost::thread_group&, uint8_t );
            
            bool checkCfg(void);
            bool checkData(void);
            const std::vector<std::shared_ptr<ObjectCfg>>& getObjects(void) const { return objects; };
        
            uint8_t basis, fillpix_method, output_data_type;

            uint32_t flags, patchSize, sequenceNumber;
            uint32_t klMinMode, klMaxMode, borderClip;
            uint32_t minIterations , maxIterations, nDoneMask;
            uint32_t gradient_method, getstep_method;
            uint32_t max_local_shift, nInitialModes, nModesIncrement;
            uint32_t pupilSize, nPatchesX, nPatchesY, ncal;

            double telescopeFocalLength, telescopeDiameter, arcSecsPerPixel, pixelSize;
            double reg_gamma, FTOL, EPS, svd_reg;

            std::vector<uint32_t> modes, imageNumbers, darkNumbers;
            std::vector<uint32_t> subImagePosX,subImagePosY;
            std::vector<double> stokesWeights;

        private:

            void checkParts( void );
            void preProcess( boost::asio::io_service&, boost::thread_group& );
            void applyLocalOffsets( PatchData::Ptr );
            void runMain( WorkSpace& );
            void postProcess( boost::asio::io_service&, boost::thread_group& );


            std::string imageDataDir;
            std::string programDataDir;
            std::string time_obs, date_obs;
            std::string pupilFile;
            std::vector<std::string> outputFiles;

            std::vector<std::shared_ptr<ObjectCfg>> objects;

            redux::util::Array<double> pupil;
            redux::util::Array<float> imageStack;

            std::map<size_t,PatchData::Ptr> patches;
        
            friend class ObjectCfg;
            friend class ChannelCfg;


        };

        const size_t dummy = MomfbdJob::jobType;       // this will trigger the registration of MomfbdJob in Job::jobMap

        /*! @} */

    }   // momfbd

}   // redux

#endif  // REDUX_MOMFBD_MOMFBDJOB_HPP
