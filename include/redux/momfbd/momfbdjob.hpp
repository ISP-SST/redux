#ifndef REDUX_MOMFBD_MOMFBDJOB_HPP
#define REDUX_MOMFBD_MOMFBDJOB_HPP

#include "redux/momfbd/object.hpp"
#include "redux/job.hpp"
#include "redux/serialization.hpp"
#include "redux/util/array.hpp"

#include <map>

#include <boost/serialization/vector.hpp>

namespace redux {

    namespace momfbd {
        
        /*! @defgroup momfbd MOMFBD
         *  @{
         */
        
        extern const std::map<std::string, int> getstepMap;
        extern const std::map<std::string, int> gradientMap;
        extern const std::map<std::string, int> fillpixMap;

        class Channel;
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

            void parseProperties( po::variables_map& vm, bpt::ptree& tree );
            bpt::ptree getPropertyTree( bpt::ptree* root=nullptr );

            size_t size(void) const;
            char* pack(char*) const;
            const char* unpack(const char*, bool);
        
        private:

            void* prePart( void );
            void postPart( void* );
            void* runPreJob( void );
            void runPostJob( void* );

            uint32_t preProcess( void );
            uint32_t postProcess( void );
            uint32_t runJob( void );
            

            char basis,fillpix_method, output_data_type;
            
            uint32_t flags, nPoints, sequenceNumber;
            uint32_t klMinMode, klMaxMode, borderClip;
            uint32_t minIterations , maxIterations, nDoneMask;    
            uint32_t gradient_method, getstep_method;
            uint32_t max_local_shift, mstart, mstep;
            uint32_t pupilSize, nsx, nsy, ncal;
            std::vector<uint32_t> modes, imageNumbers, darkNumbers;
            std::vector<uint32_t> xl, xh, yl, yh;

            double telescopeFocalLength, telescopeDiameter, arcSecsPerPixel, pixelSize;
            double reg_gamma, FTOL, EPS, svd_reg;
            std::vector<double> stokesWeights;
            
            std::string imageDataDir;
            std::string programDataDir;
            std::string time_obs, date_obs;
            std::vector<std::string> outputFiles;

            std::vector<std::shared_ptr<Object>> objects;

            redux::util::Array<double> pupil;

            // serialize data to store/send it in compact form
            template <typename Archive>
            void serialize( Archive& ar, const unsigned int version ) {

                ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( Job );
                ar & basis & modes;
                ar & imageNumbers & darkNumbers & imageDataDir;
                ar & outputFiles & nPoints & sequenceNumber;
                ar & reg_gamma & flags;

                ar & klMinMode & klMaxMode & borderClip;
                ar & telescopeFocalLength & telescopeDiameter & arcSecsPerPixel;
                ar & pixelSize & minIterations & maxIterations & nDoneMask;
                ar & FTOL & EPS & svd_reg;
                ar & programDataDir & stokesWeights;

                ar & gradient_method & getstep_method & fillpix_method;
                ar & time_obs & date_obs & max_local_shift;
                ar & mstart & mstep;

                ar & pupilSize & output_data_type;
                ar & nsx & nsy;
                ar & ncal & xl & xh & yl & yh;

                ar & objects;

                // TODO
                //     redux::util::Array<double> pupil;

            }

            friend class boost::serialization::access;
            friend class Object;
            friend class Channel;

        };

        const size_t dummy = MomfbdJob::jobType;       // this will trigger the registration of MomfbdJob in Job::jobMap

        /*! @} */

    }   // momfbd

}   // redux

#endif  // REDUX_MOMFBD_MOMFBDJOB_HPP
