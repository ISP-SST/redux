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

        extern const std::map<std::string, int> getstepMap;
        extern const std::map<std::string, int> gradientMap;
        extern const std::map<std::string, int> fillpixMap;

        class Channel;
        class MomfbdJob : public Job {

        public:
            static size_t jobType;
            static int getFromMap( std::string str, const std::map<std::string, int>& m );
            static void maybeOverride( bool value, uint32_t& set, uint32_t flag );

            MomfbdJob( void );
            ~MomfbdJob( void );

            void parseProperties( po::variables_map& vm, bpt::ptree& tree );
            bpt::ptree getProperties( bpt::ptree* );

        private:

            uint8_t basis;
            std::vector<uint32_t> modes, imageNumbers, darkNumbers;
            std::string imageDataDir;
            std::vector<std::string> outputFiles;
            int nPoints, sequenceNumber;
            double reg_gamma;
            uint32_t flags;

            void* prePart( void );
            void postPart( void* );
            void* runPreJob( void );
            void runPostJob( void* );

            uint32_t preProcess( void );
            uint32_t postProcess( void );
            uint32_t runJob( void );

            std::vector<std::shared_ptr<Object>> objects;

            int klMinMode, klMaxMode, borderClip;
            double telescopeFocalLength, telescopeDiameter, arcSecsPerPixel, pixelSize;
            int minIterations , maxIterations, nDoneMask;
            double FTOL, EPS, svd_reg;
            std::string programDataDir;
            std::vector<double> stokesWeights;

            int gradient_method, getstep_method;
            uint8_t fillpix_method;
            std::string time_obs, date_obs;
            int max_local_shift, mstart, mstep;

            redux::util::Array<double> pupil;
            int pupilSize;
            uint8_t output_data_type, fp_method;
            int nsx, nsy, ncal;

            std::vector<uint32_t> xl, xh, yl, yh;

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

                ar & pupilSize & output_data_type & fp_method;
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

    }   // momfbd

}   // redux

#endif  // REDUX_MOMFBD_MOMFBDJOB_HPP
