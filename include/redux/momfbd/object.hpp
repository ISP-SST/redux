#ifndef REDUX_MOMFBD_OBJECT_HPP
#define REDUX_MOMFBD_OBJECT_HPP

#include "redux/momfbd/config.hpp"
#include "redux/momfbd/channel.hpp"
#include "redux/momfbd/data.hpp"

#include "redux/util/array.hpp"
#include "redux/work.hpp"

#include <boost/asio.hpp>
#include <boost/property_tree/ptree.hpp>
namespace bpt = boost::property_tree;


namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */

        class MomfbdJob;
        
        /*! @brief Class containing the object-specific configuration a MOMFBD job
         * 
         */
        class Object : public ObjectCfg {

        public:

            Object( MomfbdJob&, uint16_t id=0 );
            ~Object();

            void parsePropertyTree( bpt::ptree& );
            bpt::ptree getPropertyTree( bpt::ptree& );

            uint16_t id(void) const { return ID; };
            size_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            
            uint32_t nImages(void) const;
            
            const std::vector<std::shared_ptr<Channel>>& getChannels(void) const { return channels; };
            const MomfbdJob& getJob(void) const { return myJob; };
        
            MomfbdJob& myJob;
            std::vector<std::shared_ptr<Channel>> channels;
            
            /*************   Processing on slave   ***************/
            /*************         Methods         ***************/
            void initProcessing( const Solver& );                        //!< Called once per job, to set up storage etc.
            void initPatch(ObjectData&);                                //!< Called once per patch, normalize, fit plane etc.
            void getResults(ObjectData&);                               //!< Called on patch-completion, gather up what is going back to the master.
            void initPQ(void);
            void addAllFT(void);
            void addRegGamma(double);
            void addToFT(const redux::image::FourierTransform&);
            void addDiffToFT( const redux::util::Array<complex_t>& ft, const redux::util::Array<complex_t>& oldft );
            void addToPQ(const redux::image::FourierTransform&, const redux::util::Array<complex_t>&);
            void addDiffToPQ(const redux::image::FourierTransform&, const redux::util::Array<complex_t>&, const redux::util::Array<complex_t>&);
            void addAllPQ(void);
            void fitAvgPlane(ObjectData&);
            void calcMetric(void);
            double metric(void) const { return currentMetric; };
            /*****************************************************/
            
            void dump( std::string tag );
        //private:

            
            bool checkCfg(void);
            bool checkData(void);

            void initCache(void);
            void loadData(boost::asio::io_service&, redux::util::Array<PatchData::Ptr>&);
            void writeAna(const redux::util::Array<PatchData::Ptr>&);
            void writeFits(const redux::util::Array<PatchData::Ptr>&);
            void writeMomfbd(const redux::util::Array<PatchData::Ptr>&);
            void writeResults(const redux::util::Array<PatchData::Ptr>&);
            void storePatches( WorkInProgress&, boost::asio::io_service&, uint8_t );
            
            Point16 getImageSize(void);
            
            /*************   Processing on slave   ***************/
            /*************     Local variables     ***************/
            std::mutex mtx;
            redux::util::Array<double>  ftSum;                              //! Sum of the norm of all images' fourier-transforms
            redux::util::Array<double> Q;
            redux::util::Array<complex_t> P;
            redux::util::Array<float> fittedPlane;
            redux::image::Pupil pupil;
            ModeSet modes;                                                  //!< modes used in this object
            double currentMetric;
            double reg_gamma;
            double frequencyCutoff, pupilRadiusInPixels;
            /*****************************************************/

            uint16_t ID;
            double objMaxMean;
            Point16 imgSize;
            mutable uint32_t nObjectImages;
            bpx::ptime startT, endT;
            friend class MomfbdJob;
            friend class Channel;

        };


        /*! @} */
        
    }

}

#endif  // REDUX_MOMFBD_OBJECT_HPP
