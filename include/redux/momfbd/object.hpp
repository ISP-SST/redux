#ifndef REDUX_MOMFBD_OBJECT_HPP
#define REDUX_MOMFBD_OBJECT_HPP

#include "redux/momfbd/config.hpp"
#include "redux/momfbd/channel.hpp"
#include "redux/momfbd/cache.hpp"
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
            
            size_t nImages(void) const;
            
            const std::vector<std::shared_ptr<Channel>>& getChannels(void) const { return channels; };
            const MomfbdJob& getJob(void) const { return myJob; };
        
            MomfbdJob& myJob;
            std::vector<std::shared_ptr<Channel>> channels;
            
            /*************   Processing on slave   ***************/
            /*************         Methods         ***************/
            void initProcessing(WorkSpace::Ptr);
            void initPatch(ObjectData&);
            void getResults(ObjectData&);
            void initPQ(void);
            void addAllFT(void);
            void addToFT(const redux::image::FourierTransform&, double);
            void addDiffToFT( const redux::util::Array<complex_t>& ft, const redux::util::Array<complex_t>& oldft, double );
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

            void init(void);
            void initCache(void);
            void cleanup(void);
            void loadData(boost::asio::io_service&);
            void preprocessData(boost::asio::io_service&);
            void writeAna(const redux::util::Array<PatchData::Ptr>&);
            void writeFits(const redux::util::Array<PatchData::Ptr>&);
            void writeMomfbd(const redux::util::Array<PatchData::Ptr>&);
            void writeResults(const redux::util::Array<PatchData::Ptr>&);
            void normalize(boost::asio::io_service&);
            void prepareStorage(void);
            void storePatches( WorkInProgress&, boost::asio::io_service&, uint8_t );
            
            size_t sizeOfPatch(uint32_t) const;
            
            Point16 getImageSize(void);
            
            /*************   Processing on slave   ***************/
            /*************     Local variables     ***************/
            std::mutex mtx;
            redux::util::Array<double>  ftSum;                //! Sum of the norm of all images' fourier-transforms
            redux::util::Array<double> Q;
            redux::util::Array<complex_t> P;
            redux::util::Array<float> fittedPlane;
            std::set<size_t> pupilIndices, otfIndices;                   //!< Arrays with the offsets where the pupil/otf are greater than some threshold
            double currentMetric;
            double reg_gamma;
            /*****************************************************/

            uint16_t ID;
            double objMaxMean;
            Point16 imgSize;
            uint32_t nObjectImages;
            bpx::ptime startT, endT;
            friend class MomfbdJob;
            friend class Channel;

        };


        /*! @} */
        
    }

}

#endif  // REDUX_MOMFBD_OBJECT_HPP
