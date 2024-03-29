#ifndef REDUX_MOMFBD_OBJECT_HPP
#define REDUX_MOMFBD_OBJECT_HPP

#include "redux/momfbd/config.hpp"
#include "redux/momfbd/channel.hpp"
#include "redux/momfbd/data.hpp"

#include "redux/util/array.hpp"
#include "redux/util/point.hpp"
#include "redux/util/progresswatch.hpp"
#include "redux/work.hpp"

#ifdef RDX_TRACE_JOB
#   include "redux/util/trace.hpp"
#endif

#include <boost/asio.hpp>
#include <boost/property_tree/ptree.hpp>
namespace bpt = boost::property_tree;


namespace redux {
    
    namespace logging {
        class Logger;
    }


    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */

        class MomfbdJob;
        
        /*! @brief Class containing the object-specific configuration a MOMFBD job
         * 
         */
        class Object : public ObjectCfg
#ifdef RDX_TRACE_JOB
#ifndef RDX_TRACE_PROC
    , public redux::util::TraceObject<Object>
#endif
#endif
        {

        public:
            typedef std::shared_ptr<Object> Ptr;

            explicit Object( MomfbdJob&, uint16_t id=0 );
            Object( const Object&, uint16_t id=0, int tid=-1 );
            ~Object();

            void parsePropertyTree( bpt::ptree&, redux::logging::Logger& );
            bpt::ptree getPropertyTree( bpt::ptree&, bool showAll=false );

            uint16_t id(void) const { return ID; };
            size_t size(void) const override;
            uint64_t pack(char*) const override;
            uint64_t unpack(const char*, bool) override;
            
            void cleanup(void);
            uint32_t nImages(bool reCalc=false);
            
            void addChannel( Channel::Ptr obj ) { channels.push_back(obj); }
            Channel::Ptr addChannel( void ) {
                Channel::Ptr chan( new Channel( *this, myJob, channels.size() ) );
                addChannel(chan);
                return chan;
            }
            const std::vector<Channel::Ptr>& getChannels(void) const { return channels; };
            const MomfbdJob& getJob(void) const { return myJob; };
        
            MomfbdJob& myJob;
            logging::Logger& logger;
            std::vector<Channel::Ptr> channels;
            
            /*************   Processing on slave   ***************/
            /*************         Methods         ***************/
            void initProcessing( Solver& );                        //!< Called once per job, to set up storage etc.
            void initPatch(void);                                        //!< Called once per patch, normalize, fit plane etc.
            void getInit(ObjectData&, double* alpha);                    //!< Called on patch-initialization, copy starting values *if any(.
            void initPQ(void);
            void addAllFT(void);
            void addRegGamma(double);
            void addToFT(const complex_t*);
            void addDiffToFT( const complex_t* newFT, const complex_t* oldFT );
            void addToPQ(const complex_t* pp, const double* qq);
            void addDiffToPQ(const redux::image::FourierTransform&, const redux::util::Array<complex_t>&, const redux::util::Array<complex_t>&);
            void addAllPQ(void);
            void calcHelpers(void);
            void fitAvgPlane( redux::util::Array<float>& plane, const std::vector<uint32_t>& wf );
            void fitAvgPlane(void) { fitAvgPlane( fittedPlane, waveFrontList ); };
            void calcMetric(void);
            inline double metric(void) const { return weight*currentMetric; };
            /*****************************************************/
            
            void restorePatch( ObjectData&, const std::vector<uint32_t>& wf );
            void restorePatch( ObjectData& od ) { restorePatch( od, waveFrontList ); };
            
            inline std::string idString( void ) const { return std::to_string(ID); }
            void dump( std::string tag );
        //private:

            bool checkImageScale( double&, double&, double& );
            bool checkCfg(void);
            bool checkData(bool verbose=false);

            void initObject(void);
            void reInitialize(boost::asio::io_service&, redux::util::ProgressWatch& pw, bool reset=false);
            void loadData(boost::asio::io_service&, redux::util::Array<PatchData::Ptr>&);
            void loadInit(boost::asio::io_service&, redux::util::Array<PatchData::Ptr>&);
            size_t getResultSize( void );
            void maybeInitializeStorage( void );          
            void getStorage( PatchData&, std::shared_ptr<ObjectData> );          
            void doMozaic(float**& img, size_t& imgRows, size_t& imgCols, const redux::util::Array<PatchData::Ptr>&);
            void writeAna(const redux::util::Array<PatchData::Ptr>&);
            void writeFits(const redux::util::Array<PatchData::Ptr>&);
            void writeMomfbd(const redux::util::Array<PatchData::Ptr>&);
            void writeResults(boost::asio::io_service&, const redux::util::Array<PatchData::Ptr>&);
            void writeResults(redux::util::Array<PatchData::Ptr>&);
            
            size_t estimateOutputSizeANA(void);
            size_t estimateOutputSizeFITS(void);
            size_t estimateOutputSizeMOMFBD(void);
            size_t estimateOutputSize(void);
            
            redux::util::Point16 getImageSize( bool force=false );
            
            /*************   Processing on slave   ***************/
            /*************     Local variables     ***************/
            std::mutex mtx;
            std::atomic<int> imgShifted;
            redux::util::ProgressWatch progWatch;
            std::shared_ptr<double>  ftSum;                                 //! Sum of the norm of all images' fourier-transforms
            std::shared_ptr<double> Q,PS,QS;
            std::shared_ptr<complex_t> P,PQ;
            redux::util::Array<float> fittedPlane;
            std::shared_ptr<redux::image::Pupil> pupil;
            std::shared_ptr<ModeSet> modes;                                 //!< modes used in this object
            double currentMetric;
            double reg_gamma;
            double frequencyCutoff, pupilRadiusInPixels;
            size_t patchSize2, otfSize, otfSize2;
            std::string cacheFile;
            redux::util::Array<float> results;
            /*****************************************************/

            uint16_t ID;
            int traceID;
            double normalizeTo;
            redux::util::Point16 imgSize;
            mutable uint32_t nObjectImages;
            bpx::ptime startT, endT;
            friend struct ObjectData;
            friend struct WaveFront;
            friend class MomfbdJob;
            friend class Channel;

        };

        /*! @} */
        
    }

}

#endif  // REDUX_MOMFBD_OBJECT_HPP
