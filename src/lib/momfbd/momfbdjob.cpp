#include "redux/momfbd/momfbdjob.hpp"

#include "redux/momfbd/util.hpp"

#include "redux/logger.hpp"
#include "redux/util/bitoperations.hpp"
#include "redux/util/stringutil.hpp"

#include <thread>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
using boost::algorithm::iequals;

using namespace redux::momfbd;
using namespace redux::momfbd::util;
using namespace redux::util;
using namespace redux;
using namespace std;

#define lg Logger::mlg
namespace {

    const string thisChannel = "momfbdjob";
    static Job* createMomfbdJob( void ) {
        return new MomfbdJob();
    }
}
size_t MomfbdJob::jobType = Job::registerJob( "momfbd", createMomfbdJob );


MomfbdJob::MomfbdJob( void ) {
    info.typeString = "momfbd";
}


MomfbdJob::~MomfbdJob( void ) {
    
}


uint64_t MomfbdJob::unpackParts( const char* ptr, WorkInProgress& wip, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count(0);
    if( wip.nParts ) {
        wip.parts.resize(1);
        wip.parts[0].reset(new PatchData(*this));
        count += wip.parts[0]->unpack( ptr+count, swap_endian );
        if( wip.nParts > 1 ) {
            globalData.reset( new GlobalData(*this) );
            count += globalData->unpack( ptr+count, swap_endian );
        }
    }
    return count;
}


void MomfbdJob::parsePropertyTree( bpo::variables_map& vm, bpt::ptree& tree ) {

    Job::parsePropertyTree( vm, tree );
    //LOG_DEBUG << "MomfbdJob::parsePropertyTree()";
    
    // possibly override cfg-entries with command-line arguments
    if( vm.count( "simx" ) ) tree.put( "SIM_X", vm["simx"].as<string>() );
    if( vm.count( "simy" ) ) tree.put( "SIM_Y", vm["simy"].as<string>() );
    if( vm.count( "imgn" ) ) tree.put( "IMAGE_NUM", vm["imgn"].as<string>() );
    if( vm.count( "output-file" ) ) tree.put( "output-file", vm["output-file"].as<string>() );
    if( vm.count( "force" ) ) tree.put( "OVERWRITE", true );
    if( vm.count( "swap" ) ) tree.put( "SWAP", true );

    GlobalCfg::parseProperties(tree);

    uint16_t nObj(0);
    for( auto & property : tree ) {
        if( iequals( property.first, "OBJECT" ) ) {
            Object* tmpObj = new Object( *this, nObj );
            tmpObj->parsePropertyTree( property.second );
            if( nObj < outputFiles.size() ) {
                tmpObj->outputFileName = outputFiles[nObj];
            }
            objects.push_back( shared_ptr<Object>( tmpObj ) );
            nObj++;
        }
    }
    if( outputFiles.size() > objects.size() ) {
        LOG_WARN << outputFiles.size() << " output file names specified but only " << objects.size() << " objects found.";
    }
    //LOG_DEBUG << "MomfbdJob::parsePropertyTree() done.";

}


bpt::ptree MomfbdJob::getPropertyTree( bpt::ptree* root ) {

    bpt::ptree tree = Job::getPropertyTree();         // get Job-properties

    GlobalCfg::getProperties(tree);

    for( auto& obj: objects ) {
        obj->getPropertyTree( tree );
    }

    if( root ) {
        root->push_back( bpt::ptree::value_type( "momfbd", tree ) );
    }

    return tree;
}


uint64_t MomfbdJob::size( void ) const {
    uint64_t sz = Job::size();
    sz += GlobalCfg::size();
    sz += sizeof(uint16_t);           // objects.size()
    for( const auto& obj: objects ) {
        sz += obj->size();
    }
    return sz;
}


uint64_t MomfbdJob::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = Job::pack( ptr );
    count += GlobalCfg::pack( ptr+count );
    count += pack( ptr+count, (uint16_t)objects.size() );
    for( const auto& obj: objects ) {
        count += obj->pack( ptr+count );
    }
    
    return count;
    
}


uint64_t MomfbdJob::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = Job::unpack( ptr, swap_endian );
    count += GlobalCfg::unpack( ptr+count, swap_endian );
    uint16_t tmp;
    count += unpack( ptr+count, tmp, swap_endian );
    objects.resize( tmp );
    for( auto& obj: objects ) {
        obj.reset(new Object(*this));
        count += obj->unpack( ptr+count, swap_endian );
    }
    return count;
}


size_t MomfbdJob::nImages(void) const {
    size_t nTotalImages(0);
    for( const auto& obj : objects) {
        nTotalImages += obj->nImages();
    }
    return nTotalImages;
}


void MomfbdJob::checkParts( void ) {

    uint8_t mask = 0;
    for( auto & patch : patches ) {
        /*if( (patch.second->step & JSTEP_ERR) && (patch.second->nRetries<info.maxPartRetries)) {
            patch.second->nRetries++;
            patch.second->step &= ~JSTEP_ERR;
        }*/
        mask |= patch->step;
    }

    if( mask & JSTEP_ERR ) {    // TODO: handle failed parts.

    }

    if( countBits( mask ) == 1 ) {  // if all parts have the same "step", set the whole job to that step.
        info.step.store( mask );
    }

}


bool MomfbdJob::getWork( WorkInProgress& wip, uint8_t nThreads ) {

    bool ret(false);
    uint8_t step = info.step.load();
    wip.parts.clear();

     // run pre-/postprocessing if local
    if( ((step == JSTEP_PREPROCESS)||(step == JSTEP_POSTPROCESS)) && !wip.connection ) {
        ret = true;
    }
    
    if ( step == JSTEP_QUEUED ) {                       // preprocessing ready -> start
        info.step = step = JSTEP_RUNNING;
    }

    if ( !ret && (step == JSTEP_RUNNING) ) {                      // running
        unique_lock<mutex> lock( jobMutex );
//         size_t nParts = wip.peer->status.nThreads;
//         if( info.nThreads ) nParts = std::min( wip.peer->status.nThreads, info.nThreads );
       /* if(!wip.connection) {   // local worker, check if there are results to write
            for( auto & patch : patches ) {
                if( patch->step & JSTEP_POSTPROCESS ) {
                    LOG_DEBUG << "getWork(): PP-patch   step = " << bitString(patch->step);
                    wip.parts.push_back( patch );
                }
            }
            if( wip.parts.size() ) {
                LOG_DEBUG << "getWork(): nPP = " << wip.parts.size();
                ret = true;
            }
        }*/
        if(!ret && wip.connection) {
            for( auto & patch : patches ) {
                if( patch->step == JSTEP_QUEUED ) {
                    patch->step = JSTEP_RUNNING;
                    wip.parts.push_back( patch );
                    if(wip.previousJob.get() != this) {     // First time for this slave -> include global data
                        wip.parts.push_back( globalData );
                    }
                    ret = true;
                    break;// only 1 part at a time for MomfbdJob
                }
            }
        }
    }
    
    //  LOG_DEBUG << "getWork(): step = " << (int)step << " conn = " << (bool)wip.connection;
    if( ret ) {
        unique_lock<mutex> lock( jobMutex );
        checkParts();
    }
    wip.nParts = wip.parts.size();
    return ret;
}


void MomfbdJob::ungetWork( WorkInProgress& wip ) {
    unique_lock<mutex> lock( jobMutex );
    for( auto& part : wip.parts ) {
        part->step = JSTEP_QUEUED;
    }
    wip.parts.clear();
}


void MomfbdJob::returnResults( WorkInProgress& wip ) {
    unique_lock<mutex> lock( jobMutex );
    checkParts();
    for( auto& part : wip.parts ) {
        auto patch = static_pointer_cast<PatchData>( part );
        patch->step = JSTEP_POSTPROCESS;
        patches(patch->index.y,patch->index.x) = patch;
    }
    checkParts();
}


bool MomfbdJob::run( WorkInProgress& wip, boost::asio::io_service& service, uint8_t maxThreads ) {
    
    uint8_t jobStep = info.step.load();
    uint8_t patchStep = 0;
    if(wip.parts.size() && wip.parts[0]) patchStep = wip.parts[0]->step;
    if( jobStep == JSTEP_PREPROCESS ) {
        preProcess(service);                           // preprocess on master: load, flatfield, split in patches
    }
    else if( jobStep == JSTEP_RUNNING || jobStep == JSTEP_QUEUED ) {
        uint8_t nThreads = std::min( maxThreads, info.maxThreads);
        if( patchStep == JSTEP_POSTPROCESS ) {      // patch-wise post-processing, do we need any for momfbd?  the filtering etc. is done on the slaves.
        } else {                                    // main processing
            if(!globalData) {
                 globalData.reset( new GlobalData(*this) );
            }
            if( !solver ) {
                solver.reset( new Solver(*this) );            // Initialize, allocations, etc.
                solver->init();
            }
            for( auto& part : wip.parts ) {      // momfbd jobs will only get 1 part at a time, this is just to keep things generic.
                // Run main processing
                solver->run(static_pointer_cast<PatchData>(part), service, nThreads);
                // Get results
                //it = proc->result;
            }
        }
    }
    else if( jobStep == JSTEP_POSTPROCESS ) {
        postProcess(service);                          // postprocess on master, collect results, save...
    }
    else {
        LOG << "MomfbdJob::run()  unrecognized step = " << ( int )info.step.load();
        info.step.store( JSTEP_ERR );
    }
    return false;
    
}


void MomfbdJob::preProcess( boost::asio::io_service& service ) {

    // TODO: start logging (to file)

    LOG_TRACE << "MomfbdJob::preProcess()";
    
    if ( !checkData() ) {
        LOG_ERR << "MomfbdJob::preProcess(): sanity check failed.";
        info.step.store( JSTEP_ERR );
        info.state.store( JSTATE_IDLE );
        return;
    }
    
    Point16 imageSizes;
    for( auto& obj : objects ) {
        Point16 tmp = obj->getImageSize();
        if(imageSizes == 0) {
            imageSizes = tmp;
        } else if( tmp != imageSizes ) {    // TBD: allow for different patchsizes (i.e. pixelsize/ccd-size) for different objects/channels.
            throw std::logic_error("The clipped images have different sizes for the different objects, please verify the ALIGN_CLIP values.");
        }
    }
    //service.post( std::bind( &MomfbdJob::initCache, this) );    // TBD: should cache initialization be parallelized?
    std::thread initcache(std::bind( &MomfbdJob::initCache, this));     // run in background while loading data.
    
    uint16_t halfPatchSize = patchSize/2;
    uint16_t totalOverlap = minimumOverlap+patchSize/4;     // from MvN: always overlap 25% + 16 pixels.
    // TODO: do split per channel instead, to allow for different image-scales and/or hardware
    // NOTE:  subImagePosX/Y are kept 1-based, so ffset by 1 during cut-out.
    LOG << boost::format("MomfbdJob::preProcess(): halfPatchSize=%d  overlap=%d") % halfPatchSize % totalOverlap;
    if( subImagePosX.empty() ) { // x-coordinate of patch-centre
        subImagePosX = segment<uint16_t>(halfPatchSize+1,imageSizes.x-halfPatchSize+1,patchSize,totalOverlap);
        LOG << "MomfbdJob::preProcess(): Generated patch positions  " << printArray(subImagePosX,"X");
    }
    if( subImagePosY.empty() ) { // y-coordinate of patch-centre
        subImagePosY = segment<uint16_t>(halfPatchSize+1,imageSizes.y-halfPatchSize+1,patchSize,totalOverlap);
        LOG << "MomfbdJob::preProcess(): Generated patch positions  " << printArray(subImagePosY,"Y");
    }
 
    if( subImagePosX.empty() || subImagePosY.empty() ) {
        LOG_ERR << "MomfbdJob::preProcess(): No patches specified or generated, can't continue.";
        info.step.store( JSTEP_ERR );
        info.state.store( JSTATE_IDLE );
        return;
    }

    for( uint16_t& pos : subImagePosY ) {
        uint16_t adjustedPos = std::min<uint16_t>(std::max<uint16_t>(halfPatchSize+1,pos),imageSizes.y-halfPatchSize+1);       // stay inside borders
        if( adjustedPos != pos ) {
            LOG_WARN << "MomfbdJob::preProcess() y-position of patch is too close to the border, adjusting: " << pos << " -> " << adjustedPos;
            pos = adjustedPos;
        }
    }

    for( uint16_t& pos : subImagePosX ) {
        uint16_t adjustedPos = std::min<uint16_t>(std::max<uint16_t>(halfPatchSize+1,pos),imageSizes.x-halfPatchSize+1);       // stay inside borders
        if( adjustedPos != pos ) {
            LOG_WARN << "MomfbdJob::preProcess() x-position of patch is too close to the border, adjusting: " << pos << " -> " << adjustedPos;
            pos = adjustedPos;
        }
    }
    uint64_t count(0);
    patches.resize(subImagePosY.size(),subImagePosX.size());
    Point16 ps(patchSize,patchSize);
    for( uint y=0; y<subImagePosY.size(); ++y ) {
        for( uint x=0; x<subImagePosX.size(); ++x ) {
            PatchData::Ptr patch( new PatchData(*this, y, x ) );
            patch->setPath(to_string(info.id));
            patch->step = JSTEP_QUEUED;
            patch->roi.first = Point16(subImagePosY[y]-halfPatchSize-1, subImagePosX[x]-halfPatchSize-1);   // subImagePosX/Y is 1-based
            patch->roi.last = patch->roi.first+ps-1;
            patch->id = ++count;
            patches(y,x) = std::move(patch);
        }
    }


    for( auto& obj : objects ) {
        obj->loadData(service, patches);
    }
    
    for( auto& obj : objects ) {
        for (auto& ch : obj->channels) {
            ch->unloadCalib();
        }
    }
    
    initcache.join();           // wait for background jobs.
    bool writeFailed(false);
    uint32_t nTotalImages(0);
    for( auto& obj : objects ) {
        for( auto& ch : obj->channels ) {
            writeFailed |= ch->patchWriteFail.get();            // FIXME only used as a synch-point atm. It should deal with errors eventually.
        }
        nTotalImages += obj->nImages();
    }
    if(writeFailed) LOG_ERR<< "MomfbdJob::preProcess()  Hmmmmm, that's odd, I haven't implemented any error-reporting yet.";

    LOG_DETAIL << "MomfbdJob::preProcess()  Done.  nPatches = " << patches.nElements() << "   nObjects = " << objects.size() << "   nImages = " << nTotalImages;

    info.step.store( JSTEP_QUEUED );

}

void MomfbdJob::initCache(void) {
    LOG_DETAIL << "MomfbdJob::initCache()";
    if(!globalData) {       // create GlobalData if they don't exist
        globalData.reset(new GlobalData(*this));
    }
    
    globalData->constraints.init();
 
    for( auto& obj: objects ) {
        obj->initCache();
    }

    LOG_DETAIL << "MomfbdJob::initCache()  Done.";
}
 

void MomfbdJob::storePatches( WorkInProgress& wip, boost::asio::io_service& service, uint8_t nThreads) {
    LOG << "MomfbdJob::storePatches()";
    for( auto& obj: objects ) {
        obj->storePatches(wip, service, nThreads);
    }
    
    for( auto& part: wip.parts ) {
        //auto patch = static_pointer_cast<PatchData>(part);
        part->step = JSTEP_COMPLETED;
    }

    
}


void MomfbdJob::postProcess( boost::asio::io_service& service ) {

    LOG << "MomfbdJob::postProcess()";
    for( auto& patch: patches ) {
        patch->cacheLoad(true);        // load and erase cache-file.
    }
    for( auto& obj : objects ) {
        service.post( std::bind( &Object::writeResults, obj.get(), patches ) );
    }
    runThreadsAndWait(service, 1); //objects.size());  TODO: fix multithreaded write
    
    patches.clear();
    
    info.step.store( JSTEP_COMPLETED );
    info.state.store( 0 );

}


bool MomfbdJob::check(void) {
    bool ret(false);
    unique_lock<mutex> lock(jobMutex);
    int val = info.step;
    switch (val) {
        case 0:                 ret = checkCfg(); if(ret) info.step = JSTEP_SUBMIT; break;
        case JSTEP_SUBMIT:      ret = checkData(); if(ret) info.step = JSTEP_PREPROCESS; break;
        case JSTEP_PREPROCESS: ;                  // no checks at these steps, just fall through and return true
        case JSTEP_QUEUED: ;
        case JSTEP_RUNNING: ;
        case JSTEP_POSTPROCESS: ;
        case JSTEP_COMPLETED: ret = true; break;
        default: LOG_ERR << "check(): No check defined for step = " << (int)info.step << " (" << stepString(info.step) << ")";
    }
    return ret;
}


bool MomfbdJob::checkCfg(void) {
    
    if( (runFlags&RF_FLATFIELD) && (runFlags&RF_CALIBRATE) ) {
        LOG_ERR << "Both FLATFIELD and CALIBRATE mode requested";
        return false;
    }
    if( objects.empty() ) return false;     // need at least 1 object
    
    for( auto& obj: objects ) {
        if( !obj->checkCfg() ) return false;
    }
    
    return true;
}


bool MomfbdJob::checkData(void) {

    for( auto& obj: objects ) {
        if( !obj->checkData() ) return false;
    }
    
    return true;
}
        
        
const MomfbdJob& MomfbdJob::operator=(const GlobalCfg& rhs) {
    GlobalCfg::operator=(rhs);
    return *this;
}


