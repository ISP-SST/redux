#include "redux/momfbd/momfbdjob.hpp"

#include "redux/momfbd/util.hpp"

#include "redux/logging/logger.hpp"
#include "redux/file/fileana.hpp"
#include "redux/util/bitoperations.hpp"
#include "redux/util/cache.hpp"
#include "redux/util/fileutil.hpp"
#include "redux/util/stringutil.hpp"

#include <thread>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
using boost::algorithm::iequals;

using namespace redux::logging;
using namespace redux::momfbd;
using namespace redux::momfbd::util;
using namespace redux::file;
using namespace redux::util;
using namespace redux;
using namespace std;


size_t MomfbdJob::jobType = Job::registerJob( "momfbd", momfbd::MomfbdJob::create );


int MomfbdJob::staticInit(void) {
    
    counts[ StepID(jobType,JSTEP_CHECKING) ]    = CountT(1,5);
    counts[ StepID(jobType,JSTEP_CHECKED) ]     = CountT();
    counts[ StepID(jobType,JSTEP_PREPROCESS) ]  = CountT(1,1);
    counts[ StepID(jobType,JSTEP_QUEUED) ]      = CountT(0,2);
    counts[ StepID(jobType,JSTEP_RUNNING) ]     = CountT(1,5);
    counts[ StepID(jobType,JSTEP_DONE) ]        = CountT();
    counts[ StepID(jobType,JSTEP_VERIFY) ]      = CountT(1,1);
    counts[ StepID(jobType,JSTEP_VERIFIED) ]    = CountT();
    counts[ StepID(jobType,JSTEP_POSTPROCESS) ] = CountT(1,1);
    counts[ StepID(jobType,JSTEP_WRITING) ]     = CountT(1,1);

    return 0;
}


MomfbdJob::MomfbdJob( void ) : waveFronts(*this), cfgChecked(false), dataChecked(false) {
    
    static int si RDX_UNUSED = staticInit();
    info.typeString = "momfbd";

}


MomfbdJob::~MomfbdJob( void ) {

    cleanup();
    Job::moveTo( this, Job::JSTEP_NONE );
    
}


uint64_t MomfbdJob::unpackParts( const char* ptr, WorkInProgress::Ptr wip, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count(0);
    if( wip->nParts ) {
        try {
            PatchData* tmpPD = new PatchData(*this);
            wip->parts.assign(1,Part::Ptr(tmpPD));        // if size > 1 the old "wip" has a globalData appended to it, we don't need it anymore.
            count += tmpPD->unpack( ptr+count, swap_endian );
            if( wip->nParts > 1 ) {
                globalData.reset( new GlobalData(*this) );
                count += globalData->unpack( ptr+count, swap_endian );
            }
        } catch( const std::exception& e ) {
            wip->parts.clear();
            throw;
        }
    }
    return count;
}


void MomfbdJob::parsePropertyTree( bpo::variables_map& vm, bpt::ptree& tree, redux::logging::Logger& logger ) {

    Job::parsePropertyTree( vm, tree, logger );
    
    // possibly override cfg-entries with command-line arguments
    if( vm.count( "simxy" ) ) tree.put( "SIM_XY", vm["simxy"].as<string>() );
    if( vm.count( "simx" ) ) tree.put( "SIM_X", vm["simx"].as<string>() );
    if( vm.count( "simy" ) ) tree.put( "SIM_Y", vm["simy"].as<string>() );
    if( vm.count( "imgn" ) ) tree.put( "IMAGE_NUM", vm["imgn"].as<string>() );
    if( vm.count( "output-file" ) ) tree.put( "OUTPUT_FILES", vm["output-file"].as<string>() );
    if( vm.count( "init" ) ) tree.put( "INIT_FILES", vm["init"].as<string>() );
    if( vm.count( "force" ) ) tree.put( "OVERWRITE", true );
    if( vm.count( "no-swap" ) ) tree.put( "NOSWAP", true );
    if( vm.count( "trace" ) ) tree.put( "TRACE", true );
    if( vm.count( "no-trace" ) ) tree.erase( "TRACE" );
    if( vm.count( "old-ns" ) ) tree.put( "OLD_NS", true );

    GlobalCfg::parseProperties(tree, logger);
    uint16_t nObj(0);
    for( auto & property : tree ) {
        if( iequals( property.first, "OBJECT" ) ) {
            Object* tmpObj = new Object( *this, nObj );
            tmpObj->parsePropertyTree( property.second, logger );
            if( nObj < outputFiles.size() ) {
                tmpObj->outputFileName = outputFiles[nObj];
            }
            if( nObj < initFiles.size() ) {
                tmpObj->initFile = initFiles[nObj];
            } else if( initFiles.size() ) {
                tmpObj->initFile = "OUTPUT";
            }
            if( iequals( tmpObj->initFile, "OUTPUT" ) ) tmpObj->initFile = tmpObj->outputFileName;
            objects.push_back( shared_ptr<Object>( tmpObj ) );
            nObj++;
        }
    }

    for( uint16_t i=0; i<nObj; ++i ) {
        if( objects[i] ) {
            if( i < outputFiles.size() ) {
                objects[i]->outputFileName = outputFiles[i];
            }
            if( i < initFiles.size() ) {
                objects[i]->initFile = initFiles[i];
            } else if( initFiles.size() ) {
                objects[i]->initFile = "OUTPUT";
            }
            if( iequals( objects[i]->initFile, "OUTPUT" ) ) objects[i]->initFile = objects[i]->outputFileName;
        }
    }

    if( outputFiles.size() > objects.size() ) {
        LOG_WARN << outputFiles.size() << " output file names specified but only " << objects.size() << " objects found." << ende;
    }
    
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
    sz += roi.size();
    sz += 2*sizeof(uint16_t) + 2;           // objects.size() + trace_objects.size() + cfgChecked & dataChecked
    for( const auto& obj: objects ) {
        sz += obj->size();
    }
    for( const auto& tobj: trace_objects ) {
        sz += tobj->size();
    }
    return sz;
}


uint64_t MomfbdJob::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = Job::pack( ptr );
    count += GlobalCfg::pack( ptr+count );
    count += roi.pack( ptr+count );
    count += pack( ptr+count, cfgChecked );
    count += pack( ptr+count, dataChecked );
    count += pack( ptr+count, (uint16_t)objects.size() );
    for( const auto& obj: objects ) {
        count += obj->pack( ptr+count );
    }
    count += pack( ptr+count, (uint16_t)trace_objects.size() );
    for( const auto& tobj: trace_objects ) {
        count += tobj->pack( ptr+count );
    }
    
    return count;
    
}


uint64_t MomfbdJob::unpack( const char* ptr, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count = Job::unpack( ptr, swap_endian );
    count += GlobalCfg::unpack( ptr+count, swap_endian );
    count += roi.unpack( ptr+count, swap_endian);
    count += unpack( ptr+count, cfgChecked );
    count += unpack( ptr+count, dataChecked );
    uint16_t tmp;
    count += unpack( ptr+count, tmp, swap_endian );
    objects.resize( tmp );
    for( auto& obj: objects ) {
        obj.reset(new Object(*this));
        count += obj->unpack( ptr+count, swap_endian );
    }
    count += unpack( ptr+count, tmp, swap_endian );
    trace_objects.resize( tmp );
    for( auto& tobj: trace_objects ) {
        tobj.reset(new Object(*this));
        count += tobj->unpack( ptr+count, swap_endian );
    }
    return count;
}


size_t MomfbdJob::nImages(void) const {
    size_t nTotalImages(0);
    for( const auto& obj: objects ) {
        nTotalImages += obj->nImages();
    }
    return nTotalImages;
}


void MomfbdJob::checkParts( void ) {

    uint16_t mask = 0;

    auto lock = getLock( true );
    if( !lock.owns_lock() ) return;

    for( auto & patch : patches ) {
        if( (patch->step == JSTEP_RUNNING) && (patch.use_count() == 1) ) {  // no reference in any existing WIP, so it is in limbo.
            patch->step = JSTEP_ERR;
        }
        if( patch->step & JSTEP_ERR ) {
            //if( patch->nRetries++ < info.maxPartRetries ) {
                patch->step = JSTEP_QUEUED;
            //}
        }
        mask |= patch->step;
    }
    lock.unlock();

    if( mask & JSTEP_ERR ) {
        moveTo( this, JSTEP_ERR );
        info.state.store( JSTATE_ERR );
        progWatch.clear();
        updateProgressString();
    } else if( countBits( mask ) == 1 ) {  // if all parts are "done", set the whole job to done.
         if(mask == JSTEP_DONE) moveTo( this, mask );
    }
    
}


bool MomfbdJob::getWork( WorkInProgress::Ptr wip, uint16_t nThreads, const map<Job::StepID,Job::CountT>& nActive ) {

    //uint16_t partMask = checkParts();
    
    bool ret(false);
    uint16_t step = info.step.load();

    if( (step == JSTEP_COMPLETED) || (step == JSTEP_ERR) ) {
        return false;
    }

    if( wip->isRemote ) {
        
        if( (step != JSTEP_QUEUED) && (step != JSTEP_RUNNING) ) {
            return false;
        }

        auto lock = getLock(false);
        if( info.step == JSTEP_QUEUED ) {                            // preprocessing ready -> start
            auto glock = Job::getGlobalLock();
            const CountT& limits = counts[StepID(jobType,JSTEP_RUNNING)];
            if( limits.active >= limits.max ) return ret;
            glock.unlock();
            moveTo( this, JSTEP_RUNNING );
            updateProgressString();
            progWatch.clear();
            progWatch.set( patches.nElements() );
            progWatch.setTicker( std::bind( &MomfbdJob::updateProgressString, this) );
            progWatch.setHandler([this](){
                moveTo( this, JSTEP_DONE );
                updateProgressString();
            });
            startLog();
        }
        if( info.step == JSTEP_RUNNING ) {                      // running
            for( auto & patch : patches ) {
                if( patch && (patch->step == JSTEP_QUEUED) ) {
                    patch->step = JSTEP_RUNNING;
                    patch->load();
                    wip->parts.push_back( patch );
                    auto pJob = wip->previousJob.lock();
                    if( pJob.get() != this ) {     // First time for this slave -> include global data
                        wip->parts.push_back( globalData );
                    }
                    ret = true;
                    break;// only 1 part at a time for MomfbdJob
                }
            }
        }
        wip->nParts = wip->parts.size();
        return ret;
    } else {    // local processing, i.e. on master.
        uint16_t nextStep = getNextStep(step);
        switch( step ) {
            case JSTEP_NONE:                            // Minor jobs
            case JSTEP_SUBMIT:   break;
            case JSTEP_CHECKED:  ret = true; break;     // These steps are heavy/async, return true so the worker goes to run()
            case JSTEP_DONE:     ret = true; break;
            case JSTEP_VERIFIED: ret = true; break;
            case JSTEP_CHECKING:        // When these steps are active, the manager should just wait.
            case JSTEP_PREPROCESS: 
            case JSTEP_QUEUED:          // The job will go from queued -> running when first patch is requested by a slave
            case JSTEP_RUNNING: 
            case JSTEP_VERIFY: 
            case JSTEP_POSTPROCESS: 
            case JSTEP_WRITING: 
            default:
                check();
                return false;
        }
        
        if( ret && !nThreads ) {    // don't start heavy/async jobs if nThreads=0.
            return false;
        }

        auto lock = getGlobalLock();
        const CountT& limits = counts[StepID(jobType,nextStep)];
        if( limits.active >= limits.max ) {   // Not allowed to start until some jobs are done with this step
            return false;
        } else lock.unlock();
        
        if( step == JSTEP_SUBMIT ) {    // No check done yet. It will now be run asynchronously, so return from here to avoid a race with moveTo() below.
            startLog();
            check();
            return false;
        }

        if( step == JSTEP_CHECKED ) {
            // we are also restricted by nQueued.
            const CountT& limits2 = counts[StepID(jobType,JSTEP_QUEUED)];
            if( limits2.active >= limits2.max ) {
                return false;
            }
            startLog();
            info.startedTime = boost::posix_time::second_clock::local_time();
        }
       
        if( (step == JSTEP_WRITING) && !checkWriting() ) return false;

        moveTo( this, nextStep );
        updateProgressString();
        
    }
    
    return ret;
}


void MomfbdJob::ungetWork( WorkInProgress::Ptr wip ) {
    auto lock = getLock();
    for( auto& part : wip->parts ) {
        part->step = JSTEP_QUEUED;
    }
    
}


void MomfbdJob::failWork( WorkInProgress::Ptr wip ) {
    auto lock = getLock();
    for( auto& part : wip->parts ) {
        part->step = JSTEP_ERR;
    }
}


void MomfbdJob::returnResults( WorkInProgress::Ptr wip ) {
    
    if( !wip->isRemote ) return;        // NOTE the manager doesn't return any results, but it seems to deadlock job-starting for a while
    
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    boost::posix_time::time_duration elapsed = (now - wip->workStarted);
    vector<PatchData::Ptr> pV;
    
    {
        auto lock = getLock();
        for( auto& part : wip->parts ) {
            auto tmpPatch = static_pointer_cast<PatchData>( part );
            PatchData::Ptr patch = patches( tmpPatch->index.y, tmpPatch->index.x );
            patch->copyResults(*tmpPatch);         // copies the returned results without overwriting other variables.
            patch->step = JSTEP_POSTPROCESS;
            pV.push_back(patch);
        }
        
        info.maxProcessingTime = max<uint32_t>( info.maxProcessingTime, elapsed.total_seconds() );
        if( info.maxProcessingTime ) {
            uint32_t newTimeout = info.maxProcessingTime * (20.0 - 15.0*progWatch.progress());    // TBD: start with large margin, then shrink to ~5*maxProcessingTime ?
            LOG_TRACE << "returnResults(): Adjusting job-timeout from " << info.timeout << " to " << newTimeout << " seconds." << ende;
            info.timeout = newTimeout;    // TBD: fixed timeout or 10 * maxProcTime ?
        }
    }
    
    for( auto p : pV ) {
        auto lock = getLock();
        ++progWatch;
    }

}


void MomfbdJob::cleanup(void) {
    
    auto lock = getLock();
    
    progWatch.clear();
    solver.reset();

    for( auto &o: objects ) {
        o->cleanup();
    }

    objects.clear();
    patches.clear();
    globalData.reset();
    
}


bool MomfbdJob::mayBeDeleted(void) {
    uint16_t jobStep = info.step.load();
    return ( jobStep != JSTEP_PREPROCESS && jobStep != JSTEP_POSTPROCESS );     // these are asynchronous steps and job can't be safely deleted until it finishes
}



size_t MomfbdJob::memUsage(void) {
  
  return 1;
  
}


size_t MomfbdJob::diskUsage(void) {
  if( cachePath.empty() ) return 0;
  bfs::path cp = bfs::path(Cache::get().path()) / bfs::path(cachePath);
  return getDirSize( cp.string() );
  
}


size_t MomfbdJob::procUsage(void) {
  
  return 3;
  
}
        

bool MomfbdJob::run( WorkInProgress::Ptr wip, boost::asio::io_service& service, uint16_t maxThreads ) {
    

    uint16_t jobStep = info.step.load();
    uint16_t patchStep = 0;
    uint16_t nThreads = std::min( maxThreads, info.maxThreads );
    
    logger.setContext( "job "+to_string(info.id) );
    
    if(wip->parts.size() && wip->parts[0]) patchStep = wip->parts[0]->step;

    if( jobStep == JSTEP_PREPROCESS ) {
        progWatch.set(0,0);
        progWatch.setTicker( std::bind( &MomfbdJob::updateProgressString, this) );
        preProcess(service, nThreads);                           // preprocess on master: load, flatfield, split in patches
    } else if( jobStep == JSTEP_RUNNING || jobStep == JSTEP_QUEUED ) {
        if( patchStep == JSTEP_POSTPROCESS ) {      // patch-wise post-processing, do we need any for momfbd?  the filtering etc. is done on the slaves.
        } else {                                    // main processing
            if( !globalData ) {

                LOG_ERR << "MomfbdJob::run()  Generating globalData, this should have been received from the master!!" << ende;
                globalData.reset( new GlobalData(*this) );
            }

            if( !solver ) solver.reset( new Solver(*this, service, nThreads) );
            for( auto& part : wip->parts ) {      // momfbd jobs will only get 1 part at a time, this is just to keep things generic.
                logger.setContext( "job "+to_string(info.id)+":"+to_string(part->id) );
                // Run main processing
                auto data = static_pointer_cast<PatchData>(part);
                StopWatch timer;
                data->initPatch();
                solver->run( data );
                generateTraceData( data );
                part->nThreads = nThreads;
                part->runtime_wall = timer.getSeconds();
                part->runtime_cpu = timer.getCPUSeconds();
                //solver->run_new(static_pointer_cast<PatchData>(part));
                wip->hasResults = true;
            }
        }
    } else if( jobStep == JSTEP_VERIFY ) {
        progWatch.set(0,0);
        progWatch.setTicker( std::bind( &MomfbdJob::updateProgressString, this) );
//        if( !solver ) solver.reset( new Solver(*this, service, nThreads) );            // Initialize, allocations, etc.
        service.post( std::bind( &MomfbdJob::verifyPatches, this) );                   // check for failed patches
    } else if( jobStep == JSTEP_POSTPROCESS ) {
//         progWatch.set( patches.nElements() );
//         progWatch.setHandler( std::bind( &MomfbdJob::postProcess, this, std::ref(service), nThreads) );
//         progWatch.setTicker( std::bind( &MomfbdJob::updateProgressString, this) );
//         if( !solver ) solver.reset( new Solver(*this, service, nThreads) );
//         service.post( std::bind( &MomfbdJob::loadPatchResults, this, std::ref(service), nThreads) );
// 
        writeOutput(service);
    } else {
        LOG << "MomfbdJob::run()  unrecognized step = " << ( int )info.step.load() << ende;
        moveTo( this, JSTEP_ERR );
        updateProgressString();
    }
    
    return false;
    
}


void MomfbdJob::setLogChannel(std::string channel) {
    Job::setLogChannel(channel);
    GlobalCfg::setLogChannel(channel);
    waveFronts.setLogChannel(channel);
    for( auto& obj : objects ) {
        obj->setLogChannel(channel);
        for (auto& ch : obj->channels) {
            ch->setLogChannel(channel);
        }
    }
}


bool MomfbdJob::checkPatchPositions(void) {

    shared_ptr<Object> ref = getObject(0);
    const vector<int16_t>& clip = ref->channels[0]->alignClip;
    Point16 imageSizes = ref->getImageSize();
    if( clip.size() == 4 ) {
        roi = Region16( min(clip[3],clip[2]), min(clip[1],clip[0]), max(clip[3],clip[2]), max(clip[1],clip[0]) );
    } else roi = Region16( imageSizes.y-1, imageSizes.x-1 );
    
    uint16_t halfPatchSize = patchSize/2;
    Region16 posLimits = roi+1;     // local, 1-based copy of ROI
    posLimits.shrink( halfPatchSize );
    posLimits -= roi.first;         // patch coordinates are given relative to the ROI
    
    bool xOutside(false);
    bool yOutside(false);
    if( subImagePosXY.empty() ) {
        uint16_t totalOverlap = minimumOverlap+patchSize/4;     // from MvN: always overlap 25% + 16 pixels.
        // TODO: do split per channel instead, to allow for different image-scales and/or hardware
        // NOTE:  subImagePosX/Y are kept 1-based, so offset by 1 during cut-out.
        LOG << boost::format("MomfbdJob::checkPatchPositions(): halfPatchSize=%d  overlap=%d") % halfPatchSize % totalOverlap << ende;
        if( subImagePosX.empty() ) { // x-coordinate of patch-centre
            subImagePosX = segment<uint16_t>( posLimits.first.x, posLimits.last.x, patchSize, totalOverlap );
            LOG << "MomfbdJob::checkPatchPositions(): Generated patch positions  " << printArray(subImagePosX,"X") << ende;
        }
        if( subImagePosY.empty() ) { // y-coordinate of patch-centre
            subImagePosY = segment<uint16_t>( posLimits.first.y, posLimits.last.y, patchSize, totalOverlap );
            LOG << "MomfbdJob::checkPatchPositions(): Generated patch positions  " << printArray(subImagePosY,"Y") << ende;
        }

        if( subImagePosX.empty() || subImagePosY.empty() ) {
            LOG_ERR << "MomfbdJob::checkPatchPositions(): No patches specified or generated, can't continue." << ende;
            moveTo( this, JSTEP_ERR );
            info.state.store( JSTATE_IDLE );
            return false;
        }

        uint16_t minP(-1),maxP(0);
        for( uint16_t pos : subImagePosX ) {
            if( pos < posLimits.first.x ) {
                pos = posLimits.first.x;
                xOutside = true;
            }
            if( pos > posLimits.last.x ) {
                pos = posLimits.last.x;
                xOutside = true;
            }
            if( pos < minP ) minP = pos;
            if( pos > maxP ) maxP = pos;
        }

        if( xOutside ) {
            string msg = "SIM_X is too close to the border, this will lead to strong artefacts.\n\t" + printArray(subImagePosX,"SIM_X");
            msg += "\n\tThe patches should lie in the rectangle: "+(string)posLimits;
            auto tmpX = segment<uint16_t>( minP, maxP, patchSize, totalOverlap );
            msg += "\n\tIf the cfg-file was generated by the crispred pipeline, it is recommended you re-run the prepmomfbd"
                   "\n\tcommand with extraclip large enough to force the patches inside.";
            msg += "\n\tTo autogenerate locations, you can add --simx (without locations) to the rdx_sub command.\n\t"
                   "You can also pass the locations as --simx=\"" + printArray(tmpX,"") + "\" to rdx_sub.";
            LOG_ERR << msg << ende;
        }
        minP = -1;
        maxP = 0;
        for( uint16_t pos : subImagePosY ) {
            if( pos < posLimits.first.y ) {
                pos = posLimits.first.y;
                yOutside = true;
            }
            if( pos > posLimits.last.y ) {
                pos = posLimits.last.y;
                yOutside = true;
            }
            if( pos < minP ) minP = pos;
            if( pos > maxP ) maxP = pos;
        }
        
        if( yOutside ) {
            string msg = "SIM_Y is too close to the border, this will lead to strong artefacts.\n\t" + printArray(subImagePosY,"SIM_Y");
            msg += "\n\tThe patches should lie in the rectangle: "+(string)posLimits;
            auto tmpY = segment<uint16_t>( minP, maxP, patchSize, totalOverlap );
            msg += "\n\tIf the cfg-file was generated by the crispred pipeline, it is recommended you re-run the prepmomfbd"
                   "\n\tcommand with extraclip large enough to force the patches inside.";
            msg += "\n\tTo autogenerate locations, you can add --simy (without locations) to the rdx_sub command.\n\t"
                   "You can also pass the locations as --simy=\"" + printArray(tmpY,"") + " to rdx_sub.";
            LOG_ERR << msg << ende;
        }

        if( xOutside || yOutside ) {
            return false;
        }
        
    } else {
        subImagePosX.clear();
        subImagePosY.clear();
        for ( size_t i=0; i< subImagePosXY.size(); ++i ) {
            if( i%2 ) subImagePosY.push_back( subImagePosXY[i] );
            else subImagePosX.push_back( subImagePosXY[i] );
        }
        // TODO check for the XY-version
    }
    
    return true;
    
}


void MomfbdJob::unloadCalib( boost::asio::io_service& service ) {

    for( shared_ptr<Object>& obj : objects ) {
        for( shared_ptr<Channel>& ch : obj->channels ) {
            ch->unloadCalib();
        }
    }

    if( runFlags & RF_FLATFIELD ) {
        LOG << "MomfbdJob #" << info.id << " ("  << info.name << ") flat-fielding completed." << ende;
        moveTo( this, JSTEP_COMPLETED );
    } else {
        LOG << "MomfbdJob #" << info.id << " ("  << info.name << ") pre-processed and queued:"
            << "  nObjects = " << objects.size() << "  nImages = " << nImages()
            << "  nPatches = " << patches.nElements() << ende;
        moveTo( this, JSTEP_QUEUED );
        stopLog();
        //info.progressString = "Q";
        //info.step.store( JSTEP_DONE );
        //check();
    }
    info.state.store( JSTATE_IDLE );
    updateProgressString();
    
}


void MomfbdJob::preProcess( boost::asio::io_service& service, uint16_t nThreads ) {

    LOG_TRACE << "MomfbdJob #" << info.id << " ("  << info.name << ") pre-processing..." << ende;

    generateTraceObjects();
    
    uint32_t nTotalImages = 0;
    uint32_t nTotalChannels(0);
    Point16 imageSizes;
    for( shared_ptr<Object>& obj : objects ) {
        Point16 tmp = obj->getImageSize();
        if( tmp == 0 ) {
            throw std::logic_error("An object returned imgSize=0, this probably means that the input data could not be read.");
        }
        nTotalImages += obj->nImages( true );
        nTotalChannels += obj->channels.size();
        if(imageSizes == 0) {
            imageSizes = tmp;
        } else if( tmp != imageSizes ) {    // TBD: allow for different patchsizes (i.e. pixelsize/ccd-size) for different objects/channels.
            throw std::logic_error("The clipped images have different sizes for the different objects, please verify the ALIGN_CLIP values.");
        }
    }
    //service.post( std::bind( &MomfbdJob::initCache, this) );    // TBD: should cache initialization be parallelized?

    progWatch.setTarget( nTotalImages );
    progWatch.setHandler( std::bind( &MomfbdJob::unloadCalib, this, std::ref(service)) );

    if( !(runFlags&RF_FLATFIELD) ) {    // Skip this if we are only doing flatfielding
        
        uint16_t halfPatchSize = patchSize/2;
        bool hasXY = !subImagePosXY.empty();

        unsigned int nPatchesX = subImagePosX.size();
        unsigned int nPatchesY = subImagePosY.size();
        if ( hasXY ) nPatchesY = 1;         // since the .momfbd format expects a rectangular grid of points, we store
                                            // the results as a 1 x nPatches array.
        
        progWatch.increaseTarget( 1 );

        uint64_t count(0);
        patches.resize( nPatchesY, nPatchesX );
        Point16 ps( patchSize, patchSize );
        bool use_swap = !(runFlags&RF_NOSWAP);
        if( use_swap ) {
            cachePath = to_string(Cache::pid()) +"_"+ to_string( info.id );
        }
        for( unsigned int y=0; y<nPatchesY; ++y ) {
            for( unsigned int x=0; x<nPatchesX; ++x ) {
                PatchData::Ptr patch( new PatchData(*this, y, x ) );
                if( use_swap ) patch->setPath(cachePath);
                patch->step = JSTEP_QUEUED;
                if( hasXY ) {
                    patch->position = Point16( subImagePosY[x]-1, subImagePosX[x]-1 );   // subImagePosX/Y is 1-based
                } else {
                    patch->position = Point16( subImagePosY[y]-1, subImagePosX[x]-1 );   // subImagePosX/Y is 1-based
                }
                patch->roi.first = patch->position - halfPatchSize;
                patch->roi.last = patch->roi.first+ps-1;
                patch->id = ++count;
                patches(y,x) = std::move(patch);
            }
        }   // end nPatchesY

        cachePath = cleanPath(Cache::get().path() + "/" + to_string(Cache::pid()) +"_"+ to_string( info.id )) + "/";   // this is used in Channel::loadData()
        bfs::path tmpP(cachePath);
        if( !bfs::exists(tmpP) ) {
            bfs::create_directories( tmpP );
        }

        if( use_swap ) {
            for( shared_ptr<Object>& obj: objects ) {
                obj->cacheFile = cachePath + "results_" + to_string(obj->ID);
            }
            for( shared_ptr<Object>& obj: trace_objects ) {
                obj->cacheFile = cachePath + "trace_" + to_string(obj->ID) + "_" + uIntsToString( obj->waveFrontList );
            }
        }
        service.post( [this](){
            initCache();
            ++progWatch;
        } );
        
    }       // end RF_FLATFIELD

    for( shared_ptr<Object>& obj : objects ) {
        obj->loadData( service, nThreads, patches );
    }
    
    waveFronts.loadInit( service, patches );
    
}

void MomfbdJob::initCache(void) {
    
    if( !globalData ) {
        globalData.reset(new GlobalData(*this));
    }
    globalData->constraints.init();
    
    for( shared_ptr<Object>& obj: objects ) {
        obj->initCache();
    }
    
}


void MomfbdJob::clearPatches(void) {

    auto lock = getLock();
    uint64_t nTotalThreads(0);
    double cpu_sum(0.0);
    for( const PatchData::Ptr& p: patches ) {
        nTotalThreads += p->nThreads;
        cpu_sum += p->runtime_cpu;
        p->clear();
    }
    size_t nPatches = patches.nElements();
    patches.clear();
    moveTo( this, JSTEP_COMPLETED );
    
    boost::posix_time::time_duration elapsed = getElapsed( JSTEP_CHECKING, JSTEP_CHECKED );
    boost::posix_time::time_duration total;
    string timings = "\nTiming details:\n";
    if( ! elapsed.is_not_a_date_time() ) {
        timings += "       Checking: " + to_simple_string(elapsed) +"\n";
        total += elapsed;
    }
    elapsed = getElapsed( JSTEP_PREPROCESS, JSTEP_QUEUED );
    if( ! elapsed.is_not_a_date_time() ) {
        timings += "  PreProcessing: " + to_simple_string(elapsed) +"\n";
        total += elapsed;
    }
    elapsed = getElapsed( JSTEP_RUNNING, JSTEP_DONE );
    if( ! elapsed.is_not_a_date_time() ) {
        timings += "        Running: " + to_simple_string(elapsed) +"\n";
        total += elapsed;
    }
    elapsed = getElapsed( JSTEP_WRITING, JSTEP_COMPLETED );
    if( ! elapsed.is_not_a_date_time() ) {
        timings += "        Writing: " + to_simple_string(elapsed) +"\n";
        total += elapsed;
    }
    if( ! total.is_not_a_date_time() ) {
        timings += "          TOTAL: " + to_simple_string(total) +"\n";
    }
    
    LOG << nPatches << " patches, processed in: " << to_simple_string(total)
        << ".  CPU-time: " << (cpu_sum/3600) << " hours." << ende;
    
    LOG_DETAIL << timings << ende;
    
}

 
void MomfbdJob::verifyPatches( void ) {

    moveTo( this, JSTEP_VERIFIED );
    return; // TODO fix patch verification & posTprocessing

    int nPatchesX = patches.dimSize(1);
    int nPatchesY = patches.dimSize(0);
    int nFailedPatches(0);

    size_t nModes = modeNumbers.size();
    size_t nAlpha = nModes*nImages();
    
    unique_ptr<double[]> bestFit( new double[ nAlpha ] );
    double* bPtr = bestFit.get();
    unique_ptr<double[]> tmp( new double[ nAlpha ] );
    //double* tmpPtr = tmp.get();
    unique_ptr<double[]> tmp2( new double[ nAlpha ] );
    //double* tmp2Ptr = tmp2.get();
static int bla(1);    
    for( int y=0; y<nPatchesY; ++y ) {
        for( int x=0; x<nPatchesX; ++x ) {
    cout << "verifyPatches  " << __LINE__ << "  x=" << x << "  y=" << y << endl;
            PatchData::Ptr patch = patches(y,x);
patch->dump( "ver"+to_string(bla) );
            patch->initPatch();
            patch->loadAlpha( solver->alpha.get() );
cout << printArray( solver->alpha.get(), 20, "alphaA" ) << endl;

            solver->shiftAndInit( true );
            solver->applyAlpha();
            double origMetric = solver->metric();
            
            memset( bPtr, 0, nAlpha*sizeof(double) );
            solver->shiftAndInit( bPtr, true );
            solver->applyAlpha( bPtr );
            double bestMetric = solver->metric();
            /*int lastY = min(nPatchesY-1,y+1);
            int lastX = min(nPatchesX-1,x+1);
            for( int yy=max(0,y-1); yy<=lastY; ++yy ) {
                for( int xx=max(0,x-1); xx<=lastX; ++xx ) {
                    patches( yy, xx )->loadAlpha( tmpPtr );
                    for( unsigned int m=0; m<nModes; ++m ) {
                        memcpy( tmp2Ptr, bPtr, nAlpha*sizeof(double) );
                        for( unsigned int a=m; a<nAlpha; a+=nModes ) {
                            tmp2Ptr[a] = tmpPtr[a];
                        }
                        solver->shiftAndInit( tmp2Ptr, true );
                        solver->applyAlpha( tmp2Ptr );
                        double thisMetric = solver->metric();
                        if( thisMetric < bestMetric ){
cout << __LINE__ << "  Patch"<< patch->index << " mode: " << m << "  improved by (" << yy << "," << xx << ")" << endl;
                            memcpy( bPtr, tmp2Ptr, nAlpha*sizeof(double) );
                            bestMetric = thisMetric;
                        }
                    }
                }
            }*/
            
            if( bla /*bestMetric < origMetric*/ ) {
cout << __LINE__ << "  Patch"<< patch->index << " metric improved: " << origMetric<< " -> " << bestMetric << endl;
cout << printArray( solver->alpha.get(), 20, "alphaB" ) << endl;
                patch->storeAlpha( solver->alpha.get() );
                patch->step = JSTEP_QUEUED;
                ++nFailedPatches;
                bla = 0;
            }

        }   // end nPatchesX
    }   // end nPatchesY
    
    if( nFailedPatches > 0 ) {
        LOG << "verifyPatches: " << nFailedPatches << " patches failed, returning them to queue with new initial values." << ende;
        progWatch.set( patches.nElements() );
        progWatch.setHandler([this](){
            moveTo( this, JSTEP_DONE );
            updateProgressString();
        });
        progWatch.step( patches.nElements()-nFailedPatches );
        moveTo( this, JSTEP_RUNNING );
    } else {                                // All ok, proceed to post-processing step.
        LOG_DETAIL << "verifyPatches: all ok." << ende;
        moveTo( this, JSTEP_POSTPROCESS );
    }
    
    updateProgressString();
    
}


void MomfbdJob::writeOutput( boost::asio::io_service& service ) {
    
    moveTo( this, JSTEP_WRITING );
    
    progWatch.clear();
    progWatch.set( objects.size()+trace_objects.size() );
    progWatch.setTicker( std::bind( &MomfbdJob::updateProgressString, this) );
    progWatch.setHandler( std::bind( &MomfbdJob::clearPatches, this) );
    
    updateProgressString();
    
    if( saveMask&SF_SAVE_METRIC ) {
        //progWatch.increaseTarget();
        //service.post( [this](){
            bfs::path fn = bfs::path(info.outputDir) / bfs::path(info.name+"_metrics.f0");
            vector<size_t> dims = patches.dimensions();
            Array<float> metrics( dims );
            size_t maxIterations(0);
            for( auto& patch: patches ) {
                metrics( patch->index.y, patch->index.x ) = patch->finalMetric;
                maxIterations = max( maxIterations, patch->metrics.size() );
            }
            LOG << "Writing final metrics to file: " << fn << ende;
            Ana::write( fn.string(), metrics );
            if( maxIterations ) {
                fn = bfs::path(info.outputDir) / bfs::path(info.name+"_metric_progress.f0");
                dims.push_back(maxIterations);
                metrics.resize(dims);
                metrics.zero();
                for( auto& patch: patches ) {
                    copy( patch->metrics.begin(), patch->metrics.end(), metrics.ptr( patch->index.y, patch->index.x, 0 ) );
                }
                LOG << "Writing metric progression to file: " << fn << ende;
                Ana::write( fn.string(), metrics );
            }
        //    ++progWatch;
        //});

    }
    
    for( auto obj : objects ) {
        service.post( [this,obj,&service](){
            obj->writeResults( patches );
        });
    }
    for( auto tobj : trace_objects ) {
        service.post( [this,tobj,&service](){
            tobj->writeResults( patches );
        });
    }

}


void MomfbdJob::loadPatchResults( boost::asio::io_service& service, uint16_t nThreads ) {

    for( auto& patch: patches ) {
        service.post( [this,patch](){
            ++progWatch;
        });
    }
    
}


int MomfbdJob::getReferenceObject( void ) {
    
    set<uint32_t> allWaveFronts;
    for( shared_ptr<Object>& obj: objects ) {
        if( obj ) allWaveFronts.insert( obj->waveFrontList.begin(), obj->waveFrontList.end() );
    }
    
    int ret = -1;
    vector<shared_ptr<Object>> tmpObjs( objects );
    tmpObjs.erase( std::remove_if(tmpObjs.begin(), tmpObjs.end(),
                                        [&allWaveFronts]( const shared_ptr<Object>& o ){
                                            if( !o ) return true;
                                            set<uint32_t> objWaveFronts( o->waveFrontList.begin(), o->waveFrontList.end() ); 
                                            for( auto wf: allWaveFronts ) {
                                                if( objWaveFronts.count(wf) == 0 ) return true;
                                            }
                                            return false;
                                    }), tmpObjs.end() );
    
    if( tmpObjs.size() >= 1 ) {     // if only one object contains all wave-fronts, this should be the desired reference
        ret = tmpObjs[0]->ID;
        if( tmpObjs.size() > 1 ) {
            LOG_WARN << "More than one object contained all wavefronts." << ende;
        }
        LOG_DETAIL << "Selecting object #" << ret << " as reference." << ende;
    } else if( tmpObjs.size() > 1 ) {
        LOG_DEBUG << "No suitable object found to use as reference." << ende;
    }
    
    return ret;

}


void MomfbdJob::generateTraceObjects( void ) {

    if( !trace ) return;
    
    int refID = -1;
    set<int> idSet;

    for( shared_ptr<Object>& ref_obj: objects ) {
        idSet.insert(ref_obj->ID);
        if ( ref_obj && ref_obj->traceObject ) {
            refID = ref_obj->ID;
            break;
        }
    }
    
    if( refID < 0 ) refID = getReferenceObject();
    if( refID < 0 ) return;

    shared_ptr<Object> ref_obj = objects[refID];
    size_t found =  ref_obj->outputFileName.find_first_of("_.");
    string refTag;
    if( found != string::npos ) {
        refTag = ref_obj->outputFileName.substr( 0, found );
    }
    set<uint32_t> refWaveFronts( ref_obj->waveFrontList.begin(), ref_obj->waveFrontList.end() );
    vector<string> refFiles;
    for( auto& c: ref_obj->channels ) {
        c->getFileNames(refFiles);
    }
    set<string> refFileSet;
    for( auto& fn: refFiles ) {
        bfs::path resolved = bfs::canonical(fn);
        refFileSet.insert( resolved.string() );
    }
    
    map< set<uint32_t>, vector<string> > wfSets;        // TODO better generation of output filenames for trace-objects?
                                                        // For now, we are using the same names as the old-style trace objects.
    // convert old-style trace-objects
    for( shared_ptr<Object>& obj: objects ) {
        if ( obj && (obj->ID != refID) && (obj->weight == 0) ) {
            set<uint32_t> thisWf( obj->waveFrontList.begin(), obj->waveFrontList.end() );
            for( auto& wf: thisWf ) {
                if( refWaveFronts.count(wf) == 0 ) {    // not present in reference -> can't be traced.
                    goto loopend;
                }
            }
            vector<string> objFiles;
            for( auto& c: obj->channels ) {
                c->getFileNames(objFiles);
            }
            set<string> objFileSet;
            for( auto& fn: objFiles ) {
                bfs::path resolved = bfs::canonical(fn);
                objFileSet.insert( resolved.string() );
            }
            for( auto& fn: objFileSet ) {
                if( refFileSet.count(fn) == 0 ) {    // not present in reference -> can't be traced.
                    goto loopend;
                }
            }
            auto ret = wfSets.emplace( thisWf, vector<string>(1,obj->outputFileName) );
            if( !ret.second ) {      // not inserted, so this WF-set already existed
                ret.first->second.push_back( obj->outputFileName );
            }
            obj->traceID = refID;
            obj->channels.clear();
            trace_objects.push_back( std::move(obj) );
        }
loopend: ;
    }

    objects.erase( std::remove_if(objects.begin(), objects.end(),
                                        []( const shared_ptr<Object>& o ){
                                            if( !o ) return true;
                                            return false;
                                    }), objects.end() );

    size_t nOld = trace_objects.size();

    uint16_t idCnt(0);
    while( idSet.count(idCnt) ) idCnt++;
    for( shared_ptr<Object>& obj: objects ) {
        if ( obj && (obj->ID != refID) ) {
            set<uint32_t> thisWf( obj->waveFrontList.begin(), obj->waveFrontList.end() );
            auto ret = wfSets.emplace( thisWf, vector<string>(1,obj->outputFileName) );
            if( ret.second ) {      // inserted, so this WF combination has not been generated yet.
                shared_ptr<Object> tmpObj( new Object( *ref_obj, idCnt, refID ) );
                tmpObj->outputFileName = obj->outputFileName;
                tmpObj->waveFrontList = obj->waveFrontList;
                tmpObj->nObjectImages = tmpObj->waveFrontList.size();
                trace_objects.push_back( tmpObj );
            } else {                // already in list, append outputFileName
                ret.first->second.push_back(obj->outputFileName);
            }
        }
    }
    size_t nNew = trace_objects.size() - nOld;

    if( !refTag.empty() ) {
        for( shared_ptr<Object>& obj: trace_objects ) {
            size_t found =  obj->outputFileName.find_first_of("_.");
            string objTag;
            if( found != string::npos ) {
                objTag = obj->outputFileName.substr( 0, found );
            }
            if( !objTag.empty() && (objTag != refTag) ) {
                obj->outputFileName = replace_n( obj->outputFileName, objTag, refTag );
            }
        }
    }
    
    if( nNew || nOld ) {
        LOG_DETAIL << "Generated trace-objects for reference #" << refID << ende;
        if( nNew ) LOG_DETAIL << "Generated " << nNew << " trace-objects for reference object " << refID << ende;
        if( nOld ) LOG_DETAIL << "Converted " << nOld << " old-style trace objects." << ende;
    }

    
}


void MomfbdJob::generateTraceData( PatchData::Ptr patch ) {
    
    for( shared_ptr<Object>& tobj: trace_objects ) {
        shared_ptr<Object> o = getObject( tobj->traceID );
        shared_ptr<ObjectData> od = make_shared<ObjectData>();
        od->myObject = o;
        if( runFlags & RF_FIT_PLANE ) {
            tobj->fittedPlane.resize( patchSize, patchSize );
            o->fitAvgPlane( tobj->fittedPlane, tobj->waveFrontList );
        }
        o->restorePatch( *od, tobj->waveFrontList );
        patch->trace_data.push_back( od );
    }
    
}


void MomfbdJob::postProcess( boost::asio::io_service& service, uint16_t nThreads ) {

    LOG_TRACE << "MomfbdJob::postProcess()" << ende;
 

    progWatch.set( 1 );
    progWatch.setHandler( std::bind( &MomfbdJob::writeOutput, this, std::ref(service)) );

//     int nPatchesX = patches.dimSize(1);
//     int nPatchesY = patches.dimSize(0);
//     
//     if( trace_objects.size() ) {
//         //ProgressWatch pw;
//         for( int y=0; y<nPatchesY; ++y ) {
//             for( int x=0; x<nPatchesX; ++x ) {
//                 PatchData::Ptr patch = patches(y,x);
//                 patch->cacheLoad();
//                 patch->load();
//                 patch->initPatch();
//                 patch->loadAlpha( solver->alpha.get() );
//                 solver->shiftAndInit( true );
//                 solver->applyAlpha();
//                 //pw.set( tmpObjects.size() );
//                 for( shared_ptr<Object>& tobj: trace_objects ) {
//                 //    service.post([this,&pw,&objIdMap,patch,tobj](){
//                         shared_ptr<ObjectData> tmpOD( make_shared<ObjectData>(tobj) );
//                         tmpOD->channels = patch->objects[ objIdMap[tobj->ID] ]->channels;
//                         objects[ objIdMap[tobj->ID] ]->fitAvgPlane(tobj->fittedPlane, tobj->waveFrontList );
//                         tobj->restorePatch( *tmpOD );
//                         patch->objects.push_back( tmpOD );
//                  //       ++pw;
//                 //    });
//                 }
//                 //pw.wait();
//             }
//         }
//     }
    ++progWatch;

    
}


void MomfbdJob::updateProgressString(void) {

    memset( info.progressString, 0, RDX_JOB_PROGSTRING_LENGTH );
    float prog = (100.0*progWatch.progress());
    if( !isnormal(prog) ) prog = 0;
    switch( info.step ) {
        case JSTEP_CHECKING:    snprintf( info.progressString, RDX_JOB_PROGSTRING_LENGTH, "Checking" ); break;
        case JSTEP_PREPROCESS:  snprintf( info.progressString, RDX_JOB_PROGSTRING_LENGTH, "Pre: (%03.1f%%)", prog ); break;
        case JSTEP_VERIFY:      snprintf( info.progressString, RDX_JOB_PROGSTRING_LENGTH, "Verifying: (%03.1f%%)", prog ); break;
        case JSTEP_POSTPROCESS: ;
        case JSTEP_WRITING:     snprintf( info.progressString, RDX_JOB_PROGSTRING_LENGTH, "Writing: (%03.1f%%)", prog ); break;
        case JSTEP_QUEUED:      snprintf( info.progressString, RDX_JOB_PROGSTRING_LENGTH, "Q" ); break;
        case JSTEP_RUNNING:     snprintf( info.progressString, RDX_JOB_PROGSTRING_LENGTH, "(%03.1f%%)", prog ); break;
        case JSTEP_COMPLETED:   snprintf( info.progressString, RDX_JOB_PROGSTRING_LENGTH, "(completed)" ); break;
        case JSTEP_ERR:         snprintf( info.progressString, RDX_JOB_PROGSTRING_LENGTH, "(failed)" ); break;
        default: ;
    }


}


bool MomfbdJob::active(void) {
    switch (info.step) {
        //case JSTEP_PREPROCESS: ;
        case JSTEP_QUEUED:  return true;
        //case JSTEP_RUNNING: return true;
        //case JSTEP_POSTPROCESS: return true;
        default: return false;
    }
}


bool MomfbdJob::check(void) {
  
    bool ret(false);

    switch (info.step) {
        case JSTEP_NONE: {   // rdx_sub will check with a default-value (from constructor) of JSTEP_NONE
            moveTo( this, JSTEP_CHECKING );
            updateProgressString();
            ret =  (cfgChecked || checkCfg());
            ret &= (dataChecked || checkData(false));
            if(ret) {
                moveTo( this, JSTEP_CHECKED );
                updateProgressString();
            }
            stopLog();
            break;
        }
        case JSTEP_SUBMIT: {    // When checking is delegated to the manager.
            moveTo( this, JSTEP_CHECKING );
            updateProgressString();
            std::thread([this]() {
                auto lock = getLock();
                bool ret =  (cfgChecked || checkCfg());
                ret &= (dataChecked || checkData(false));
                if(ret) {
                    moveTo( this, JSTEP_CHECKED );
                    updateProgressString();
                } else {
                    moveTo( this, JSTEP_ERR );
                    updateProgressString();
                }
                stopLog();
            }).detach();
            break;
        }
        case JSTEP_RUNNING: {    // this check should find orphan parts etc.
            std::thread( std::bind(&MomfbdJob::checkParts,this) ).detach();
            break;
        }
        case JSTEP_CHECKING:
        case JSTEP_CHECKED:                  // no (manager-)checks at these steps, just fall through and return false
        case JSTEP_PREPROCESS:
        case JSTEP_QUEUED:
        case JSTEP_DONE:
        case JSTEP_VERIFY:
        case JSTEP_VERIFIED:
        case JSTEP_POSTPROCESS:
        case JSTEP_WRITING: break; 
        default: LOG_ERR << "MomfbdJob::check(): No check defined for step = " << (int)info.step << ende;
    }
    
    return ret;
    
}


bool MomfbdJob::checkCacheUsage(void) {
 
    if( runFlags&RF_FLATFIELD ) {    // No cache usage if we are only doing flatfielding
        return true;
    }
    
    size_t estimatedNeededCache(0);
    size_t nPatches = subImagePosXY.size();
    if( !nPatches ) {
        nPatches = subImagePosX.size()*subImagePosY.size();
    }
    for( shared_ptr<Object>& o: objects ) {
        size_t tmp = o->patchSize + 2*o->maxLocalShift;
        estimatedNeededCache += tmp*tmp*o->nImages();
    }
   
    estimatedNeededCache *= nPatches;
    
    bfs::path cachePath( Cache::get().path() );
    boost::system::error_code ec;
    bfs::space_info si = bfs::space(cachePath,ec);
    if( ec ) {
        LOG_WARN << "Failed to check available space for path" << cachePath << ": " << ec.message()
                 << "\n\tProceeding anyway!" << ende;
    } else {
        double diskFree = static_cast<double>(si.available)/si.capacity;
        double diskFreeGB = static_cast<double>(si.available)/(1<<30);
        double neededGB = static_cast<double>(estimatedNeededCache)/(1<<30);
        if( si.available < estimatedNeededCache ) {
            LOG_ERR << "Only " << (int)diskFreeGB << "Gb (" << (int)(diskFree*100)
                     << "%) free space on the cache drive (" << cachePath << ")."
                     << " This job needs ~" << (int)neededGB << "Gb free!!." << ende;
            return false;
        }
    }
    
    return true;
}


bool MomfbdJob::checkOutputUsage(void) {
    size_t estimatedOutputSize(0);
    for( shared_ptr<Object>& o: objects ) {
        estimatedOutputSize += o->estimateOutputSize();
    }
    return true;
}


bool MomfbdJob::checkCfg(void) {

    if( cfgChecked ) return true;

    if( patchSize%4 ) {
        LOG_ERR << "The patch-size has to be a multiple of 4." << ende;
        return false;
    }
    
    if( (runFlags&RF_FLATFIELD) && (runFlags&RF_CALIBRATE) ) {
        LOG_ERR << "Both FLATFIELD and CALIBRATE mode requested" << ende;
        return false;
    }

    if( subImagePosXY.size() && (subImagePosXY.size()%2) ) {
        LOG_ERR << "SIM_XY has to have an even number of entries, since it's  (x,y)-pairs." << ende;
        return false;
    }

    if( !checkPatchPositions() ) return false;

    if( subImagePosXY.empty() && (subImagePosX.empty() || subImagePosY.empty()) ) {
        LOG_ERR << "Patch X and/or Y positions are not specified (and autogeneration failed)." << ende;
        return false;
    }
    
    if( objects.empty() ) {
        LOG_ERR << "The configuration file contains no objects." << ende;
        return false;
    }
    
    for( size_t i=0; i<objects.size(); ++i ) {
        if( objects[i] ) {
            if( objects[i]->channels.empty() ) {
                LOG_ERR << "Object " << i << " has no channel specified." << ende;
                return false;
            }
            if( !objects[i]->nImages() ) {
                LOG_ERR << "Object " << i << " has no images." << ende;
                return false;
            }
            if( !objects[i]->checkCfg() ) return false;
        } else {
            LOG_ERR << "Object #" << i << " is invalid." << ende;
            return false;
        }
    }

    boost::system::error_code ec;
    bfs::path outDir( info.outputDir );
    
    if( !outDir.empty() && !bfs::exists(outDir) ){
        if( !bfs::create_directories(outDir,ec) ){
            LOG_FATAL << boost::format( "output directory %s not writable: %s" ) % outDir % ec.message() << ende;
            return false;
        } else LOG_TRACE << boost::format( "created output directory %s" ) % outDir << ende;
    }

    bfs::space_info si = bfs::space(outDir,ec);
    if( ec ) {
        LOG_WARN << "Failed to check available space for path" << outDir << ": " << ec.message()
                 << "\n\tProceeding anyway!" << ende;
    } else {
        double diskFree = static_cast<double>(si.available)/si.capacity;
        double diskFreeGB = static_cast<double>(si.available)/(1<<30);
        if (false) if( diskFree < 0.05 || diskFreeGB < 100 ) {
            LOG_WARN << "Only " << (int)diskFreeGB << "Gb (" << (int)(diskFree*100)
                     << "%) free space for the output data (" << outDir << ")."
                     << "\n\tYour jobs will fail if you run out of space!!" << ende;
            return false;
        }
    }

    cfgChecked = true;
    return true;

}


bool MomfbdJob::checkData(bool verbose) {
    
    if( dataChecked ) return true;

    set<uint32_t> wfs;
    for( shared_ptr<Object>& obj: objects ) {
        if( !obj->checkData(verbose) ) return false;
        wfs.insert( obj->waveFrontList.begin(), obj->waveFrontList.end() );
    }
    
    waveFrontList.resize( wfs.size() );
    std::copy( wfs.begin(), wfs.end(), waveFrontList.begin() );
    
    dataChecked = true;
    return true;
}
        
        
bool MomfbdJob::checkPre(void) {
    
    if( !globalData || !globalData->verify() ) return false;

    for( auto& obj: objects ) {
        if( !obj->progWatch.verify() ) {
            return false;
        }
        for( auto& ch : obj->channels ) {
            if( !ch->progWatch.verify() ) return false;
        }
    }
    
    return true;
    
}


bool MomfbdJob::checkPost(void) {
    
    return true;
    
}


bool MomfbdJob::checkWriting(void) {
    
    if( patches.nElements() || !progWatch.verify() ) {
        return false;
    }
    
    return true;
}


uint16_t MomfbdJob::getNextStep( uint16_t step ) const {
    
    if( step == JSTEP_NONE ) step = info.step;
    switch( step ) {
        case JSTEP_NONE:
        case JSTEP_SUBMIT:      return JSTEP_CHECKING;
        case JSTEP_CHECKING:    return JSTEP_CHECKED;
        case JSTEP_CHECKED:     return JSTEP_PREPROCESS;
        case JSTEP_PREPROCESS:  return JSTEP_QUEUED;
        case JSTEP_QUEUED:      return JSTEP_RUNNING;
        case JSTEP_RUNNING:     return JSTEP_DONE;
        case JSTEP_DONE:        return JSTEP_VERIFY;
        case JSTEP_VERIFY:      return JSTEP_VERIFIED;
        case JSTEP_VERIFIED:    return JSTEP_POSTPROCESS;
        case JSTEP_POSTPROCESS: return JSTEP_POSTPROCESS;
        case JSTEP_WRITING:     return JSTEP_COMPLETED;
        default: return Job::getNextStep( step );
    }
    
}


const shared_ptr<Object> MomfbdJob::getObject( uint16_t id ) const {
    for( auto& o: objects ) {
        if( o && (o->ID == id) ) return o;
    }
    for( auto& o: trace_objects ) {
        if( o && (o->ID == id) ) return o;
    }
    return shared_ptr<Object>();
}


const vector<shared_ptr<Channel>>& MomfbdJob::getChannels( uint16_t objID ) const {
    shared_ptr<Object> o = getObject(objID);
    if( o && (o->traceID >= 0) ) o = getObject( o->traceID );       // for trace-objects we return the channel-list for the reference object.
    if( !o ) throw out_of_range("invalid object-id: " + to_string(objID) );
    return o->getChannels();
}


const MomfbdJob& MomfbdJob::operator=(const GlobalCfg& rhs) {
    GlobalCfg::operator=(rhs);
    return *this;
}


