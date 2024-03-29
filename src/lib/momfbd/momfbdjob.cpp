#include "redux/momfbd/momfbdjob.hpp"

#include "redux/momfbd/util.hpp"

#ifdef DEBUG_
#   define TRACE_THREADS
#endif

#include "redux/logging/logger.hpp"
#include "redux/file/fileana.hpp"
#include "redux/util/bitoperations.hpp"
#include "redux/util/cache.hpp"
#include "redux/util/fileutil.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/util/trace.hpp"

#include <thread>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/info_parser.hpp>
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
    THREAD_MARK;
    MomfbdJob::cleanup();
    Job::moveTo( this, Job::JSTEP_NONE );
    THREAD_MARK
}


uint64_t MomfbdJob::packParts( char* ptr, WorkInProgress::Ptr wip  ) const {
    
    using redux::util::pack;
    if( !wip ) throw job_error( info.name + ": Can not pack parts without a valid WIP instance."  );
    if( wip->parts.empty() ) throw job_error( info.name + ": Can not pack an empty list of parts."  );
    THREAD_MARK
    wip->parts.resize(1);
    THREAD_MARK
    uint64_t count = wip->parts[0]->pack( ptr );
    THREAD_MARK
    if( wip->jobID != info.id ) {           // First part from this job, so we need to send the global info
        THREAD_MARK
        if( globalData ) {
            THREAD_MARK
            count += globalData->pack( ptr+count );
            THREAD_MARK
            wip->parts.push_back( globalData );
            THREAD_MARK
        } else throw job_error( info.name + ": No globalData defined."  );
    }
    wip->nParts = wip->parts.size();
    THREAD_MARK
    return count;
    
}


uint64_t MomfbdJob::unpackParts( const char* ptr, WorkInProgress::Ptr wip, bool swap_endian ) {
    
    using redux::util::unpack;
    if( !wip ) throw job_error( info.name + ": Can not unpack parts without a valid WIP instance."  );

    wip->parts.clear();
    uint64_t count(0);
    if( wip->nParts ) {
        PatchData* tmpPD = new PatchData(*this);
        count += tmpPD->unpack( ptr+count, swap_endian );
        if( wip->nParts > 1 ) {
            globalData.reset( new GlobalData(*this) );
            count += globalData->unpack( ptr+count, swap_endian );
        }
        wip->parts.push_back( Part::Ptr(tmpPD) );
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


bpt::ptree MomfbdJob::getPropertyTree( bpt::ptree* root, bool showAll ) {

    bpt::ptree tree = Job::getPropertyTree( nullptr, showAll );         // get Job-properties

    GlobalCfg::getProperties(tree);

    for( auto& obj: objects ) {
        obj->getPropertyTree( tree, showAll );
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
    
    if( packed.packedSize ) {
        memcpy( ptr, packed.data.get(), packed.packedSize );
        return packed.packedSize;
    }
    
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


void MomfbdJob::prePack( bool force ) {
    
    if( packed.packedSize && !force ) {
        return;
    }
    packed.size = size();
    if( !packed.size ) {
        return;
    }
    packed.data = rdx_get_shared<char>(packed.size);
    packed.packedSize = pack( packed.data.get() );
    
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

    THREAD_MARK
    auto lock = getLock( true );
    if( !lock.owns_lock() ) {
        THREAD_UNMARK
        return;
    }
    THREAD_MARK
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
    THREAD_MARK
    if( mask & JSTEP_ERR ) {
        moveTo( this, JSTEP_ERR );
        info.state.store( JSTATE_ERR );
        progWatch.clear();
        updateProgressString();
    } else if( countBits( mask ) == 1 ) {  // if all parts are "done", set the whole job to done.
         if(mask == JSTEP_DONE) moveTo( this, mask );
    }
    THREAD_UNMARK
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

        auto lock = getLock(true);
        if( lock && (info.step == JSTEP_QUEUED) ) {                            // preprocessing ready -> start
            THREAD_MARK
            auto glock = Job::getGlobalLock();
            const CountT& limits = counts[StepID(jobType,JSTEP_RUNNING)];
            if( limits.active >= limits.max ) return ret;
            glock.unlock();
            THREAD_MARK
            moveTo( this, JSTEP_RUNNING );
            updateProgressString();
            progWatch.clear();
            progWatch.set( patches.nElements() );
            progWatch.setTicker( std::bind( &MomfbdJob::updateProgressString, this) );
            progWatch.setHandler([this](){
                THREAD_MARK
                moveTo( this, JSTEP_DONE );
                updateProgressString();
                THREAD_UNMARK
            });
            startLog();
            THREAD_MARK
            prePack();
            THREAD_MARK
        }
        if( lock && (info.step == JSTEP_RUNNING) ) {                      // running
            THREAD_MARK
            wip->resetParts();
            for( auto & patch : patches ) {
                if( patch && (patch->step == JSTEP_QUEUED) ) {
                    patch->step = JSTEP_RUNNING;
                    THREAD_MARK
                    wip->parts.push_back( patch );
                    if( globalData && (wip->jobID != info.id) ) {     // First time for this slave -> include global data
                        wip->parts.push_back( globalData );
                    }
                    ret = true;
                    break;// only 1 part at a time for MomfbdJob
                }
            }
            THREAD_MARK
        }
        wip->nParts = wip->parts.size();
        return ret;
    } else {    // local processing, i.e. on master.
        THREAD_MARK
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
        THREAD_MARK
        if( ret && !nThreads ) {    // don't start heavy/async jobs if nThreads=0.
            return false;
        }
        THREAD_MARK
        auto lock = getGlobalLock();
        const CountT& limits = counts[StepID(jobType,nextStep)];
        if( limits.active >= limits.max ) {   // Not allowed to start until some jobs are done with this step
            return false;
        } else lock.unlock();
        THREAD_MARK
        if( step == JSTEP_SUBMIT ) {    // No check done yet. It will now be run asynchronously, so return from here to avoid a race with moveTo() below.
            startLog();
            check();
            return false;
        }
        THREAD_MARK
        if( step == JSTEP_CHECKED ) {
            // we are also restricted by nQueued.
            const CountT& limits2 = counts[StepID(jobType,JSTEP_QUEUED)];
            if( limits2.active >= limits2.max ) {
                return false;
            }
            startLog();
            info.startedTime = boost::posix_time::second_clock::universal_time();
        }
        THREAD_MARK
        if( (step == JSTEP_WRITING) && !checkWriting() ) return false;
        THREAD_MARK
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
    
    boost::posix_time::ptime now = boost::posix_time::second_clock::universal_time();
    boost::posix_time::time_duration elapsed = (now - wip->workStarted);
    
    {
        THREAD_MARK
        for( auto& part : wip->parts ) {
            auto tmpPatch = static_pointer_cast<PatchData>( part );
            PatchData::Ptr patch = patches( tmpPatch->index.y, tmpPatch->index.x );
            patch->copyResults(*tmpPatch);         // copies the returned results without overwriting other variables.
            patch->step = JSTEP_POSTPROCESS;
            ++progWatch;
        }
        
        info.maxProcessingTime = max<uint32_t>( info.maxProcessingTime, elapsed.total_seconds() );
        if( info.maxProcessingTime ) {
            uint32_t newTimeout = info.maxProcessingTime * (20.0 - 15.0*progWatch.progress());    // TBD: start with large margin, then shrink to ~5*maxProcessingTime ?
            LOG_TRACE << "returnResults(): Adjusting job-timeout from " << info.timeout << " to " << newTimeout << " seconds." << ende;
            info.timeout = newTimeout;    // TBD: fixed timeout or 10 * maxProcTime ?
        }
    }

    THREAD_MARK

}


void MomfbdJob::cleanup(void) {
    
    THREAD_MARK;
    clearPatches();
    
    THREAD_MARK;
    auto lock = getLock();
    
    pool.interrupt_all();
    workLoop.reset();
    ioService.stop();
    pool.join_all();
    
    progWatch.clear();
    solver.reset();

    THREAD_MARK
    for( auto &o: objects ) {
        o->cleanup();
    }

    objects.clear();
    globalData.reset();
    
    THREAD_MARK
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
        

bool MomfbdJob::run( WorkInProgress::Ptr wip, uint16_t maxThreads ) {
    
    THREAD_MARK
    uint16_t jobStep = info.step.load();
    uint16_t patchStep = 0;
    uint16_t nThreads = std::min( maxThreads, info.maxThreads );
    
    if( !workLoop ) {
        ioService.reset();
        workLoop.reset( new boost::asio::io_service::work(ioService) );
        addThread(maxThreads);
    }

    logger.setContext( "job "+to_string(info.id) );
    
    if(wip->parts.size() && wip->parts[0]) patchStep = wip->parts[0]->step;

    if( jobStep == JSTEP_PREPROCESS ) {
        THREAD_MARK
        progWatch.set(0,0);
        progWatch.setTicker( std::bind( &MomfbdJob::updateProgressString, this) );
        preProcess();                           // preprocess on master: load, flatfield, split in patches
        THREAD_MARK
    } else if( jobStep == JSTEP_RUNNING || jobStep == JSTEP_QUEUED ) {
        if( patchStep == JSTEP_POSTPROCESS ) {      // patch-wise post-processing, do we need any for momfbd?  the filtering etc. is done on the slaves.
        } else {                                    // main processing
            if( !globalData ) {

                LOG_ERR << "MomfbdJob::run()  Generating globalData, this should have been received from the master!!" << ende;
                globalData.reset( new GlobalData(*this) );
            }

            if( !solver ) solver.reset( new Solver(*this, ioService, maxThreads) );
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
        THREAD_MARK
        verifyPatches();
        updateProgressString();
    } else if( jobStep == JSTEP_POSTPROCESS ) {
        THREAD_MARK
        writeOutput();
    } else {
        LOG << "MomfbdJob::run()  unrecognized step = " << ( int )info.step.load() << ende;
        moveTo( this, JSTEP_ERR );
        updateProgressString();
    }
    
    THREAD_UNMARK
    return false;
    
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
    
    Point roi_span = roi.span();
    if( roi_span.x <= patchSize ) {     // only 1 patch possible 
        uint16_t pos =  roi.first.x + halfPatchSize + 1;
        posLimits.first.x = posLimits.last.x = pos;
    }
    if( roi_span.y <= patchSize ) {     // only 1 patch possible 
        uint16_t pos =  roi.first.y + halfPatchSize + 1;
        posLimits.first.y = posLimits.last.y = pos;
    }

    if( subImagePosXY.empty() ) {
        
        bool xOutside(false);
        bool yOutside(false);
        
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


void MomfbdJob::unloadCalib( void ) {

    THREAD_MARK
    for( shared_ptr<Object>& obj : objects ) {
        for( shared_ptr<Channel>& ch : obj->channels ) {
            ch->unloadCalib();
        }
    }
    THREAD_MARK
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
    THREAD_UNMARK
    
}


void MomfbdJob::preProcess( void ) {

    LOG_TRACE << "MomfbdJob #" << info.id << " ("  << info.name << ") pre-processing..." << ende;

    THREAD_MARK;
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
    //ioService.post( std::bind( &MomfbdJob::initCache, this) );    // TBD: should cache initialization be parallelized?

    progWatch.setTarget( nTotalImages );
    progWatch.setHandler( std::bind( &MomfbdJob::unloadCalib, this ) );

    THREAD_MARK
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

        ioService.post( [this](){
            THREAD_MARK
            initCache();
            THREAD_MARK
            ++progWatch;
            THREAD_UNMARK
        } );
        
    }       // end RF_FLATFIELD

    THREAD_MARK
    for( shared_ptr<Object>& obj : objects ) {
        obj->loadData( ioService, patches );
    }
    
    THREAD_MARK
    waveFronts.loadInit( ioService, patches );
    
}

void MomfbdJob::initCache(void) {
    
    THREAD_MARK;
    if( !globalData ) {
        globalData.reset(new GlobalData(*this));
    }
    
    THREAD_MARK;
    for( shared_ptr<Object>& obj: objects ) {
        if( obj ) {
            obj->initObject();
            nModes = obj->modes->dimSize(0);     // FIXME: Now we are assuming all objects have the same number of modes, and that #2/3 are tilts, this should be checked!!
        }
    }
    globalData->constraints.init();
    
    THREAD_MARK;
    generateTraceObjects();

    if( !(runFlags&RF_NOSWAP) ) {
        for( shared_ptr<Object>& obj: objects ) {
            obj->cacheFile = cachePath + "results_" + to_string(obj->ID);
        }
        for( shared_ptr<Object>& obj: trace_objects ) {
            obj->cacheFile = cachePath + "trace_" + to_string(obj->ID) + "_" + uIntsToString( obj->waveFrontList );
        }
    }
    
    THREAD_MARK;
    globalData->prePack();

}


void MomfbdJob::clearPatches(void) {

    THREAD_MARK
    auto lock = getLock();
    
    THREAD_MARK;
    for( const PatchData::Ptr& p: patches ) {
        THREAD_MARK;
        if( p ) p->clear();
    }
    
    THREAD_MARK;
    patches.clear();

    THREAD_MARK;
    
}

 
void MomfbdJob::verifyPatches( void ) {

    THREAD_MARK
    moveTo( this, JSTEP_VERIFIED );
    return; // TODO fix patch verification & posTprocessing

    int nPatchesX = patches.dimSize(1);
    int nPatchesY = patches.dimSize(0);
    int nFailedPatches(0);

    size_t nAlpha = nModes*nImages();
    
    shared_ptr<double> bestFit = rdx_get_shared<double>(nAlpha);
    double* bPtr = bestFit.get();
    shared_ptr<double> tmp = rdx_get_shared<double>(nAlpha);
    //double* tmpPtr = tmp.get();
    shared_ptr<double> tmp2 = rdx_get_shared<double>(nAlpha);
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
            THREAD_MARK
            moveTo( this, JSTEP_DONE );
            updateProgressString();
            THREAD_MARK
        });
        progWatch.step( patches.nElements()-nFailedPatches );
        moveTo( this, JSTEP_RUNNING );
    } else {                                // All ok, proceed to post-processing step.
        LOG_DETAIL << "verifyPatches: all ok." << ende;
        moveTo( this, JSTEP_POSTPROCESS );
    }
    
    updateProgressString();
    
}


void MomfbdJob::writeOutput( void ) {
    
    THREAD_MARK
    moveTo( this, JSTEP_WRITING );
    
    progWatch.clear();
    progWatch.set( objects.size()+trace_objects.size() );
    progWatch.setTicker( std::bind( &MomfbdJob::updateProgressString, this) );

    progWatch.setHandler( [this](){
        
        THREAD_MARK;
        size_t nPatches = patches.nElements();
        uint64_t nTotalThreads(0);
        double cpu_sum(0.0);
        for( const PatchData::Ptr& p: patches ) {
            nTotalThreads += p->nThreads;
            cpu_sum += p->runtime_cpu;
        }
        
        clearPatches();
        
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
        
        THREAD_MARK;

        moveTo( this, JSTEP_COMPLETED );
        updateProgressString();
        
    });
    
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
    
    THREAD_MARK
    for( auto obj : objects ) {
        ioService.post( [this,obj](){
            THREAD_MARK
            obj->writeResults( patches );
            THREAD_UNMARK
        });
    }
    THREAD_MARK
    for( auto tobj : trace_objects ) {
        ioService.post( [this,tobj](){
            THREAD_MARK
            tobj->writeResults( patches );
            THREAD_UNMARK
        });
    }
    THREAD_MARK

}


void MomfbdJob::loadPatchResults( void ) {

    for( auto& patch: patches ) {
        ioService.post( [this,patch](){
            THREAD_MARK
            ++progWatch;
            THREAD_UNMARK
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
    size_t last_slash =  ref_obj->outputFileName.find_last_of("/");     // so we don't detect "._" in dirnames
    size_t found =  ref_obj->outputFileName.find_first_of("_.",last_slash);
    string refTag;
    if( last_slash == string::npos ) {
        last_slash = 0;
    }
    if( found != string::npos ) {
        refTag = ref_obj->outputFileName.substr( last_slash, found );
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
            if( thisWf == refWaveFronts ) {    // matches the total reference -> no need to generate a subset
                continue;
            }
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
            last_slash =  obj->outputFileName.find_last_of("/");     // so we don't detect "._" in dirnames
            size_t separator_found =  obj->outputFileName.find_first_of("_.",last_slash);
            string objTag;
            if( last_slash == string::npos ) {
                last_slash = 0;
            }
            if( separator_found != string::npos ) {
                objTag = obj->outputFileName.substr( last_slash, separator_found );
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


void MomfbdJob::postProcess( void ) {

    LOG_TRACE << "MomfbdJob::postProcess()" << ende;
 

    progWatch.set( 1 );
    progWatch.setHandler( std::bind( &MomfbdJob::writeOutput, this) );

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


void MomfbdJob::dumpConfig(void) {

    if( cachePath.empty() ) {
        return;
    }
    string cfgName = cachePath + "cfg";
    
    try {
        ofstream cfgOut( cfgName, ofstream::trunc );
        bpt::ptree cfgTree;
        getPropertyTree( &cfgTree );
        bpt::write_info( cfgOut, cfgTree );
    } catch ( exception& e ) {
        LOG_ERR << "Error when trying to write config \"" << cfgName << "\": " << e.what() << ende;
    }
}


void MomfbdJob::updateProgressString(void) {

    memset( info.progressString, 0, RDX_JOB_PROGSTRING_LENGTH );
    float prog = (100.0*progWatch.progress());
    if( !isnormal(prog) ) prog = 0;
    switch( static_cast<int>(info.step) ) {
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


void MomfbdJob::updateStatus(void) {

    dumpConfig();
  
}


bool MomfbdJob::active(void) {
    switch( static_cast<int>(info.step) ) {
        //case JSTEP_PREPROCESS: ;
        case JSTEP_QUEUED:  return true;
        //case JSTEP_RUNNING: return true;
        //case JSTEP_POSTPROCESS: return true;
        default: return false;
    }
}


bool MomfbdJob::check(void) {
  
    bool ret(false);
    THREAD_MARK;
    switch( static_cast<int>(info.step) ) {
        case JSTEP_NONE: {   // rdx_sub will check with a default-value (from constructor) of JSTEP_NONE
            moveTo( this, JSTEP_CHECKING );
            updateProgressString();
            ret = (cfgChecked || checkCfg());
            if(ret) ret &= (dataChecked || checkData(false));
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
                THREAD_MARK
                auto lock = getLock();
                bool all_ok =  (cfgChecked || checkCfg());
                all_ok &= (dataChecked || checkData(false));
                if( all_ok ) {
                    moveTo( this, JSTEP_CHECKED );
                    updateProgressString();
                } else {
                    moveTo( this, JSTEP_ERR );
                    updateProgressString();
                }
                stopLog();
//                THREAD_UNMARK;
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
    THREAD_MARK
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
            LOG_ERR << "Object #" << i << " is null." << ende;
            return false;
        }
    }
    
    if( modeList.size() < 2 ) {
        // FIXME: sanity check: either file(s) specified or MODES, same number of modes has to exist in all objects AND the tilts has to be the first 2 modes. 
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
        //if( nModes && nModes != obj->modes.dimSize(0) ) return false;     // FIXME: implement testing nModes
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


Point16 MomfbdJob::getSmallestImageSize( void ) {
    Point16 ret;
    ret = numeric_limits<uint16_t>::max();
    for( auto& o: objects ) {
        if( o ) ret = o->getImageSize( true ).min(ret);
    }
    return ret;
}


const MomfbdJob& MomfbdJob::operator=(const GlobalCfg& rhs) {
    GlobalCfg::operator=(rhs);
    return *this;
}


