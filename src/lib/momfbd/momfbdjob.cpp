#include "redux/momfbd/momfbdjob.hpp"

#include "redux/momfbd/util.hpp"

#include "redux/logging/logger.hpp"
#include "redux/file/fileana.hpp"
#include "redux/util/bitoperations.hpp"
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


atomic<int> MomfbdJob::nActivePre(0);
atomic<int> MomfbdJob::nActivePost(0);
std::map<uint16_t,uint16_t> MomfbdJob::maxActive = { {JSTEP_PREPROCESS,1},
                                                     {JSTEP_QUEUED,2},
                                                     {JSTEP_VERIFY,1},
                                                     {JSTEP_POSTPROCESS,1} };

MomfbdJob::MomfbdJob( void ) {
    
    info.typeString = "momfbd";

}


MomfbdJob::~MomfbdJob( void ) {

    cleanup();
    
}


// size_t MomfbdJob::getTypeID(void) {
//     return momfbd::MomfbdJobDummy;
// }
// 

uint64_t MomfbdJob::unpackParts( const char* ptr, WorkInProgress::Ptr wip, bool swap_endian ) {
    using redux::util::unpack;
    uint64_t count(0);
    if( wip->nParts ) {
        try {
            wip->parts.resize(1);        // if size > 1 the old "wip" has a globalData appended to it, we don't need it anymore.
            PatchData* tmpPD = new PatchData(*this);
            wip->parts[0].reset(tmpPD);
            count += tmpPD->unpack( ptr+count, swap_endian );
            if( wip->nParts > 1 ) {
                LOG_DEBUG << "Unpacking globalData for job #" << info.id << ende;
                globalData.reset( new GlobalData(*this) );
                count += globalData->unpack( ptr+count, swap_endian );
            }
        } catch( const std::exception& e ) {
            LOG_ERR << "Unpacking of parts failed: " << e.what() << ende;
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
    if( vm.count( "swap" ) ) tree.put( "SWAP", true );

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

    //outputFiles.clear();
    //initFiles.clear();
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


uint16_t MomfbdJob::checkParts( void ) {

    uint16_t mask = 0;
    map<uint16_t,uint16_t> counts;
    auto lock = getLock();
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
        counts[patch->step]++;
    }
    lock.unlock();

    if( mask & JSTEP_ERR ) {
        info.step.store( JSTEP_ERR );
        info.state.store( JSTATE_ERR );
        progWatch.clear();
        info.progressString = "error";
    } else if( countBits( mask ) == 1 ) {  // if all parts have the same "step", set the whole job to that step.
        uint16_t tmp = mask;
        if( !info.step.compare_exchange_strong( tmp, tmp ) ) {  // TODO make this neater
            info.step.store( mask );
        }

    }
    
    return mask;
    
}


bool MomfbdJob::getWork( WorkInProgress::Ptr wip, uint16_t nThreads, const map<uint16_t,uint16_t>& nActive ) {

    //uint16_t partMask = checkParts();
    
    bool ret(false);
    uint16_t step = info.step.load();
    uint16_t origStep = step;
    if( (step == JSTEP_COMPLETED) || (step == JSTEP_ERR) ) {
        return false;
    }

    if( wip->isRemote ) {
        auto lock = getLock( true );
        if( lock.owns_lock() ) {
            if( step == JSTEP_QUEUED ) {                       // preprocessing ready -> start
                progWatch.set(patches.nElements());
                progWatch.setHandler([this](){
                    info.step = JSTEP_DONE;
                    updateProgressString();
                });
                info.step = step = JSTEP_RUNNING;
                updateProgressString();
            }

            if( step == JSTEP_RUNNING ) {                      // running
                 for( auto & patch : patches ) {
                    if( patch && (patch->step == JSTEP_QUEUED) ) {
                        //LOG_DETAIL << "Starting patch: #" << patch->id << "   step=" << (int)patch->step << "  ptr = " << hexString(patch.get()) << ende;
                        patch->step = JSTEP_RUNNING;
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
        }
        wip->nParts = wip->parts.size();
        return ret;
        
    } else {    // local processing, i.e. on master.
        
        if( step == JSTEP_DONE ) {
            return true;
        }

        if( step == JSTEP_SUBMIT ) {
            startLog();
            if( checkData() ) {
                step = info.step = JSTEP_CHECKED;
            } else {
                step = info.step = JSTEP_ERR;
            }
            stopLog();
        }
        
//         if( (step == JSTEP_PREPROCESS) && checkPre() ) {
//             step = info.step = JSTEP_QUEUED;
//         }
//         
        if( (step == JSTEP_VERIFIED) && checkPost() ) {
            ret = true;
        }
        
        if( (step == JSTEP_WRITING) && checkWriting() ) {
            step = info.step = JSTEP_COMPLETED;
        }

        // check against maximum allowed active
        if( (maxActive.count(step) == 0) || (nActive.count(step) > 0) || (nActive.at(step) <= maxActive.at(step)) ) {
            
            if( step == JSTEP_CHECKED ) {
                size_t nPrepQueued(0);
                if( nActive.count(JSTEP_PREPROCESS) ) nPrepQueued += nActive.at(JSTEP_PREPROCESS);
                if( nPrepQueued >= maxActive[JSTEP_PREPROCESS] ) return false;
                if( nActive.count(JSTEP_QUEUED) ) nPrepQueued += nActive.at(JSTEP_QUEUED);
                if( nPrepQueued < maxActive[JSTEP_QUEUED] ) {
                    startLog();
                    info.step = JSTEP_PREPROCESS;
                    info.startedTime = boost::posix_time::second_clock::local_time();
                    ret = true;
                }
            }
        
        }
        
    }
    
    if( info.step != origStep ) {
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
        p->cacheStore(true);
        auto lock = getLock();
        ++progWatch;
    }

}


void MomfbdJob::cleanup(void) {
    
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
                data->initPatch();
                solver->run( data );
                
                //solver->run_new(static_pointer_cast<PatchData>(part));
                wip->hasResults = true;
            }
        }
    } else if( jobStep == JSTEP_DONE ) {
        info.step.store( JSTEP_VERIFY );
        progWatch.set(0,0);
        progWatch.setTicker( std::bind( &MomfbdJob::updateProgressString, this) );
        if( !solver ) solver.reset( new Solver(*this, service, nThreads) );            // Initialize, allocations, etc.
        service.post( std::bind( &MomfbdJob::verifyPatches, this) );                   // check for failed patches
    } else if( jobStep == JSTEP_VERIFIED ) {
        info.step.store( JSTEP_POSTPROCESS );
//         progWatch.set( patches.nElements() );
//         progWatch.setHandler( std::bind( &MomfbdJob::postProcess, this, std::ref(service), nThreads) );
//         progWatch.setTicker( std::bind( &MomfbdJob::updateProgressString, this) );
//         if( !solver ) solver.reset( new Solver(*this, service, nThreads) );
//         service.post( std::bind( &MomfbdJob::loadPatchResults, this, std::ref(service), nThreads) );
// 
        writeOutput(service);
    } else {
        LOG << "MomfbdJob::run()  unrecognized step = " << ( int )info.step.load() << ende;
        info.step.store( JSTEP_ERR );
        updateProgressString();
    }
    
    return false;
    
}


void MomfbdJob::setLogChannel(std::string channel) {
    Job::setLogChannel(channel);
    GlobalCfg::setLogChannel(channel);
    for( auto& obj : objects ) {
        obj->setLogChannel(channel);
        for (auto& ch : obj->channels) {
            ch->setLogChannel(channel);
        }
    }
}


bool MomfbdJob::checkPatchPositions(void) {
    
    const vector<int16_t>& clip = objects[0]->channels[0]->alignClip;
    Point16 imageSizes = objects[0]->getImageSize();
    if( clip.size() == 4 ) {
        roi = Region16( abs(clip[3]-clip[2]), abs(clip[1]-clip[0]) );
    } else roi = Region16( imageSizes.y-1, imageSizes.x-1 );
    
    uint16_t halfPatchSize = patchSize/2;
    Region16 posLimits = roi+1;     // local, 1-based copy of ROI
    posLimits.shrink( halfPatchSize );

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
            info.step.store( JSTEP_ERR );
            info.state.store( JSTATE_IDLE );
            return false;
        }

        uint16_t minP(-1),maxP(0);
        for( uint16_t pos : subImagePosX ) {
            if( pos < posLimits.first.x || pos > posLimits.last.x ) {
                pos = std::min<uint16_t>(std::max<uint16_t>(posLimits.first.x,pos),posLimits.last.x);
                xOutside = true;
            }
            if( pos < minP ) minP = pos;
            if( pos > maxP ) maxP = pos;
        }

        if( xOutside ) {
            string msg = "SIM_X is too close to the border, this will lead to strong artefacts.\n\t" + printArray(subImagePosX,"SIM_X");
            auto tmpX = segment<uint16_t>( minP, maxP, patchSize, totalOverlap );
            msg += "\n\tIf the cfg-file was generated by the crispred pipeline, it is recommended you re-run the prepmomfbd"
                   "\n\tcommand with extraclip large enough to force the patches inside.";
            msg += "\n\tTo autogenerate locations, you can add --simx (without locations) to the rdx_sub command.\n\t"
                   "You can also pass the locations as --simx=\"" + printArray(tmpX,"") + "\" to rdx_sub.";
            //msg += "\n\tTo FORCE the use of your current patch-positions, add the --no-check flag to the rdx_sub command. (NOT recommended)";
            LOG_ERR << msg << ende;
        }
        minP = -1;
        maxP = 0;
        for( uint16_t pos : subImagePosY ) {
            if( pos < posLimits.first.y || pos > posLimits.last.y ) {
                pos = std::min<uint16_t>(std::max<uint16_t>(posLimits.first.y,pos),posLimits.last.y);
                yOutside = true;
            }
            if( pos < minP ) minP = pos;
            if( pos > maxP ) maxP = pos;
        }
        
        if( yOutside ) {
            string msg = "SIM_Y is too close to the border, this will lead to strong artefacts.\n\t" + printArray(subImagePosY,"SIM_Y");
            auto tmpY = segment<uint16_t>( minP, maxP, patchSize, totalOverlap );
            msg += "\n\tIf the cfg-file was generated by the crispred pipeline, it is recommended you re-run the prepmomfbd"
                   "\n\tcommand with extraclip large enough to force the patches inside.";
            msg += "\n\tTo autogenerate locations, you can add --simy (without locations) to the rdx_sub command.\n\t"
                   "You can also pass the locations as --simy=\"" + printArray(tmpY,"") + " to rdx_sub.";
            //msg += "\n\tTo FORCE the use of your current patch-positions, add the --no-check flag to the rdx_sub command. (NOT recommended)";
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
        info.step.store( JSTEP_COMPLETED );
    } else {
        LOG << "MomfbdJob #" << info.id << " ("  << info.name << ") pre-processed and queued:"
            << "  nPatches = " << patches.nElements() << "  nObjects = " << objects.size() << ende;
        info.step.store( JSTEP_QUEUED );
        updateProgressString();
        //info.progressString = "Q";
        //info.step.store( JSTEP_DONE );
        //check();
    }
    info.state.store( JSTATE_IDLE );
    nActivePre--;
    
}


void MomfbdJob::preProcess( boost::asio::io_service& service, uint16_t nThreads ) {

    LOG_TRACE << "MomfbdJob #" << info.id << " ("  << info.name << ") pre-processing..." << ende;

    uint32_t nTotalImages(0);
    uint32_t nTotalChannels(0);
    Point16 imageSizes;
    for( shared_ptr<Object>& obj : objects ) {
        Point16 tmp = obj->getImageSize();
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

//    progWatch.setTicker([this](){
//        cout << "Tick: "  << progWatch.progressString() << endl;
//    });
    

    if( !(runFlags&RF_FLATFIELD) ) {    // Skip this if we are only doing flatfielding
        
        uint16_t halfPatchSize = patchSize/2;
        bool hasXY = !subImagePosXY.empty();

        unsigned int nPatchesX = subImagePosX.size();
        unsigned int nPatchesY = subImagePosY.size();
        if ( hasXY ) nPatchesY = 1;         // since the .momfbd format expects a rectangular grid of points, we store
                                            // the results as a 1 x nPatches array.
        
        progWatch.increaseTarget( nPatchesX*nPatchesY*nTotalChannels );

        uint64_t count(0);
        patches.resize( nPatchesY, nPatchesX );
        Point16 ps( patchSize, patchSize );
        cachePath = to_string(Cache::pid()) +"_"+ to_string( info.id );
        for( unsigned int y=0; y<nPatchesY; ++y ) {
            for( unsigned int x=0; x<nPatchesX; ++x ) {
                PatchData::Ptr patch( new PatchData(*this, y, x ) );
                patch->setPath(cachePath);
                patch->step = JSTEP_QUEUED;
                if( hasXY ) {
                    patch->position = Point16( subImagePosY[x]-1, subImagePosX[x]-1 );   // subImagePosX/Y is 1-based
                } else {
                    patch->position = Point16( subImagePosY[y]-1, subImagePosX[x]-1 );   // subImagePosX/Y is 1-based
                }
                patch->roi.first = patch->position - halfPatchSize;
                patch->roi.last = patch->roi.first+ps-1;
                patch->id = ++count;
//                 for( auto& obj: objects ) {
//                     for (auto& ch: obj->channels) {
//                         service.post( [this,patch,obj,ch](){
//                             auto chData = patch->objects[obj->ID]->channels[ch->ID];
//                             //ch->adjustCutout( *chData, patch->roi );
//                             //++progWatch;
//                         } );
//                         
//                     }
//                 }
                patches(y,x) = std::move(patch);
            }
        }   // end nPatchesY

        service.post( [this](){
            initCache();
        } );
        
    }       // end RF_FLATFIELD

    for( shared_ptr<Object>& obj : objects ) {
        obj->loadData( service, nThreads, patches );
    }
    
}

void MomfbdJob::initCache(void) {
    
    auto lock = getLock();

    if( !globalData ) {
        globalData.reset(new GlobalData(*this));
    }
    globalData->constraints.init();
//globalData->constraints.dump();
    for( auto& obj: objects ) {
        obj->initCache();
    }
    
}


void MomfbdJob::clearPatches(void) {

    nActivePost--;
    patches.clear();
    info.step = JSTEP_COMPLETED;

}

 
void MomfbdJob::verifyPatches( void ) {

info.step = JSTEP_VERIFIED;
return;

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
            patch->cacheLoad();
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
            info.step = JSTEP_DONE;
            updateProgressString();
        });
        progWatch.step( patches.nElements()-nFailedPatches );
        info.step = JSTEP_RUNNING;
    } else {                                // All ok, proceed to post-processing step.
        LOG_DETAIL << "verifyPatches: all ok." << ende;
        info.step = JSTEP_POSTPROCESS;
    }
    
    updateProgressString();
    
}


void MomfbdJob::writeOutput( boost::asio::io_service& service ) {
    
    auto jlock = getLock();
    info.step = JSTEP_WRITING;
    
    progWatch.clear();
    progWatch.set( objects.size() );
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
//         service.post( [this,obj,&service](){
//             obj->writeResults( service, patches );
//         });
        service.post( [this,obj,&service](){
            obj->writeResults( patches );
        });
    }

}


void MomfbdJob::loadPatchResults( boost::asio::io_service& service, uint16_t nThreads ) {

    for( auto& patch: patches ) {
        service.post( [this,patch](){
            patch->cacheLoad(false);
            ++progWatch;
        });
    }
    
}


void MomfbdJob::postProcess( boost::asio::io_service& service, uint16_t nThreads ) {

    LOG_DEBUG << "MomfbdJob::postProcess()" << ende;
 

    progWatch.set( patches.nElements() );
    progWatch.setHandler( std::bind( &MomfbdJob::writeOutput, this, std::ref(service)) );
//     progWatch.setHandler([this,&service](){      // triggered when all patches are loaded from the cache-files.
//         progWatch.set( objects.size() );
//         if( saveMask&SF_SAVE_METRIC ) {
//             progWatch.increaseTarget();
//             service.post( [this](){
//                 bfs::path fn = bfs::path(info.outputDir) / bfs::path(info.name+"_metrics.f0");
//                 Array<float> metrics( patches.dimensions() );
//                 for( auto& patch: patches ) {
//                     metrics( patch->index.y, patch->index.x ) = patch->finalMetric; 
//                 }
//                 LOG << "Writing metrics to file: " << fn << ende;
//                 Ana::write( fn.string(), metrics );
//                 ++progWatch;
//             });
// 
//         }
//     
//         for( auto obj : objects ) {
//             service.post( [this,obj](){
//                 obj->writeResults( patches );
//                 ++progWatch;
//             });
//         }
//         
//         progWatch.setHandler([this](){      // triggered when all object results (and possible the metric data) are written
//             patches.clear();
//             info.step = JSTEP_COMPLETED;
//         });
// 
// 
//     });

    for( auto& patch: patches ) {
        service.post( [this,patch](){
            patch->cacheLoad(false); //true);        // load and erase cache-file.
            ++progWatch;
        });
    }
    
}


void MomfbdJob::updateProgressString(void) {
    
    //cout << "updateProgressString0: " << endl;
    //static int blark = sleep(1);
    //cout << "updateProgressString1: " << hexString(this) << endl;
    //sleep(1);
    //cout << "updateProgressString2: " << progWatch.progressString() << endl;
    //sleep(1);
    
    //info.progressString = "-";
    //info.progressString = progWatch.dump();
    //return;
    switch( info.step ) {
        case JSTEP_PREPROCESS:  info.progressString = "Preprocessing"; break; //"(P:" + progWatch.progressString() + ")"; break;
        case JSTEP_VERIFY:      info.progressString = "V"; break;
        case JSTEP_POSTPROCESS: ;
        case JSTEP_WRITING:     info.progressString = "Writing"; break; // + progWatch.progressString() + ")"; break;
        case JSTEP_QUEUED:      info.progressString = "Q"; break;
        case JSTEP_RUNNING:     info.progressString = "(" + progWatch.progressString() + ")"; break;
        case JSTEP_COMPLETED:   info.progressString = "(completed)"; break;
        case JSTEP_ERR: info.progressString = "(failed)"; break;
        default: info.progressString = "";
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
    uint16_t step = info.step;
           bool skipCheck = ((info.flags&Job::CHECKED) || (info.flags&Job::NOCHECK));
    switch (step) {
        case JSTEP_NONE: {
             ret = skipCheck || (checkCfg() && checkData(false));// && checkCacheUsage();
            if(ret) info.step = JSTEP_SUBMIT;
            break;
        }
        case JSTEP_SUBMIT: {
            ret = skipCheck || checkData(true); // && checkCacheUsage();
            if(ret) info.step = JSTEP_CHECKED;
            break;
        }
        case JSTEP_RUNNING:     ret = true; checkParts(); break;    // this check should find orphan parts etc.
        case JSTEP_VERIFY:      return false;
//         case JSTEP_PREPROCESS:  return checkPre();
//         case JSTEP_WRITING:     return checkWriting();
//         case JSTEP_CHECKED: ;                  // no checks at these steps, just fall through and return true
//         case JSTEP_POSTPROCESS: ;
//         case JSTEP_QUEUED: ;
//         case JSTEP_COMPLETED: ret = true; break;
//         case JSTEP_ERR: ret = false; break;
        default: LOG_ERR << "check(): No check defined for step = " << (int)info.step << ende;
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

    return true;

}


bool MomfbdJob::checkData(bool verbose) {
    
    for( shared_ptr<Object>& obj: objects ) {
        if( !obj->checkData(verbose) ) return false;
    }
    
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
        
        
const MomfbdJob& MomfbdJob::operator=(const GlobalCfg& rhs) {
    GlobalCfg::operator=(rhs);
    return *this;
}


