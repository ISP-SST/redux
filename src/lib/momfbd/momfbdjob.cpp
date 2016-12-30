#include "redux/momfbd/momfbdjob.hpp"

#include "redux/momfbd/util.hpp"

#include "redux/logging/logger.hpp"
#include "redux/file/fileana.hpp"
#include "redux/util/bitoperations.hpp"
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


void MomfbdJob::parsePropertyTree( bpo::variables_map& vm, bpt::ptree& tree ) {

    Job::parsePropertyTree( vm, tree );
    //LOG_DEBUG << "MomfbdJob::parsePropertyTree()" << ende;
    
    // possibly override cfg-entries with command-line arguments
    if( vm.count( "simx" ) ) tree.put( "SIM_X", vm["simx"].as<string>() );
    if( vm.count( "simy" ) ) tree.put( "SIM_Y", vm["simy"].as<string>() );
    if( vm.count( "imgn" ) ) tree.put( "IMAGE_NUM", vm["imgn"].as<string>() );
    if( vm.count( "output-file" ) ) tree.put( "OUTPUT_FILES", vm["output-file"].as<string>() );
    if( vm.count( "init" ) ) tree.put( "INIT_FILES", vm["init"].as<string>() );
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
    //LOG_DEBUG << "MomfbdJob::parsePropertyTree() done." << ende;

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
cout << "mask = " << bitString(mask) << printArray(counts,"\ncounts") << endl; 
    if( mask & JSTEP_ERR ) {
        info.step.store( JSTEP_ERR );
        info.state.store( JSTATE_ERR );
        progWatch.clear();
        info.progressString = "blaha";
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
    if( (step == JSTEP_COMPLETED) || (step == JSTEP_ERR) ) {
        return false;
    }

    if( wip->isRemote ) {
        auto lock = getLock();
        if ( step == JSTEP_QUEUED ) {                       // preprocessing ready -> start
            progWatch.set(patches.nElements());
            progWatch.setHandler([this](){
                info.step = JSTEP_DONE;
                updateProgressString();
            });
            info.step = step = JSTEP_RUNNING;
            info.progress[0] = 0;
            info.progress[1] = patches.nElements();
            updateProgressString();
        }

        if ( step == JSTEP_RUNNING ) {                      // running
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

        wip->nParts = wip->parts.size();
        return ret;
        
    } else {    // local processing, i.e. on master.
        
        if( step == JSTEP_DONE ) {
            return true;
        }
            
        if( (step == JSTEP_SUBMIT) && checkData() ) {
            info.step = JSTEP_CHECKED;
        }
        
        if( (step == JSTEP_PREPROCESS) && checkPre() ) {
            info.step = JSTEP_QUEUED;
        }
        
        if( (step == JSTEP_POSTPROCESS) && checkPost() ) {
            ret = true;
        }
        
        if( (step == JSTEP_WRITING) && checkWriting() ) {
            info.step = JSTEP_COMPLETED;
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
                    info.progress[0] = 0;
                    info.startedTime = boost::posix_time::second_clock::local_time();
                    ret = true;
                }
            }
        
        }
        
    }
    
    if( info.step != step ) {
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
    
    auto lock = getLock();
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();

    boost::posix_time::time_duration elapsed = (now - wip->workStarted);

    for( auto& part : wip->parts ) {
        auto tmpPatch = static_pointer_cast<PatchData>( part );
        PatchData::Ptr patch = patches( tmpPatch->index.y, tmpPatch->index.x );
        patch->copyResults(*tmpPatch);         // copies the returned results without overwriting other variables.
        patch->step = JSTEP_POSTPROCESS;
        patch->cacheStore(true);
        ++progWatch;
        info.progress[0]++;
    }

    info.maxProcessingTime = max<uint32_t>( info.maxProcessingTime, elapsed.total_seconds() );
    if( info.maxProcessingTime ) {
        uint32_t newTimeout = info.maxProcessingTime * (20 - info.progress[0]*15.0/info.progress[1]);    // TBD: start with large margin, then shrink to ~5*maxProcessingTime ?
        LOG_TRACE << "returnResults(): Adjusting job-timeout from " << info.timeout << " to " << newTimeout << " seconds." << ende;
        info.timeout = newTimeout;    // TBD: fixed timeout or 10 * maxProcTime ?
    }
    
}


void MomfbdJob::cleanup(void) {
    
    progWatch.clear();
    for( auto &o: objects ) {
        o->cleanup();
    }

    objects.clear();
    patches.clear();
    globalData.reset();
    solver.reset();
    
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
            if( !solver ) {
                LOG_TRACE << "Initializing new Solver for job #" << info.id << ende;
                solver.reset( new Solver(*this, service, nThreads) );            // Initialize, allocations, etc.
                solver->init();
            }
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
        progWatch.set(0,0);
        progWatch.setTicker( std::bind( &MomfbdJob::updateProgressString, this) );
        verifyPatches( service, nThreads );                          // check for failed patches
    } else if( jobStep == JSTEP_POSTPROCESS ) {
        progWatch.set(0,0);
        progWatch.setTicker( std::bind( &MomfbdJob::updateProgressString, this) );
//         progWatch.setTicker([this](){
//             auto lock = getLock();
//             info.progressString = "(W:" + progWatch.progressString() + ")";
//         });
        postProcess(service, nThreads);                          // postprocess on master, collect results, save...
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


void MomfbdJob::unloadData( boost::asio::io_service& service ) {

    for( auto& obj : objects ) {
        for( auto& ch : obj->channels ) {
            service.post( bind(&Channel::unloadData, ch.get()) );
        }
    }

    if( runFlags & RF_FLATFIELD ) {
        LOG << "MomfbdJob #" << info.id << " ("  << info.name << ") flat-fielding completed." << ende;
        info.step.store( JSTEP_COMPLETED );
    } else {
        LOG << "MomfbdJob #" << info.id << " ("  << info.name << ") pre-processed and queued:"
            << "  nPatches = " << patches.nElements() << "  nObjects = " << objects.size() << ende;
        //info.step.store( JSTEP_QUEUED );
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
    for( auto& obj : objects ) {
        Point16 tmp = obj->getImageSize();
        nTotalImages += obj->nImages();
        nTotalChannels += obj->channels.size();
        if(imageSizes == 0) {
            imageSizes = tmp;
        } else if( tmp != imageSizes ) {    // TBD: allow for different patchsizes (i.e. pixelsize/ccd-size) for different objects/channels.
            throw std::logic_error("The clipped images have different sizes for the different objects, please verify the ALIGN_CLIP values.");
        }
    }
    //service.post( std::bind( &MomfbdJob::initCache, this) );    // TBD: should cache initialization be parallelized?
    if( nTotalChannels ) {
        const vector<int16_t>& clip = objects[0]->channels[0]->alignClip;
        if( clip.size() == 4 ) {
            roi = Region16( clip[2]-1, clip[0]-1, clip[3]-1, clip[1]-1 );
        }
    }

    progWatch.setTarget( nTotalImages );
    progWatch.setHandler( std::bind( &MomfbdJob::unloadData, this, std::ref(service)) );

//    progWatch.setTicker([this](){
//        cout << "Tick: "  << progWatch.progressString() << endl;
//    });
    

    if( !(runFlags&RF_FLATFIELD) ) {    // Skip this if we are only doing flatfielding
        
        uint16_t halfPatchSize = patchSize/2;
        uint16_t totalOverlap = minimumOverlap+patchSize/4;     // from MvN: always overlap 25% + 16 pixels.
        // TODO: do split per channel instead, to allow for different image-scales and/or hardware
        // NOTE:  subImagePosX/Y are kept 1-based, so ffset by 1 during cut-out.
        LOG << boost::format("MomfbdJob::preProcess(): halfPatchSize=%d  overlap=%d") % halfPatchSize % totalOverlap << ende;
        if( subImagePosX.empty() ) { // x-coordinate of patch-centre
            subImagePosX = segment<uint16_t>(halfPatchSize+1,imageSizes.x-halfPatchSize+1,patchSize,totalOverlap);
            LOG << "MomfbdJob::preProcess(): Generated patch positions  " << printArray(subImagePosX,"X") << ende;
        }
        if( subImagePosY.empty() ) { // y-coordinate of patch-centre
            subImagePosY = segment<uint16_t>(halfPatchSize+1,imageSizes.y-halfPatchSize+1,patchSize,totalOverlap);
            LOG << "MomfbdJob::preProcess(): Generated patch positions  " << printArray(subImagePosY,"Y") << ende;
        }
     
        if( subImagePosX.empty() || subImagePosY.empty() ) {
            LOG_ERR << "MomfbdJob::preProcess(): No patches specified or generated, can't continue." << ende;
            info.step.store( JSTEP_ERR );
            info.state.store( JSTATE_IDLE );
            return;
        }

        for( uint16_t& pos : subImagePosY ) {
            uint16_t adjustedPos = std::min<uint16_t>(std::max<uint16_t>(halfPatchSize+1,pos),imageSizes.y-halfPatchSize+1);       // stay inside borders
            if( adjustedPos != pos ) {
                LOG_WARN << "MomfbdJob::preProcess() y-position of patch is too close to the border, adjusting: " << pos << " -> " << adjustedPos << ende;
                pos = adjustedPos;
            }
        }

        for( uint16_t& pos : subImagePosX ) {
            uint16_t adjustedPos = std::min<uint16_t>(std::max<uint16_t>(halfPatchSize+1,pos),imageSizes.x-halfPatchSize+1);       // stay inside borders
            if( adjustedPos != pos ) {
                LOG_WARN << "MomfbdJob::preProcess() x-position of patch is too close to the border, adjusting: " << pos << " -> " << adjustedPos << ende;
                pos = adjustedPos;
            }
        }
       
        unsigned int nPatchesX = subImagePosX.size();
        unsigned int nPatchesY = subImagePosY.size();
        
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
                patch->position = Point16( subImagePosY[y]-1, subImagePosX[x]-1 );   // subImagePosX/Y is 1-based
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

    for( auto& obj : objects ) {
        obj->loadData( service, nThreads, patches );
    }
    
}

void MomfbdJob::initCache(void) {
    
    auto lock = getLock();

    if( !globalData ) {
        globalData.reset(new GlobalData(*this));
        globalData->constraints.init();
    }
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

 
void MomfbdJob::verifyPatches( boost::asio::io_service& service, uint16_t nThreads ) {

    //LOG << "Verifying patches for job #" << info.id << ende;
info.step = JSTEP_POSTPROCESS;
return;
    if( !solver ) {
        //LOG_TRACE << "Initializing new Solver for job #" << info.id << ende;
        solver.reset( new Solver(*this, service, nThreads) );            // Initialize, allocations, etc.
        solver->init();
    }

    unsigned int nPatchesX = patches.dimSize(1);
    unsigned int nPatchesY = patches.dimSize(0);

    cout << "verifyPatches  " << __LINE__ << "  nPatchesX=" << nPatchesX << "  nPatchesY=" << nPatchesY << endl;
    //progWatch.increaseTarget( nPatchesX*nPatchesY );

    for( unsigned int y=0; y<nPatchesY; ++y ) {
        for( unsigned int x=0; x<nPatchesX; ++x ) {
    cout << "verifyPatches  " << __LINE__ << "  x=" << x << "  y=" << y << endl;
            PatchData::Ptr patch = patches(y,x);
            patch->cacheLoad();
            patch->initPatch();
            solver->loadInit( patch, solver->alpha );
            //memset( solver->alpha, 0, solver->nParameters*sizeof(double) );
            cout << __LINE__ << printArray(solver->alpha, 7,"  alpha1") << endl;
            solver->shiftAndInit( solver->alpha, true );
            cout << __LINE__ << printArray(solver->alpha, 7,"  alpha2") << endl;
            solver->applyAlpha( solver->alpha );
            double mm = solver->metric2();
            //solver->dump( "vp_"+(string)patch->index );
            //LOG << "Patch " << patch->index << "   metric = " << patch->finalMetric << ende;
            cout << "Patch " << patch->index << "   metric = " << patch->finalMetric << "  mm= " << mm << endl;
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
//            patches(y,x) = std::move(patch);
        }   // end nPatchesX
    }   // end nPatchesY
    

    
/*    progWatch.clear();
    progWatch.set( objects.size() + (saveMask&SF_SAVE_METRIC) );

    progWatch.setHandler( std::bind( &MomfbdJob::clearPatches, this) );
    if( saveMask&SF_SAVE_METRIC ) {
        service.post( [this](){
            bfs::path fn = bfs::path(info.outputDir) / bfs::path(info.name+"_metrics.f0");
            Array<float> metrics( patches.dimensions() );
            for( auto& patch: patches ) {
                metrics( patch->index.y, patch->index.x ) = patch->finalMetric; 
            }
            LOG << "Writing metrics to file: " << fn << ende;
            Ana::write( fn.string(), metrics );
            ++progWatch;
        });

    }
    
    for( auto obj : objects ) {
        service.post( [this,obj](){
            obj->writeResults( patches );
            ++progWatch;
        });
    }
    
    
            for( auto& part : wip->parts ) {      // momfbd jobs will only get 1 part at a time, this is just to keep things generic.
                logger.setContext( "job "+to_string(info.id)+":"+to_string(part->id) );
                // Run main processing
                auto data = static_pointer_cast<PatchData>(part);
                data->initPatch();
                solver->run( data );
                //solver->run_new(static_pointer_cast<PatchData>(part));
                wip->hasResults = true;
            }

    
    
    
    
    
*/

    info.step = JSTEP_POSTPROCESS;
        
}


void MomfbdJob::writeOutput( boost::asio::io_service& service ) {
    
    
    auto jlock = getLock();
    
    unsigned int nPatchesX = patches.dimSize(1);
    unsigned int nPatchesY = patches.dimSize(0);
    int nFailedPatches(0);
    //progWatch.increaseTarget( nPatchesX*nPatchesY );
static int dummy(0);
    for( unsigned int y=0; y<nPatchesY; ++y ) {
        for( unsigned int x=0; x<nPatchesX; ++x ) {
            PatchData::Ptr patch = patches(y,x);
            patch->cacheLoad();
            if( dummy && y == 0 && x == 0 ) {
                patch->step = JSTEP_QUEUED;
                ++nFailedPatches;
                dummy = 0;
            }
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
        }
    }   // end nPatchesY

    if( nFailedPatches > 0 ) {
        LOG << "writeOutput: " << nFailedPatches << " patches failed, returning them to queue." << ende;
        progWatch.set( patches.nElements() );
//         progWatch.setTicker([this](){
//             info.progressString = "(" + progWatch.progressString() + ")";
//         });
        progWatch.setHandler([this](){
            info.step = JSTEP_POSTPROCESS;
        });
        progWatch.step( patches.nElements()-nFailedPatches );
        info.step = JSTEP_RUNNING;
        updateProgressString();
        return;
    }
    

    progWatch.clear();
    progWatch.set( objects.size() );

    progWatch.setHandler( std::bind( &MomfbdJob::clearPatches, this) );
    if( saveMask&SF_SAVE_METRIC ) {
        progWatch.increaseTarget();
        service.post( [this](){
            bfs::path fn = bfs::path(info.outputDir) / bfs::path(info.name+"_metrics.f0");
            Array<float> metrics( patches.dimensions() );
            for( auto& patch: patches ) {
                metrics( patch->index.y, patch->index.x ) = patch->finalMetric; 
            }
            LOG << "Writing metrics to file: " << fn << ende;
            Ana::write( fn.string(), metrics );
            ++progWatch;
        });

    }
    
    for( auto obj : objects ) {
        service.post( [this,obj,&service](){
            obj->writeResults( service, patches );
        });
    }
        
}


void MomfbdJob::postProcess( boost::asio::io_service& service, uint16_t nThreads ) {

    LOG_DEBUG << "MomfbdJob::postProcess()" << ende;
    info.progress[0] = 0;
    info.progress[1] = objects.size();

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

    info.step = JSTEP_WRITING;

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
        case JSTEP_PREPROCESS:  info.progressString = "P"; break; //"(P:" + progWatch.progressString() + ")"; break;
        case JSTEP_POSTPROCESS: ;
        case JSTEP_WRITING:     info.progressString = "(W:"; break;// + progWatch.progressString() + ")"; break;
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
    switch (step) {
        case JSTEP_NONE:        ret = checkCfg(); if(ret) info.step = JSTEP_SUBMIT; break;
        case JSTEP_SUBMIT:      ret = checkData();  if(ret) info.step = JSTEP_CHECKED; break;
//         case JSTEP_PREPROCESS:  return checkPre();
//         case JSTEP_WRITING:     return checkWriting();
//         case JSTEP_CHECKED: ;                  // no checks at these steps, just fall through and return true
//         case JSTEP_POSTPROCESS: ;
//         case JSTEP_QUEUED: ;
//         case JSTEP_RUNNING: ;
//         case JSTEP_COMPLETED: ret = true; break;
//         case JSTEP_ERR: ret = false; break;
        default: LOG_ERR << "check(): No check defined for step = " << (int)info.step << ende;
    }
    return ret;
}


bool MomfbdJob::checkCfg(void) {
    
    if( (runFlags&RF_FLATFIELD) && (runFlags&RF_CALIBRATE) ) {
        LOG_ERR << "Both FLATFIELD and CALIBRATE mode requested" << ende;
        return false;
    }

    if( objects.empty() ) return false;     // need at least 1 object
    
    for( auto& obj: objects ) {
        if( !obj->checkCfg() ) return false;
    }

    return true;
}


bool MomfbdJob::checkData(void) {
    
    info.progress[1] = 0;
    for( auto& obj: objects ) {
        if( !obj->checkData() ) return false;
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


