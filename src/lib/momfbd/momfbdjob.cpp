#include "redux/momfbd/momfbdjob.hpp"

#include "redux/momfbd/util.hpp"
#include "redux/momfbd/workspace.hpp"

#include "redux/logger.hpp"
#include "redux/util/bitoperations.hpp"
#include "redux/util/stringutil.hpp"

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
        if(!wip.parts[0]) wip.parts[0].reset(new PatchData(*this));
        count += wip.parts[0]->unpack( ptr+count, swap_endian );
        if( wip.nParts > 1 ) {
            globalData.reset( new GlobalData(*this) );
            count += globalData->unpack( ptr+count, swap_endian );
        }
    }
    return count;
}
/*
    for( shared_ptr<Object>& obj: objects ) {
        ObjectData::Ptr objData(new ObjectData(obj));
        patch->data.push_back(objData);
        objData->init(patch->index.y,patch->index.x);
        //obj->applyLocalOffsets(objData);
    }
*/
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
    for( auto & it : tree ) {
        if( iequals( it.first, "OBJECT" ) ) {
            Object* tmpObj = new Object( *this, nObj );
            tmpObj->parsePropertyTree( it.second );
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

    for( shared_ptr<Object>& obj: objects ) {
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
    for( const shared_ptr<Object>& obj: objects ) {
        sz += obj->size();
    }
    return sz;
}


uint64_t MomfbdJob::pack( char* ptr ) const {
    using redux::util::pack;
    uint64_t count = Job::pack( ptr );
    count += GlobalCfg::pack( ptr+count );
    count += pack( ptr+count, (uint16_t)objects.size() );
    for( const shared_ptr<Object>& obj: objects ) {
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
    for( shared_ptr<Object>& obj: objects ) {
        obj.reset(new Object(*this));
        count += obj->unpack( ptr+count, swap_endian );
    }
    return count;
}


size_t MomfbdJob::nImages(void) const {
    size_t nTotalImages(0);
    for( const shared_ptr<Object>& obj : objects) {
        nTotalImages += obj->nImages();
    }
    return nTotalImages;
}


void MomfbdJob::checkParts( void ) {

    uint8_t mask = 0;
    for( auto & it : patches ) {
        /*if( (it.second->step & JSTEP_ERR) && (it.second->nRetries<info.maxPartRetries)) {
            it.second->nRetries++;
            it.second->step &= ~JSTEP_ERR;
        }*/
        mask |= it->step;
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
            for( auto & it : patches ) {
                if( it->step & JSTEP_POSTPROCESS ) {
                    LOG_DEBUG << "getWork(): PP-patch   step = " << bitString(it->step);
                    wip.parts.push_back( it );
                }
            }
            if( wip.parts.size() ) {
                LOG_DEBUG << "getWork(): nPP = " << wip.parts.size();
                ret = true;
            }
        }*/
        if(!ret /*&& wip.connection*/) {
            for( auto & it : patches ) {
                if( it->step == JSTEP_QUEUED ) {
                    it->step = JSTEP_RUNNING;
                    it->cacheLoad(true);
                    wip.parts.push_back( it );
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
    for( Part::Ptr& it : wip.parts ) {
        it->step = JSTEP_QUEUED;
    }
    wip.parts.clear();
}


#include "redux/file/fileana.hpp"

void MomfbdJob::returnResults( WorkInProgress& wip ) {
    unique_lock<mutex> lock( jobMutex );
    checkParts();
    for( Part::Ptr& it : wip.parts ) {
        PatchData::Ptr patch = static_pointer_cast<PatchData>( it );
        patch->step = JSTEP_POSTPROCESS;
        patches(patch->index.y,patch->index.x) = patch;
        patch->cacheStore(true);
        //patches[it->id]->result = patch->result;
    //redux::file::Ana::write( "patch_" + to_string(patch->index.x) + "_" + to_string(patch->index.y) + ".f0", patch->images );
    }
    wip.parts.clear();
    checkParts();
}


void MomfbdJob::init(void) {
    for( shared_ptr<Object>& obj: objects ) {
        obj->init();
    }
}


void MomfbdJob::cleanup(void) {
    for( shared_ptr<Object>& obj: objects ) {
        obj->cleanup();
    }
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
        if( patchStep == JSTEP_POSTPROCESS ) {      // store results
            //storePatches(wip, service, nThreads);
        } else {                                    // main processing
            if(!globalData) {
                 globalData.reset( new GlobalData(*this) );
            }
            if( !proc ) {
                proc.reset( new WorkSpace(*this) );            // Initialize, allocations, etc.
                proc->init();
            }
            for( Part::Ptr& it : wip.parts ) {      // momfbd jobs will only get 1 part at a time, this is just to keep things generic.
                // Run main processing
                proc->run(static_pointer_cast<PatchData>(it), service, nThreads);
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
    for( shared_ptr<Object>& obj : objects ) {
        Point16 tmp = obj->getImageSize();
        if(imageSizes.x == 0) {
            imageSizes = tmp;
        } else if( tmp != imageSizes ) {    // TBD: allow for different patchsizes (i.e. pixelsize/ccd-size) for different objects/channels.
            throw std::logic_error("The clipped images have different sizes for the different objects, please verify the ALIGN_CLIP values.");
        }
    }
    service.post( std::bind( &MomfbdJob::initCache, this) );    // TBD: should cache initialization be parallelized?
    //runThreadsAndWait(service, 1); //info.maxThreads);
    
    // load shared files synchronously (dark,gain,psf,offset...)
    /*if( pupil.nDimensions() < 2 && pupilFile != "" ) {
        service.post( std::bind( util::loadPupil, pupilFile, std::ref(pupil), 0 ) );
    }*/
    for( shared_ptr<Object>& obj : objects ) {
        obj->loadData(service);
    }

    //info.maxThreads = 12;
    runThreadsAndWait(service, info.maxThreads);

    // Done loading files -> start the preprocessing (flatfielding etc.)
    for( shared_ptr<Object>& obj : objects ) {
        obj->preprocessData(service);
    }
    runThreadsAndWait(service, info.maxThreads);
    
    // Done pre-processing -> normalize within each object
    for( shared_ptr<Object>& obj : objects) {
        obj->normalize(service);
    }
    
    // Done normalizing -> collect images to master-stack
    
    /*imageStack.resize(nTotalImages,imageSizes.y,imageSizes.x);
    for( shared_ptr<Object>& obj : objects ) {
        obj->collectImages(imageStack);
    }
    redux::file::Ana::write( "masterstack.f0", imageStack );
    */
    
    // Done normalizing -> split in patches
    
    //int patchSeparation = 3 * patchSize / 4 - minimumOverlap; // target separation
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
            patch->setPath(to_string(info.id)+"/patch_"+(string)patch->index);      // will cache PatchData in sysDataDir/jobID/patch_(?,?)
            patch->step = JSTEP_QUEUED;
            patch->roi.first = Point16(subImagePosY[y]-halfPatchSize-1, subImagePosX[x]-halfPatchSize-1);   // subImagePosX/Y is 1-based
            patch->roi.last = patch->roi.first+ps-1;
            patch->id = ++count;
            service.post( std::bind( &PatchData::getData, patch.get() ) );
            patches(y,x) = patch;
        }
    }

    //service.post( std::bind( &MomfbdJob::initCache, this) );    // TBD: should cache initialization be parallelized?
    
    LOG_DETAIL << "MomfbdJob::preProcess()  nPatches = " << patches.nElements();
    runThreadsAndWait(service, info.maxThreads);

    info.step.store( JSTEP_QUEUED );

    LOG_DETAIL << "MomfbdJob::preProcess()  Done.";
}

void MomfbdJob::initCache(void) {
    LOG_DETAIL << "MomfbdJob::initCache()";
    if(!globalData) {       // create GlobalData if they don't exist
        globalData.reset(new GlobalData(*this));
    }
    
    globalData->constraints.init();
 
    for( shared_ptr<Object>& obj: objects ) {
        obj->initCache();
    }

    LOG_DETAIL << "MomfbdJob::initCache()  Done.";
}

void MomfbdJob::initPatchData( PatchData::Ptr patch ) {
    
    //LOG << "MomfbdJob::applyLocalOffsets() #" << patch->id;
//     size_t totalPatchSize(0);
//     for( auto & it : objects ) {
//         totalPatchSize += it->sizeOfPatch(patch->nPixels());
//     }
// 
//     patch->dataSize = totalPatchSize;
//     patch->data = sharedArray<char>(totalPatchSize);
//     char* ptr = patch->data.get();
//     uint64_t count(0);
//     patch->objects.resize(objects.size());
//     auto it = patch->objects.begin();
//     for( shared_ptr<Object>& obj: objects ) {
//         it->channels.resize(obj->channels.size());
//         auto it2 = it->channels.begin();
//         for( const shared_ptr<Channel>& ch : obj->channels ) {
//             ch->getPatchData(*it2++,patch->index.y,patch->index.x);
//         }
// //         ObjectData::Ptr objData(new ObjectData(obj));
// //         patch->objects.push_back(objData);
// //         it->init(patch->index.y,patch->index.x);
//         it++;
//         //obj->applyLocalOffsets(objData);
//     }
//     
//     if(count != totalPatchSize) {
//         LOG_WARN << "Estimation of patch data-size was wrong:  est = " << totalPatchSize << "  real = " << ptrdiff_t(ptr-patch->data.get());
//     }
    // TODO: compress and store in swapfile
 //   LOG_TRACE << "MomfbdJob::applyLocalOffsets() #" << patch->id << "  All done !!";
}
 

void MomfbdJob::storePatches( WorkInProgress& wip, boost::asio::io_service& service, uint8_t nThreads) {
    LOG << "MomfbdJob::storePatches()";
    for( shared_ptr<Object>& obj: objects ) {
        obj->storePatches(wip, service, nThreads);
    }
    
    for( Part::Ptr& part: wip.parts ) {
        //PatchData::Ptr patch = static_pointer_cast<PatchData>(part);
        part->step = JSTEP_COMPLETED;
    }

    
}


void MomfbdJob::postProcess( boost::asio::io_service& service ) {

    LOG << "MomfbdJob::postProcess()";
    for( shared_ptr<Object>& obj : objects ) {
        service.post( std::bind( &Object::writeResults, obj.get(), patches ) );
    }
    runThreadsAndWait(service, 1); //objects.size());  TODO: fix multithreaded write
    
//     auto image = sharedArray<int16_t>( ySize, xSize );
//     int16_t** img = image.get();
//
//     int64_t minPID, maxPID, minID, maxID, minSID, maxSID;
//     minPID = minID = minSID = UINT32_MAX;
//     maxPID = maxID = maxSID = 0;
//     for( auto & it : jobParts ) {
//
//         auto ptr = static_pointer_cast<DebugPart>( it.second );
//
//         uint32_t sizeX = ptr->xPixelH - ptr->xPixelL + 1;
//         uint32_t sizeY = ptr->yPixelH - ptr->yPixelL + 1;
//
//         auto blaha = reshapeArray( ptr->result.ptr( 0 ), sizeY, sizeX );
//         auto res = blaha.get();
//
//         for( uint32_t ix = 0; ix < sizeX; ++ix ) {
//             for( uint32_t iy = 0; iy < sizeY; ++iy ) {
//                 int64_t tmp = res[iy][ix];
//                 if( tmp < 0 ) {
//                     continue;      // to skip the contour for the normalization
//                 }
//                 if( ix < iy ) {
//                     if( tmp > maxSID ) maxSID = tmp;
//                     if( tmp < minSID ) minSID = tmp;
//                 }
//                 else if( ix > ( sizeY - iy ) ) {
//                     if( tmp > maxID ) maxID = tmp;
//                     if( tmp < minID ) minID = tmp;
//                 }
//                 else {
//                     if( tmp > maxPID ) maxPID = tmp;
//                     if( tmp < minPID ) minPID = tmp;
//                 }
//             }
//         }
//     }
//
//     for( auto & it : jobParts ) {
//
//         auto ptr = static_pointer_cast<DebugPart>( it.second );
//
//         uint32_t sizeX = ptr->xPixelH - ptr->xPixelL + 1;
//         uint32_t sizeY = ptr->yPixelH - ptr->yPixelL + 1;
//
//         auto blaha = reshapeArray( ptr->result.ptr( 0 ), sizeY, sizeX );
//         auto res = blaha.get();
//
//         for( uint32_t ix = 0; ix < sizeX; ++ix ) {
//             for( uint32_t iy = 0; iy < sizeY; ++iy ) {
//                 size_t tmp = res[iy][ix];
//
//                 if( tmp < 0 ) {
//                     img[ptr->yPixelL + iy][ptr->xPixelL + ix] = 0;
//                     continue;
//                 }
//
//                 if( ix < iy ) {
//                     if( maxSID == minSID ) img[ptr->yPixelL + iy][ptr->xPixelL + ix] = 0;
//                     else img[ptr->yPixelL + iy][ptr->xPixelL + ix] = ( tmp - minSID + 1 ) * 1.0 / ( maxSID - minSID + 1 ) * INT16_MAX;
//                 }
//                 else if( ix > ( sizeY - iy ) ) {
//                     if( maxID == minID ) img[ptr->yPixelL + iy][ptr->xPixelL + ix] = 0;
//                     else img[ptr->yPixelL + iy][ptr->xPixelL + ix] = ( tmp - minID + 1 ) * 1.0 / ( maxID - minID + 1 ) * INT16_MAX;
//                 }
//                 else {
//                     if( maxPID == minPID ) img[ptr->yPixelL + iy][ptr->xPixelL + ix] = 0;
//                     else img[ptr->yPixelL + iy][ptr->xPixelL + ix] = ( tmp - minPID + 1 ) * 1.0 / ( maxPID - minPID + 1 ) * INT16_MAX;
//                 }
//
//             }
//         }
//
//     }
//
//
//     Ana::Ptr hdr( new Ana() );
//
//     hdr->m_ExtendedHeader = "DebugJob";
//     hdr->m_Header.datyp = Ana::ANA_WORD;
//
//     hdr->m_Header.ndim = 2;
//     hdr->m_Header.dim[0] = xSize;
//     hdr->m_Header.dim[1] = ySize;
//
//     std::ofstream file( "debugjob_output.f0" );
//
//     Ana::write( file, reinterpret_cast<char*>( *img ), hdr );

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
    
    for( shared_ptr<Object>& obj: objects ) {
        if( !obj->checkCfg() ) return false;
    }
    
    return true;
}


bool MomfbdJob::checkData(void) {

    for( shared_ptr<Object>& obj: objects ) {
        if( !obj->checkData() ) return false;
    }
    
    return true;
}
        
        
const MomfbdJob& MomfbdJob::operator=(const GlobalCfg& rhs) {
    GlobalCfg::operator=(rhs);
    return *this;
}


