#include "redux/momfbd/momfbdjob.hpp"

#include "redux/momfbd/defines.hpp"
#include "redux/momfbd/util.hpp"

#include "redux/translators.hpp"
#include "redux/file/fileio.hpp"
#include "redux/constants.hpp"
#include "redux/logger.hpp"
#include "redux/util/bitoperations.hpp"

#include <boost/algorithm/string.hpp>

using namespace redux::momfbd;
using namespace redux::file;
using namespace redux::util;
using namespace redux;
using namespace std;
using boost::algorithm::iequals;


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


uint64_t MomfbdJob::unpackParts( const char* ptr, std::vector<Part::Ptr>& parts, bool swap_endian ) {
    using redux::util::unpack;
    size_t nParts;
    uint64_t count = unpack( ptr, nParts, swap_endian );
    parts.resize( nParts );
    for( auto & it : parts ) {
        it.reset( new PatchData );
        count += it->unpack( ptr+count, swap_endian );
    }
    return count;
}


void MomfbdJob::parsePropertyTree( po::variables_map& vm, bpt::ptree& tree ) {

    Job::parsePropertyTree( vm, tree );
    LOG_DEBUG << "MomfbdJob::parsePropertyTree()";
    
    // possibly override cfg-entries with command-line arguments
    if( vm.count( "simx" ) ) tree.put( "SIM_X", vm["simx"].as<string>() );
    if( vm.count( "simy" ) ) tree.put( "SIM_Y", vm["simy"].as<string>() );
    if( vm.count( "imgn" ) ) tree.put( "IMAGE_NUM", vm["imgn"].as<string>() );
    if( vm.count( "output-file" ) ) tree.put( "output-file", vm["output-file"].as<string>() );
    if( vm.count( "force" ) ) tree.put( "OVERWRITE", true );
    if( vm.count( "swap" ) ) tree.put( "SWAP", true );

    GlobalCfg::parseProperties(tree);

    size_t nObj( 0 );
    for( auto & it : tree ) {
        if( iequals( it.first, "OBJECT" ) ) {
            Object* tmpObj = new Object( *this );
            tmpObj->parsePropertyTree( it.second );
            if( nObj < outputFiles.size() ) {
                tmpObj->outputFileName = outputFiles[nObj++];
            }
            objects.push_back( shared_ptr<Object>( tmpObj ) );
        }
    }
    if( outputFiles.size() > objects.size() ) {
        LOG_WARN << outputFiles.size() << " output file names specified but only " << objects.size() << " objects found.";
    }
    LOG_DEBUG << "MomfbdJob::parsePropertyTree() done.";

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
    sz += pupil.size();
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
    count += pupil.pack( ptr+count);
    
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
    count += pupil.unpack( ptr+count, swap_endian );
    return count;
}


void MomfbdJob::checkParts( void ) {

    uint8_t mask = 0;
    for( auto & it : patches ) {
        /*if( (it.second->step & JSTEP_ERR) && (it.second->nRetries<info.maxPartRetries)) {
            it.second->nRetries++;
            it.second->step &= ~JSTEP_ERR;
        }*/
        mask |= it.second->step;
    }

    if( mask & JSTEP_ERR ) {    // TODO: handle failed parts.

    }

    LOG << "checkParts(): mask = " << bitString(mask);
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
        if(!wip.connection) {   // local worker, check if there are results to write
            LOG_DEBUG << "getWork(): LOCAL " << (bool)wip.connection;
            for( auto & it : patches ) {
                //    LOG_DEBUG << "getWork(): patch " << it.second->id <<  "  step = " << bitString(it.second->step);
                if( it.second->step & JSTEP_POSTPROCESS ) {
                    LOG_DEBUG << "getWork(): PP-patch   step = " << bitString(it.second->step);
                    wip.parts.push_back( it.second );
                }
            }
            if( wip.parts.size() ) {
                LOG_DEBUG << "getWork(): nPP = " << wip.parts.size();
                ret = true;
            }
        }
        if(!ret && wip.connection) {
            for( auto & it : patches ) {
                if( it.second->step == JSTEP_QUEUED ) {
                    it.second->step = JSTEP_RUNNING;
                    wip.parts.push_back( it.second );
                    ret = true;
                    break;// only 1 part at a time for MomfbdJob
                }
            }
        }
    }
    
        LOG_DEBUG << "getWork(): step = " << (int)step << " conn = " << (bool)wip.connection;
    if( ret ) {
        unique_lock<mutex> lock( jobMutex );
        checkParts();
    }
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
        auto patch = static_pointer_cast<PatchData>( it );
        patches[it->id]->step = patch->step;
        //patches[it->id]->result = patch->result;
    redux::file::Ana::write( "patch_" + to_string(patch->index.x) + "_" + to_string(patch->index.y) + ".f0", patch->images );
    }
    wip.parts.clear();
    checkParts();
}


void MomfbdJob::init(void) {
    
    //redux::image::KL_cfg* coeff = nullptr;
    if( klMinMode || klMaxMode ) {
        LOG << "calculating Karhunen-Loeve coefficients";
        //coeff = legacy::klConfig(klMinMode,klMaxMode);
    }
    
    for( shared_ptr<Object>& obj: objects ) {
        obj->init(); //coeff);
    }
}


void MomfbdJob::cleanup(void) {
    for( shared_ptr<Object>& obj: objects ) {
        obj->cleanup();
    }
}


bool MomfbdJob::run( WorkInProgress& wip, boost::asio::io_service& service, uint8_t maxThreads ) {
    
    uint8_t jobStep = info.step.load();
    uint8_t patchStep = wip.parts.size() ? wip.parts[0]->step : 0;
    if( jobStep == JSTEP_PREPROCESS ) {
        preProcess(service);                           // preprocess on master: load, flatfield, split in patches
    }
    else if( jobStep == JSTEP_RUNNING || jobStep == JSTEP_QUEUED ) {
        size_t nThreads = std::min( maxThreads, info.maxThreads);
        if( patchStep == JSTEP_POSTPROCESS ) {      // store results
            storePatches(wip, service, nThreads);
        } else {                                    // main processing
            for( Part::Ptr& it : wip.parts ) {      // momfbd jobs will only get 1 part at a time, this is just to keep things generic.
                //LOG_DETAIL << "Configuring slave";
                WorkSpace ws( *this, static_pointer_cast<PatchData>(it) );
                ws.init(service);            // Global setup, allocations, etc.
                runThreadsAndWait(service, nThreads);
                // Run
                service.post( std::bind( &MomfbdJob::runMain, this, boost::ref( ws ) ) );
                runThreadsAndWait(service, nThreads);
                // Get/return results and cleanup.
                ws.collectResults();
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

    // load shared files synchronously (dark,gain,psf,offset...)
    if( pupil.nDimensions() < 2 && pupilFile != "" ) {
        service.post( std::bind( util::loadPupil, pupilFile, std::ref(pupil), 0 ) );
    }
    for( shared_ptr<Object>& obj : objects ) {
        obj->loadData(service);
    }

    info.maxThreads = 12;
    runThreadsAndWait(service, info.maxThreads);

    // Done loading files -> start the preprocessing (flatfielding etc.)
    
    Point16 imageSizes;
    for( shared_ptr<Object>& obj : objects ) {
        Point16 tmp = obj->clipImages();
        if(imageSizes.x == 0) {
            imageSizes = tmp;
        } else if( tmp != imageSizes ) {
            throw std::logic_error("The clipped images have different sizes for the different objects, please verify the ALIGN_CLIP values.");
        }
        obj->preprocessData(service);
    }
    runThreadsAndWait(service, info.maxThreads);
    
    // Done pre-processing -> normalize within each object
    
    size_t nTotalImages(0);
    for( shared_ptr<Object>& obj : objects) {
        obj->normalize(service);
        nTotalImages += obj->nImages(nTotalImages);
    }
    runThreadsAndWait(service, info.maxThreads);
    
    // Done normalizing -> collect images to master-stack
    
    imageStack.resize(nTotalImages,imageSizes.y,imageSizes.x);
    for( shared_ptr<Object>& obj : objects ) {
        obj->collectImages(imageStack);
    }

    redux::file::Ana::write( "masterstack.f0", imageStack );
    
    // Done normalizing -> split in patches
    
    int minimumOverlap = 16;                // desired width of blending zone in pixels
    //int patchSeparation = 3 * patchSize / 4 - minimumOverlap; // target separation
    if( subImagePosX.empty() ) {
        // TODO: do we need to handle single patches ?? (this only supports nPatches >= 2)
        // TODO: verify michiel's splitting method
        int firstPos = maxLocalShift + patchSize/2;
        int lastPos = imageSizes.x - firstPos - 1;
        nPatchesX = 2;
        double separation = (lastPos-firstPos)/static_cast<double>(nPatchesX-1);
        double overlap = std::max(patchSize-separation,0.0);
        while(overlap < minimumOverlap) {
            ++nPatchesX;
            separation = (lastPos-firstPos)/static_cast<double>(nPatchesX-1);
            overlap = std::max(patchSize-separation,0.0);
        }
        for( size_t i = 0; i < nPatchesX; ++i ) {
            subImagePosX.push_back(static_cast<uint32_t>( i*separation + firstPos ) );
        }
        LOG << "MomfbdJob::preProcess(): Generated patch positions  " << printArray(subImagePosX,"x-pos");
    }

    if( subImagePosY.empty() ) {
        int firstPos = maxLocalShift + patchSize/2;
        int lastPos = imageSizes.y - firstPos - 1;
        nPatchesY = 2;
        double separation = (lastPos-firstPos)/static_cast<double>(nPatchesY-1);
        double overlap = std::max(patchSize-separation,0.0);
        while(overlap < minimumOverlap) {
            ++nPatchesY;
            separation = (lastPos-firstPos)/static_cast<double>(nPatchesY-1);
            overlap = std::max(patchSize-separation,0.0);
        }
        for( size_t i = 0; i < nPatchesY; ++i ) {
            subImagePosY.push_back(static_cast<uint32_t>( i*separation + firstPos ) );
        }
        LOG << "MomfbdJob::preProcess(): Generated patch positions  " << printArray(subImagePosY,"y-pos");
    }
 
    if( subImagePosX.empty() || subImagePosY.empty() ) {
        LOG_ERR << "MomfbdJob::preProcess(): No patches specified or generated, can't continue.";
        info.step.store( JSTEP_ERR );
        info.state.store( JSTATE_IDLE );
        return;
    }

    for( shared_ptr<Object>& obj : objects ) {
        service.post( std::bind( &Object::prepareStorage, obj.get() ) );
    }
    runThreadsAndWait(service, 1); //objects.size());  TODO: fix multithreaded write
    
    size_t count = 0;
    uint16_t yid=0;
    uint16_t halfBlockSize = patchSize/2 + maxLocalShift;
    for( const uint16_t posY : subImagePosY ) {       // y-coordinate of patch-centre
        uint16_t xid=0;
        uint16_t trimmedPosY = std::min(std::max(halfBlockSize,posY),uint16_t(imageSizes.y-halfBlockSize));       // stay inside borders
        if( trimmedPosY != posY ) LOG_WARN << "MomfbdJob::preProcess() y-position of patch was outside the image area and was trimmed: " << posY << " -> " << trimmedPosY;
        for( const uint16_t posX : subImagePosX ) {   // x-coordinate of patch-centre
            uint16_t trimmedPosX = std::min(std::max(halfBlockSize,posX),uint16_t(imageSizes.x-halfBlockSize));   // stay inside borders
            if( trimmedPosX != posX ) LOG_WARN << "MomfbdJob::preProcess() x-position of patch was outside the image area and was trimmed: " << posX << " -> " << trimmedPosX;
            PatchData::Ptr patch( new PatchData() );
            patch->step = JSTEP_QUEUED;
            patch->first.x = maxLocalShift;
            patch->first.y = maxLocalShift;
            patch->last.x = patch->first.x+patchSize-1;
            patch->last.y = patch->first.y+patchSize-1;
            patch->id = ++count;
            patch->images = Array<float>(imageStack,0,nTotalImages-1,trimmedPosY-halfBlockSize, trimmedPosY+halfBlockSize-1, trimmedPosX-halfBlockSize, trimmedPosX+halfBlockSize-1);
            patch->setIndex( yid, xid++ );
            service.post( std::bind( &MomfbdJob::applyLocalOffsets, this, patch ) );
            patches.insert( make_pair( patch->id, patch ) );
        }
        yid++;
    }

    LOG_DETAIL << "MomfbdJob::preProcess()  nPatches = " << patches.size();

    runThreadsAndWait(service, 1); //info.maxThreads);  TODO: fix multithreaded localoffstets

    info.step.store( JSTEP_QUEUED );

    LOG_DETAIL << "MomfbdJob::preProcess()  Done.";
}


void MomfbdJob::applyLocalOffsets( PatchData::Ptr patch ) {
    
//    LOG_TRACE << "MomfbdJob::applyLocalOffsets() #" << patch->id;
//     size_t totalPatchSize(0);
//     for( auto & it : objects ) {
//         totalPatchSize += it->sizeOfPatch(patch->nPixels());
//     }
// 
//     patch->dataSize = totalPatchSize;
//     patch->data = sharedArray<char>(totalPatchSize);
//     char* ptr = patch->data.get();
//     uint64_t count(0);
    for( shared_ptr<Object>& obj: objects ) {
        obj->applyLocalOffsets(patch);
    }
//     
//     if(count != totalPatchSize) {
//         LOG_WARN << "Estimation of patch data-size was wrong:  est = " << totalPatchSize << "  real = " << ptrdiff_t(ptr-patch->data.get());
//     }
    // TODO: compress and store in swapfile
 //   LOG_TRACE << "MomfbdJob::applyLocalOffsets() #" << patch->id << "  All done !!";
}
 

void MomfbdJob::runMain( WorkSpace& ws ) {



    LOG << "MomfbdJob::runMain()  patch#" << ws.data->id << " (" << ws.data->index.x << "," << ws.data->index.y << ")  x1=" << ws.data->first.x << "  y1=" << ws.data->first.y;
usleep(100000);
    //     // temporaries, to avoid cache collisions.
//     uint32_t sizeX = pptr->xPixelH - pptr->xPixelL + 1;
//     uint32_t sizeY = pptr->yPixelH - pptr->yPixelL + 1;
//     double stepX = ( pptr->endX - pptr->beginX ) / ( sizeX - 1 );
//     double stepY = ( pptr->endY - pptr->beginY ) / ( sizeY - 1 );
//     double beginX = pptr->beginX;
//     double beginY = pptr->beginY;
//
//     size_t id = pptr->id;
//     size_t sid = pptr->sortedID;
//     uint32_t max_iters = maxIterations;
//
//     int32_t pid = getpid();
//
//     auto tmp = sharedArray<int64_t>( sizeY, sizeX );
//     auto ptr = tmp.get();
//     double x, y;
//     for( uint32_t ix = 0; ix < sizeX; ++ix ) {
//         x = beginX + ix * stepX;
//         for( uint32_t iy = 0; iy < sizeY; ++iy ) {
//             y = beginY + iy * stepY;
//
//             ptr[iy][ix] = mandelbrot( complex<double>( x, y ), max_iters );
//
//             if( ptr[iy][ix] < 0 ) continue;
//
//             if( ix < iy ) {                                 // top-left triangle showing the real part-ID (should increase upwards and to the right)
//                 ptr[iy][ix] = sid;
//             }
//             else if( ix > ( sizeY - iy ) ) {                // right triangle: the unsorted part-ID (=processing order)
//                 ptr[iy][ix] = id;
//             }
//             else  {                                         // bottom left triangle: pid, to distinguish parts processed on different machines or instances.
//                 ptr[iy][ix] = pid;
//             }
//         }
//     }
//
//     pptr->result.reset( sizeY, sizeX );
//     memcpy( pptr->result.ptr(), tmp.get()[0], sizeY * sizeX * sizeof( int64_t ) );

    //sleep(1);
    ws.data->step = JSTEP_POSTPROCESS;

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
    info.state.store( JSTATE_IDLE );

}


bool MomfbdJob::check(void) {
    bool ret(false);
    unique_lock<mutex> lock(jobMutex);
    switch (info.step) {
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
    
    if( !GlobalCfg::check() ) return false;
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


