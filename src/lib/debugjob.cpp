#include "redux/debugjob.hpp"

#include "redux/translators.hpp"
#include "redux/file/fileana.hpp"
#include "redux/logger.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/bitoperations.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"

#include <algorithm>
#include <thread>

using namespace redux::file;
using namespace redux::util;
using namespace redux;
using namespace std;

#define lg Logger::mlg
namespace {

    const string thisChannel = "debugjob";
    static Job* createDebugJob( void ) {
        return new DebugJob();
    }

    int64_t mandelbrot( const complex<double>& p, int maxIterations = 100 ) {

        complex<double> Z( 0, 0 );
        int32_t i( 0 );
        for( i = 1; i <= maxIterations; ++i ) {
            Z = Z * Z + p;
            if( std::norm( Z ) > 2 ) {
                break;
            }
        }
        if( i > maxIterations || ( log( i ) < log( maxIterations ) / 2 ) ) return 0;

        return -i;

    }



}
size_t DebugJob::jobType = Job::registerJob( "DebugJob", createDebugJob );


uint64_t DebugJob::unpackParts( const char* ptr, WorkInProgress& wip, bool swap_endian ) {

    using redux::util::unpack;
    size_t nParts;
    uint64_t count = unpack( ptr, nParts, swap_endian );
    wip.parts.resize( nParts );
    for( auto& it : wip.parts ) {
        if(!it) it.reset( new DebugPart() );
        count += it->unpack( ptr+count, swap_endian );
    }
    return count;
}


DebugJob::DebugJob( void ) : maxIterations( 1000 ), gamma( 1 ), xSize( 1920 ), ySize( 1080 ), coordinates { -1.9, 1.9, -0.9, 0.9 } {
    info.typeString = "debugjob";
}


DebugJob::~DebugJob( void ) {

}


void DebugJob::parsePropertyTree( bpo::variables_map& vm, bpt::ptree& tree ) {
    
    Job::parsePropertyTree( vm, tree );
    
    maxIterations = tree.get<uint32_t>( "MAX_ITERATIONS", 1000 );
    patchSize = tree.get<uint32_t>( "PATCH_SIZE", 200 );
    gamma = tree.get<double>( "GAMMA", 1.0 );
    
    vector<uint32_t> tmp = tree.get<vector<uint32_t>>( "IMAGE_SIZE", {1920, 1080} );
    if( tmp.size() == 2 ) {
        xSize = tmp[0];
        ySize = tmp[1];
    }
    
    vector<double> tmpD = tree.get<vector<double>>( "COORDINATES", { -1.9, 1.9, -0.9, 0.9 } );
    if( tmpD.size() != 4 ) {
        tmpD = { -1.9, 1.9, -0.9, 0.9 };
    }
    
    for( size_t i = 0; i < 4; ++i ) coordinates[i] = tmpD[i];
    
}


bpt::ptree DebugJob::getPropertyTree( bpt::ptree* root ) {

    bpt::ptree tree = Job::getPropertyTree();         // get common properties

    tree.put( "MAX_ITERATIONS", maxIterations );
    tree.put( "PATCH_SIZE", patchSize );
    tree.put( "GAMMA", gamma );
    vector<double> tmpD( 4 );
    for( size_t i = 0; i < 4; ++i ) tmpD[i] = coordinates[i];
    tree.put( "COORDINATES", tmpD );
    vector<uint32_t> tmp = { xSize, ySize };
    tree.put( "IMAGE_SIZE", tmp );

    if( root ) {
        root->push_back( bpt::ptree::value_type( "debugjob", tree ) );
    }

    return tree;
}


size_t DebugJob::size( void ) const {
    size_t sz = Job::size();
    sz += 4 * sizeof( uint32_t ) + 5 * sizeof( double );
    return sz;
}


uint64_t DebugJob::pack( char* ptr ) const {

    using redux::util::pack;

    uint64_t count = Job::pack( ptr );
    count += pack( ptr+count, maxIterations );
    count += pack( ptr+count, patchSize );
    count += pack( ptr+count, gamma );
    count += pack( ptr+count, xSize );
    count += pack( ptr+count, ySize );
    count += pack( ptr+count, coordinates, 4 );

#ifdef DEBUG_
    LOG_TRACE << "DebugJob::pack():  count=" << count << " sz=" << size();
#endif
    
    return count;

}


uint64_t DebugJob::unpack( const char* ptr, bool swap_endian ) {

    using redux::util::unpack;

    uint64_t count = Job::unpack( ptr, swap_endian );
    count += unpack( ptr+count, maxIterations, swap_endian );
    count += unpack( ptr+count, patchSize, swap_endian );
    count += unpack( ptr+count, gamma, swap_endian );
    count += unpack( ptr+count, xSize, swap_endian );
    count += unpack( ptr+count, ySize, swap_endian );
    count += unpack( ptr+count, coordinates, 4, swap_endian );

#ifdef DEBUG_
    LOG_TRACE << "DebugJob::unpack():  swap_endian=" << swap_endian << "  count=" << count << " sz=" << size();
#endif

    return count;

}


void DebugJob::checkParts( void ) {
    
    uint8_t mask = 0;
    for( auto & it : jobParts ) {
        /*if( it.second->step & JSTEP_ERR && (it.second->nRetries<info.maxPartRetries)) {    // TODO: handle failed parts.
            it.second->nRetries++;
            it.second->step &= ~JSTEP_ERR;
        }*/
        mask |= it.second->step;
    }

#ifdef DEBUG_
    LOG_TRACE << "DebugJob::checkParts():   nParts=" << jobParts.size() << "   state=" << bitString(mask);
#endif
    
    if( mask & JSTEP_ERR ) {    // TODO: handle failed parts.

    }

    if( countBits( mask ) == 1 ) {  // if all parts have the same "step", set the whole job to that step.
        info.step.store( mask );
    }

}

bool DebugJob::check(void) {
    bool ret(false);
    unique_lock<mutex> lock(jobMutex);
    int val = info.step;
    switch (val) {
        case 0:                 ret = true; if(ret) info.step = JSTEP_SUBMIT; break;
        case JSTEP_SUBMIT:      ret = true; if(ret) info.step = JSTEP_PREPROCESS; break;
        case JSTEP_PREPROCESS: ;                  // no checks at these steps, just fall through and return true
        case JSTEP_QUEUED: ;
        case JSTEP_RUNNING: ;
        case JSTEP_POSTPROCESS: ;
        case JSTEP_COMPLETED: ret = true; break;
        default: LOG_ERR << "check(): No check defined for step = " << (int)info.step << " (" << stepString(info.step) << ")";
    }
    return ret;
}

bool DebugJob::getWork( WorkInProgress& wip, uint8_t nThreads ) {

#ifdef DEBUG_
    LOG_TRACE << "DebugJob::getWork("<<(int)nThreads<<")";
#endif
    
    uint8_t step = info.step.load();
    wip.parts.clear();
     // run pre-/postprocessing if local
    if( ((step == JSTEP_PREPROCESS)||(step == JSTEP_POSTPROCESS)) && !wip.connection ) {
        return true;
    }
    
    if ( step == JSTEP_QUEUED ) {                       // starting processing
        info.step = step = JSTEP_RUNNING;
    }
    
    if( step == JSTEP_RUNNING ) {
        unique_lock<mutex> lock( jobMutex );
        size_t nParts = std::min( nThreads, info.maxThreads) * 2;
        for( auto & it : jobParts ) {
            if( it.second->step == JSTEP_QUEUED ) {
                it.second->step = JSTEP_RUNNING;
                wip.parts.push_back( it.second );
                info.state.store( JSTATE_ACTIVE );
                if( wip.parts.size() == nParts ) break;
            }
        }
        checkParts();
    }
    return wip.parts.size();
}


void DebugJob::ungetWork( WorkInProgress& wip ) {
    unique_lock<mutex> lock( jobMutex );
    for( auto & it : wip.parts ) {
        it->step = JSTEP_QUEUED;
    }
    wip.parts.clear();
}


void DebugJob::returnResults( WorkInProgress& wip ) {
    unique_lock<mutex> lock( jobMutex );
    checkParts();
    for( auto & it : wip.parts ) {
        auto dpart = static_pointer_cast<DebugPart>( it );
        jobParts[it->id]->step = dpart->step;
        jobParts[it->id]->result = dpart->result;
    }
    wip.parts.clear();
    checkParts();
}


bool DebugJob::run( WorkInProgress& wip, boost::asio::io_service& service, uint8_t maxThreads ) {

#ifdef DEBUG_
    LOG_TRACE << "DebugJob::run("<<(int)maxThreads<<") ";
#endif

    uint8_t step = info.step.load();
    if( step < JSTEP_SUBMIT ) {
        info.step.store( JSTEP_SUBMIT );        // do nothing before submitting
        return true;                            // run again
    }
    else if( step == JSTEP_PREPROCESS ) {
        preProcess();                           // preprocess on master, split job in parts
    }
    else if( step == JSTEP_RUNNING || step == JSTEP_QUEUED ) {          // main processing
        service.reset();
        for( auto & it : wip.parts ) {
            service.post( boost::bind( &DebugJob::runMain, this, boost::ref( it ) ) );
        }
        
        boost::thread_group pool;
        size_t nThreads = std::min( maxThreads, info.maxThreads)*2;
        for( size_t t = 0; t < nThreads; ++t ) {
            pool.create_thread( boost::bind( &boost::asio::io_service::run, &service ) );
        }

        LOG_DEBUG << "DebugJob::run ThreadCount = "<<(int)pool.size();
        pool.join_all();

    }
    else if( step == JSTEP_POSTPROCESS ) {
        postProcess();                          // postprocess on master, collect results, save...
    }
    else {
        LOG << "DebugJob::run()  unrecognized step = " << ( int )info.step.load();
        info.step.store( JSTEP_ERR );
    }
    return false;
}


void DebugJob::preProcess( void ) {

#ifdef DEBUG_
    LOG_TRACE << "DebugJob::preProcess()";
#endif

    if( xSize < 2 || ySize < 2 ) return;

    double stepX = ( coordinates[1] - coordinates[0] ) / ( xSize - 1 );
    double stepY = ( coordinates[3] - coordinates[2] ) / ( ySize - 1 );
    size_t count = 0;
    unique_lock<mutex> lock;
    vector<size_t> indices;
    vector<PartPtr> pts;

    for( uint32_t i = 0; i < xSize; i += patchSize ) {
        size_t lX = std::min( i + patchSize - 1, xSize - 1 );
        double x = coordinates[0] + i * stepX;
        for( uint32_t j = 0; j < ySize; j += patchSize ) {
            size_t lY = std::min( j + patchSize - 1, ySize - 1 );
            double y = coordinates[2] + j * stepY;
            PartPtr part( new DebugPart() );
            part->id = ++count;
            part->step = JSTEP_QUEUED;
            part->sortedID = part->id;
            part->xPixelL = i; part->xPixelH = lX;
            part->yPixelL = j; part->yPixelH = lY;
            part->beginX = x; part->endX = coordinates[0] + lX * stepX;
            part->beginY = y; part->endY = coordinates[2] + lY * stepY;
            pts.push_back( part );
            indices.push_back( count );
        }

    }

    std::random_shuffle( indices.begin(), indices.end() );
    count = 0;
    for( auto & it : pts ) {
        it->id = indices[count++];
        jobParts.insert( pair<size_t, PartPtr>( it->id, it ) );
    }
    info.step.store( JSTEP_QUEUED );

}


void DebugJob::runMain( Part::Ptr& part ) {

#ifdef DEBUG_
    LOG_TRACE << "DebugJob::runMain()";
#endif
    
    auto pptr = static_pointer_cast<DebugPart>( part );
    //part->step = JSTEP_POSTPROCESS;
    //return;
    // temporaries, to avoid cache collisions.
    uint32_t sizeX = pptr->xPixelH - pptr->xPixelL + 1;
    uint32_t sizeY = pptr->yPixelH - pptr->yPixelL + 1;
    double stepX = ( pptr->endX - pptr->beginX ) / ( sizeX - 1 );
    double stepY = ( pptr->endY - pptr->beginY ) / ( sizeY - 1 );
    double beginX = pptr->beginX;
    double beginY = pptr->beginY;

    size_t id = pptr->id;
    size_t sid = pptr->sortedID;
    uint32_t max_iters = maxIterations;

    int32_t pid = getpid();

//     auto tmp = sharedArray<int64_t>( sizeY, sizeX );
//     auto ptr = tmp.get();
    pptr->result.resize( sizeY, sizeX );
    auto tmp = pptr->result.get(sizeY, sizeX);
    auto ptr = tmp.get();

    for( uint32_t ix = 0; ix < sizeX; ++ix ) {
        double x = beginX + ix * stepX;
        for( uint32_t iy = 0; iy < sizeY; ++iy ) {
            double y = beginY + iy * stepY;

            ptr[iy][ix] = mandelbrot( complex<double>( x, y ), max_iters );

            if( ptr[iy][ix] < 0 ) continue;

            if( ix < iy ) {                                 // top-left triangle showing the real part-ID (should increase upwards and to the right)
                ptr[iy][ix] = sid;
            }
            else if( ix > ( sizeY - iy ) ) {                // right triangle: the unsorted part-ID (=processing order)
                ptr[iy][ix] = id;
            }
            else  {                                         // bottom left triangle: pid, to distinguish parts processed on different machines or instances.
                ptr[iy][ix] = pid;
            }
        }
    }

    //pptr->result.resize( sizeY, sizeX );
    //memcpy( pptr->result.ptr(), tmp.get()[0], sizeY * sizeX * sizeof( int64_t ) );

    part->step = JSTEP_POSTPROCESS;

}



void DebugJob::postProcess( void ) {

#ifdef DEBUG_
    LOG_TRACE << "DebugJob::postProcess()";
#endif
    //info.step.store( JSTEP_COMPLETED );
    //info.state.store( JSTATE_IDLE );
    //return;
    auto image = sharedArray<int16_t>( ySize, xSize );
    int16_t** img = image.get();

    int64_t minPID, maxPID, minID, maxID, minSID, maxSID;
    minPID = minID = minSID = UINT32_MAX;
    maxPID = maxID = maxSID = 0;
    for( auto & it : jobParts ) {

        auto ptr = static_pointer_cast<DebugPart>( it.second );

        uint32_t sizeX = ptr->xPixelH - ptr->xPixelL + 1;
        uint32_t sizeY = ptr->yPixelH - ptr->yPixelL + 1;

        auto blaha = reshapeArray( ptr->result.ptr( 0 ), sizeY, sizeX );
        auto res = blaha.get();

        for( uint32_t ix = 0; ix < sizeX; ++ix ) {
            for( uint32_t iy = 0; iy < sizeY; ++iy ) {
                int64_t tmp = res[iy][ix];
                if( tmp < 0 ) {
                    continue;      // to skip the contour for the normalization
                }
                if( ix < iy ) {
                    if( tmp > maxSID ) maxSID = tmp;
                    if( tmp < minSID ) minSID = tmp;
                }
                else if( ix > ( sizeY - iy ) ) {
                    if( tmp > maxID ) maxID = tmp;
                    if( tmp < minID ) minID = tmp;
                }
                else {
                    if( tmp > maxPID ) maxPID = tmp;
                    if( tmp < minPID ) minPID = tmp;
                }
            }
        }
    }

    for( auto & it : jobParts ) {

        auto ptr = static_pointer_cast<DebugPart>( it.second );

        uint32_t sizeX = ptr->xPixelH - ptr->xPixelL + 1;
        uint32_t sizeY = ptr->yPixelH - ptr->yPixelL + 1;

        auto blaha = reshapeArray( ptr->result.ptr( 0 ), sizeY, sizeX );
        auto res = blaha.get();

        for( uint32_t ix = 0; ix < sizeX; ++ix ) {
            for( uint32_t iy = 0; iy < sizeY; ++iy ) {
                int64_t tmp = res[iy][ix];

                if( tmp < 0 ) {
                    img[ptr->yPixelL + iy][ptr->xPixelL + ix] = 0;
                    continue;
                }

                if( ix < iy ) {
                    if( maxSID == minSID ) img[ptr->yPixelL + iy][ptr->xPixelL + ix] = 0;
                    else img[ptr->yPixelL + iy][ptr->xPixelL + ix] = ( tmp - minSID + 1 ) * 1.0 / ( maxSID - minSID + 1 ) * INT16_MAX;
                }
                else if( ix > ( sizeY - iy ) ) {
                    if( maxID == minID ) img[ptr->yPixelL + iy][ptr->xPixelL + ix] = 0;
                    else img[ptr->yPixelL + iy][ptr->xPixelL + ix] = ( tmp - minID + 1 ) * 1.0 / ( maxID - minID + 1 ) * INT16_MAX;
                }
                else {
                    if( maxPID == minPID ) img[ptr->yPixelL + iy][ptr->xPixelL + ix] = 0;
                    else img[ptr->yPixelL + iy][ptr->xPixelL + ix] = ( tmp - minPID + 1 ) * 1.0 / ( maxPID - minPID + 1 ) * INT16_MAX;
                }

            }
        }

    }


    Ana::Ptr hdr( new Ana() );

    hdr->m_ExtendedHeader = "DebugJob";
    hdr->m_Header.datyp = Ana::ANA_WORD;

    hdr->m_Header.ndim = 2;
    hdr->m_Header.dim[0] = xSize;
    hdr->m_Header.dim[1] = ySize;

    Ana::write( "debugjob_output.f0", reinterpret_cast<char*>( *img ), hdr );

    info.step.store( JSTEP_COMPLETED );
    info.state.store( JSTATE_IDLE );

}


size_t DebugJob::DebugPart::size( void ) const {
    size_t sz = Part::size();
    sz += 4 * sizeof( uint32_t ) + 4 * sizeof( double ) + result.size() + sizeof( size_t );
    return sz;
}


uint64_t DebugJob::DebugPart::pack( char* ptr ) const {

    using redux::util::pack;

    uint64_t count = Part::pack( ptr );
    count += pack( ptr+count, xPixelL );
    count += pack( ptr+count, xPixelH );
    count += pack( ptr+count, yPixelL );
    count += pack( ptr+count, yPixelH );
    count += pack( ptr+count, beginX );
    count += pack( ptr+count, endX );
    count += pack( ptr+count, beginY );
    count += pack( ptr+count, endY );
    count += pack( ptr+count, sortedID );
    count += result.pack( ptr+count );

    return count;
    
}


uint64_t DebugJob::DebugPart::unpack( const char* ptr, bool swap_endian ) {

    using redux::util::unpack;

    uint64_t count = Part::unpack( ptr, swap_endian );
    count += unpack( ptr+count, xPixelL, swap_endian );
    count += unpack( ptr+count, xPixelH, swap_endian );
    count += unpack( ptr+count, yPixelL, swap_endian );
    count += unpack( ptr+count, yPixelH, swap_endian );
    count += unpack( ptr+count, beginX, swap_endian );
    count += unpack( ptr+count, endX, swap_endian );
    count += unpack( ptr+count, beginY, swap_endian );
    count += unpack( ptr+count, endY, swap_endian );
    count += unpack( ptr+count, sortedID, swap_endian );
    count += result.unpack( ptr+count, swap_endian );

    return count;

}
