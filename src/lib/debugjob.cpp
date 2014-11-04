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
size_t DebugJob::jobType = Job::registerJob( "Debug", createDebugJob );


const char* DebugJob::unpackParts( const char* ptr, std::vector<Part::Ptr>& parts, bool swap_endian ) {

    using redux::util::unpack;
    size_t nParts;
    ptr = unpack( ptr, nParts, swap_endian );
    parts.resize( nParts );
    for( auto & it : parts ) {
        it.reset( new DebugPart() );
        ptr = it->unpack( ptr, swap_endian );
    }
    return ptr;
}


DebugJob::DebugJob( void ) : maxIterations( 1000 ), gamma( 1 ), xSize( 1920 ), ySize( 1080 ), coordinates { -1.9, 1.9, -0.9, 0.9 } {
    info.typeString = "Debug";

}


DebugJob::~DebugJob( void ) {

}


void DebugJob::parseProperties( po::variables_map& vm, bpt::ptree& tree ) {
    maxIterations = tree.get<uint32_t>( "MAX_ITERATIONS", 1000 );
    patchSize = tree.get<uint32_t>( "PATCH_SIZE", 200 );
    gamma = tree.get<double>( "GAMMA", 1.0 );
    vector<uint32_t> tmp = tree.get<vector<uint32_t>>( "IMAGE_SIZE", {1920, 1080} );
    cout << printArray( tmp, "parsed sizes" ) << endl;
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
        root->push_back( bpt::ptree::value_type( "debug", tree ) );
    }

    return tree;
}


size_t DebugJob::size( void ) const {
    size_t sz = Job::size();
    sz += 4 * sizeof( uint32_t ) + 5 * sizeof( double );
    return sz;
}


char* DebugJob::pack( char* ptr ) const {

    using redux::util::pack;

    ptr = Job::pack( ptr );
    ptr = pack( ptr, maxIterations );
    ptr = pack( ptr, patchSize );
    ptr = pack( ptr, gamma );
    ptr = pack( ptr, xSize );
    ptr = pack( ptr, ySize );
    ptr = pack( ptr, coordinates, 4 );

    return ptr;

}


const char* DebugJob::unpack( const char* ptr, bool swap_endian ) {

    using redux::util::unpack;

    ptr = Job::unpack( ptr, swap_endian );
    ptr = unpack( ptr, maxIterations, swap_endian );
    ptr = unpack( ptr, patchSize, swap_endian );
    ptr = unpack( ptr, gamma, swap_endian );
    ptr = unpack( ptr, xSize, swap_endian );
    ptr = unpack( ptr, ySize, swap_endian );
    ptr = unpack( ptr, coordinates, 4, swap_endian );

    return ptr;

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

    if( mask & JSTEP_ERR ) {    // TODO: handle failed parts.

    }

    if( countBits( mask ) == 1 ) {  // if all parts have the same "step", set the whole job to that step.
        info.step.store( mask );
    }

}

size_t DebugJob::getParts( WorkInProgress& wip ) {

    uint8_t step = info.step.load();
    wip.parts.clear();
    if( step == JSTEP_QUEUED || step == JSTEP_RUNNING ) {
        unique_lock<mutex> lock( jobMutex );
        size_t nParts = wip.peer->stat.nThreads;
        if( info.nThreads ) nParts = std::min( wip.peer->stat.nThreads, info.nThreads );
        for( auto & it : jobParts ) {
            if( it.second->step == JSTEP_QUEUED ) {
                it.second->step = JSTEP_RUNNING;
                wip.parts.push_back( it.second );
                info.step.store( JSTEP_RUNNING );
                info.state.store( JSTATE_ACTIVE );
                if( wip.parts.size() == nParts ) break;
            }
        }
        checkParts();
    }
    return wip.parts.size();
}


void DebugJob::ungetParts( WorkInProgress& wip ) {
    unique_lock<mutex> lock( jobMutex );
    for( auto & it : wip.parts ) {
        it->step = JSTEP_QUEUED;
    }
    wip.parts.clear();
}

void DebugJob::returnParts( WorkInProgress& wip ) {
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


bool DebugJob::run( WorkInProgress& wip, boost::asio::io_service& service, boost::thread_group& pool ) {

    uint8_t step = info.step.load();
    if( step < JSTEP_SUBMIT ) {
        info.step.store( JSTEP_SUBMIT );        // do nothing before submitting
        return true;                            // run again
    }
    else if( step == JSTEP_RECEIVED ) {
        preProcess();                           // preprocess on master, split job in parts
    }
    else if( step == JSTEP_RUNNING || step == JSTEP_QUEUED ) {          // main processing
        size_t nThreads = wip.peer->stat.nThreads;
        if( info.nThreads ) nThreads = std::min( wip.peer->stat.nThreads, info.nThreads );
        service.reset();

        for( auto & it : wip.parts ) {
            service.post( boost::bind( &DebugJob::runMain, this, boost::ref( it ) ) );
        }
        for( size_t t = 0; t < nThreads * 2; ++t ) {
            pool.create_thread( boost::bind( &boost::asio::io_service::run, &service ) );
        }

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

    auto pptr = static_pointer_cast<DebugPart>( part );

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

    auto tmp = sharedArray<int64_t>( sizeY, sizeX );
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

    pptr->result.resize( sizeY, sizeX );
    memcpy( pptr->result.ptr(), tmp.get()[0], sizeY * sizeX * sizeof( int64_t ) );

    part->step = JSTEP_POSTPROCESS;

}



void DebugJob::postProcess( void ) {

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
                size_t tmp = res[iy][ix];

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


char* DebugJob::DebugPart::pack( char* ptr ) const {

    using redux::util::pack;

    ptr = Part::pack( ptr );
    ptr = pack( ptr, xPixelL );
    ptr = pack( ptr, xPixelH );
    ptr = pack( ptr, yPixelL );
    ptr = pack( ptr, yPixelH );
    ptr = pack( ptr, beginX );
    ptr = pack( ptr, endX );
    ptr = pack( ptr, beginY );
    ptr = pack( ptr, endY );
    ptr = pack( ptr, sortedID );
    ptr = result.pack( ptr );

    return ptr;
}


const char* DebugJob::DebugPart::unpack( const char* ptr, bool swap_endian ) {

    using redux::util::unpack;

    ptr = Part::unpack( ptr, swap_endian );
    ptr = unpack( ptr, xPixelL, swap_endian );
    ptr = unpack( ptr, xPixelH, swap_endian );
    ptr = unpack( ptr, yPixelL, swap_endian );
    ptr = unpack( ptr, yPixelH, swap_endian );
    ptr = unpack( ptr, beginX, swap_endian );
    ptr = unpack( ptr, endX, swap_endian );
    ptr = unpack( ptr, beginY, swap_endian );
    ptr = unpack( ptr, endY, swap_endian );
    ptr = unpack( ptr, sortedID, swap_endian );
    ptr = result.unpack( ptr, swap_endian );

    return ptr;

}
