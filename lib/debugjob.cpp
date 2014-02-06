#include "redux/debugjob.hpp"

#include "redux/logger.hpp"

#include <boost/serialization/export.hpp>

BOOST_CLASS_EXPORT( redux::DebugJob )

using namespace redux;
using namespace std;

#define lg Logger::mlg
namespace {
    const string thisChannel = "debugjob";
    static Job* createDebugJob( void ) {
        return new DebugJob();
    }
}
size_t DebugJob::jobType = Job::registerJob( "DebugJob", createDebugJob );


DebugJob::DebugJob( void )  {
    LOG_DEBUG << "DebugJob::DebugJob()";
}


DebugJob::~DebugJob( void ) {
    LOG_DEBUG << "DebugJob::~DebugJob()";
}


void DebugJob::parseProperties( po::variables_map& vm, bpt::ptree& tree ) {
    LOG_DETAIL << "DebugJob::parseProperties()";
}

uint32_t DebugJob::preProcess( void ) {
    LOG_DETAIL << "DebugJob::preProcess()";
    return 0;
}


uint32_t DebugJob::postProcess( void ) {
    LOG_DETAIL << "DebugJob::postProcess()";
    return 0;
}


uint32_t DebugJob::runJob( void ) {
    LOG_DETAIL << "DebugJob::runJob()";
    return 0;
}


void* DebugJob::prePart( void ) {
    LOG_DETAIL << "DebugJob::prePart()";
    return nullptr;
}


void DebugJob::postPart( void* ) {
    LOG_DETAIL << "DebugJob::postPart()";
}


void* DebugJob::runPreJob( void ) {
    LOG_DETAIL << "DebugJob::runPreJob()";
    return nullptr;
}


void DebugJob::runPostJob( void* ) {
    LOG_DETAIL << "DebugJob::runPostJob()";
}
