#include <redux/logging/logger.hpp>

#include <redux/logging/logtofile.hpp>
#include <redux/logging/logtostream.hpp>
#include <redux/logging/logtonetwork.hpp>
#include <redux/network/protocol.hpp>
#include <redux/util/datautil.hpp>
#include <redux/util/stringutil.hpp>

#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

using namespace redux::logging;
using namespace redux::network;
using namespace redux::util;
using namespace std;

uint8_t Logger::defaultLevelMask = LOG_UPTO(LOG_LEVEL_NORMAL);
thread_local LogItem Logger::threadItem;

int Logger::getDefaultLevel( void ) {
    if( defaultLevelMask == 0 ) return 0;
    int cnt(1);
    uint8_t tmp = defaultLevelMask;
    while( (tmp>>=1) ) cnt++;
    return cnt;
}


pair<string, string> Logger::customParser( const string& s ) { // custom parser to handle multiple -q/-v flags (e.g. -vvvv)

    if( s.find( "-v" ) == 0 || s.find( "--verbose" ) == 0 ) { //
        int count = std::count( s.begin(), s.end(), 'v' );
        while( count-- ) defaultLevelMask = (defaultLevelMask<<1)+1;
    }
    else if( s.find( "-q" ) == 0 || s.find( "--quiet" ) == 0 ) { //
        int count = std::count( s.begin(), s.end(), 'q' );
        while( count-- ) defaultLevelMask >>= 1;
    }
    return make_pair( string(), string() );                 // no need to return anything, we handle the verbosity directly.
}


string Logger::environmentMap( const string &envName ) {

    static map<string, string> vmap;
    if( vmap.empty() ) {
//        vmap["RDX_LOGFILE"] = "log-file";
        vmap["RDX_VERBOSITY"] = "verbosity";
    }
    map<string, string>::const_iterator ci = vmap.find( envName );
    if( ci == vmap.end() ) {
        return "";
    } else {
        return ci->second;
    }
}

bpo::options_description Logger::getOptions( const string& application_name ) {

    bpo::options_description logging( "Logging Options" );
    logging.add_options()
    ( "verbosity", bpo::value< int >(), "Specify verbosity level (0-8, 0 means no output)."
      " The environment variable RDX_VERBOSITY will be used as default if it exists." )
    ( "verbose,v", bpo::value<vector<string>>()->implicit_value( vector<string>( 1, "1" ), "" )
      ->composing(), "More output. (ignored if --verbosity is specified)" )
    ( "quiet,q", bpo::value<vector<string>>()->implicit_value( vector<string>( 1, "-1" ), "" )
      ->composing(), "Less output. (ignored if --verbosity is specified)" )

    ( "log-file,L", bpo::value< vector<string> >()->implicit_value( vector<string>( 1, "" ), "" )
      ->composing(),
      "Print output to file."
      /*" The environment variable RDX_LOGFILE will be used as default if it exists."*/ )
    ( "log-stdout,d", "Debug mode. Will write output from all channels to stdout."
      " --log-file can not be used together with this option." )
    ;

    return logging;
}



Logger::Logger( bpo::variables_map& vm ) : LogOutput( defaultLevelMask, 1 ) {

    if( vm.count( "verbosity" ) > 0 ) {         // if --verbosity N is specified, use it.
        defaultLevelMask = LOG_UPTO(vm["verbosity"].as<int>());
    }
    
    mask = defaultLevelMask;
    
    if( vm.count( "log-stdout" ) ) {
        addStream( cout, mask );
    }
    else if( vm.count( "log-file" ) ) {
        vector<string> logfiles = vm["log-file"].as<vector<string>>();
        bool hasDefault = false;
        for( auto & filename : logfiles ) {
            if( filename == "" ) {
                if( hasDefault ) {
                    continue;
                }
                filename = vm["appname"].as<string>() + ".log";
                hasDefault = true;
            }
            addFile( filename, mask, false );
        }
    }


}



Logger::Logger(void) : LogOutput(defaultLevelMask,1) {

}


Logger::~Logger() {

    this->flushAll();
    
    // clear them now so that the flushBuffer call from the LogOutput destructor does not attempt to access deleted items.
    connections.clear();
    outputs.clear();

}


void Logger::append( LogItem &i ) {
    
    if( mask && !(i.entry.getMask() & mask) ) {
        return;
    }

    LogItemPtr tmpItem( new LogItem() );
    tmpItem->setLogger( this );
    tmpItem->entry = i.entry;
    tmpItem->context = i.context;
    
    addItem( tmpItem );

}


void Logger::flushBuffer( void ) {

    unique_lock<mutex> lock( queueMutex );
    vector<LogItemPtr> tmpQueue( itemQueue.begin(), itemQueue.end() );
    itemQueue.clear();
    itemCount = 0;
    lock.unlock();

    unique_lock<mutex> lock2( outputMutex );
    for( auto &it: outputs ) {
        it.second->addItems( tmpQueue );
    }
    
}


void Logger::flushAll( void ) {

    flushBuffer();

    unique_lock<mutex> lock( outputMutex );
    for( auto &it: outputs ) {
        it.second->flushBuffer();
    }

    
}


void Logger::addLogger( Logger& out ) {

    if( &out == this ) return;                  // avoid infinite loop
    out.removeOutput( hexString(this) );        // avoid infinite loop
    
    string name = hexString(&out);
    unique_lock<mutex> lock( outputMutex );
    OutputMap::iterator it = outputs.find( name );
    if( it == outputs.end() ) {
        std::shared_ptr<LogOutput> output( &out, [](LogOutput* p){} );
        outputs.insert(make_pair(name,output));

    }
    
}


void Logger::addStream( ostream& strm, uint8_t m, unsigned int flushPeriod ) {
    
    if( m == 0 ) {
        m = getMask();
    }
    string name = hexString(&strm);
    unique_lock<mutex> lock( outputMutex );
    OutputMap::iterator it = outputs.find( name );
    if( it == outputs.end() ) {
        std::shared_ptr<LogOutput> output( new LogToStream( strm, m, flushPeriod) );
        outputs.insert(make_pair(name,output));
    }
    
}


void Logger::addFile( const std::string &filename, uint8_t m, bool replace, unsigned int flushPeriod ) {
    
    if( m == 0 ) {
        m = getMask();
    }
    bfs::path tmpPath = cleanPath( filename );
    string name = tmpPath.string().c_str();
    unique_lock<mutex> lock( outputMutex );
    OutputMap::iterator it = outputs.find( name );
    if( it == outputs.end() ) {
        std::shared_ptr<LogOutput> output( new LogToFile( name, m, replace, flushPeriod) );
        outputs.insert(make_pair(name,output));

    }

}


void Logger::addNetwork( const TcpConnection::Ptr conn, uint32_t id, uint8_t m, unsigned int flushPeriod ) {
    
    if( !conn ) {
        return;
    }
    
    if( m == 0 ) {
        m = getMask();
    }
    
    try {
        auto test RDX_UNUSED = conn->socket().remote_endpoint();  // check if endpoint exists
        string name = hexString( conn.get() );
        unique_lock<mutex> lock( outputMutex );
        OutputMap::iterator it = outputs.find( name );
        if( it == outputs.end() ) {
            Command cmd;
            *conn << CMD_LOG_CONNECT << id;
            *conn >> cmd;
            if( cmd != CMD_OK ) {
                getItem(LOG_MASK_ERROR) << "Failed to connect to logging server (reply: " << cmd << ")" << ende;
            } else {
                std::shared_ptr<LogOutput> output( new LogToNetwork( conn, id, m, flushPeriod) );
                outputs.insert(make_pair( name, output ));
            }
        }
    } catch ( std::exception& e ) {
        getItem(LOG_MASK_ERROR) << "Logger::addNetwork exception: " << e.what() << ende;
        throw;
    }
    
}


void Logger::removeOutput( const string& name ) {
    
    unique_lock<mutex> lock( outputMutex );
    OutputMap::iterator it = outputs.find( name );
    if( it != outputs.end() ) {
        it->second->flushBuffer();
        outputs.erase(it);
    }
}


void Logger::removeAllOutputs( void ) {

    unique_lock<mutex> lockq( queueMutex );
    unique_lock<mutex> lock( outputMutex );
    for( auto &op: outputs ) {
        op.second->flushBuffer();
    }
    outputs.clear();

}


void Logger::addConnection( TcpConnection::Ptr conn, network::Host::Ptr host ) {
    
    conn->setCallback( bind( &Logger::netReceive, this, std::placeholders::_1 ) );
    conn->setErrorCallback( bind( &Logger::removeConnection, this, std::placeholders::_1 ) );
    unique_lock<mutex> lock( outputMutex );
    connections.insert( make_pair(conn, host) );
    
}


void Logger::removeConnection( TcpConnection::Ptr conn ) {
    
    unique_lock<mutex> lock( outputMutex );
    try {
        auto it = connections.find( conn );
        if( it != connections.end() ) {
            connections.erase( it );
        }
        if( conn ) {
            conn->setErrorCallback(nullptr);
            conn->setCallback(nullptr);
            conn->socket().close();
            //conn->idle();
        }
    } catch( std::exception& e ) {
        getItem(LOG_MASK_ERROR) << "Exception caught while removing a connection: " << e.what() << ende; 
    }
    
}


void Logger::netReceive( TcpConnection::Ptr conn ) {
    
    Command cmd = CMD_ERR;
    try {
        *conn >> cmd;
        if( cmd != CMD_PUT_LOG ) {
            throw std::runtime_error("Unexpected input.");
        }
    } catch( const std::exception& e ) {      // disconnected or wrong first byte -> remove connection and return.  TODO narrower catch
        //cout << "netReceive() exception: " << e.what() << endl;
        removeConnection(conn);
        return;
    } catch( ... ) {      // disconnected or wrong first byte -> remove connection and return.  TODO narrower catch
        //cout << "netReceive() uncaught exception. " << hexString(conn.get()) << endl;
        removeConnection(conn);
        return;
    }

    
    unique_lock<mutex> lock( outputMutex );
    auto it = connections.find( conn );
    string hostname = "client";
    if( it != connections.end() ) {
        if( it->second ) {
            it->second->touch();
            hostname = it->second->info.name;
            size_t pos = hostname.find_first_of(". ");
            if( pos != string::npos) {
                hostname.erase( pos );
            }
        }
    }
    lock.unlock();
    
    try {
        auto test RDX_UNUSED = conn->socket().remote_endpoint();  // check if endpoint exists
        if( !conn->socket().is_open() ) {
            throw runtime_error("Connection closed.");
        }

        size_t blockSize;
        shared_ptr<char> buf = conn->receiveBlock( blockSize );               // reply
        *conn << CMD_OK;

        if( blockSize ) {
            vector<LogItemPtr> tmpQueue;
            uint64_t count(0);
            char* ptr = buf.get();
            while( count < blockSize ) {
                LogItemPtr tmpItem( new LogItem() );
                tmpItem->setLogger( this );
                if( tmpItem->context.empty() ) {
                    tmpItem->context = hostname;
                }
                count += tmpItem->unpack( ptr+count, conn->getSwapEndian() );
                tmpQueue.push_back( tmpItem );
            }
            addItems(tmpQueue);
        }

    } catch ( const std::exception& e ) {
        getItem(LOG_MASK_WARNING) << "Exception caught while receiving log messages from " << hostname
            << ": " << e.what() << ende; 
        removeConnection(conn);
        return;
    }
    
    conn->idle();

}

