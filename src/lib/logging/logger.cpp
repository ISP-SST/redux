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


bpo::options_description Logger::getOptions( const string& application_name ) {

    bpo::options_description logging( "Logging Options" );
    logging.add_options()
    ( "verbosity", bpo::value< int >(), "Specify verbosity level" )
    ( "verbose,v", bpo::value<vector<string>>()->implicit_value( vector<string>( 1, "1" ), "" )
      ->composing(), "More output. (ignored if --verbosity is specified)" )
    ( "quiet,q", bpo::value<vector<string>>()->implicit_value( vector<string>( 1, "-1" ), "" )
      ->composing(), "Less output. (ignored if --verbosity is specified)" )

    ( "log-file,L", bpo::value< vector<string> >()->default_value( vector<string>( 1, "" ), "" )
      ->composing(),
      "Print logevents to file. If no log-file name is"
      " specified the log will be written to ./<name>.log,"
      " where <name> is the name given to this application instance." )
    ( "log-stdout,d", "Debug mode. Will write output from all channels to stdout. --log-file can not be used together with this option." )
    ;

    return logging;
}



Logger::Logger( bpo::variables_map& vm ) : LogOutput(defaultLevelMask,1) {

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


void Logger::addNetwork( const TcpConnection::Ptr& conn, uint32_t id, uint8_t m, unsigned int flushPeriod ) {
    
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
                conn->socket().close();
            } else {
                std::shared_ptr<LogOutput> output( new LogToNetwork( conn, id, m, flushPeriod) );
                outputs.insert(make_pair( name, output ));
            }
        }
    } catch ( ... ) { }
    
}


void Logger::removeOutput( const string& name ) {
    
    unique_lock<mutex> lock( outputMutex );
    OutputMap::iterator it = outputs.find( name );
    if( it != outputs.end() ) {
        outputs.erase (it);
    }
}


void Logger::removeAllOutputs( void ) {

    unique_lock<mutex> lock( outputMutex );
    outputs.clear();

}


void Logger::addConnection( TcpConnection::Ptr& conn, network::Host::Ptr& host ) {
    
    TcpConnection::callback oldCallback = conn->getCallback();
    conn->setCallback( bind( &Logger::netReceive, this, std::placeholders::_1 ) );
    unique_lock<mutex> lock( outputMutex );
    auto cinfo = make_pair( host, oldCallback );
    connections.insert( make_pair(conn, cinfo) );
    
}


void Logger::removeConnection( TcpConnection::Ptr& conn ) {
    
    TcpConnection::callback oldCallback = conn->getCallback();
    unique_lock<mutex> lock( outputMutex );
    auto it = connections.find( conn );
    if( it != connections.end() ) {
        conn->setCallback( it->second.second );
        connections.erase( it );
    }
    
}


void Logger::netReceive( TcpConnection::Ptr conn ) {
    
    
    auto it = connections.find( conn );
    string hostname = "client";
    if( it != connections.end() ) {
        if( it->second.first ) {
            hostname = it->second.first->info.name;
            size_t pos = hostname.find_first_of(". ");
            if( pos != string::npos) {
                hostname.erase( pos );
            }
        }
    }

    try {
        auto test RDX_UNUSED = conn->socket().remote_endpoint();  // check if endpoint exists
        if( !conn->socket().is_open() ) {
            throw exception();
        }

        size_t blockSize;
        shared_ptr<char> buf = conn->receiveBlock( blockSize );               // reply

        if( blockSize ) {
            vector<LogItemPtr> tmpQueue;
            uint64_t count(0);
            char* ptr = buf.get();
            while( count < blockSize ) {
                LogItemPtr tmpItem( new LogItem() );
                tmpItem->setLogger( this );
                tmpItem->context.scope = hostname;
                count += tmpItem->unpack( ptr+count, conn->getSwapEndian() );
                tmpQueue.push_back( tmpItem );
            }
            addItems(tmpQueue);
        }

    } catch ( ... ) {
        removeConnection(conn);
        return;
    }
    
    conn->idle();

}

