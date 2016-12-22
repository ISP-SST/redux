#include "redux/application.hpp"
#include "redux/logging/logger.hpp"
#include "redux/debugjob.hpp"
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/network/tcpconnection.hpp"
#include "redux/network/host.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/stringutil.hpp"

#include <sstream>
#include <thread>

#include <boost/bind.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/program_options.hpp>

using namespace redux::logging;
using namespace redux::util;
using namespace redux::network;
using namespace redux;

using namespace std;


namespace {

    const string logChannel = "rctl";

    // define options specific to this binary
    bpo::options_description getOptions( void ) {

        bpo::options_description options( "Program Options" );
        options.add_options()
        ( "master,m", bpo::value<string>()->default_value( "localhost" ),
          "Hostname/IP of a master to connect to."
          " The environment variable RDX_MASTER can be used to override the default value." )
        ( "port,p", bpo::value<string>()->default_value( "30000" ),
          "Port to use when connecting to a master."
          " The environment variable RDX_PORT can be used to override the default value." )
        ( "force,f", "Bypass safeguards" )
        ( "kill,k", "Send exit command to Server." )
        ;

        return options;
    }

    // define environment variables to use as defaults if the corresponding command-line option is not specified
    string environmentMap( const string &envName ) {

        static map<string, string> vmap;
        if( vmap.empty() ) {
            vmap["RDX_VERBOSITY"] = "verbosity";  // For debugging this might be convenient.
            vmap["RDX_MASTER"] = "master";        // If it exists, it will override the default value (localhost) above
            vmap["RDX_PORT"] = "port";            // If it exists, it will override the default value (30000) above
        }
        map<string, string>::const_iterator ci = vmap.find( envName );
        if( ci == vmap.end() ) {
            return "";
        }
        else {
            return ci->second;
        }
    }
    
    bool swap_endian;
    
}


void checkReply( TcpConnection::Ptr conn, Logger& logger ) {


    uint8_t cmd = CMD_ERR;
    boost::asio::read( conn->socket(), boost::asio::buffer(&cmd,1) );

    try {
        
        switch( cmd ) {
            case CMD_OK: LOG_DEBUG << "Command successful." << ende; break;
            case CMD_NOTICE: {          // reply will have some messages appended
                uint64_t blockSize, received;
                received = boost::asio::read( conn->socket(), boost::asio::buffer( &blockSize, sizeof(uint64_t) ) );
                if( received == sizeof(uint64_t) ) {
                    if( swap_endian ) swapEndian( blockSize );
                    if( blockSize ) {
                        shared_ptr<char> buf( new char[blockSize+1], []( char* p ){ delete[] p; } );
                        char* ptr = buf.get();
                        memset( ptr, 0, blockSize+1 );
                        received = boost::asio::read( conn->socket(), boost::asio::buffer( ptr, blockSize ) );
                        if( received ) {
                            vector<string> messages;
                            unpack( ptr, messages, swap_endian );
                            if( !messages.empty() ) {
                                string msgText = "\nServer reply:";
                                for( auto& msg: messages ) {
                                    msgText += "\n\t" + msg;
                                }
                                LOG << msgText << ende;
                            }
                        } else LOG_ERR << "Server replied with an empty block." << ende;
                    }
                } else LOG_ERR << "Failed to receive server reply." << ende;
                break;
            }
            default: LOG_ERR << "Command failed: " << bitString(cmd) << ende;
        }

    } catch ( exception &e ) {
        LOG_ERR << "Error while checking server reply: " << e.what() << ende;
    }

}


void killServer( TcpConnection::Ptr conn, Logger& logger) {
   
    uint8_t cmd = CMD_DIE;
    boost::asio::write( conn->socket(), boost::asio::buffer(&cmd,1) );
    
    checkReply( conn, logger );
   
}


int main (int argc, char *argv[]) {
    
    bpo::variables_map vm;
    bpo::options_description programOptions = getOptions();

    try {
        bpo::options_description& allOptions = Application::parseCmdLine( argc, argv, vm, &programOptions );

        // load matched environment variables according to the environmentMap() above.
        bpo::store( bpo::parse_environment( allOptions, environmentMap ), vm );
#if BOOST_VERSION > 104800  // TODO check which version notify appears in
        vm.notify();
#endif
    }
    catch( const exception &e ) {
        cerr << "Error parsing commandline: " << e.what() << endl;// << programOptions << endl;
        return EXIT_FAILURE;
    }


    try {
        
        vm.erase("log-file");       // always log to cout for rctl
        vm.insert( std::make_pair("log-stdout", bpo::variable_value()) );
        
        Logger logger( vm );
        
        boost::asio::io_service ioservice;
        auto conn = TcpConnection::newPtr(ioservice);
        conn->connect( vm["master"].as<string>(), vm["port"].as<string>() );

        if( conn->socket().is_open() ) {
                    
            Host::HostInfo me, master;
            uint8_t cmd = CMD_CONNECT;
            boost::asio::write(conn->socket(),boost::asio::buffer(&cmd,1));
            boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));

            if( cmd == CMD_AUTH ) {
                // implement
            }
            if( cmd == CMD_CFG ) {  // handshake requested
                if( vm.count("force") ) me.peerType = Host::TP_MASTER;
                *conn << me;
                *conn >> master;
                boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));       // ok or err
            }
            
            if( cmd == CMD_OK ) {
                
                swap_endian = (me.littleEndian != master.littleEndian);
                
                if( vm.count("kill") ) killServer( conn, logger );
                
            } else {
                LOG_ERR << "Failed to connect to reduxd at " << vm["master"].as<string>() << ":" << vm["port"].as<string>() << ende;
            }
            
        }

    }
    catch( const exception &e ) {
        cerr << "Uncaught exception (fatal): " << e.what() << endl;
    }

    return EXIT_SUCCESS;

}

