#include "redux/application.hpp"
#include "redux/logger.hpp"
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
namespace bpo = boost::program_options;
namespace bpt = boost::property_tree;

using namespace redux::util;
using namespace redux::network;
using namespace redux;

using namespace std;


#define lg Logger::lg
namespace {

    const string thisChannel = "jsub";

    // define options specific to this binary
    bpo::options_description getOptions( void ) {

        bpo::options_description options( "Program Options" );
        options.add_options()
        ( "master,m", bpo::value<string>()->default_value( "localhost" ),
          "Hostname/IP of a master to connect to."
          " The environment variable REDUX_MASTER can be used to override the default value." )
        ( "port,p", bpo::value<string>()->default_value( "30000" ),
          "Port to use when connecting to a master."
          " The environment variable REDUX_PORT can be used to override the default value." )
        ( "force,f", "Overwrite output file if exists" )
        ( "swap,s", "swap mode: write compressed data to swap file instead of keeping it in memory (useful for large problems)" )
        ( "config,c", bpo::value<string>()->default_value( "momfbd.cfg" ), "Configuration file to process." )
        ( "simx", bpo::value<string>(), "x coordinate[s] of subimages to restore" )
        ( "simy", bpo::value<string>(), "y coordinate[s] of subimages to restore" )
        ( "imgn,n", bpo::value<string>(), "Image numbers" )
        ( "sequence", bpo::value<string>(), "sequence number to insert in filename template." )
        ( "print,P", "(debug) print the parsed configuration to console and exit without uploading." )
        ( "no-check", "Don't verify the configuration." )
        ( "output-file,o", "Comma separated list of output file base names."
          "File names are applied to the objects in the order they are found in the config file."
          "If insufficient names are provided, a default name will be created for the remaining objects."
          "Excess names will generate a warning but are ignored otherwise. The names are base names only,"
          "an appropriate suffix will be attached (.fits/.f0)." )
        ;

        return options;
    }

    // define environment variables to use as defaults if the corresponding command-line option is not specified
    string environmentMap( const string &envName ) {

        static map<string, string> vmap;
        if( vmap.empty() ) {
            vmap["REDUX_VERBOSITY"] = "verbosity";  // For debugging this might be convenient.
            vmap["REDUX_MASTER"] = "master";        // If it exists, it will override the default value (localhost) above
            vmap["REDUX_PORT"] = "port";            // If it exists, it will override the default value (30000) above
        }
        map<string, string>::const_iterator ci = vmap.find( envName );
        if( ci == vmap.end() ) {
            return "";
        }
        else {
            return ci->second;
        }
    }
}

void uploadJobs(TcpConnection::Ptr conn, vector<Job::JobPtr>& jobs) {
    
    Host::HostInfo me, master;
    uint8_t cmd = CMD_CONNECT;
    boost::asio::write(conn->socket(),boost::asio::buffer(&cmd,1));
    boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));

    if( cmd == CMD_AUTH ) {
        // implement
    }
    if( cmd == CMD_CFG ) {  // handshake requested
        *conn << me;
        *conn >> master;
        boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));       // ok or err
    }
    if( cmd != CMD_OK ) {
        LOG_ERR << "Handshake with server failed.";
        return;
    }
   
    // all ok, upload jobs.
    size_t bufferSize = 0;
    for( auto &job: jobs ) {
        bufferSize += job->size();
    }
    
    auto buf = sharedArray<char>( bufferSize+sizeof(size_t)+1 );
    char* ptr = buf.get();
    uint64_t packedBytes = pack(ptr,CMD_ADD_JOB);
    packedBytes += pack(ptr+packedBytes,bufferSize);

    bufferSize += sizeof(size_t)+1;
    for( auto &job: jobs ) {
        packedBytes += job->pack(ptr+packedBytes);
    }
    if( packedBytes != bufferSize ) {
        LOG_ERR << "Packing of jobs failed:  bufferSize = " << bufferSize << "  packedBytes = " << packedBytes;
        return;
    }

    conn->writeAndCheck( buf, bufferSize );
    boost::asio::write(conn->socket(),boost::asio::buffer(buf.get(),bufferSize));
    boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));

    bool swap_endian = (me.littleEndian != master.littleEndian);
    uint64_t count, received;
    ptr = buf.get();

    if( cmd == CMD_OK ) {
        received = boost::asio::read( conn->socket(), boost::asio::buffer( ptr, sizeof(uint64_t) ) );
        if( received == sizeof(uint64_t) ) {
            unpack(ptr, count, swap_endian);
            size_t thisSize = count*sizeof(uint64_t);
            received = boost::asio::read( conn->socket(), boost::asio::buffer( ptr, thisSize ) );
            if( received == thisSize ) {
                if( count ) LOG << "Upload of " << count << " job(s) completed successfully. " << printArray(reinterpret_cast<size_t*>(buf.get()),count,"IDs");
            } else {
                LOG_ERR << "Failed to read job IDs.  received=" << received << " thisSize=" << thisSize;
            }
        } else LOG_ERR << "Failed to read number of job IDs.";
    } else {
        LOG_ERR << "Failure while sending jobs  (server reply = " << (int)cmd << "   " << bitString(cmd) << ")";
    }
    
    ptr = buf.get();
    received = boost::asio::read( conn->socket(), boost::asio::buffer( ptr, sizeof(uint64_t) ) );
    if( received == sizeof(uint64_t) ) {
        unpack(ptr, count, swap_endian);
        if( count ) {
            received = boost::asio::read( conn->socket(), boost::asio::buffer( ptr, count ) );    // hopefully the buffer allocated for the cfg-data is big enough for any messages.
            vector<string> messages;
            unpack(ptr, messages,swap_endian);
            if( !messages.empty() ) {
                string msgText = "Server messages:";
                for( auto& msg: messages ) {
                    msgText += "\n\t" + msg;
                }
                LOG << msgText;
            }
        }
    }


}

bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    bool found = false;
    while(start_pos != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        found = true;
        start_pos = str.find(from);
    }
    return found;
}


string filterOldCfg(string filename) {
    
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (!in) {
        throw (errno);
    }

    std::string text;
    in.seekg (0, std::ios::end);
    text.reserve (in.tellg());
    in.seekg (0, std::ios::beg);
    std::copy ( (std::istreambuf_iterator<char> (in)), std::istreambuf_iterator<char>(), std::back_inserter (text));
    in.close();

    bool found = replace(text, "channel{", "channel {");
    found |= replace(text, "object{", "object {");
    found |= replace(text, "=", " ");
    
    if (found) { // old cfg file, wrap inside momfbd { ... }
        text = "momfbd { \n" + text + "\n}";
    }

    return text;
}


int main (int argc, char *argv[]) {


    bpo::variables_map vm;
    bpo::options_description programOptions = getOptions();

    try {
        bpo::options_description& allOptions = Application::parseCmdLine( argc, argv, vm, &programOptions );

        // load matched environment variables according to the environmentMap() above.
        bpo::store( bpo::parse_environment( allOptions, environmentMap ), vm );
        vm.notify();
    }
    catch( const exception &e ) {
        cerr << "Error parsing commandline: " << e.what() << endl;// << programOptions << endl;
        return EXIT_FAILURE;
    }


    try {
        Logger logger( vm );
        bpt::ptree momfbd;
        std::stringstream filteredCfg(filterOldCfg(vm["config"].as<string>()));
        bpt::read_info (filteredCfg , momfbd);
        bool check = (vm.count ("no-check") + vm.count ("print")) == 0;
        vector<Job::JobPtr> jobs = Job::parseTree (vm, momfbd, check);

        if (vm.count ("print")) {       // dump configuration to console and exit
            bpt::ptree dump;
            for( auto & job : jobs ) {
                job->getPropertyTree( &dump );
            }
            bpt::write_info( cout<<endl, dump );
            return EXIT_SUCCESS;
        }

        boost::asio::io_service ioservice;
        TcpConnection::Ptr conn = TcpConnection::newPtr(ioservice);
        conn->connect( vm["master"].as<string>(), vm["port"].as<string>() );

        if( conn->socket().is_open() ) {
            std::vector<std::shared_ptr<std::thread> > threads;
            for ( int i=0; i<5; ++i) {
                shared_ptr<thread> t( new thread( boost::bind( &boost::asio::io_service::run, &ioservice ) ) );
                threads.push_back( t );
            }
            shared_ptr<thread> t( new thread( boost::bind( uploadJobs, conn, jobs) ) );
            threads.push_back( t );
            //thread t( boost::bind( &boost::asio::io_service::run, &ioservice ) );
            //thread tt( boost::bind( &boost::asio::io_service::run, &ioservice ) );
            //uploadJobs(conn, jobs);
            //ioservice.run();
            //t.join();
            //tt.join();
            for( auto & t : threads ) {
                t->join();
            }
        }

    }
    catch( const exception &e ) {
        cerr << "Uncaught exception (fatal): " << e.what() << endl;
    }

    return EXIT_SUCCESS;

}

