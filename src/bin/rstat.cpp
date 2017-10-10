#include "redux/application.hpp"
#include "redux/job.hpp"
#include "redux/network/tcpconnection.hpp"
#include "redux/network/host.hpp"
#include "redux/util/endian.hpp"

#include <boost/date_time/posix_time/time_formatters.hpp>
#include <boost/program_options.hpp>

using namespace redux::util;
using namespace redux::network;
using namespace redux;

using namespace std;
using namespace boost::posix_time;


namespace {

    const string logChannel = "rdx_stat";

    // define options specific to this binary
    bpo::options_description getOptions( void ) {

        bpo::options_description options( "Program Options" );
        options.add_options()
        ( "master,m", bpo::value<string>()->default_value( "localhost" ),
          "Hostname/IP of a master to connect to."
          " The environment variable RDX_HOST can be used to override the default value." )
        ( "port,p", bpo::value<string>()->default_value( "30000" ), "Port to use, either when connecting to the master." )
        ( "jobs,j", bpo::value<int>()->implicit_value( 0 ), "Get joblist" )
        ( "slaves,s", bpo::value<int>()->implicit_value( 0 ), "Get slavelist" )
        ( "count,c", bpo::value<int>()->implicit_value( 0 ), "List only n first slaves/jobs, and the last." )
        ( "time,t", bpo::value<int>()->implicit_value( 1 ), "Loop and display list every (n) seconds" )
        ( "runtime,r", bpo::value<int>()->implicit_value( 1 ), "Sort slaves according to runtime." )
        ;

        return options;
    }

    // define environment variables to use as defaults if the corresponding command-line option is not specified
    string environmentMap( const string &envName ) {

        static map<string, string> vmap;
        if( vmap.empty() ) {
            vmap["RDX_VERBOSITY"] = "verbosity";  // For debugging this might be convenient.
            vmap["RDX_HOST"] = "master";        // If it exists, it will override the default value (localhost) above
            vmap["RDX_PORT"] = "port";            // This means the environment variable RDX_PORT will override the
                                                  // default value of 30000 specified above.
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

void printJobList( TcpConnection::Ptr conn, int nj ) {

    uint8_t cmd = CMD_JSTAT;
    boost::asio::write(conn->socket(),boost::asio::buffer(&cmd,1));

    uint64_t blockSize;
    shared_ptr<char> buf = conn->receiveBlock( blockSize );

    if ( !blockSize ) return;

    const char* ptr = buf.get();
    uint64_t count(0);
    try {
        typedef std::shared_ptr<Job::Info> InfoPtr;
        vector<InfoPtr> infos;
        size_t maxLength[2] = {0,0};       // jobName & user
        //cout << Job::Info::printHeader() << endl;
        while( count < blockSize ) {
            InfoPtr info(new Job::Info);
            count += info->unpack(ptr+count,conn->getSwapEndian());
            maxLength[0] = max(maxLength[0],info->name.size());
            maxLength[1] = max(maxLength[1],info->user.size()+info->host.size());
            infos.push_back(info);
        }
        std::sort( infos.begin(), infos.end(), [](const InfoPtr& a, const InfoPtr& b) {
            if(a->step != b->step) return (a->step > b->step);
            if(a->priority != b->priority) return (a->priority > b->priority);
            return (a->id < b->id);
            
        } );
        int nJobs = infos.size();
        bool addComma(false);
        if ( nj != 0 ) {
            if ( nj < nJobs ) addComma = true;
            nJobs = min(nj,nJobs);
        }
        cout << alignRight("#", 4) << alignRight("ID", 5) << alignCenter("type", 10) << alignCenter("started", 20); // + alignCenter("started", 20);
        cout << alignCenter("name", maxLength[0]+4) + alignCenter("user", maxLength[1]+5) + alignCenter("priority", 8) + alignCenter("state", 8);
        cout << endl;
        ptime now = boost::posix_time::second_clock::local_time();
        for( int i=0; i<nJobs; ++i ) {
            if( i == nJobs-1 ) {
                i = infos.size()-1;
                if (addComma) cout << "   :" << endl;
            }
                
            string info = alignRight(std::to_string(infos[i]->id), 5) + alignCenter(infos[i]->typeString, 10);
            string startedString;
            if ( infos[i]->startedTime < now ) startedString = to_iso_extended_string(infos[i]->startedTime);
            else startedString = to_iso_extended_string(infos[i]->submitTime);
            info += alignCenter(startedString, 20);
            info += alignCenter(infos[i]->name, maxLength[0]+4) + alignLeft(infos[i]->user + "@" + infos[i]->host, maxLength[1]+5);
            info += alignCenter(std::to_string(infos[i]->priority), 8);
            info += alignCenter(Job::stateTag(infos[i]->state), 3) + alignLeft(infos[i]->progressString, 15);
            cout << alignRight(to_string(i+1),4) << info << endl;
        }
    } catch ( const exception& e) {
        cerr << "printJobList: Exception caught while parsing block: " << e.what() << endl;
    }
    if( count != blockSize ) {
        cerr << "printJobList: Parsing of datablock failed,  count = " << count << "   blockSize = " << blockSize << "  bytes." << endl;
    }

}


namespace {
    
    enum{ by_name=1, by_runtime };
    
    bool host_by_runtime( const Host::Ptr& a, const Host::Ptr& b ){
        return (a->status.lastActive < b->status.lastActive);
    }
    
}


void printPeerList( TcpConnection::Ptr conn, int ns, int sorting ) {

    uint8_t cmd = CMD_PSTAT;
    boost::asio::write(conn->socket(),boost::asio::buffer(&cmd,1));

    uint64_t blockSize;
    shared_ptr<char> buf = conn->receiveBlock( blockSize );

    if ( !blockSize ) return;

    const char* ptr = buf.get();
    uint64_t count(0);
    try {
        vector<Host::Ptr> hosts;
        cout << "    " << Host::printHeader() << endl;
        while( count < blockSize ) {
            Host::Ptr peer(new Host);
            count += peer->unpack( ptr+count, conn->getSwapEndian() );
            hosts.push_back(peer);
        }
        if( sorting == by_runtime ) std::sort( hosts.begin()+1, hosts.end(), host_by_runtime );
        int nHosts = hosts.size();
        bool addComma(false);
        if ( ns != 0 ) {
            if ( ns < nHosts ) addComma = true;
            nHosts = min(ns,nHosts);
        }
        for( int i=0; i<nHosts; ++i ) {
            if ( i == 0 ) {
                cout << "    ";
            } else if( i == nHosts-1 ) {
                i = hosts.size()-1;
                if( addComma ) cout << "   :" << endl;
            }
            if ( i ) cout << alignRight(to_string(i),4);
            cout << hosts[i]->print() << endl;
        }
    } catch ( const exception& e) {
        cerr << "printPeerList: Exception caught while parsing block: " << e.what() << endl;
    }
    if( count != blockSize ) {
        cerr << "printPeerList: Parsing of datablock failed, count = " << count << "  blockSize = " << blockSize << "  bytes." << endl;
    }
}

int main( int argc, char *argv[] ) {

    bpo::variables_map vm;
    bpo::options_description programOptions = getOptions();

    bpo::options_description& allOptions = Application::parseCmdLine( argc, argv, vm, &programOptions );

    // load matched environment variables according to the getOptionName() above.
    bpo::store( bpo::parse_environment( allOptions, environmentMap ), vm );

#if BOOST_VERSION > 104800  // TODO check which version notify appears in
    vm.notify();
#endif

    bool loop = vm.count( "time" );
    int sorting = by_name;
    try {
        boost::asio::io_service ioservice;
        auto conn = TcpConnection::newPtr( ioservice );
        conn->connect( vm["master"].as<string>(), vm["port"].as<string>() );

        if( conn->socket().is_open() ) {
            Host::HostInfo me;
            Host master;
            uint8_t cmd = CMD_CONNECT;
            boost::asio::write(conn->socket(),boost::asio::buffer(&cmd,1));
            boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));
            if( cmd == CMD_AUTH ) {
                // implement
            }
            if( cmd == CMD_CFG ) {  // handshake requested
                *conn << me;
                *conn >> master.info;
                boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));       // ok or err
            }
            if( cmd != CMD_OK ) {
                cerr << "Handshake with server failed." << endl;
                return EXIT_FAILURE;
            }
            int maxItems = 0;
            int hasFlags = vm.count( "jobs" ) + vm.count( "slaves" );
            if( vm.count( "count" ) ) maxItems = vm["count"].as<int>();
            int maxJobs = maxItems;
            if( vm.count( "jobs" ) ) maxJobs = vm["jobs"].as<int>();
            if( vm.count( "jobs" ) || !hasFlags ) printJobList( conn, maxJobs );
            
            if( hasFlags != 1 ) cout << endl;   // separator
            
            int maxSlaves = maxItems;
            if( vm.count( "runtime" ) ) sorting = by_runtime;
            if( vm.count( "slaves" ) ) maxSlaves = vm["slaves"].as<int>();
            if( vm.count( "slaves" ) || !hasFlags ) printPeerList( conn, maxSlaves, sorting );
            while( loop ) {
                sleep(vm["time"].as<int>());
                if( vm.count( "jobs" ) || !hasFlags ) printJobList( conn, maxJobs );
                if( hasFlags != 1 ) cout << endl;
                if( vm.count( "slaves" ) || !hasFlags ) printPeerList( conn, maxSlaves, sorting );
            }
        } else {
            cout << "Connection failed: " << vm["master"].as<string>() << ":" << vm["port"].as<string>() << endl;
        }
    }
    catch( const exception &e ) {
        cerr << "Uncaught exception (fatal): " << e.what() << endl;
    }
    return EXIT_SUCCESS;

}

