#include "redux/application.hpp"
#include "redux/job.hpp"
#include "redux/logger.hpp"
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

#define lg Logger::lg
namespace {

    const string logChannel = "jstat";

    // define options specific to this binary
    bpo::options_description getOptions( void ) {

        bpo::options_description options( "Program Options" );
        options.add_options()
        ( "master,m", bpo::value<string>()->default_value( "localhost" ), "Hostname/IP of the master" )
        ( "port,p", bpo::value<string>()->default_value( "30000" ), "Port to use, either when connecting to the master." )
        ( "jobs,j", "Get joblist" )
        ( "slaves,s", "Get slavelist" )
        ( "count,c", bpo::value<int>()->implicit_value( 0 ), "List only n first slaves/jobs, and the last." )
        ( "time,t", bpo::value<int>()->implicit_value( 1 ), "Loop and display list every (n) seconds" )
        ;

        return options;
    }

    // define environment variables to use as defaults if the corresponding command-line option is not specified
    string environmentMap( const string &envName ) {

        static map<string, string> vmap;
        if( vmap.empty() ) {
            vmap["RDX_VERBOSITY"] = "verbosity";  // For debugging this might be convenient.
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
        LOG_ERR << "printJobList: Exception caught while parsing block: " << e.what();
    }
    if( count != blockSize ) {
        LOG_ERR << "printJobList: Parsing of datablock failed,  count = " << count << "   blockSize = " << blockSize << "  bytes.";
    }

}

void printPeerList( TcpConnection::Ptr conn, int ns ) {

    uint8_t cmd = CMD_PSTAT;
    boost::asio::write(conn->socket(),boost::asio::buffer(&cmd,1));

    uint64_t blockSize;
    shared_ptr<char> buf = conn->receiveBlock( blockSize );

    if ( !blockSize ) return;

    const char* ptr = buf.get();
    uint64_t count(0);
    int hostCount(0);
    try {
        Host peer;
        cout << "    " << peer.printHeader() << endl;
        while( count < blockSize ) {
            count += peer.unpack(ptr+count,conn->getSwapEndian());
            if ( hostCount == 0 ) cout << "    " << peer.print() << endl;
            else if ( ns==0 || hostCount < ns) cout << alignRight(to_string(hostCount),4) << peer.print() << endl;
            hostCount++;
        }
        if (ns!=0 && hostCount > ns) cout << "   :" << endl << alignRight(to_string(hostCount),4) << peer.print() << endl;
    } catch ( const exception& e) {
        LOG_ERR << "printPeerList: Exception caught while parsing block: " << e.what();
    }
    if( count != blockSize ) {
        LOG_ERR << "printPeerList: Parsing of datablock failed, count = " << count << "  blockSize = " << blockSize << "  bytes.";
    }
}

int main( int argc, char *argv[] ) {

    bpo::variables_map vm;
    bpo::options_description programOptions = getOptions();

    bpo::options_description& allOptions = Application::parseCmdLine( argc, argv, vm, &programOptions );

    // load matched environment variables according to the getOptionName() above.
    bpo::store( bpo::parse_environment( allOptions, environmentMap ), vm );
    vm.notify();

    if( vm.count( "help" ) ) {
        cout << allOptions << endl;
        return EXIT_SUCCESS;
    }

    bool loop = vm.count( "time" );
    try {
        Logger logger( vm );
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
                LOG_ERR << "Handshake with server failed.";
                return EXIT_FAILURE;
            }
            int c = 0;
            int js = vm.count( "jobs" ) + vm.count( "slaves" );
            if( vm.count( "count" ) ) c = vm["count"].as<int>();
            if( vm.count( "jobs" ) || !js ) printJobList( conn, c );
            if( js != 1 ) cout << endl;
            if( vm.count( "slaves" ) || !js ) printPeerList( conn, c );
            while( loop ) {
                sleep(vm["time"].as<int>());
                if( vm.count( "jobs" ) || !js ) printJobList( conn, c );
                if( js != 1 ) cout << endl;
                if( vm.count( "slaves" ) || !js ) printPeerList( conn, c );
            }
        }
    }
    catch( const exception &e ) {
        cerr << "Uncaught exception (fatal): " << e.what() << endl;
    }
    return EXIT_SUCCESS;

}

