#include "redux/application.hpp"
#include "redux/file/fileio.hpp"
#include "redux/logging/logger.hpp"
#include "redux/debugjob.hpp"
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/network/tcpconnection.hpp"
#include "redux/network/host.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/util/stringutil.hpp"

#include <sstream>
#include <thread>

#include <boost/bind/bind.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/program_options.hpp>

using namespace redux::file;
using namespace redux::logging;
using namespace redux::util;
using namespace redux::network;
using namespace redux;

using namespace std;


namespace {

    // define options specific to this binary
    bpo::options_description getOptions( void ) {

        bpo::options_description options( "Program Options" );
        options.add_options()
        ( "master,m", bpo::value<string>()->default_value( "localhost" ),
          "Hostname/IP of a master to connect to."
          " The environment variable RDX_HOST can be used to override the default value." )
        ( "port,p", bpo::value<string>()->default_value( "30000" ),
          "Port to use when connecting to a master."
          " The environment variable RDX_PORT can be used to override the default value." )
        ( "priority", bpo::value<int>()->default_value( 10 ), "Job priority" )
        ( "reg_alpha", bpo::value<float>(), "REG_ALPHA override" )
        ( "force,f", "Overwrite output file if exists" )
        ( "no-swap", "disable swap mode: keeping temporary files in memory." )
        ( "config,c", bpo::value< vector<string> >()->multitoken(), "Configuration file(s) to process." )
        ( "name", bpo::value<string>(), "Name to use for the supplied configurations." )
        ( "user", bpo::value<string>(), "Username to use for the supplied configurations (default is to use current user)." )
        ( "simxy", bpo::value<string>(), "(x,y) coordinate[s] of subimages to restore" )
        ( "simx", bpo::value<string>()->implicit_value(""), "x coordinate[s] of subimages to restore" )
        ( "simy", bpo::value<string>()->implicit_value(""), "y coordinate[s] of subimages to restore" )
        ( "imgn,n", bpo::value<string>(), "Image numbers" )
//        ( "sequence", bpo::value<string>(), "sequence number to insert in filename template." )
        ( "print,P", "(debug) print the parsed configuration to console and exit without uploading." )
        ( "no-check", "Don't verify the configuration." )
        ( "trace", "Generate trace objects." )
        ( "no-trace", "Don't generate trace objects." )
        ( "old-ns", "Don't use the new way of generating the basis for the constraint nullspace." )
        ( "output-dir,O", bpo::value<string>(), "Output directory. If left blank, the current directory is used.")
        ( "output-file,o", bpo::value<string>(), "Output file base names." )
        ( "init", bpo::value<string>()->implicit_value(""), "File with initial values for alpha. If no argument is provided, the output is used.")
        ;

        return options;
    }

    // define environment variables to use as defaults if the corresponding command-line option is not specified
    string environmentMap( const string &envName ) {

        static map<string, string> vmap;
        if( vmap.empty() ) {
            vmap["RDX_VERBOSITY"] = "verbosity";  // For debugging this might be convenient.
            vmap["RDX_HOST"] = "master";        // If it exists, it will override the default value (localhost) above
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
}


void uploadJobs(TcpConnection::Ptr conn, vector<Job::JobPtr>& jobs, int prio, Logger& logger) {
    
    const Host& me = Host::myInfo();
    Host::HostInfo master;
    uint8_t cmd = CMD_CONNECT;
    try {
        boost::asio::write(conn->socket(),boost::asio::buffer(&cmd,1));
        boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));

        if( cmd == CMD_AUTH ) {
            // implement
        }
        if( cmd == CMD_CFG ) {  // handshake requested
            *conn << me.info;
            *conn >> master;
            boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));       // ok or err
        }
        if( cmd != CMD_OK ) {
            LOG_ERR << "Handshake with server failed." << ende;
            return;
        }
       
        // all ok, upload jobs.
        uint64_t jobsSize(0);
        for( auto &job: jobs ) {
            if( job ) {
                job->info.priority = prio;
                jobsSize += job->size();
            }
        }
        
        shared_ptr<char> buf( new char[jobsSize+sizeof(uint64_t)+1], []( char* p ){ delete[] p; } );
        char* ptr = buf.get()+sizeof(uint64_t)+1;
        uint64_t packedBytes(0);
        for( auto &job: jobs ) {
            if( job ) {
                packedBytes += job->pack(ptr+packedBytes);
            }
        }
        
        uint64_t totalSize = packedBytes+sizeof(uint64_t)+1;                    // cmd & blockSize will be sent before block
        ptr = buf.get();
        pack( ptr, CMD_ADD_JOB);
        pack( ptr+1, packedBytes );
        
        //conn->asyncWrite( buf, packedBytes );
        conn->syncWrite( buf.get(), totalSize );
        //boost::asio::write(conn->socket(),boost::asio::buffer(buf.get(),packedBytes));

        boost::asio::read(conn->socket(),boost::asio::buffer(&cmd,1));

        bool swap_endian = (me.info.littleEndian != master.littleEndian);
        ptr = buf.get();

        if( cmd == CMD_OK ) {
            uint64_t received = boost::asio::read( conn->socket(), boost::asio::buffer( ptr, sizeof(uint64_t) ) );
            if( received == sizeof(uint64_t) ) {
                uint64_t count;
                unpack(ptr, count, swap_endian);
                size_t thisSize = count*sizeof(uint64_t);
                received = boost::asio::read( conn->socket(), boost::asio::buffer( ptr, thisSize ) );
                if( received == thisSize ) {
                    if( count ) LOG << "Upload of " << count << " job(s) completed successfully. " << printArray(reinterpret_cast<size_t*>(buf.get()),count,"IDs") << ende;
                } else {
                    LOG_ERR << "Failed to read job IDs.  received=" << received << " thisSize=" << thisSize << ende;
                }
            } else LOG_ERR << "Failed to read number of job IDs." << ende;
        } else {
            LOG_ERR << "Failure while sending jobs  (server reply = " << (int)cmd << "   " << bitString(cmd) << ")" << ende;
        }
        
        vector<string> messages;
        *conn >> messages;
        if( !messages.empty() ) {
            string msgText = "Server messages:";
            for( auto& msg: messages ) {
                msgText += "\n\t" + msg;
            }
            LOG_WARN << "Error uploading jobs: " << msgText << ende;
        }

    }
    catch( const exception &e ) {
        LOG_ERR << "Error uploading jobs: " << e.what() << ende;
    }
    
    logger.flushAll();

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


string filterOldCfg( const string& filename, const string& jobname, const string& logfile, const string& outputDir ) {
    
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

    replace(text, " ", "");     // spaces not allowed in cfg-values (e.g. DIVERSITY=0.00 mm)
    bool found = replace(text, "channel{", "channel {");
    found |= replace(text, "object{", "object {");
    
    if (found) { // old cfg file, wrap inside momfbd { ... }
        bool newLine(true);
        for( auto& c: text ) {
            if( newLine && c == '=' ) {
                newLine = false;
                c = ' ';
            } else if( c == '\n' ) {
                newLine = true;
            }
        }
        text = "momfbd { \n" + text;
        if( *(text.rbegin()) != '\n') text += "\n";
        text += "NAME " + jobname + "\n";
        text += "LOGFILE " + logfile + "\n";
        text += "OUTPUT_DIR " + outputDir + "\n";
        if( *(text.rbegin()) != '\n') text += "\n";
        text += "}";
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
#if BOOST_VERSION > 104800  // TODO check which version notify appears in
        vm.notify();
#endif
    }
    catch( const exception &e ) {
        cerr << "Error parsing commandline: " << e.what() << endl;// << programOptions << endl;
        return EXIT_FAILURE;
    }


    try {
        
        string globalLog;
        if( vm.count ("log-file") ) {
            vector<string> logFiles = vm["log-file"].as<vector<string>>();
            if( logFiles.size() ) {
                globalLog = logFiles[0];
            }
            if( logFiles.size() > 1 ) {
                cerr << "Only 1 log-file supported at the moment. Using: " << globalLog << endl;
            }
        }

        vm.erase("log-file");       // always log to cout for rsub
        vm.insert( std::make_pair("log-stdout", bpo::variable_value()) );
        
        Logger logger( vm );
        bpt::ptree momfbd;
        
        if( !vm.count ("config") ) {
            LOG_FATAL << "No configuration file supplied." << ende;
            return 0;
        }

        vector<string> files = vm["config"].as<vector<string>>();

        string globalName;
        if( vm.count ("name") ) {
            globalName = vm["name"].as<string>();
        }
        
        bfs::path outputDir = bfs::current_path();
        if( vm.count ("output-dir") ) {
            bfs::path tmpPath = vm["output-dir"].as<string>();
            if( isRelative( tmpPath ) && !outputDir.empty() ) {
                outputDir = outputDir / tmpPath;
            }
        }

        stringstream filteredCfg;
        for( auto it: files ) {
            if( ! bfs::is_regular_file( it ) ) {
                LOG_WARN << "No such file: " << it << ende;
                continue;
            }
            string bn = bfs::basename(it);
            string jobName = globalName.empty() ? bn : globalName;
            string logFile = globalLog.empty() ? bn + ".log" : globalLog;
            string tmpS = filterOldCfg(it, jobName, logFile, outputDir.string());
            filteredCfg.write(tmpS.c_str(),tmpS.size());
        }

        bpt::read_info (filteredCfg , momfbd);
        bool check = (vm.count ("no-check") == 0 && !vm.count ("print"));
        vector<Job::JobPtr> jobs;
        try {
            jobs = Job::parseTree (vm, momfbd, logger, check);
        } catch ( const exception& e ) {
            LOG_ERR << "Error while parsing cfg file(s):\n" << e.what() << ende;
        }
        
        if( jobs.empty() ) {
            LOG_WARN << "No jobs to upload." << ende;
            return EXIT_SUCCESS;
        }
        
        if( vm.count ("reg_alpha") ) {       // FIXME: this is just while testing...
            for( auto & job : jobs ) {
                static_pointer_cast<momfbd::MomfbdJob>(job)->reg_alpha = vm["reg_alpha"].as<float>();
            }
        }
        
        if( vm.count ("user") ) {
            string user = vm["user"].as<string>();
            Host::myInfo().info.user = user;
            for( auto & job : jobs ) {
                job->info.user = user;
            }
        }

        if( vm.count ("print") ) {       // dump configuration to console and exit
            bpt::ptree dump;
            bool showAll = vm.count ("verbose");
            for( auto & job : jobs ) {
                job->getPropertyTree( &dump, showAll );
            }
            bpt::write_info( cout<<endl, dump );
            return EXIT_SUCCESS;
        }
        
        boost::asio::io_service ioservice;
        auto conn = TcpConnection::newPtr(ioservice);
        conn->connect( vm["master"].as<string>(), vm["port"].as<string>() );

        if( conn->socket().is_open() ) {
            std::vector<std::shared_ptr<std::thread> > threads;
            for ( int i=0; i<5; ++i) {
                shared_ptr<thread> t( new thread( boost::bind( &boost::asio::io_service::run, &ioservice ) ) );
                threads.push_back( t );
            }
            int priority = vm["priority"].as<int>();
            shared_ptr<thread> t( new thread( boost::bind( uploadJobs, conn, jobs, priority, std::ref(logger)) ) );
            threads.push_back( t );
            //thread t( boost::bind( &boost::asio::io_service::run, &ioservice ) );
            //thread tt( boost::bind( &boost::asio::io_service::run, &ioservice ) );
            //uploadJobs(conn, jobs, logger);
            //ioservice.run();
            //t.join();
            //tt.join();
            for( auto & it : threads ) {
                it->join();
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

