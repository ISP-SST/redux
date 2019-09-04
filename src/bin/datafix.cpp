#define FIX_FITS_201908

#include "redux/file/filefits.hpp"
#include "redux/util/array.hpp"
#include "redux/util/ricecompress.hpp"
#include "redux/util/recursepath.hpp"
#include "redux/util/stringutil.hpp"

#include "fitswriter.hpp"

#include <iostream>
#include <mutex>
#include <string>
#include <thread>

#include <sys/types.h>
#include <sys/stat.h>

#include <boost/asio.hpp>
#include <boost/program_options.hpp>
#include <boost/thread/thread.hpp>
namespace bpo = boost::program_options;

using namespace redux;
using namespace redux::file;
using namespace redux::util;
using namespace std;

/*

*** Verify it's a broken FITS!
*** Verify that running on good file generates identical outfile

*/

namespace {
    
    boost::asio::io_service ioService;
    boost::thread_group pool;
    
    int nThreads = std::thread::hardware_concurrency();
    bool doIt(false);
    bool running(false);
    bool verbose(false);
    
    bfs::path out_base;        // output directory, either supplied or generated from in_base
    bfs::path in_base;         // supplied in-path, the output name will be the subpath under this directory appended to the out_base
    
    mutex gmtx;
    atomic<int> nCurrentThreads(0);
    atomic<int> activeCount(0);
    struct acWrapper {
        acWrapper( atomic<int>& c ) : cnt(c) { cnt++; };
        ~acWrapper() { cnt--; };
        atomic<int>& cnt;
    };
    
    void threadLoop( void ) {
        acWrapper aw(nCurrentThreads);
        while( running ) {
            try {
                boost::this_thread::interruption_point();
                ioService.run();
                std::this_thread::sleep_for( std::chrono::milliseconds(1000) );
            } catch( const boost::thread_interrupted& ) {
                lock_guard<mutex> lock(gmtx);
                if( verbose ) cout << "Thread interrupted. " << endl;
                break;
            } catch( exception& e ) {
                lock_guard<mutex> lock(gmtx);
                cerr << "Exception in thread: " << e.what() << endl;
            } catch( ... ) {
                lock_guard<mutex> lock(gmtx);
                cerr << "Unhandled exception in thread." << endl;
            }
        }
    }
    
    void fix_perms( const bfs::path in_file, const bfs::path out_file ) {
        if( bfs::exists(in_file) && bfs::exists(out_file) ) {
            if( access( out_file.string().c_str(), W_OK) ) {
                perror( out_file.string().c_str() );
                return;
            }
            lock_guard<mutex> lock(gmtx);
            struct stat in_stat;
            struct stat out_stat;
            if( stat( in_file.string().c_str(), &in_stat ) < 0 ) {
                perror( in_file.string().c_str() );
                return;
            }
            if( stat( out_file.string().c_str(), &out_stat ) < 0 ) {
                perror( out_file.string().c_str() );
                return;
            }
            
            if( (out_stat.st_uid != in_stat.st_uid) || (out_stat.st_gid != in_stat.st_gid) ) {
                if( chown( out_file.string().c_str(), in_stat.st_uid, in_stat.st_gid ) < 0 ) {
                    perror( out_file.string().c_str() );
                    return;
                }
            }
            if( (out_stat.st_mode != in_stat.st_mode) ) {
                if( chmod( out_file.string().c_str(), in_stat.st_mode ) < 0 ) {
                    perror( out_file.string().c_str() );
                    return;
                }
            }
            struct timespec file_times[] = { in_stat.st_atim, in_stat.st_mtim };
            if( utimensat( AT_FDCWD, out_file.string().c_str(), file_times, 0 ) < 0 ) {
                perror( out_file.string().c_str() );
                return;
            }
        }
    }
    
    bfs::path make_out_path( const bfs::path& p ) {
        bfs::path subpath = bfs::relative( p ,in_base );
        return (out_base / subpath);
    }
    
    void fix_file( const bfs::path in_file, const bfs::path out_file ) {
        
        acWrapper aw(activeCount);

        if( doIt ) {
            
            Array<int16_t> data;
            shared_ptr<Fits> hdr;
            try {
                Fits::read( in_file.string(), data, hdr );
            } catch( std::exception& e ) {
                lock_guard<mutex> lock(gmtx);
                cerr << "Failed to read file: " << in_file << "  reason: " << e.what() << endl;
                cerr << "    Perhaps this is not a corrupt file?  It will be SKIPPED!" << endl;
                return;
            }
            
            if( verbose ) {
                lock_guard<mutex> lock(gmtx);
                cout << "Fixing file: " << "  " << in_file << " -> " << out_file << endl;
            }
            FitsWriter fw( out_file.string(), 1, true );
            if( hdr ) {
                fw.save_meta( hdr->primaryHDU.cards );
            }
            fw.save( data, hdr );
            fw.wait();
            
            fix_perms( in_file, out_file );
            
            bfs::path d1 = in_file;
            bfs::path d2 = out_file;

            d1 = d1.parent_path();
            d2 = d2.parent_path();

            while( !d2.empty() && (d2 != out_base) ) {      // just to make sure all permissions/times are fixed
                fix_perms( d1, d2 );
                d2 = d2.parent_path();
            }
            
        } else if( verbose ) {
            lock_guard<mutex> lock(gmtx);
            cout << "NOT fixing file: " << "  " << in_file << " -> " << out_file << endl;
        }

        std::this_thread::sleep_for( std::chrono::seconds(2) );
    }

    
    bool process_path( const bfs::path& p ) {
        
        if( !bfs::exists( p ) ) {       // bailout if the path does not exist
            return false;
        }
        
        acWrapper aw(activeCount);
        bfs::path out_file = make_out_path( p );
        if( bfs::is_directory( p ) ) {  // Just keep recursing the directory, but don't do any operations on the directory
            return true;
        }
        
        // Here we process files.
        if( p.extension() == ".fits" ) {
            ioService.post( boost::bind( fix_file, p, out_file ) );
        }
        return false;
        
    }
    
    // define options
    bpo::options_description getOptions( void ) {

        bpo::options_description options( "datafix options" );
        options.add_options()
            ( "for-real,f", "Actually do it, if this flag is not specified it will just be a dry-run." )
            ( "in-dir,i", bpo::value<string>()->required(), "Input-path to recursively process." )
            ( "out-dir,o", bpo::value<string>(), "Input-path where to store the results."
                " If left blank, the Input-path will be used, but with \"_FIXED\" appended." )
            ( "threads,t", bpo::value<uint16_t>(), "Max number of threads to use.")
            ( "verbose,v", "Print info. (default is to just print errors/warnings.")
        ;

        return options;
    }


}


int main( int argc, char *argv[] ) {

    try {

        bpo::variables_map vm;
        bpo::options_description programOptions = getOptions();

        bfs::path tmpPath = bfs::path(string(argv[0])).filename();
        
        bpo::command_line_parser parser( argc, argv );
        parser.options( programOptions );
        bpo::store( parser.run(), vm );

        vm.notify();

        if( vm.count("for-real") ) {
            doIt = true;
        }
        
        if( vm.count("verbose") ) {
            verbose = true;
        }
        
        if( vm.count("threads") ) {
            if( vm["in-dir"].as<int>() > 0 ) {
                nThreads = vm["in-dir"].as<int>();
            }
        }
        
        in_base = vm["in-dir"].as<string>();
        if( !bfs::exists(in_base) ) {
            cerr << "Input-path does not exist: " << in_base << endl;
        }
        in_base = bfs::canonical( in_base, "/" );               // verifies the path exists
        

        if( vm.count("out-dir") ) {
            out_base = vm["out-dir"].as<string>();
        } else {        // construct from in_base
            out_base = in_base.parent_path();
            out_base /= in_base.leaf().string()+"_FIXED";
        }
        
        if( doIt ) {
            maybeCreateDir(out_base);
            out_base = bfs::canonical( out_base, "/" );     // verifies the path exists
        }
        
        running = true;
        boost::asio::io_service::work work(ioService);
        for( int i=0; i < nThreads; ++i ) {
            pool.create_thread( std::bind( &threadLoop ) );
        }
        
        std::this_thread::sleep_for( std::chrono::seconds(3) );
        RecursePath( in_base, process_path );

        ioService.post( [&](void){
            while( running ) {
                std::this_thread::sleep_for( std::chrono::seconds(1) );
                int ac = activeCount;
                if( !ac ) {
                    break;
                }
            }
            running = false;
            ioService.stop();
        });
        
        pool.join_all();

    } catch( const exception& e ) {
        cerr << "Failed to parse command-line. Reason: " << e.what() << endl;
    }

    
    return EXIT_SUCCESS;

}

