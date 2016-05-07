#ifndef REDUX_APPLICATION_HPP
#define REDUX_APPLICATION_HPP

#include "redux/logger.hpp"

#include <string>

#include <boost/noncopyable.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

namespace bpt = boost::property_tree;
namespace bpo = boost::program_options;

namespace redux {

    class Application : private boost::noncopyable {
        
    public:
        enum RunMode { LOOP=0, EXIT, RESET };

        struct KillException : public std::exception {
            KillException( void ) {}
            virtual const char *what( void ) const throw() {
                return "application killed";
            }
        };

        struct ResetException : public std::exception {
            ResetException( void ) {}
            virtual const char *what( void ) const throw() {
                return "application reset";
            }
        };

        Application( bpo::variables_map& vm, RunMode=EXIT );
        virtual ~Application( void );

        static void getOptions( bpo::options_description& options, const std::string& );
        static bpo::options_description& getOptionsDescription( const std::string& name = "application" );

        typedef void ( parserFunction )( bpo::options_description&, bpo::variables_map& );

        static std::pair<std::string, std::string> appCmdParser( const std::string& s );
        static bpo::options_description& parseCmdLine( int argc, const char* const argv[], bpo::variables_map& vm,
                                                      bpo::options_description* programOptions = nullptr,
                                                      bpo::positional_options_description * positionalOptions = nullptr,
                                                      parserFunction customParser = nullptr );

        static void checkGeneralOptions( bpo::options_description& desc, bpo::variables_map& vm );

        int run( void );                    //!< application entry-point, basically just a loop that calls \c doWork()
        
        virtual void reset( void ) { runMode = RESET; };
        virtual void stop( void ) { runMode = EXIT; };

        std::string getName( void ) const;
        static std::string executableName;

    protected:
        /*! Application main method, default application returns false immediately (i.e. exits).
         *   Overload doWork with something interesting. 
         */
        virtual bool doWork( void ) { return false; };
        
        volatile RunMode runMode;
        
        int returnValue;

        std::string applicationName;
        std::string settingsFile;

        bpt::ptree propTree;
        redux::Logger logger;
        
        
    };


}

#endif // REDUX_APPLICATION_HPP
