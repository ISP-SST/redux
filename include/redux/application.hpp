#ifndef REDUX_APPLICATION_HPP
#define REDUX_APPLICATION_HPP

#include <string>

#include <boost/noncopyable.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace redux {

    class Application : private boost::noncopyable {

    public:
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

        Application( po::variables_map& vm );
        virtual ~Application( void );

        static void getOptions( po::options_description& options, const std::string& );
        static po::options_description& getOptionsDescription( const std::string& name = "application" );

        typedef void ( parserFunction )( po::options_description&, po::variables_map& );

        static std::pair<std::string, std::string> customParser( const std::string& s );
        static po::options_description& parseCmdLine( int argc, const char* const argv[], po::variables_map& vm,
                                                      po::options_description* programOptions = nullptr, parserFunction parser = nullptr );

        static void checkGeneralOptions( po::options_description& desc, po::variables_map& vm );

        int run( void );


        static std::string executableName;

    private:

        void reset( void );
        void kill( void );

        bool dispatch( void );

        std::string getName( void ) const;

        std::string m_Name;
        std::string m_IniFile;

        volatile bool m_ShouldStop;
        volatile bool m_ShouldRestart;

    };


}

#endif // REDUX_APPLICATION_HPP
