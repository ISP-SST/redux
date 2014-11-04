#include "redux/application.hpp"

#include "redux/gui/reduxmw.hpp"
#include "redux/gui/reduxwidget.hpp"
#include "redux/debugjob.hpp"
#include "redux/momfbd/momfbdjob.hpp"

#include <boost/program_options.hpp>
namespace bpo = boost::program_options;

using namespace redux::gui;
using namespace redux;
using namespace std;

#define lg Logger::lg
namespace {

    const string thisChannel = "redux-gui";

    // define options specific to this binary
    bpo::options_description getOptions( void ) {

        bpo::options_description options( "REDUXd Options" );
        options.add_options()
        ( "master,m", po::value<string>()->default_value( "" ), "name or ip of master." )
        ( "port,p", bpo::value<uint16_t>()->default_value( 30000 ), "Port to listen on."
          " The environment variable REDUX_PORT will be used instead of the default if this option"
          " is not specified." )
        ;

        return options;
    }

    // define environment variables to use as defaults if the corresponding command-line option is not specified
    string environmentMap( const string &envName ) {

        static map<string, string> vmap;
        if( vmap.empty() ) {
            vmap["REDUX_VERBOSITY"] = "verbosity";  // For debugging this might be convenient.
            vmap["REDUX_PORT"] = "port";            // This means the environment variable REDUX_PORT will override the
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


int main( int argc, char *argv[] ) {

    bpo::variables_map vm;
    bpo::options_description programOptions = getOptions();

    bpo::options_description& allOptions = Application::parseCmdLine( argc, argv, vm, &programOptions );

    // load matched environment variables according to the getOptionName() above.
    bpo::store( bpo::parse_environment( allOptions, environmentMap ), vm );
    vm.notify();

    try {

        QApplication app( argc, argv );

        QCoreApplication::setOrganizationName( "redux" );
        QCoreApplication::setOrganizationDomain( "redux" );
        QCoreApplication::setApplicationName( "redux-gui" );
        QCoreApplication::setApplicationVersion( "alpha" );

        ReduxMW mw;
        mw.setWindowTitle( "redux-gui" );
        ReduxWidget* w = new ReduxWidget(&mw);
        mw.setCentralWidget ( w );
        mw.show();
        return app.exec();

    }
    catch( const exception &e ) {
        LOG_ERR << "Uncaught exception (fatal): " << e.what();
    }

    return EXIT_SUCCESS;

}

