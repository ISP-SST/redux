
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/logging/logger.hpp"
#include "redux/util/arrayutil.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/test/unit_test.hpp>


using namespace redux::logging;
using namespace redux::momfbd;
using namespace redux::util;
 
using namespace std;
namespace bpt = boost::property_tree;


namespace testsuite {

    namespace momfbd {
        
        namespace {
//            MomfbdJob mjob;
            /*MomfbdJob mjob;     // unit-wide job for testing.
            
            
            int logInit(void) {
                redux::logging::Logger& lg = mjob.getLogger();
                lg.setLevel( 5 );
                lg.setContext( "testsuite" );
                lg.addStream( cout, 0, 1 );
                return 1;
            }

            
            int initOnce(void) {
                static int init_done = logInit();
                return init_done;
            }*/
        }

        void cfgTest( void ) {
            return;
            //int dummy RDX_UNUSED = initOnce();
            MomfbdJob mjob;
            redux::logging::Logger& logger = mjob.getLogger();
            logger.setLevel( 8 );
            logger.setContext( "testsuite" );
            logger.addStream( cout, 0, 1 );


            bpt::ptree tree;

            
            //return;
            ChannelCfg ccfg;
            {
                // Set some non-default values for ChannelCfg
                ccfg.noiseFudge = ccfg.weight = ccfg.rotationAngle = 138;
                ccfg.alignClip = {13,34,14,35};
                ccfg.borderClip = 139;
                ccfg.incomplete = true;
                ccfg.subImagePosX = {88,99,111,122};
                ccfg.subImagePosY = {89,98,112,121};
                ccfg.imageDataDir = "path/to/data/";
                ccfg.imageTemplate = "imgname_with_number_%07d.ext";
                ccfg.darkTemplate = "darkname_with_number_%07d.ext";
                ccfg.gainFile = "gainFile.ext";
                ccfg.responseFile = "responseFile.ext";
                ccfg.backgainFile = "backgainFile.ext";
                ccfg.psfFile = "psfFile.ext";
                ccfg.mmFile = "mmFile.ext";
                ccfg.mmRow = ccfg.imageNumberOffset = 141;
                ccfg.mmWidth = 4;
                ccfg.xOffsetFile = "xOffsetFile.ext";
                ccfg.yOffsetFile = "yOffsetFile.ext";
                ccfg.fileNumbers = {23,45};
                ccfg.waveFrontList = {2,4,6};
                ccfg.darkNumbers = {145,88};
                ccfg.stokesWeights = {0.99,0.11,0.12,0.13};
                // export settings as config file and parse/compare
                tree.clear();
                ccfg.getProperties(tree);
                stringstream cfg;
                bpt::write_info( cfg, tree );   // test that the ptree is valid
                bpt::read_info( cfg, tree );
                ChannelCfg tmp;                 // default values
                BOOST_CHECK( !(ccfg == tmp) );
                tmp.parseProperties( tree, mjob.getLogger() );
                BOOST_CHECK( ccfg == tmp );
                // pack/unpack and compare
                auto buf = sharedArray<char>( ccfg.size() );
                char* ptr = buf.get();
                uint64_t count = ccfg.pack( ptr );
                BOOST_CHECK_EQUAL( count, ccfg.size() );
                tmp = ChannelCfg();             // reset default values and test unpack
                count = tmp.unpack( ptr, false );
                BOOST_CHECK_EQUAL( count, ccfg.size() );
                BOOST_CHECK_EQUAL( count, tmp.size() );
                BOOST_CHECK( ccfg == tmp );
                tmp.borderClip = 15;       // modify
                BOOST_CHECK( !(ccfg == tmp) );
                tmp = ccfg;                     // assign
                BOOST_CHECK( ccfg == tmp );
            }
return;
            ObjectCfg ocfg;
            {
                // Set some non-default values for ObjectCfg
                ocfg.telescopeF = ocfg.arcSecsPerPixel = ocfg.pixelSize = 137;
                ocfg.maxLocalShift = ocfg.minimumOverlap = 139;
                ocfg.patchSize = ocfg.pupilPixels = 140;
                ocfg.saveMask = ocfg.wavelength = 58;
                ocfg.pupilFile = "pupilfile.ext";
                ocfg.modeFile = "modefile.ext";
                ocfg.outputFileName = "filename.ext";
                // export settings as config file and parse/compare
                tree.clear();
                ocfg.getProperties(tree);
                stringstream cfg;
                bpt::write_info( cfg, tree );   // test that the ptree is valid
                bpt::read_info( cfg, tree );
                ObjectCfg tmp;                  // default values
                BOOST_CHECK( !(ocfg == tmp) );
                tmp.parseProperties( tree, mjob.getLogger() );
                BOOST_CHECK( ocfg == tmp );
                // pack/unpack and compare
                auto buf = sharedArray<char>( ocfg.size() );
                char* ptr = buf.get();
                uint64_t count = ocfg.pack( ptr );
                BOOST_CHECK_EQUAL( count, ocfg.size() );
                tmp = ObjectCfg();              // reset default values and test unpack
                count = tmp.unpack( ptr, false );
                BOOST_CHECK_EQUAL( count, ocfg.size() );
                BOOST_CHECK_EQUAL( count, tmp.size() );
                BOOST_CHECK( ocfg == tmp );
                BOOST_CHECK( ChannelCfg() == tmp );
                tmp.patchSize = 15;              // modify
                BOOST_CHECK( !(ocfg == tmp) );
                tmp = ocfg;                     // assign
                BOOST_CHECK( ocfg == tmp );
                tmp.noiseFudge = 15;            // modify some ChannelCfg parameter
                BOOST_CHECK( !(ocfg == tmp) );
                tmp = ChannelCfg();             // assign default ChannelCfg
                BOOST_CHECK( ocfg == tmp );
            }

            GlobalCfg gcfg;
            {
                // Set some non-default values for GlobalCfg
                gcfg.runFlags = 4095;
                gcfg.modeBasis = KARHUNEN_LOEVE;
                gcfg.klMinMode = gcfg.klMaxMode = gcfg.klCutoff = gcfg.nInitialModes = gcfg.nModeIncrement = 118;
                gcfg.telescopeD = gcfg.minIterations = gcfg.maxIterations = gcfg.targetIterations = 119;
                gcfg.fillpixMethod = gcfg.getstepMethod = 3;
                gcfg.gradientMethod = 2;
                gcfg.badPixelThreshold = gcfg.FTOL = gcfg.EPS = gcfg.reg_alpha = gcfg.sequenceNumber = 117;
                gcfg.outputFileType = gcfg.outputDataType = 1;
                gcfg.observationTime = "time";
                gcfg.observationDate = "date";
                gcfg.tmpDataDir = "/datadir/";
                gcfg.outputFiles = {"file1","file2"};
                // export settings as config file and parse/compare
                tree.clear();
                gcfg.getProperties(tree);
                stringstream cfg;
                bpt::write_info( cfg, tree );   // test that the ptree is valid
                bpt::read_info( cfg, tree );
                GlobalCfg tmp;                  // default values
                BOOST_CHECK( !(gcfg == tmp) );
                tmp.parseProperties( tree, mjob.getLogger() );
                BOOST_CHECK( gcfg == tmp );
                // pack/unpack and compare
                auto buf = sharedArray<char>( gcfg.size() );
                char* ptr = buf.get();
                uint64_t count = gcfg.pack( ptr );
                BOOST_CHECK_EQUAL( count, gcfg.size() );
                tmp = GlobalCfg();              // reset default values and test unpack
                count = tmp.unpack( ptr, false );
                BOOST_CHECK_EQUAL( count, gcfg.size() );
                BOOST_CHECK_EQUAL( count, tmp.size() );
                BOOST_CHECK( gcfg == tmp );
                BOOST_CHECK( ObjectCfg() == tmp );
                BOOST_CHECK( ChannelCfg() == tmp );
                tmp.telescopeD = 15;             // modify
                BOOST_CHECK( !(gcfg == tmp) );
                tmp = gcfg;                     // assign
                tmp.saveMask = 15;              // modify ObjectCfg
                BOOST_CHECK( !(gcfg == tmp) );
                tmp = ObjectCfg();              // assign default ObjectCfg
                BOOST_CHECK( gcfg == tmp );
                tmp.noiseFudge = 15;            // modify ChannelCfg
                BOOST_CHECK( !(gcfg == tmp) );
                tmp = ChannelCfg();             // assign default ChannelCfg
                BOOST_CHECK( gcfg == tmp );
            }

            BOOST_CHECK( !(ccfg == ocfg) );     // using ChannelCfg::operator==() ->  should not be equal
            ocfg = ccfg;
            BOOST_CHECK( ccfg == ocfg );        // now they are equal

            BOOST_CHECK( !(ocfg == gcfg) );     // using ObjectCfg::operator==()  ->  should not be equal
            gcfg = ocfg;
            BOOST_CHECK( ocfg == gcfg );        // now they are equal

        //    MomfbdJob mjob;
            mjob = gcfg;                        // assign cfg to job-class
            BOOST_CHECK( gcfg == mjob );        // using GlobalCfg::operator==()  ->  should be equal

        //         gcfg.getProperties(tree);
        //         bpt::write_info( cout<<endl, tree );


            return;
        //         bpo::variables_map vm;
        //         bpt::ptree tree;

        //         boost::any v(0);
        //         bpo::variable_value vv(v,false);
        //         vv.value() = 8;
        //         vm.insert( std::make_pair(std::string("verbosity"), vv ) );
        //         //vm.at("verbosity") = bpo::variable_value(0);
        //        redux::Logger logger( vm );
        //        logger.addNullLog();
        // logger.addStream( cerr, redux::Logger::getDefaultMask() );
            /*
                    stringstream cfg;
                    cfg << "object { }";

                    bpt::read_info( cfg, tree );
                    mjob.parsePropertyTree( vm, tree );

                    bpt::write_info( cout<<endl, tree );*/

        }


        using namespace boost::unit_test;
        void add_config_tests( test_suite* ts ) {

            ts->add( BOOST_TEST_CASE_NAME( &cfgTest, "Config" ) );

        }

    }

}
