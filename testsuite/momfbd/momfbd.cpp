// #define private public  //TODO How to fix private access for unit tests ??!!

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/logging/logger.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>


using namespace redux::logging;
using namespace redux::image;
using namespace redux::momfbd;
using namespace redux::util;
using namespace std;

#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;

namespace testsuite {

    namespace momfbd {
        
        void test_structure( void ) {
            
            MomfbdJob job;
            // Logging
            redux::logging::Logger& logger = job.getLogger();
            logger.addStream( cout, 2, 1 );
            logger.setLevel( 0 );
            logger.setContext( "testsuite" );
            
            float EPS = 1E-6;
            // Instances 
            Object::Ptr obj = job.addObject();
            BOOST_REQUIRE_MESSAGE( obj, "Got null Object, can't continue." );
            Channel::Ptr chan = obj->addChannel();
            BOOST_REQUIRE_MESSAGE( chan, "Got null Channel, can't continue." );
            PatchData::Ptr patch( new PatchData(job, 0, 0) );
            ObjectData::Ptr od = patch->getObjectData(0);
            ChannelData::Ptr cd = patch->getChannelData(0,0);
            BOOST_REQUIRE_MESSAGE( od, "Got null ObjectData, can't continue." );
            BOOST_REQUIRE_MESSAGE( cd, "Got null ChannelData, can't continue." );


            {   // test MomfbdJob
                // image size should be 0 by default
                BOOST_CHECK( job.getSmallestImageSize() == 0 );
                chan->setImageSize(123);    // set something
                BOOST_CHECK( job.getSmallestImageSize() == 123 );
                chan->setImageSize(0);      // reset to 0 and see that it clears
                BOOST_CHECK( job.getSmallestImageSize() == 0 );
                
            }
            
            {   // test Object
                
            }
            
            {   // test Channel
                
            }
            
            
            { // test adjustCutout
                
                // verify throw if called with faulty input
                chan->setImageSize(0);
                BOOST_CHECK_THROW( chan->adjustCutout( *cd, nullptr ), std::runtime_error );
                BOOST_CHECK_THROW( chan->adjustCutout( *cd, patch ), std::runtime_error );
                
                obj->patchSize = 1;         // We're just testing the center position here, so keep it trivial.
                obj->maxLocalShift = 0;
                
                chan->setImageSize(1000);       // use image size = (1000,1000)
                patch->position = 503;      // place patch near midpoint
                
                chan->alignMap.clear();
                chan->alignClip.clear();
                
                // test that the corners map as expected
                patch->position = 0;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutPosition, PointI(0,0) );
                chan->setClip({ 999, 0, 999, 0 });
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutPosition, PointI(999,999) );

                // patch near midpoint
                chan->alignClip.clear();        // No clip
                patch->position = 503;
                chan->adjustCutout( *cd, patch );
                PointF diff = cd->exactPatchPosition - patch->position;
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );    // should be exact

                // without flips
                chan->setClip({ 0, 999, 0, 999 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(503,503);     // compare to expected
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                chan->setClip({ 1, 999, 2, 999 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(505,504);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );

                // flip + shift
                chan->setClip({ 1, 998, 997, 2 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(494,504);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                BOOST_CHECK_SMALL( cd->residualOffset.maxAbs(), EPS );
                chan->setClip({ 998, 1, 2, 997 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(505,495);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                BOOST_CHECK_SMALL( cd->residualOffset.maxAbs(), EPS );
                chan->setClip({ 998, 1, 997, 2 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(494,495);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                BOOST_CHECK_SMALL( cd->residualOffset.maxAbs(), EPS );
                chan->setClip({ 997, 1, 998, 2 });              // test with odd clip-sizes
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(495,494);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                BOOST_CHECK_SMALL( cd->residualOffset.maxAbs(), EPS );
                
                // simple maps  (integer shifts and mirrorings)
                chan->setImageSize(1000);       // map needs imgSize to be set
                chan->setClip({});
                chan->setMap({ 1, 0, 0, 0, 1, 0, 0, 0, 1 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(503,503);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                BOOST_CHECK_SMALL( cd->residualOffset.maxAbs(), EPS );
                chan->setMap({ -1, 0, 990, 0, 1, 0, 0, 0, 1 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(503,487);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                BOOST_CHECK_SMALL( cd->residualOffset.maxAbs(), EPS );
                chan->setMap({ 1, 0, 0, 0, -1, 1009, 0, 0, 1 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(506,503);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                BOOST_CHECK_SMALL( cd->residualOffset.maxAbs(), EPS );
                
                // decimal maps
                EPS = 1E-4;
                chan->setMap({ 1, 0, 5.6, 0, 1, 3.8, 0, 0, 1 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointF(506.8,508.6);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                diff = cd->residualOffset - PointF(-0.2,-0.4);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                chan->setMap({ -1, 0, 993.3, 0, 1, 1.2, 0, 0, 1 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointF(504.2,490.3);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                diff = cd->residualOffset - PointF(0.2,0.3);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                chan->setMap({ 1, 0, 0, 0, -1, 997.33, 0, 0, 1 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointF(494.33,503);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                diff = cd->residualOffset - PointF(0.33,0.0);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );

                // offset files
                chan->setMap({});
                chan->setClip({ 0, 999, 0, 999 });
                Image<int16_t> xoff(1000,1000), yoff(1000,1000);
                
                // initial test, check that arrays are copied as references.
                BOOST_CHECK_EQUAL( chan->getOffsetAt( PointI(50,50) ), PointF(0,0) );
                chan->setOffsets( xoff, yoff );
                xoff = -33;
                yoff = 145;
                diff = chan->getOffsetAt( PointI(50,50) ) - PointF(1.45,-0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointF(504.45,502.67);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                diff = cd->residualOffset - PointF(0.45,-0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                chan->setClip({ 999, 0, 0, 999 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointF(504.45,495.67);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                diff = cd->residualOffset - PointF(0.45,-0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                chan->setClip({ 999, 0, 999, 0 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointF(497.45,495.67);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                diff = cd->residualOffset - PointF(0.45,-0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );

                // cutoutRegion
                // only mirroring & shifts
                xoff = 0; yoff = 0;
                obj->patchSize = 1;
                obj->maxLocalShift = 0;
                chan->setClip({});
                patch->position = 0;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(0,0,0,0) );
                patch->position = 999;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(999,999,999,999) );
                chan->setClip({1,998,2,997});
                patch->position = 0;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(2,1,2,1) );
                patch->position = 995;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(997,996,997,996) );
                chan->setClip({998,1,997,2});
                patch->position = 0;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(997,998,997,998) );
                patch->position = 995;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(2,3,2,3) );

                // add offset files with residuals, verify shifts
                xoff = 133; yoff = 255;
                obj->patchSize = 20;
                obj->maxLocalShift = 10;
                patch->position = 20;
                chan->setClip({});
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(3,1,42,40) );
                diff = cd->residualOffset - PointF(-0.45, 0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );

                patch->position = 977;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(960,958,998,997) ); // yh will be reduced by 1 by soft cutout limit
                diff = cd->residualOffset - PointF(-0.45, 0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                chan->setClip({1,998,2,997});
                patch->position = 20;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(5,2,44,41) );
                diff = cd->residualOffset - PointF(-0.45, 0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                patch->position = 973;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(958,955,996,994) );
                diff = cd->residualOffset - PointF(-0.45, 0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                chan->setClip({997,1,995,2});
                patch->position = 23;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(955,955,994,994) );
                diff = cd->residualOffset - PointF(-0.45, 0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                patch->position = 973;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(5,5,44,44) );
                diff = cd->residualOffset - PointF(-0.45, 0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                
                // test clipping near the edges
//logger.setLevel( 8 );
//LOG_MARK



//LOG_MARK
//logger.setLevel( 2 );

                // clip only. With clip, the patch-position is interpreted in clipped coordinates, i.e. (501,501) => 
//                 chan->setClip({ 11, 990, 990, 11 });  // cut 10px on each side, flip X
//                 chan->adjustCutout( *cd, patch );
//                 diff = cd->patchPosition - PointI(511,494);
//                 BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                //BOOST_CHECK_NE( cd, nullptr );
            }
        }
        
        void test_aux( void ) {
            
        }
        
        void add_config_tests( test_suite* ts );    // defined in config.cpp
        void add_data_tests( test_suite* ts );      // defined in data.cpp

        void add_tests( test_suite* ts ) {

            ts->add( BOOST_TEST_CASE_NAME( &test_structure, "Test the overall structure (classes etc.)"  ) );
            
            add_config_tests( ts );
            add_data_tests( ts );

        }

    }

}
