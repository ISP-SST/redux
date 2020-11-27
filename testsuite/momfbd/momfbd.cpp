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
            PointF diff;
            
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
                
                obj->patchSize = 2;         // We're just testing the center position here, so keep it trivial.
                obj->maxLocalShift = 0;
                
                chan->setImageSize(100);       // use image size = (100,100)
                
                chan->alignMap.clear();
                chan->setClip({});
                
                // test that the corners map as expected
                patch->position = 1;            // 0-based coordinates, not 1-based as cfg-file.
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutPosition, PointI(1,1) );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(0,0,1,1) );
                BOOST_CHECK_EQUAL( cd->patchStart, PointI(0,0) );
                chan->setMap({ -1, 0, 99, 0, -1, 99, 0, 0, 1 });
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutPosition, PointI(99,99) );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(98,98,99,99) );
                BOOST_CHECK_EQUAL( cd->patchStart, PointI(0,0) );
                chan->setMap({});
                chan->setClip({ 99, 0, 99, 0 });
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutPosition, PointI(99,99) );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(98,98,99,99) );
                BOOST_CHECK_EQUAL( cd->patchStart, PointI(0,0) );

                // same, but add maxLocalShift
                obj->maxLocalShift = 1;
                patch->position = 2;
                chan->setClip({});
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutPosition, PointI(2,2) );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(0,0,3,3) );
                BOOST_CHECK_EQUAL( cd->patchStart, PointI(1,1) );
                chan->setMap({ -1, 0, 99, 0, -1, 99, 0, 0, 1 });
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutPosition, PointI(98,98) );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(96,96,99,99) );
                BOOST_CHECK_EQUAL( cd->patchStart, PointI(1,1) );
                chan->setMap({});
                chan->setClip({ 99, 0, 99, 0 });
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutPosition, PointI(98,98) );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(96,96,99,99) );
                BOOST_CHECK_EQUAL( cd->patchStart, PointI(1,1) );
                
                obj->maxLocalShift = 0;

                // centered
                chan->setClip({});              // No clip
                patch->position = 51;
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - patch->position;
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );    // should be exact
                BOOST_CHECK_EQUAL( cd->cutoutPosition, PointI(51,51) );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(50,50,51,51) );
                BOOST_CHECK_EQUAL( cd->patchStart, PointI(0,0) );

                // without flips
                chan->setClip({ 0, 99, 0, 99 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(51,51);     // compare to expected
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(50,50,51,51) );
                BOOST_CHECK_EQUAL( cd->patchStart, PointI(0,0) );
                chan->setClip({ 1, 99, 2, 99 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(53,52);
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(52,51,53,52) );
                BOOST_CHECK_EQUAL( cd->patchStart, PointI(0,0) );
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );

                // flip + shift
                chan->setClip({ 1, 98, 97, 2 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(47,52);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                BOOST_CHECK_SMALL( cd->residualOffset.maxAbs(), EPS );
                chan->setClip({ 98, 1, 2, 97 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(53,48);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                BOOST_CHECK_SMALL( cd->residualOffset.maxAbs(), EPS );
                chan->setClip({ 98, 1, 97, 2 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(47,48);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                BOOST_CHECK_SMALL( cd->residualOffset.maxAbs(), EPS );
                chan->setClip({ 97, 1, 98, 2 });              // test with odd clip-sizes
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(48,47);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                BOOST_CHECK_SMALL( cd->residualOffset.maxAbs(), EPS );

                // simple maps  (integer shifts and mirrorings)
                patch->position = 51;
                chan->setImageSize(100);       // map needs imgSize to be set
                chan->setClip({});
                chan->setMap({ 1, 0, 0, 0, 1, 0, 0, 0, 1 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(51,51);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                BOOST_CHECK_SMALL( cd->residualOffset.maxAbs(), EPS );
                chan->setMap({ -1, 0, 89, 0, 1, 0, 0, 0, 1 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(51,39);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                BOOST_CHECK_SMALL( cd->residualOffset.maxAbs(), EPS );
                chan->setMap({ 1, 0, 0, 0, -1, 108, 0, 0, 1 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointI(58,51);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                BOOST_CHECK_SMALL( cd->residualOffset.maxAbs(), EPS );

                // decimal maps
                EPS = 1E-4;
                chan->setMap({ 1, 0, 5.6, 0, 1, 3.8, 0, 0, 1 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointF(54.8,56.6);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                diff = cd->residualOffset - PointF(-0.2,-0.4);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                chan->setMap({ -1, 0, 93.3, 0, 1, 1.2, 0, 0, 1 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointF(52.2,43.3);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                diff = cd->residualOffset - PointF(0.2,0.3);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                chan->setMap({ 1, 0, 0, 0, -1, 97.33, 0, 0, 1 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointF(47.33,51);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                diff = cd->residualOffset - PointF(0.33,0.0);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );

                // offset files
                chan->setMap({});
                chan->setClip({ 0, 99, 0, 99 });
                patch->position = 51;
                Image<int16_t> xoff(100,100), yoff(100,100);
                
                // initial test, check that arrays are copied as references.
                BOOST_CHECK_EQUAL( chan->getOffsetAt( PointI(50,50) ), PointF(0,0) );
                chan->setOffsets( xoff, yoff );
                xoff = -33;
                yoff = 145;
                diff = chan->getOffsetAt( PointI(50,50) ) - PointF(1.45,-0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointF(52.45,50.67);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                diff = cd->residualOffset - PointF(0.45,-0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                chan->setClip({ 99, 0, 0, 99 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointF(52.45,48.67);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                diff = cd->residualOffset - PointF(0.45,-0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                chan->setClip({ 99, 0, 99, 0 });
                chan->adjustCutout( *cd, patch );
                diff = cd->exactPatchPosition - PointF(50.45,48.67);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                diff = cd->residualOffset - PointF(0.45,-0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );

                // cutoutRegion
                // only mirroring & shifts
                xoff = 0; yoff = 0;
                obj->patchSize = 2;
                obj->maxLocalShift = 0;
                chan->setClip({});
                patch->position = 1;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(0,0,1,1) );
                patch->position = 99;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(98,98,99,99) );

                chan->setClip({1,98,2,97});
                patch->position = 1;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(2,1,3,2) );
                patch->position = 95;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(96,95,97,96) );
                chan->setClip({98,1,97,2});
                patch->position = 1;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(96,97,97,98) );
                patch->position = 95;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(2,3,3,4) );

                // add maxLocalShift
                obj->maxLocalShift = 10;
                chan->setClip({});
                patch->position = 11;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(0,0,21,21) );
                BOOST_CHECK_EQUAL( cd->patchStart, PointI(10,10) );
                patch->position = 89;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(78,78,99,99) );
                BOOST_CHECK_EQUAL( cd->patchStart, PointI(10,10) );
                chan->setClip({1,98,2,97});
                patch->position = 11;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(2,1,23,22) );
                BOOST_CHECK_EQUAL( cd->patchStart, PointI(10,10) );
                patch->position = 85;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(76,75,97,96) );
                BOOST_CHECK_EQUAL( cd->patchStart, PointI(10,10) );
                chan->setClip({98,1,97,2});
                patch->position = 11;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(76,77,97,98) );
                BOOST_CHECK_EQUAL( cd->patchStart, PointI(10,10) );
                patch->position = 84;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(3,4,24,25) );
                BOOST_CHECK_EQUAL( cd->patchStart, PointI(10,10) );

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

                patch->position = 77;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(60,58,99,97) );
                diff = cd->residualOffset - PointF(-0.45, 0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                chan->setClip({1,98,2,97});
                patch->position = 20;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(5,2,44,41) );
                diff = cd->residualOffset - PointF(-0.45, 0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                patch->position = 73;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(58,55,97,94) );
                diff = cd->residualOffset - PointF(-0.45, 0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                chan->setClip({97,1,95,2});
                patch->position = 23;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(56,56,95,95) );
                diff = cd->residualOffset - PointF(-0.45, 0.33);
                BOOST_CHECK_SMALL( diff.maxAbs(), EPS );
                patch->position = 73;
                chan->adjustCutout( *cd, patch );
                BOOST_CHECK_EQUAL( cd->cutoutRegion, RegionI(6,6,45,45) );
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
