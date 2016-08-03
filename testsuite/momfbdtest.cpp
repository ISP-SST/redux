#include <boost/test/unit_test.hpp>

#include "redux/momfbd/modes.hpp"
#include "redux/momfbd/momfbdjob.hpp"

#include "redux/constants.hpp"
#include "redux/logging/logger.hpp"
#include "redux/file/fileana.hpp"
#include "redux/math/linalg.hpp"
#include "redux/util/arrayutil.hpp"
#include "redux/image/grid.hpp"
#include "redux/image/utils.hpp"
#include "redux/image/zernike.hpp"

//#include "io.h"
//#include "modes.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/program_options.hpp>
namespace bpo = boost::program_options;
namespace bpt = boost::property_tree;

using namespace redux::momfbd;
using namespace redux::math;
using namespace redux::util;
using namespace redux::image;

using namespace std;
using namespace boost::unit_test_framework;

namespace {

void cfgTest( void ) {

    bpt::ptree tree;
    
    ChannelCfg ccfg;
    {
        // Set some non-default values for ChannelCfg
        ccfg.noiseFudge = ccfg.weight = ccfg.rotationAngle = 138;
        ccfg.diversity = {10.5,32.2};
        ccfg.diversityModes = {11,33};
        ccfg.diversityTypes = {12,34};
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
        ccfg.imageNumbers = {23,45};
        ccfg.wfIndex = {2,4,6};
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
        tmp.parseProperties( tree );
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

    ObjectCfg ocfg;
    {
        // Set some non-default values for ObjectCfg
        ocfg.telescopeF = ocfg.arcSecsPerPixel = ocfg.pixelSize = 137;
        ocfg.alphaToPixels = ocfg.pixelsToAlpha = 138.33;
        ocfg.alphaToDefocus = ocfg.defocusToAlpha = 138.66;
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
        tmp.parseProperties( tree );
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
        gcfg.modeBasis = 2;
        gcfg.klMinMode = gcfg.klMaxMode = gcfg.klCutoff = gcfg.nInitialModes = gcfg.nModeIncrement = 118;
        gcfg.telescopeD = gcfg.minIterations = gcfg.maxIterations = gcfg.targetIterations = 119;
        gcfg.fillpixMethod = gcfg.getstepMethod = 3;
        gcfg.gradientMethod = 2;
        gcfg.badPixelThreshold = gcfg.FTOL = gcfg.EPS = gcfg.reg_gamma = gcfg.sequenceNumber = 117;
        gcfg.outputFileType = gcfg.outputDataType = 1;
        gcfg.modeNumbers = {2,56};
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
        tmp.parseProperties( tree );
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

    MomfbdJob mjob;
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

    /*
            stringstream cfg;
            cfg << "object { }";

            bpt::read_info( cfg, tree );
            mjob.parsePropertyTree( vm, tree );

            bpt::write_info( cout<<endl, tree );*/

}



void packTest( void ) {
    //return;
    MomfbdJob mjob;
    {   // default-constructed job
        auto buf = sharedArray<char>( mjob.size() );
        char* ptr = buf.get();
        uint64_t count = mjob.pack( ptr );
        BOOST_CHECK_EQUAL( count, mjob.size() );
        string tmpS = string( ptr );
        redux::Job::JobPtr job = redux::Job::newJob( tmpS );
        BOOST_CHECK( job );
        count = job->unpack( ptr, false );
        BOOST_CHECK_EQUAL( count, mjob.size() );
        BOOST_CHECK_EQUAL( count, job->size() );
        //             job->info.id = ++jobCounter;
        //             job->info.step.store( Job::JSTEP_PREPROCESS );
        //             job->info.name = "job_" + to_string( job->info.id );
        //             job->info.submitTime = boost::posix_time::second_clock::local_time();
        //             ids.push_back( jobCounter );
        //             ids[0]++;
        //             jobs.push_back( job );

    }

    {   // no image data
        PatchData pd(mjob),pd2(mjob);
        pd.id = 123;
        pd.index = {1,2};
        auto buf = sharedArray<char>( pd.size() );
        char* ptr = buf.get();
        uint64_t count = pd.pack( ptr );
        BOOST_CHECK_EQUAL( count, pd.size() );
        count = pd2.unpack( ptr, false );
        BOOST_CHECK_EQUAL( count, pd.size() );
        BOOST_CHECK_EQUAL( count, pd2.size() );
        BOOST_CHECK( pd == pd2 );

        // with images
        /*pd.images.resize(10,100,120);
        float cnt(0.3);
        for(auto& it: pd.images) it = (cnt += 1);
        buf = sharedArray<char>( pd.size() );
        ptr = buf.get();
        count = pd.pack( ptr );
        BOOST_CHECK_EQUAL( count, pd.size() );
        count = pd2.unpack( ptr, false );
        BOOST_CHECK_EQUAL( count, pd.size() );
        BOOST_CHECK_EQUAL( count, pd2.size() );
        BOOST_CHECK( pd == pd2 );*/

    }
    /*




            Array<T> tmp;
            count = tmp.unpack( ptr,false );
            BOOST_CHECK_EQUAL( count, array.size() );

            BOOST_CHECK( tmp == array );*/

}


void modeTest(void) {

    //return;
    //io_class bla;

    double lambda = 6.30200e-07; //2; //7.77200e-07; //1000;
    double rc = 19.4376; //80; //87.1075; //29.7;
    double angle = 0; //45*redux::PI/180.0;
    int nPixels = 44; //164;
    double cutoff = 0.001;

    int first_mode = 2;
    int last_mode = 50;
    int nModes = last_mode - first_mode + 1;

    const redux::image::Grid& grid = redux::image::Grid::get(nPixels);
    Pupil pupil(nPixels,rc);
    float** aPtr = grid.angle.get();
    Array<float> awrapper(*aPtr, nPixels, nPixels);
    redux::file::Ana::write("modetest_angle.f0", awrapper);
    awrapper.wrap(*grid.distance.get(),nPixels, nPixels);
    redux::file::Ana::write("grid_distance.f0", awrapper);
    redux::file::Ana::write("modetest_pupil.f0", pupil);
    
    cout << "PupilArea = " << pupil.area << endl;
    //PupilMode::KL_cfg* m_new_cfg = legacy::klConfig(first_mode, last_mode);
    //klmc* m_cfg = kl_cfg(first_mode, last_mode);

    //const std::map<int, Modes::KL_cfg>& kle = cache.karhunenLoeveExpansion( first_mode, last_mode );
    auto kle = Zernike::karhunenLoeveExpansion(first_mode, last_mode);
    auto kle2 = Zernike::karhunenLoeveExpansion(first_mode, last_mode + 1);


    for(int m = 0; m < nModes; ++m) {
        //cout << "old:  #" << (first_mode + m) << printArray(m_cfg[first_mode + m].c + 1, m_cfg[first_mode + m].nm, "  coeffs") << endl;
        //cout << "old:  #" << (first_mode + m) << printArray(m_cfg[first_mode + m].m + 1, m_cfg[first_mode + m].nm, "  modes") << endl;
        vector<double> coeffs;
        //for(auto &it : m_new_cfg[m].zernikeWeights ) coeffs.push_back(it.second);
        //cout << "ned:  #" << (first_mode+m) << printArray(coeffs,"  coeffs") << endl;
        //coeffs.clear();
        for(auto & weight : kle[first_mode + m]->zernikeWeights)
            coeffs.push_back(weight.second);
        //cout << "new:  #" << (first_mode + m) << printArray(coeffs, "  coeffs") << endl;
        //cout << "      v = " << m_cfg[first_mode + m].v << "    cov = " << kle[first_mode + m]->covariance << endl;
    }


    // use temporary storage z to avoid recomputing too many Zernikes
    //mde **z = new mde* [nModes];
    //memset(z, 0, nModes * sizeof(mde*));
    Array<double> newZModes(nModes,nPixels,nPixels);
    Array<double> newKLModes(nModes,nPixels,nPixels);
    Array<double> oldZModes(nModes,nPixels,nPixels);
    Array<double> oldKLModes(nModes,nPixels,nPixels);
    Array<double> newz_slice(newZModes, 0, 0, 0, nPixels-1, 0, nPixels-1);          // subarray @ first mode
    Array<double> newkl_slice(newKLModes, 0, 0, 0, nPixels-1, 0, nPixels-1);        // subarray @ first mode
    Array<double> oldz_slice(oldZModes, 0, 0, 0, nPixels-1, 0, nPixels-1);          // subarray @ first mode
    Array<double> oldkl_slice(oldKLModes, 0, 0, 0, nPixels-1, 0, nPixels-1);        // subarray @ first mode

    for(int j = first_mode; j <= last_mode; ++j) {

        PupilMode zmode(j, nPixels, rc, angle);
        PupilMode klmode(first_mode, last_mode, j, nPixels, rc, angle, cutoff);

        //pmd oldzmode(lambda, rc, nPixels, j, angle, bla);                //  pmd(lambda,r_c,nph,mn,angle,io));
        //pmd oldklmode(lambda, rc, nPixels, j, m_cfg, cutoff, z - first_mode, angle, bla);
        //Array<double> zwrapper(*(oldzmode.mode + 1) + 1, nPixels, nPixels);
        //Array<double> klwrapper(*(oldklmode.mode + 1) + 1, nPixels, nPixels);
        
        //newz_slice = zmode.copy();
        zmode.copy(newz_slice);
        newz_slice *= pupil;
        newz_slice /= lambda;
        newz_slice.shift(0,1);
        
        //newkl_slice = klmode.copy();
        klmode.copy(newkl_slice);
        newkl_slice *= pupil;
        newkl_slice /= lambda;
        newkl_slice.shift(0,1);
        
//         oldz_slice = zwrapper;
//         oldz_slice *= pupil.first;
//         oldz_slice.shift(0,1);
//         
//         oldkl_slice = klwrapper;
//         oldkl_slice *= pupil.first;
//         oldkl_slice.shift(0,1);
        
        //cout << printArray(mode,"mode") << endl;
        //redux::file::Ana::write("newzmode_" + to_string(j) + ".f0", zmode);
        //redux::file::Ana::write("newklmode_" + to_string(j) + ".f0", klmode);
        // redux::file::Ana::write("modetest_" + to_string(j) + "b.f0", wrapper);

    }
/*
npixels=44
nmodes=20
imgperrow=5
lambda=6.30200e-07
datadir='/home/tomas/build/redux/'
rdir=datadir+'src/bin/'
mdir=datadir+'src/momfbd/'
pupil=f0(datadir+'modetest_pupil.f0')
otfsupport=f0(rdir+'otfsupport.f0') 
ang=f0(datadir+'modetest_angle.f0')
;rq=f0(rdir+'Q.f0')
;mq=f0(mdir+'Q.f0')
newz=f0(datadir+'newzmodes.f0')
newkl=f0(datadir+'newklmodes.f0')
;oldz=f0(datadir+'oldzmodes.f0')
;oldkl=f0(datadir+'oldklmodes.f0')
oldz=transpose(f0(datadir+'oldzmodes.f0'),[1,0,2])
oldkl=transpose(f0(datadir+'oldklmodes.f0'),[1,0,2])
diff=newz-oldz
signs=[-1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1]
dist=f0(datadir+'grid_distance.f0')
ap=f0(datadir+'aperture.f0')
ap2=f0(datadir+'aperture2.f0')
diffa = ap-ap2
for q=0,nmodes-1 do tvscl,newz(*,*,q),npixels*(q mod imgperrow),npixels*(q/imgperrow)
for q=0,nmodes-1 do tvscl,oldz(*,*,q),npixels*(q mod imgperrow),npixels*(q/imgperrow)
for q=0,nmodes-1 do tvscl,diff(*,*,q),npixels*(q mod imgperrow),npixels*(q/imgperrow)
for q=0,nmodes-1 do tvscl,newkl(*,*,q),npixels*(q mod imgperrow),npixels*(q/imgperrow)
for q=0,nmodes-1 do tvscl,oldkl(*,*,q),npixels*(q mod imgperrow),npixels*(q/imgperrow)
;for q=0,nmodes-1 do print,min(oldz(*,*,q)),max(oldz(*,*,q)),min(newz(*,*,q)),max(newz(*,*,q))
;for q=0,nmodes-1 do print,min(oldz(*,*,q))/min(newz(*,*,q)),max(oldz(*,*,q))/max(newz(*,*,q))
for q=0,nmodes-1 do print,min(newz(*,*,q)-oldz(*,*,q)),max(newz(*,*,q)-oldz(*,*,q))
;for q=0,nmodes-1 do print,min(oldkl(*,*,q)),max(oldkl(*,*,q)),min(newkl(*,*,q)),max(newkl(*,*,q))
for q=0,nmodes-1 do print,min(newkl(*,*,q))-min(oldkl(*,*,q)),max(newkl(*,*,q))-max(oldkl(*,*,q))
for q=0,nmodes-1 do print,min(diff(*,*,q)),max(diff(*,*,q))
    
    
otfsupport=f0(rdir+'otfsupport.f0') 
mq=transpose(f0(mdir+'Q.f0'))   
rq=f0(rdir+'Q.f0')              
mp=transpose(f0(mdir+'P.f0'))
rp=f0(rdir+'P.f0')           
tvscl,alog(real_part(mp)>0.1)
print,min(rp),max(rp),mean(rp)
print,min(mp),max(mp),mean(mp)
print,min(rq),max(rq),mean(rq)
print,min(mq),max(mq),mean(mq)

rrpup=f0(rdir+'blu_pup.f0')              
mmpup=f0(mdir+'blapup.f0')              
diff=rrpup-mmpup
print,min(diff),max(diff),mean(diff)
tvscl,diff

mphisum=transpose(f0(mdir+'phisum.f0'))   
mpf=transpose(f0(mdir+'pf.f0'))
rpf=f0(rdir+'blu_pf1.f0')
diff=mpf-rpf
print,min(diff),max(diff),mean(diff)
motf=transpose(f0(mdir+'otf.f0'))   
rotf=f0(rdir+'blu_otf1.f0')   
diff=motf-rotf           
print,min(diff),max(diff),mean(diff)


rpupil=f0(rdir+'pupil.f0')
mpupil=f0(mdir+'pupil1.f0')
diff=rpupil-mpupil
diff=rrpup-mpupil
print,min(diff),max(diff),mean(diff)
tvscl,diff



rvp=f0(rdir+'vogel_p.f0')
mvp=f0(mdir+'vogel_p.f0')
diff=mvp-conj(rvp)
print,min(diff),max(diff),mean(diff)
rvq=f0(rdir+'vogel_q.f0')
mvq=f0(mdir+'vogel_q.f0')
diff=mvq-rvq
print,min(diff),max(diff),mean(diff)

rvpj=f0(rdir+'vogel_pj.f0')
mvpj=f0(mdir+'vogel_pj.f0')
diff=mvpj-rvpj
print,min(diff),max(diff),mean(diff)

rvhj=f0(rdir+'vogel_hj2.f0')
mvhj=f0(mdir+'vogel_hj2.f0')
diff=mvhj-rvhj
print,min(diff),max(diff),mean(diff)

rvotf=f0(rdir+'vogel_otf.f0')
mvotf=f0(mdir+'vogel_otf.f0')
diff=mvotf-rvotf
print,min(diff),max(diff),mean(diff)

rvglft=f0(rdir+'vogel_glft.f0')
mvglft=f0(mdir+'vogel_glft.f0')
diff=mvglft-rvglft
print,min(diff),max(diff),mean(diff)


oldpsf=f0('src/bin/oldpsf.f0')                 
oldotf=f0('src/bin/oldotf.f0')
newotf=f0('src/bin/newotf.f0')
newpsf=f0('src/bin/newpsf.f0')
tvscl,oldpsf
tvscl,newpsf

mode3=f0('src/bin/mode_3_0.f0')
mode4=f0('src/bin/mode_4_1.f0')
phimode3=f0('src/bin/phi-mode_3_0.f0')
phimode4=f0('src/bin/phi-mode_4_1.f0')
phi=f0('src/bin/phi_1.f0')
tvscl,mode3*pupil
tvscl,mode4*pupil,npixels,0
tvscl,phimode3*pupil,2*npixels,0
tvscl,phimode4*pupil,3*npixels,0
tvscl,phi*pupil,4*npixels,0
print,mean(mode3),mean(mode4),mean(phimode3),mean(phimode4),mean(phi)

pupil=f0(datadir+'modetest_pupil.f0')
zpol_1_1.f0


zpol_1_1=f0(datadir+'zpol_1_1.f0')
zpol_1_5=f0(datadir+'zpol_1_5.f0')
zpol_1_3=f0(datadir+'zpol_1_3.f0')
zpol_0_4=f0(datadir+'zpol_0_4.f0')
zpol_0_2=f0(datadir+'zpol_0_2.f0')
zpol_2_2=f0(datadir+'zpol_2_2.f0')
zpol_2_4=f0(datadir+'zpol_2_4.f0')
zpol_3_3=f0(datadir+'zpol_3_3.f0')
zpol_3_5=f0(datadir+'zpol_3_5.f0')
zpol_4_4=f0(datadir+'zpol_4_4.f0')
zpol_5_5=f0(datadir+'zpol_5_5.f0')
print,zpol_1_1
print,zpol_1_5
print,zpol_1_3
print,zpol_0_4
print,zpol_0_2
print,zpol_2_2
print,zpol_2_4
print,zpol_3_3
print,zpol_3_5
print,zpol_4_4
print,zpol_5_5







wborig=f0('src/bin/camXX.00010.7772.0014668.orig')
wborig1=f0('src/bin/camXX.00010.7772.0014668.1.orig')
wborig2=f0('src/bin/camXX.00010.7772.0014668.2.orig')
wborig3=f0('src/bin/camXX.00010.7772.0014668.3.orig')
wborig4=f0('src/bin/camXX.00010.7772.0014668.4.orig')
wborig5=f0('src/bin/camXX.00010.7772.0014668.5.orig')
tvscl,wborig
tvscl,wborig1
tvscl,wborig2
tvscl,wborig3
tvscl,wborig4
tvscl,wborig5

img=f0('src/bin/windowed_1_0.f0')
imgft=f0('src/bin/windowedft_1_0.f0')
tvscl,img

ftsum_0=f0('src/bin/ftsum_0.f0')
ftsum_1=f0('src/bin/ftsum_1.f0')
ftsum_2=f0('src/bin/ftsum_2.f0')
ftsum_3=f0('src/bin/ftsum_3.f0')
ftsum_4=f0('src/bin/ftsum_4.f0')
ftsum_5=f0('src/bin/ftsum_5.f0')
ftsum_6=f0('src/bin/ftsum_6.f0')
ftsum_7=f0('src/bin/ftsum_7.f0')
ftsum_8=f0('src/bin/ftsum_8.f0')
tvscl,alog10(ftsum_0)

Output from testsuite:
Area = 13198.1  Area_mvn  = 13121.6   r^2*PI = 13197.8   ->  Area error:   T=+0.0023%  MvN=-0.58%

*/
    redux::file::Ana::write("newzmodes.f0", newZModes);
    redux::file::Ana::write("newklmodes.f0", newKLModes);
    //redux::file::Ana::write("oldzmodes.f0", oldZModes);
    //redux::file::Ana::write("oldklmodes.f0", oldKLModes);
    //delete[] z;


    // kl_uncfg(m_cfg, first_mode, last_mode);

    nPixels = 150;
    double radius = 0.4321*nPixels;
    Array<double> aperture;
    double area = makePupil_old(aperture,nPixels,radius);
    //redux::file::Ana::write("aperture.f0", aperture);

    double** ap = makePointers(aperture.ptr(), nPixels, nPixels);
    double area_mvn = makePupil_mvn(ap,nPixels,radius);
    //redux::file::Ana::write("aperture2.f0", aperture);
    delPointers(ap);

    Pupil pup(nPixels,radius);
    //redux::file::Ana::write("aperture3.f0", pup_pair.first);

    cout << "Area = " << area << "  Area_mvn  = " << area_mvn << "  ratio  = " << (area/area_mvn) << "   r^2*PI = " << (radius*radius*redux::PI) << endl;
//     cout << "ap=f0('/home/tomas/build/redux/aperture.f0')\n"
//          << "ap2=f0('/home/tomas/build/redux/aperture2.f0')\n"
//          << "ap[" <<(nPixels/2-2) << ":" << (nPixels/2+2) << ","<<(nPixels/2-2) << ":" << (nPixels/2+2) << "]\n"
//          << "ap2[" <<(nPixels/2-1) << ":" << (nPixels/2+3) << ","<<(nPixels/2-1) << ":" << (nPixels/2+3) << "]\n"
//          << "diff = ap-ap2[1:" << nPixels << ",1:" << nPixels << "]\n"
//          << "print,min(diff),max(diff),mean(diff)" << endl;

}









//     #define ROWS 512
//     #define COLS 512
//     #include <boost/timer/timer.hpp>
//
//     void qrTest(void) {
//
//         Array<double> A(ROWS, COLS);
//         Array<double> Q(ROWS, ROWS);
//         Array<double> R(ROWS, COLS);
//         double** AA = newArray<double>(ROWS, COLS);
//         //double** Q = newArray<double>(ROWS, ROWS);
//         //double** R = newArray<double>(ROWS, COLS);
//         // set Lotkin matrix
//         // a_rc = 1 (r = 0) or 1/(r+c-1) (r > 0) (i.e. Hilbert matrix with first row set to 1)
//         for(int c = 0; c < COLS; ++c) A(0,c) = 1;
//         for(int r = 1; r < ROWS; ++r) {
//             for(int c = 0; c < COLS; ++c) {
//                 A(r,c) = 1.0 / (double)(r + c + 1);
//             }
//         }
//
//     #if 0
//         // set Frank matrix
//         // a_ij = DIM - min(i,j) + 1
//         for(int r = 1; r < ROWS; ++r) {
//             for(int c = 0; c < COLS; ++c) {
//                 A(r, c) = ROWS - std::min(r, c) + 1;
//             }
//         }
//     #endif
//         memcpy(*AA,A.get(),ROWS*COLS*sizeof(double));
//         printf("Matrix(%d,%d)\n", ROWS, COLS);
//         if(ROWS < 20 && COLS < 20) {
//             for(int r = 0; r < ROWS; ++r) {
//                 printf("%3d: ", r);
//                 for(int c = 0; c < COLS; ++c) {
//                     cout << alignRight(to_string(A(r,c)), 8) << "  ";
//                 }
//                 printf("\n");
//             }
//             printf("\n");
//         }
//         double** q;
//         {
//             boost::timer::auto_cpu_timer timer;
//             q = legacy::qr(AA, ROWS, COLS, 0);
//             printf("Q(%d,%d)\n", ROWS, ROWS);
//             if(ROWS < 20 && COLS < 20) {
//                 for(int r = 0; r < ROWS; ++r) {
//                     printf("%3d: ", r);
//                     for(int c = 0; c < ROWS; ++c) {
//                         cout << alignRight(to_string(q[r][c]), 8) << "  ";
//                     }
//                     printf("\n");
//                 }
//                 printf("\n");
//             }
//
//         }
//
//         //memcpy(*AA,*A,ROWS*COLS*sizeof(double));
//         {
//             boost::timer::auto_cpu_timer timer;
//             qr_decomp(A.get(), ROWS, COLS, Q.get(), R.get());
//             printf("Q2(%d,%d)\n", ROWS, ROWS);
//             if(ROWS < 20 && COLS < 20) {
//                 for(int r = 0; r < ROWS; ++r) {
//                     printf("%3d: ", r);
//                     for(int c = 0; c < ROWS; ++c) {
//                         cout << alignRight(to_string(Q(r,c)), 8) << "  ";
//                     }
//                     printf("\n");
//                 }
//                 printf("\n");
//             }
//         }
//
//         //cout << printArray(AA,0,ROWS-1,0,COLS-1,"R1") << endl;
//         //cout << printArray(q,0,ROWS-1,0,ROWS-1,"Q1") << endl;
//         //cout << printArray(A,0,ROWS-1,0,COLS-1,"A") << endl;
//         //cout << printArray(Q,0,ROWS-1,0,ROWS-1,"Q") << endl;
//         //cout << printArray(R,0,ROWS-1,0,COLS-1,"R") << endl;
//
//         redux::file::Ana::write("a.f0", A);
//         redux::file::Ana::write("q.f0", Q);
//         redux::file::Ana::write("r.f0", R);
//
//         Array<double> qwrapper(*q, ROWS, ROWS);
//         Array<double> awrapper(*AA, ROWS, COLS);
//
//         redux::file::Ana::write("q2.f0", qwrapper);
//         redux::file::Ana::write("a2.f0", awrapper);
//
//         //delArray(A);
//         //delArray(R);
//         //delArray(Q);
//         delArray(q);
//         delArray(AA);
//
//     }
//
//     #include "redux/image/fouriertransform.hpp"
//
//     #define XOFF 0
//     #define YOFF 0

//     void ftTest(void) {
//         using namespace redux::image;
//     srand(time(NULL));
//
//         Array<double> img(ROWS, COLS);
//         Array<double> img2(ROWS, COLS);
//         Array<redux::complex_t> cimg(ROWS, COLS);
//         Array<double> psf(ROWS, COLS);
//
//         //Array<double> delta(ROWS, COLS);
//         //delta.zero();
//         //delta(ROWS/2, COLS/2) = 107.4;
//         img.zero();
//         img2.zero();
//         cimg.zero();
//         psf.zero();
//         redux::complex_t iii(0,1);
//         for( int i=0; i<COLS; ++i) {
//             for( int j=0; j<ROWS; ++j) {
//          //       cimg(j,i) = img(j,i) = rand()%10000;
//                 //img(j,i) = 1+sin(4*redux::PI*i/(COLS-1));
//     //            cimg(j,i) = img(j,i) = cos(redux::PI*(i+2*j)/24.0)*exp(-(i-COLS/2)*(i-COLS/2)/20000.0);
//                // cimg(j,i) *= iii;
//     //             img(j,i) = 1+sin(20*redux::PI*i/(COLS-1)); // +
//     //                        4+4*sin(20*redux::PI*j/(ROWS-1)) +
//     //                        16+16*cos(14*redux::PI*(i/(COLS-1)+j/(ROWS-1)));
//     //                        64+64*sin(14*redux::PI*(i/(COLS-1)-j/(ROWS-1)));
//
//                 psf(j,i) = exp(-(i-COLS/2)*(i-COLS/2)/5.0-(j-ROWS/2)*(j-ROWS/2)/100.0);
//                 if(sqrt((i-COLS/2-XOFF)*(i-COLS/2-XOFF)+(j-ROWS/2-YOFF)*(j-ROWS/2-YOFF)) < 30) cimg(j,i) = img(j,i) = 1;
//                // if(sqrt((i-COLS/2+XOFF)*(i-COLS/2+XOFF)+(j-ROWS/2+YOFF)*(j-ROWS/2+YOFF)) < 30) img2(j,i) = 1;
//                 if(abs(i-COLS/2+XOFF)<30 && abs(j-ROWS/2+YOFF) < 30) img2(j,i) = 1;
//     //            img(j,i) = rand()%10000;
//             }
//         }
//        // img /= 9999.0;
//
//         redux::image::FourierTransform ft(img);
//         redux::image::FourierTransform cft(cimg);
//         redux::image::FourierTransform cftr(cimg,FT_REORDER);
//         redux::image::FourierTransform psfft(psf,FT_REORDER);
//
//         redux::file::Ana::write("img.f0", img);
//         redux::file::Ana::write("img2.f0", img2);
//         redux::file::Ana::write("cimg.f0", cimg);
//         redux::file::Ana::write("psf.f0", psf);
//
//         redux::file::Ana::write("ft.f0", ft);
//         redux::file::Ana::write("ftpower.f0", ft.power());
//         redux::file::Ana::write("cft.f0", cft);
//         redux::file::Ana::write("cftpower.f0", cft.power());
//         redux::file::Ana::write("cftr.f0", cftr);
//         redux::file::Ana::write("cftrpower.f0", cftr.power());
//         //img = psfft.convolve(img);
//         //ft.inv(img);
//         //redux::file::Ana::write("ftinv.f0", img);
//         redux::file::Ana::write("ftcorr.f0", ft.correlate(img2));
//         //redux::file::Ana::write("psff.f0", ft.convolve(cimg,FT_REORDER));
//         redux::file::Ana::write("fpsf.f0", psfft.convolve(img+img2));
//         //redux::file::Ana::write("fpsf.f0", ft.convolve(psf));
//         redux::file::Ana::write("cpsf.f0", cft.convolve(psf));
//
//         ft.inv(img);
//         redux::file::Ana::write("ftinv.f0", img);
//     //    ft.autocorrelate();
//     //     ft.inv(img);
//         redux::file::Ana::write("ftinv_ac.f0", img);
//
//         //for(auto &it: ft) it = it.real();
//         //ft.inv(img);
//         //redux::file::Ana::write("ftinvre.f0", img);
//
//         //redux::file::Ana::write("delta.f0", delta);
//         //FourierTransform::normalize(delta);
//         //redux::file::Ana::write("deltan.f0", delta);
//
//
//         /*
//          http://www.fmwconcepts.com/misc_tests/FFT_tests/
//          http://qsimaging.com/ccd_noise_interpret_ffts.html
//          http://www.fftw.org/fftw3_doc/Multi_002dDimensional-DFTs-of-Real-Data.html#Multi_002dDimensional-DFTs-of-Real-Data
//          http://hebb.mit.edu/courses/9.29/2002/readings/c13-2.pdf
//
//     ;qq=dindgen(5000,20000)
//     tvscl,qq
//     img=f0('img.f0')
//     img2=f0('img2.f0')
//     psf=f0('psf.f0')
//     psff=f0('psff.f0')
//     imgc=f0('imgc.f0')
//     ftinv=f0('ftinv.f0')
//     ft=f0('ft.f0')
//     cft=f0('cft.f0')
//     cftr=f0('cftr.f0')
//     ftcorr=f0('ftcorr.f0')
//     fpsf=f0('fpsf.f0')
//     cpsf=f0('cpsf.f0')
//     ftinv_ac=f0('ftinv_ac.f0')
//     cimg=f0('cimg.f0')
//     imgr=f0('img_reordered.f0')
//     ftpower=f0('ftpower.f0')
//     cftpower=f0('cftpower.f0')
//     cftrpower=f0('cftrpower.f0')
//     ;ftabs=f0('ftabs.f0')
//
//
//     margin=5
//     cols=512
//     hcols=cols/2+1
//     ;tvscl,real_part(cft)
//     ;tvscl,ALOG10(cftpower>0.00001),cols+margin,0
//     ;tvscl,ALOG10(cftrpower>0.00001),2*(cols+margin),0
//     ;tvscl,cpsf,3*(cols+margin),0
//
//
//     tvscl,img+img2
//     tvscl,ftcorr,cols+margin,0
//     tvscl,real_part(fpsf),2*(cols+margin),0
//     tvscl,cpsf,3*(cols+margin),0
//     ;print,min(cpsf),max(cpsf)
//     ;print,min(real_part(psff)),max(real_part(psff))
//
//
//     ;tvscl,ALOG10(ftpower>0.00001),cols+margin,0
//     ;tvscl,ftabs,777,0
//     ;tvscl,ALOG10(cftpower>0.00001),777,0
//     ;tvscl,ftcorr,1330+10,0
//     ;tvscl,ftconv,1850+10,0
//
//
//        u(i,j) = 1 + ...
//     2^0*(1+sin(2*pi*((i-1)/nx)*200))+...
//     2^2.*(1+sin(2*pi*((j-1)/ny)*200))+...
//     2^4.*(1+cos(2*pi*((i-1)/nx+(j-1)/ny)*141))+...
//     2^6.*(1+sin(2*pi*((i-1)/nx-(j-1)/ny)*141));
//
//         */
//
//         Array<double> ftpower = ft.power();
//         redux::file::Ana::write("ftpower.f0", ftpower);
//         Array<double> cftpower = cft.power();
//         redux::file::Ana::write("cftpower.f0", cftpower);
//
//     //     ft.reorder();
//     //     it = ftpower.begin();
//     //     for( auto& it3: ft ) *it++ = log10(norm(it3));
//     //     redux::file::Ana::write("ftpower.f0", ftpower);
//
//         //redux::file::Ana::write("ftabs2.f0", ftr2);
//
//        // it = ftpower.begin();
//         //it2 = ftr2.begin();
//        // for( auto& it3: ft ) *it++ = arg(it3);
//
//
//         //for( auto& it3: ft2 ) *it2++ = arg(it3);
//        // redux::file::Ana::write("ftarg.f0", ftpower);
//         //redux::file::Ana::write("ftarg2.f0", ftr2);
//
//         //FourierTransform::reorder(img);
//        // redux::file::Ana::write("img_reordered.f0", img);
//
//
//         //redux::file::Ana::write("ft.f0", ft);
//
//
//     }
//
//
}


namespace testsuite {

namespace momfbd {

void momfbdTest(void) {

    test_suite* ts = BOOST_TEST_SUITE("MOMFBD");

    ts->add(BOOST_TEST_CASE(&cfgTest));
    ts->add(BOOST_TEST_CASE(&packTest));
    ts->add(BOOST_TEST_CASE(&modeTest));
    //ts->add(BOOST_TEST_CASE(&qrTest));
    //ts->add(BOOST_TEST_CASE(&ftTest));


    framework::master_test_suite().add(ts);

}

}

}
