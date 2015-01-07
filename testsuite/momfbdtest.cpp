#include <boost/test/unit_test.hpp>

#include "redux/constants.hpp"
#include "redux/logger.hpp"
#include "redux/file/fileana.hpp"
#include "redux/momfbd/modes.hpp"
#include "redux/momfbd/modecache.hpp"
#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/patch.hpp"
//#include "redux/momfbd/legacy.hpp"
#include "redux/math/linalg.hpp"
#include "redux/util/arrayutil.hpp"

//#include "modes.h"
#include "redux/image/utils.hpp"

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
        stringstream cfg;

        ChannelCfg ccfg;
        {
            ccfg.arcSecsPerPixel = ccfg.pixelSize = ccfg.rotationAngle = ccfg.weight = ccfg.borderClip = ccfg.maxLocalShift = 137; // just some non-default values
            ccfg.imageDataDirl = "path/to/data/";
            ccfg.imageNumbers = {23,45};
            ccfg.darkNumbersl = {145,88};
            // export settings as config file and parse/compare
            tree.clear();
            ccfg.getProperties(tree);
            cfg.clear();
            bpt::write_info( cfg, tree );   // basically just to test that the ptree is valid
            bpt::read_info( cfg, tree );
            ChannelCfg tmp;
            BOOST_CHECK( !(ccfg == tmp) );
            tmp.parseProperties( tree );
            BOOST_CHECK( ccfg == tmp );
            // pack/unpack and compare
            auto buf = sharedArray<char>( ccfg.size() );
            char* ptr = buf.get();
            uint64_t count = ccfg.pack( ptr );
            BOOST_CHECK_EQUAL( count, ccfg.size() );
            tmp = ChannelCfg();
            count = tmp.unpack( ptr, false );
            BOOST_CHECK_EQUAL( count, ccfg.size() );
            BOOST_CHECK_EQUAL( count, tmp.size() );
            BOOST_CHECK( ccfg == tmp );
            tmp.arcSecsPerPixel = 15;       // modify
            BOOST_CHECK( !(ccfg == tmp) );
            tmp = ccfg;                     // assign
            BOOST_CHECK( ccfg == tmp );
        }
    /*      return (saveMask == rhs.saveMask) &&
           (patchSize == rhs.patchSize) &&
           (pupilSize == rhs.pupilSize) &&
           (wavelength == rhs.wavelength) &&
           (outputFileName == rhs.outputFileName) &&
           (pupilFile == rhs.pupilFile) &&
           ChannelCfg::operator==(rhs);
  */
        ObjectCfg ocfg;
        {
            ocfg.saveMask = ocfg.patchSize = ocfg.pupilSize = ocfg.wavelength = 122; // just some non-default values
            ocfg.outputFileName = "filename.ext";
            ocfg.pupilFile = "pupilfile.ext";
            // export settings as config file and parse/compare
            tree.clear();
            ocfg.getProperties(tree);
            cfg.clear();
            bpt::write_info( cfg, tree );   // basically just to test that the ptree is valid
            bpt::read_info( cfg, tree );
            ObjectCfg tmp;
            BOOST_CHECK( !(ocfg == tmp) );
            tmp.parseProperties( tree );
            BOOST_CHECK( ocfg == tmp );
            // pack/unpack and compare
            auto buf = sharedArray<char>( ocfg.size() );
            char* ptr = buf.get();
            uint64_t count = ocfg.pack( ptr );
            BOOST_CHECK_EQUAL( count, ocfg.size() );
            tmp = ObjectCfg();
            count = tmp.unpack( ptr, false );
            BOOST_CHECK_EQUAL( count, ocfg.size() );
            BOOST_CHECK_EQUAL( count, tmp.size() );
            BOOST_CHECK( ocfg == tmp );
            BOOST_CHECK( ChannelCfg() == tmp );
            tmp.patchSize = 15;              // modify
            BOOST_CHECK( !(ocfg == tmp) );
            tmp = ocfg;                     // assign
            BOOST_CHECK( ocfg == tmp );
            tmp.arcSecsPerPixel = 15;       // modify ChannelCfg
            BOOST_CHECK( !(ocfg == tmp) );
            tmp = ChannelCfg();             // assign default ChannelCfg
            BOOST_CHECK( ocfg == tmp );
        }

        GlobalCfg gcfg;
        {
            gcfg.runFlags = 4095;
            gcfg.modeBasis = 2;
            gcfg.klMinMode = gcfg.klMaxMode = gcfg.klCutoff = gcfg.nInitialModes = gcfg.nModeIncrement = 118;
            gcfg.telescopeD = gcfg.telescopeF = gcfg.minIterations = gcfg.maxIterations = 119;
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
            cfg.clear();
            bpt::write_info( cfg, tree );   // basically just to test that the ptree is valid
            bpt::read_info( cfg, tree );
            GlobalCfg tmp;
            BOOST_CHECK( !(gcfg == tmp) );
            tmp.parseProperties( tree );
            BOOST_CHECK( gcfg == tmp );
            // pack/unpack and compare
            auto buf = sharedArray<char>( gcfg.size() );
            char* ptr = buf.get();
            uint64_t count = gcfg.pack( ptr );
            BOOST_CHECK_EQUAL( count, gcfg.size() );
            tmp = GlobalCfg();
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
            tmp.arcSecsPerPixel = 15;       // modify ChannelCfg
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
            PatchData pd,pd2;
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
            pd.images.resize(10,100,120);
            float cnt(0.3);
            for(auto& it: pd.images) it = (cnt += 1);
            buf = sharedArray<char>( pd.size() );
            ptr = buf.get();
            count = pd.pack( ptr );
            BOOST_CHECK_EQUAL( count, pd.size() );
            count = pd2.unpack( ptr, false );
            BOOST_CHECK_EQUAL( count, pd.size() );
            BOOST_CHECK_EQUAL( count, pd2.size() );
            BOOST_CHECK( pd == pd2 );
            
        }
/*




        Array<T> tmp;
        count = tmp.unpack( ptr,false );
        BOOST_CHECK_EQUAL( count, array.size() );

        BOOST_CHECK( tmp == array );*/

    }

    
       
        
    void modeTest(void) {

return;
       // io_class bla;

        double lambda = 100;
        double rc = 72;
        double angle = 0;
        int nPixels = 150;
        double svd_reg = 0;

        int first_mode = 1;
        int last_mode = 25;
        int nModes = last_mode - first_mode + 1;

        static ModeCache& cache = ModeCache::getCache();
        const redux::image::Grid& grid = cache.grid(nPixels);
        float** aPtr = grid.angle.get();
        Array<float> awrapper(*aPtr, nPixels, nPixels);
        redux::file::Ana::write("modetest_angle.f0", awrapper);

        //PupilMode::KL_cfg* m_new_cfg = legacy::klConfig(first_mode, last_mode);
        //klmc* m_cfg = kl_cfg(first_mode, last_mode);

        //const std::map<int, Modes::KL_cfg>& kle = cache.karhunenLoeveExpansion( first_mode, last_mode );
        auto kle = cache.karhunenLoeveExpansion(first_mode, last_mode);
        auto kle2 = cache.karhunenLoeveExpansion(first_mode, last_mode + 1);


        for(int m = 0; m < nModes; ++m) {
            //cout << "old:  #" << (first_mode + m) << printArray(m_cfg[first_mode + m].c + 1, m_cfg[first_mode + m].nm, "  coeffs") << endl;
            vector<double> coeffs;
            //for(auto &it : m_new_cfg[m].zernikeWeights ) coeffs.push_back(it.second);
            //cout << "ned:  #" << (first_mode+m) << printArray(coeffs,"  coeffs") << endl;
            //coeffs.clear();
            for(auto & it : kle[first_mode + m].zernikeWeights) coeffs.push_back(it.second);
            //cout << "new:  #" << (first_mode + m) << printArray(coeffs, "  coeffs") << endl;
        }



        PupilMode **zz = new PupilMode* [nModes];
        memset(zz, 0, nModes * sizeof(PupilMode*));

    // use temporary storage z to avoid recomputing too many Zernikes
       /// mde **z = new mde* [nModes];
       // memset(z, 0, nModes * sizeof(mde*));

        for(int j = first_mode; j <= last_mode; ++j) {

            //PupilMode mode(j, nPixels, rc, lambda, angle);           //  PupilMode ( int modeIndex, int nPoints, double r_c = 1.0, double lambda = 1.0, double angle = 0.0 ); // Zernike
            //PupilMode mode( m_new_cfg, j, nPixels, zz-first_mode, rc, lambda, angle, svd_reg );
            PupilMode mode(first_mode, last_mode, j, nPixels, rc, lambda, angle, svd_reg);
    //cout << "loop: " << hexString(&mode) << endl;
            //pmd oldmode(lambda, rc, nPixels, j, angle, bla);                //  pmd(lambda,r_c,nph,mn,angle,io));
            //pmd oldmode(lambda, rc, nPixels, j, m_cfg, svd_reg, z - first_mode, angle, bla);
           // Array<double> wrapper(*(oldmode.mode + 1) + 1, nPixels, nPixels);

            //cout << printArray(mode,"mode") << endl;
            redux::file::Ana::write("modetest_" + to_string(j) + ".f0", mode);
           // redux::file::Ana::write("modetest_" + to_string(j) + "b.f0", wrapper);

        }

       // delete[] z;
        delete[] zz;


       // kl_uncfg(m_cfg, first_mode, last_mode);
        
        nPixels = 150;
        double radius = 0.4321*nPixels;
        Array<double> aperture;
        double area = makePupil(aperture,nPixels,radius);
        redux::file::Ana::write("aperture.f0", aperture);

        Array<double> aperture2(nPixels+1,nPixels+1);    // add extra space to make room for Michiels +1 offset to array indexing
        auto a2Ptr = aperture2.get(nPixels+1,nPixels+1);
        double area_mvn = makePupil_mvn(a2Ptr.get(),nPixels,radius);
        redux::file::Ana::write("aperture2.f0", aperture2);
        
        auto pup_pair = cache.pupil(nPixels,radius);
        redux::file::Ana::write("aperture3.f0", pup_pair.first);
        
        cout << "Area = " << area << "  Area_mvn  = " << area_mvn << "   r^2*PI = " << (radius*radius*redux::PI) << endl;
        cout << "ap=f0('aperture.f0')\n"
             << "ap2=f0('aperture2.f0')\n"
             << "ap[" <<(nPixels/2-2) << ":" << (nPixels/2+2) << ","<<(nPixels/2-2) << ":" << (nPixels/2+2) << "]\n"
             << "ap2[" <<(nPixels/2-1) << ":" << (nPixels/2+3) << ","<<(nPixels/2-1) << ":" << (nPixels/2+3) << "]\n"
             << "diff = ap-ap2[1:" << nPixels << ",1:" << nPixels << "]\n"
             << "print,min(diff),max(diff),mean(diff)" << endl;

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
