
#include "redux/file/fileana.hpp"
#include "redux/image/grid.hpp"
#include "redux/image/utils.hpp"
#include "redux/image/zernike.hpp"
#include "redux/momfbd/modes.hpp"
#include "redux/util/arrayutil.hpp"
 
#include <boost/test/unit_test.hpp>
 
using namespace redux::momfbd;
using namespace redux::util;
using namespace redux::image;

using namespace std;


namespace testsuite {

    namespace image {


        void zernikeTest(void) {

            return;
            //io_class bla;

            double lambda = 6.30200e-07; //2; //7.77200e-07; //1000;
            double rc = 19.4376; //80; //87.1075; //29.7;
            double angle = 0; //45*M_PI/180.0;
            int nPixels = 44; //164;
            double cutoff = 0.001;

            int first_mode = 2;
            int last_mode = 50;
            int nModes = last_mode - first_mode + 1;

            shared_ptr<Grid> grid = Grid::get( nPixels );
            shared_ptr<double*> a = grid->angle2D();
            shared_ptr<double*> d = grid->dist2D();
            double** aPtr = a.get();
            double** dPtr = d.get();

            Pupil pupil( nPixels, rc );
            
            Array<double> awrapper( *aPtr, nPixels, nPixels );
            redux::file::Ana::write("modetest_angle.f0", awrapper);
            awrapper.wrap( *dPtr, nPixels, nPixels );
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
            double** ap = makePointers(aperture.ptr(), nPixels, nPixels);
            double area = makePupil(ap,nPixels,radius);
            //redux::file::Ana::write("aperture.f0", aperture);

            //double area_mvn = makePupil_mvn(ap,nPixels,radius);
            //redux::file::Ana::write("aperture2.f0", aperture);
            delPointers(ap);

            Pupil pup(nPixels,radius);
            //redux::file::Ana::write("aperture3.f0", pup_pair.first);

            cout << "Area = " << area << "   r^2*PI = " << (radius*radius*M_PI) << endl;
        //     cout << "ap=f0('/home/tomas/build/redux/aperture.f0')\n"
        //          << "ap2=f0('/home/tomas/build/redux/aperture2.f0')\n"
        //          << "ap[" <<(nPixels/2-2) << ":" << (nPixels/2+2) << ","<<(nPixels/2-2) << ":" << (nPixels/2+2) << "]\n"
        //          << "ap2[" <<(nPixels/2-1) << ":" << (nPixels/2+3) << ","<<(nPixels/2-1) << ":" << (nPixels/2+3) << "]\n"
        //          << "diff = ap-ap2[1:" << nPixels << ",1:" << nPixels << "]\n"
        //          << "print,min(diff),max(diff),mean(diff)" << endl;

        }


        using namespace boost::unit_test;
        void add_zernike_tests( test_suite* ts ) {

            ts->add( BOOST_TEST_CASE_NAME( &zernikeTest, "Zernike" ) );

        }

    }

}
