#ifndef REDUX_MOMFBD_UTIL_HPP
#define REDUX_MOMFBD_UTIL_HPP

#include "redux/util/array.hpp"

#include <cmath>
#include <string>

namespace redux {

    namespace momfbd {
        
        namespace util {

            template <typename T>
            inline T sqr(const T& a) { return a * a; }
            inline double sign(double a) { return (a >= 0.0) ? 1 : -1; }
            inline double sign(double a, double b) { return (b >= 0.0) ? fabs(a) : -fabs(a); }

            int64_t extractTime(const char* header);
            void loadPupil(const std::string&, redux::util::Array<double>&, uint32_t pupilSize=0);
             
            double cf2pix(double arcsecperpix,double telescope_d);
            double pix2cf(double arcsecperpix,double telescope_d);
            double def2cf(double pd_defocus,double telescope_r);
            double cf2def(double alpha,double telescope_r);

            template <typename T>
            std::vector<T> segment(T first, T last, T segmentLength, T minimumOverlap=0) {
                int nSegments=2;
                double separation = (last-first)/static_cast<double>(nSegments-1);
                double overlap = std::max(segmentLength-separation,0.0);
                while(overlap < minimumOverlap) {
                    ++nSegments;
                    separation = (last-first)/static_cast<double>(nSegments-1);
                    overlap = std::max(segmentLength-separation,0.0);
                }
                std::vector<T> ret;
                for(size_t i = 0; i < nSegments; ++i) {
                    ret.push_back(static_cast<T>(i*separation+first) );
                }
                return ret;
            }
        }

    }

}

#endif // REDUX_MOMFBD_UTIL_HPP
