#include "redux/momfbd/wavefront.hpp"

#include "redux/util/stringutil.hpp"

using namespace redux::momfbd;
using namespace redux::util;
using namespace std;

#define lg Logger::mlg
namespace {

    const std::string thisChannel = "wavefront";

}

WaveFront::WaveFront(void) : nImages(0) {

}


void WaveFront::init(size_t pupilSize) {
    phi.resize(modes.size(),pupilSize,pupilSize);
    otf.resize(modes.size(),2*pupilSize,2*pupilSize);
}


void WaveFront::reset (void) {
    images.clear();
}


size_t WaveFront::setPointers(double* a, double* g ) {
    alpha=a;
    grad=g;
    for ( auto& it: modes ) {
        it.second.value = a++;
    }
    return modes.size();
}


void WaveFront::addImage (std::shared_ptr<SubImage> im) {
    if (im) {
        images.insert (im);
        im->setWaveFront (shared_from_this());
    }
}


void WaveFront::addWeight (uint16_t modenumber, double w) {
    modes[modenumber].weight += w;
    //cout << hexString (this) << "  WaveFront::addWeight():   modeAddr = " << hexString (&modes[modenumber]) << "   Mode: " << modenumber << "   w: " << w << "   weight: " << modes[modenumber].weight << endl;
}


void WaveFront::zeroAlphaWeights (void) {
    for (auto & it : modes) {
        it.second.weight = 0;
    }
}


void  WaveFront::setAlpha (const modeinfo_map& a) {
    modes = a;
}


void WaveFront::setAlphas(boost::asio::io_service& service) {
    for (const shared_ptr<SubImage>& it : images) {
        if (it) {
            service.post (std::bind((void(SubImage::*)(const double*))&SubImage::setAlphas, it.get(), alpha));
            //service.post (std::bind((void(SubImage::*)(const double*))&SubImage::addPhases, it.get(), alpha)); // FIXME: for testing only!
        }
    }
}


void WaveFront::setAlphasAndUpdate(boost::asio::io_service& service, bool calcVogel) {
    for (const shared_ptr<SubImage>& it : images) {
        if (it) {
            service.post ([this, it, calcVogel] {    // use a lambda to ensure these calls are sequential
                it->setAlphas(alpha);
                //it->addPhases(alpha);   // FIXME: this is for testing only
                it->update(calcVogel);
            });
        }
    }
}


void WaveFront::calcGradient (boost::asio::io_service& service, const grad_t& gradient) {
    int count(0);
    for ( auto& m: modes ) {
        service.post ([this,gradient,count,&m] {
            for (const shared_ptr<SubImage>& it : images) {
                grad[count] += gradient(*it, m.first, 1E-2); //, otf.ptr(count,0,0), phi.ptr(count,0,0) );
            }
        });
        count++;
    }
}


double WaveFront::metric(boost::asio::io_service& service) {
    for (const shared_ptr<SubImage>& it : images) {
        if (it) {
            service.post ([it,this] {    // use a lambda to ensure these two calls are sequential
                //it->addScaledPhases(alpha);
                it->setAlphas(alpha);   // FIXME: for testing only !!
                it->calcOTF();
            });
        }
    }
    return 0;
}


double WaveFront::metric(boost::asio::io_service& service, const double* dAlpha) {
    
    return 0;
}


double WaveFront::coefficientMetric (void) {

    double sum = 0.0;
    size_t cnt(0);
    for (auto& a : modes) {
        sum += a.second.weight*alpha[cnt]*alpha[cnt];
        cnt++;
    }

    return sum;

}

