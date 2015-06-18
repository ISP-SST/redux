#include "redux/momfbd/wavefront.hpp"

#include "redux/util/stringutil.hpp"

using namespace redux::momfbd;
using namespace redux::util;
using namespace std;

#define lg Logger::mlg
namespace {

    const std::string thisChannel = "wavefront";

}

void WaveFront::clear (void) {
    cout << "  WaveFront::clear()" << endl;
    images.clear();
}


void WaveFront::addImage (std::shared_ptr<SubImage> im) {
    if (im) {
        images.insert (im);
        im->setWaveFront (shared_from_this());
    }

    cout << hexString (this) << "  WaveFront::addImage()  im = " << hexString (im.get()) << "   nIm = " << images.size() << endl;
}


void WaveFront::addWeight (uint16_t modenumber, double w) {
    modes[modenumber].weight += w;
    cout << hexString (this) << "  WaveFront::addWeight():   modeAddr = " << hexString (&modes[modenumber]) << "   Mode: " << modenumber << "   w: " << w << "   weight: " << modes[modenumber].weight << endl;
}


void WaveFront::zeroAlphaWeights (void) {
    for (auto & it : modes) {
        it.second.weight = 0;
    }
}


void  WaveFront::setAlpha (const modeinfo_map& a) {
    modes = a;
}



void WaveFront::computePhases (boost::asio::io_service& service) {
    for (const shared_ptr<SubImage>& it : images) {
        if (it) {
            service.post (std::bind((void(SubImage::*)(const double*))&SubImage::addPhases, it.get(), alpha));
//             it->computePhases();
//             it->OTF();
//             it->PSF();
        }
    }
}


void WaveFront::computePhasesAndOTF (boost::asio::io_service& service) {
    for (const shared_ptr<SubImage>& it : images) {
        if (it) {
            service.post ([it,this] {    // use a lambda to ensure these two calls are sequential
                it->addPhases(alpha);
                it->calcOTF();
            });
        }
    }
}


double WaveFront::coefficientMetric (void) {

    double sum = 0.0;
    size_t cnt(0);
    for (auto& a : modes) {
        sum += a.second.weight*alpha[cnt]*alpha[cnt];
        cnt++;
    }

    cout << hexString (this) << "  WaveFront::coefficientMetric()  nAlpha = " << modes.size() << "   sum = " << sum << endl;
    return sum;

}

void WaveFront::gradientTest (void) {
    cout << hexString (this) << "  WaveFront::gradientTest()  nAlpha = " << modes.size() << "  nImages = " << images.size() << endl;
    vector<double> grad (modes.size(), 0.0);

    for (const shared_ptr<SubImage>& it : images) {
        if (it) {
            //it->gradientFiniteDifference (grad);
          //  it->gradientVogel(grad.data());
        } else cout << hexString (this) << "  WaveFront::gradientTest()  it == NULL" << endl;
    }

    cout << hexString (this) << "  WaveFront::gradientTest() " << printArray (grad, "grad") << endl;
}
