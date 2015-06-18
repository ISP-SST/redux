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
    zeroValues();
}


void WaveFront::addImage (std::shared_ptr<SubImage> im) {
    if (im) {
        images.insert (im);
        im->setWaveFront (shared_from_this());
    }

    cout << "  WaveFront::addImage()  im = " << hexString (im.get()) << "   nIm = " << images.size() << endl;
}


void WaveFront::addWeight (uint16_t modenumber, double w) {
    alpha[modenumber].weight += w;
    cout << "WaveFront::addWeight():   modeAddr = " << hexString (&alpha[modenumber]) << "   Mode: " << modenumber << "   w: " << w << "   weight: " << alpha[modenumber].weight << endl;
}


void WaveFront::zeroWeights (void) {
    for (auto & it : alpha) {
        it.second.weight = 0;
    }
}


void WaveFront::zeroValues (void) {
    for (auto & it : alpha) {
        it.second.value = 0;
    }
}


void  WaveFront::setAlpha (const alpha_map& a) {
    alpha = a;
}


void WaveFront::computePhases (boost::asio::io_service& service) {
    for (const shared_ptr<SubImage>& it : images) {
        if (it) {
            service.post (std::bind (&SubImage::computePhases, it.get()));
//             it->computePhases();
//             it->OTF();
//             it->PSF();
        }
    }
}


void WaveFront::computePhasesAndOTF (boost::asio::io_service& service) {
    for (const shared_ptr<SubImage>& it : images) {
        if (it) {
            service.post ([it] {    // use a lambda to ensure these two calls are sequential
                it->computePhases();
                it->calcOTF();
            });
        }
    }
}


double WaveFront::coefficientMetric (void) {

    double sum = 0.0;

    for (auto & a : alpha) {
        sum += a.second.metric();
    }

    cout << hexString (this) << "  WaveFront::coefficientMetric()  nAlpha = " << alpha.size() << "   sum = " << sum << endl;
    return sum;

}

void WaveFront::gradientTest (void) {
    cout << hexString (this) << "  WaveFront::gradientTest()  nAlpha = " << alpha.size() << "  nImages = " << images.size() << endl;
    vector<double> grad (alpha.size(), 0.0);

    for (const shared_ptr<SubImage>& it : images) {
        if (it) {
            cout << hexString (this) << "  WaveFront::gradientTest()  it != NULL" << endl;
            it->oldGradientDiff (grad);
        }
    }

    cout << hexString (this) << "  WaveFront::gradientTest() " << printArray (grad, "grad") << endl;
}
