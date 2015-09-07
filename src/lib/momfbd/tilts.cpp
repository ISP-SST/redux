#include "redux/momfbd/tilts.hpp"

#include "redux/momfbd/channel.hpp"
#include "redux/job.hpp"

#include "redux/util/stringutil.hpp"
#include "redux/file/fileana.hpp"

using namespace redux::momfbd;
using namespace redux::file;
using namespace redux::util;
using namespace redux;
using namespace std;

#define lg Logger::mlg
namespace {

    const std::string thisChannel = "tilts";

}


Tilts::Tilts(Channel& c, uint16_t m) : channel(c), mode(m), nFreeAlpha(1), partialGrad(nullptr), anchorChannel(false) {

}

Tilts::~Tilts() {
    delete[] partialGrad;
}


size_t Tilts::nFreeParameters(void) const {
    size_t ret = nFreeAlpha;
    for( auto& it: relativeTilts ) ret += it->nFreeParameters();
    return ret;
}


size_t Tilts::setPointers(double* a, double* g ) {
    size_t ret = nFreeAlpha;
    x_  = a;
    df_ = g;
    for( auto& it: relativeTilts ) ret += it->setPointers(a+ret, g+ret);
    return ret;
}


void Tilts::init(void) {
    phi.resize(channel.imageNumbers.size(),channel.pupilPixels,channel.pupilPixels);
    otf.resize(channel.imageNumbers.size(),2*channel.pupilPixels,2*channel.pupilPixels);
    if(!anchorChannel) {     // only for relative tilts/channels
        partialGrad = new double[channel.imageNumbers.size()];
    }
}


void Tilts::addRelativeTilt(std::shared_ptr<Tilts>& t) {
    t->refTilt = shared_from_this();
    relativeTilts.push_back(t);
}


void Tilts::applyTiltToImages(boost::asio::io_service& service) {

    for( size_t i=0; i<channel.subImages.size(); ++i ) {
        service.post ([i,this] {
            channel.subImages[i]->setAlpha( mode, x_[i] );              // anchor-channel apply tilts to images
        });
    }
    
    for( auto& it: relativeTilts ) {
        it->addTiltToImages(service, x_);
    }

}


void Tilts::addTiltToImages(boost::asio::io_service& service, double* ref) {
    for( size_t i=0; i<channel.subImages.size(); ++i ) {
        service.post ([ref,i,this] {
            channel.subImages[i]->setAlpha( mode, ref[i] + *x_ );       // non-anchor, apply "anchor-alpha + delta" to images
        });
    }
}


void Tilts::calcGradient(boost::asio::io_service& service, uint8_t nThreads, const grad_t& gradient) {
    
    if( !anchorChannel ) {   // this is a single offset value, to be paired with each tilt in refTilt
        for( size_t i=0; i<channel.subImages.size(); ++i ) {
            service.post ([this,gradient,i] {
                partialGrad[i] = gradient(*channel.subImages[i], mode, 1E-2); //, otf.ptr(i,0,0), phi.ptr(i,0,0) );
            });
        }
    } else {        // this is the anchor channel
        for( auto& it: relativeTilts ) {
            it->calcGradient(service, nThreads, gradient);
        }
        runThreadsAndWait( service, nThreads );
        for( uint i=0; i<channel.subImages.size(); ++i ) {
            service.post ([this,gradient,i] {
                df_[i] = gradient(*channel.subImages[i], mode, 1E-2); //, otf.ptr(i,0,0), phi.ptr(i,0,0) );
               // cout << "   DFI = " << df_[i] << endl;
                Array<complex_t> tmpOTF(otf,i,i,0,otf.dimSize(1)-1,0,otf.dimSize(2)-1);
                //Ana::write ("gm"+to_string(mode)+"_otf.f0", tmpOTF);
                Array<double> tmpPhi(phi,i,i,0,phi.dimSize(1)-1,0,phi.dimSize(2)-1);
                //Ana::write ("gm"+to_string(mode)+"_phi.f0", tmpPhi);

            });
        }
        runThreadsAndWait( service, nThreads );
        for(uint i=0; i<relativeTilts.size(); ++i ) {
            relativeTilts[i]->addPartials();
            for( uint j=0; j<channel.subImages.size(); ++j ) df_[j] += relativeTilts[i]->partialGrad[j];
        }
        double sum=0.0;
        double minVal=fabs(df_[0]);
        size_t minIndex=0;
        for( uint i=0; i<channel.subImages.size(); ++i ) {
            if( fabs(df_[i]) < minVal ) {
                minVal=fabs(df_[i]);
                minIndex = i;
            }
            sum += df_[i];
        }
       // df_[minIndex] -= sum;
//        cout << endl << printArray(df_,channel.subImages.size(),"rawtiltgrad", 10) << endl;
        sum /= channel.subImages.size();
        for( uint i=0; i<channel.subImages.size(); ++i ) {
            df_[i] -= sum;
        }
//        cout << endl << printArray(df_,channel.subImages.size(),"tiltgrad", 10) << endl;
    }
    
}


double Tilts::getPartial(size_t i) {
    if(i < channel.subImages.size() && !anchorChannel ) {
        return partialGrad[i];
    }
    return 0;
}


void Tilts::addPartials(void) {
    if( !anchorChannel ) {
        for(uint i=0; i<channel.subImages.size(); ++i) {
            *df_ += partialGrad[i];
        }
    }
}
