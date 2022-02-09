#include "redux/momfbd/object.hpp"

#include "redux/momfbd/momfbdjob.hpp"
#include "redux/momfbd/channel.hpp"
#include "redux/momfbd/data.hpp"
#include "redux/momfbd/util.hpp"

#ifdef DEBUG_
#   define TRACE_THREADS
#endif

#include "redux/file/fileana.hpp"
#include "redux/file/filefits.hpp"
#include "redux/file/filemomfbd.hpp"
#include "redux/image/utils.hpp"
#include "redux/image/zernike.hpp"
#include "redux/math/functions.hpp"
#include "redux/translators.hpp"
#include "redux/util/cache.hpp"
#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/util/trace.hpp"
#include "redux/constants.hpp"
#include "redux/logging/logger.hpp"
#include "redux/revision.hpp"

#include <cstdio>
#include <limits>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/info_parser.hpp>

namespace bfs = boost::filesystem;
using namespace redux::file;
using namespace redux::image;
using namespace redux::logging;
using namespace redux::math;
using namespace redux::momfbd;
using namespace redux::util;
using namespace redux;
using namespace std;
using boost::algorithm::iequals;


namespace {
    
}

Object::Object( MomfbdJob& j, uint16_t id ): ObjectCfg(j), myJob(j), logger(j.logger), currentMetric(0), reg_gamma(0),
    frequencyCutoff(0),pupilRadiusInPixels(0), patchSize2(0), otfSize(0), otfSize2(0),
    ID(id), traceID(-1), normalizeTo(0), imgSize(0), nObjectImages(0),
    startT(bpx::not_a_date_time), endT(bpx::not_a_date_time) {

}


Object::Object( const Object& rhs, uint16_t id, int tid ) : ObjectCfg(rhs), myJob(rhs.myJob), logger(rhs.logger),
    channels(rhs.channels), currentMetric(rhs.currentMetric), reg_gamma(rhs.reg_gamma),
    frequencyCutoff(rhs.frequencyCutoff), pupilRadiusInPixels(rhs.pupilRadiusInPixels),
    patchSize2(rhs.patchSize2), otfSize(rhs.otfSize), otfSize2(rhs.otfSize2),
    ID (id), traceID(tid), normalizeTo(rhs.normalizeTo), imgSize(rhs.imgSize), nObjectImages(rhs.nObjectImages),
    startT(rhs.startT), endT(rhs.endT) {

}


Object::~Object() {

    THREAD_MARK
    cleanup( );
    THREAD_MARK
}


void Object::parsePropertyTree( bpt::ptree& tree, redux::logging::Logger& logger ){

    ObjectCfg::parseProperties(tree, logger, myJob );

    uint16_t nCh(0 );
    for( auto & property : tree ){
        if( iequals( property.first, "CHANNEL" ) ){
            Channel* tmpCh = new Channel( *this, myJob, nCh++ );
            tmpCh->parsePropertyTree( property.second, logger );
            channels.push_back( shared_ptr<Channel>( tmpCh ) );
        }
    }

}


bpt::ptree Object::getPropertyTree( bpt::ptree& tree, bool showAll ){

    bpt::ptree node;

    for( auto& ch : channels ){
        ch->getPropertyTree( node, showAll );
    }

    ObjectCfg::getProperties(node, myJob);

    tree.push_back( bpt::ptree::value_type( "object", node ) );

    return node;

}


size_t Object::size(void )const {
    size_t sz = ObjectCfg::size( );
    sz += 2*sizeof(uint16_t ) + 4*sizeof(double) + sizeof(traceID);                   // channels.size( )+ ID + maxMean
    for( const auto& ch : channels ){
        sz += ch->size( );
    }
    sz += imgSize.size( );
    sz += sizeof( nObjectImages );
    return sz;
}


uint64_t Object::pack(char* ptr )const {
    using redux::util::pack;
    uint64_t count = ObjectCfg::pack(ptr );
    count += pack(ptr+count, ID );
    count += pack(ptr+count, traceID );
    count += pack(ptr+count, currentMetric );
    count += pack(ptr+count, frequencyCutoff );
    count += pack(ptr+count, pupilRadiusInPixels );
    count += pack(ptr+count, normalizeTo );
    count += imgSize.pack(ptr+count );
    count += pack(ptr+count,( uint16_t)channels.size() );
    for( const auto& ch : channels ){
        count += ch->pack(ptr+count );
    }
    count += pack( ptr + count, nObjectImages );
    if( count != size() ){
        LOG_ERR << "(" << hexString( this )<< "): Packing failed, there is a size mismatch:  count = " << count << "  sz = " << size( )<< ende;
    }
    return count;
}


uint64_t Object::unpack(const char* ptr, bool swap_endian ){
    using redux::util::unpack;

    uint64_t count = ObjectCfg::unpack(ptr, swap_endian );
    count += unpack(ptr+count, ID, swap_endian );
    count += unpack(ptr+count, traceID, swap_endian );
    count += unpack(ptr+count, currentMetric, swap_endian );
    count += unpack(ptr+count, frequencyCutoff, swap_endian );
    count += unpack(ptr+count, pupilRadiusInPixels, swap_endian );
    count += unpack(ptr+count, normalizeTo, swap_endian );
    count += imgSize.unpack(ptr+count, swap_endian );
    uint16_t tmp;
    count += unpack(ptr+count, tmp, swap_endian );
    channels.resize(tmp);
    for( auto& ch : channels ){
        ch.reset(new Channel(*this,myJob) );
        count += ch->unpack(ptr+count, swap_endian );
    }
    count += unpack( ptr + count, nObjectImages, swap_endian );
    return count;
}


void Object::cleanup(void ){

    THREAD_MARK
    channels.clear( );
    ftSum.reset();
    Q.reset();
    P.reset();
    PQ.reset();
    PS.reset();
    QS.reset();
    fittedPlane.clear( );
    pupil.reset( );
    modes.reset( );

    THREAD_MARK
    if( !cacheFile.empty() ) {
        bfs::path tmpP(cacheFile);
        if( bfs::exists(tmpP) ) {
            try {
                bfs::remove(tmpP);
            } catch( exception& e ) {
                LOG_ERR << "Object::cleanup: failed to remove cacheFile: " << cacheFile << "  reason: " << e.what() << endl;
            }
        }
    }
    THREAD_MARK

}


uint32_t Object::nImages( bool reCalc ) {
    if( reCalc ) nObjectImages = 0;
    if( nObjectImages ) return nObjectImages;
    for( const auto& ch : channels ) nObjectImages += ch->nImages( );
    return nObjectImages;
}


void Object::initProcessing( Solver& ws ){

    if( patchSize && pupilPixels ){
        patchSize2 = patchSize*patchSize;
        otfSize = 2*pupilPixels;
        otfSize2 = otfSize*otfSize;
        P = rdx_get_shared<complex_t>(otfSize2);
        Q = rdx_get_shared<double>(otfSize2);
        PQ = rdx_get_shared<complex_t>(otfSize2);
        PS = rdx_get_shared<double>(otfSize2);
        QS = rdx_get_shared<double>(otfSize2);
        ftSum = rdx_get_shared<double>(otfSize2);
        
        // basically just to shut valgrind up about "uninitialized values"
        std::fill_n( P.get(), otfSize2, complex_t(0) );
        std::fill_n( Q.get(), otfSize2, double(0) );
        std::fill_n( PQ.get(), otfSize2, complex_t(0) );
        std::fill_n( PS.get(), otfSize2, double(0) );
        std::fill_n( QS.get(), otfSize2, double(0) );
        std::fill_n( ftSum.get(), otfSize2, double(0) );
        
        if( myJob.runFlags & RF_FIT_PLANE ){
            fittedPlane.resize( patchSize, patchSize );
            fittedPlane.zero();
        }
        
        
        for( auto& ch : channels ){
            ch->initProcessing( ws );
        }
        initObject( );            // load global pupil/modes

        double mode_scale = 1.0/wavelength;
        ModeInfo info( modeFile, pupilPixels );
        if( modeFile.empty() ) {
            info = ModeInfo( myJob.klMinMode, myJob.klMaxMode, myJob.modeList, pupilPixels, pupilRadiusInPixels, rotationAngle, myJob.klCutoff );
        }
        
        modes = myJob.globalData->get( info, mode_scale, pupil, modes );    // get properly scaled modes.


        if( modes->empty() ) {
            throw std::logic_error("The mode-set is empty.");
        } 
        auto& mdims = modes->dimensions();
        if( mdims.size() != 3 || mdims[1] != myJob.pupilPixels || mdims[2] != myJob.pupilPixels ) {
            Cache::clear< std::pair<ModeInfo, double>, shared_ptr<ModeSet>>();
            throw std::logic_error("The mode-set has the wrong dimensions.");
        }
    
        auto& pdims = pupil->dimensions();
        if( pdims.size() != 2 || pdims[0] != myJob.pupilPixels || pdims[0] != myJob.pupilPixels ) {
            throw std::logic_error("The pupil has the wrong dimensions.");
        }
    
        modes->tiltMode = -1;   // FIXME: Why do we need to force a recalc? Jacobian should be sent from master.
        modes->measureJacobian( *pupil, wavelength*(0.5*frequencyCutoff)/util::pix2cf( arcSecsPerPixel, myJob.telescopeD ) );

    } else {
        LOG_ERR << "Object patchSize is 0 !!!" << ende;
    }
    
}


void Object::initPatch( void ){
    unique_lock<mutex> lock( mtx );
    reg_gamma = 0;
    std::fill_n( ftSum.get(), otfSize2, 0.0 );
}


void Object::restorePatch( ObjectData& od, const vector<uint32_t>& wf ) {

    unique_lock<mutex> lock( mtx );
    set<uint32_t> wfSet( wf.begin(), wf.end() );
    vector<shared_ptr<SubImage>> images;
    for( const shared_ptr<Channel>& ch : channels ) {
        if( ch->noRestore ) continue;
        for( size_t i=0; i<ch->waveFrontList.size(); ++i) {
            if( wfSet.count( ch->waveFrontList[i] ) && ch->subImages[i] ) {
                images.push_back( ch->subImages[i] );
            }
        }
    }
    
    if( images.empty() ) return;
    
    size_t nModes = myJob.nModes;
    
    size_t otfSize = pupilPixels<<1;
    
    FourierTransform avgObjFT(otfSize, otfSize, FULLCOMPLEX|REORDER_FT); //|NORMALIZE_FT );
    Array<complex_t> tmpC(otfSize, otfSize);
    Array<double> tmpD(otfSize, otfSize);
    avgObjFT.zero();
    tmpD.zero();
    complex_t* aoPtr = avgObjFT.get();
    double* dPtr = tmpD.get();
    double avgNoiseVariance = 0.0;
    PointD avgShift;

    uint32_t imgCount(0);
    for( shared_ptr<SubImage>& si: images ) {
        if( si ) {
            si->restore( aoPtr, dPtr );
            avgNoiseVariance += sqr(si->stats.noise);
            avgShift += si->currentShift;
            imgCount++;
        }
    }
    avgObjFT.conj( );    // This is because we re recycling addPQ in SubImage.restore(),
                        // which actually returns the complex conjugate of the deconvolution.
    avgNoiseVariance /= imgCount;
    avgShift *= 1.0/imgCount;

    avgObjFT.safeDivide( tmpD );     // TBD: non-zero cutoff by default? based on reg_gamma ?
    
    if (!(myJob.runFlags&RF_NO_FILTER)) {
        LOG_DEBUG << boost::format("Object %d: Applying Scharmer filter with frequency-cutoff = %g and noise-variance = %g") % ID % (0.9*frequencyCutoff) % avgNoiseVariance << ende;
        ScharmerFilter( aoPtr, dPtr, otfSize, otfSize, avgNoiseVariance, 0.90 * frequencyCutoff);
    }

    od.img.resize( otfSize, otfSize );
    avgObjFT.ift( od.img.get( ) );
    
    // PSF
    uint16_t nPSF = 0;
    if( saveMask & (SF_SAVE_PSF|SF_SAVE_PSF_AVG) ) {
        nPSF =( saveMask&SF_SAVE_PSF_AVG)? 1 : imgCount;
        od.psf.resize(nPSF, otfSize, otfSize );
        od.psf.zero( );
        Array<float> view( od.psf, 0, 0, 0, otfSize-1, 0, otfSize-1 );
        tmpD.zero( );
        if( nPSF > 1 ){
            for( shared_ptr<SubImage>& si: images ) {
                si->getPSF( tmpD.get( ) );
                view.assign( tmpD );
                view.shift( 0, 1 );
            }
        } else if( nPSF == 1 ){
            for( shared_ptr<SubImage>& si: images ) {
                si->addPSF( tmpD.get( ) );
            }
            view.assign( tmpD );
            od.psf *=( 1.0/imgCount );
        }
    } else {
        od.psf.clear();
    }

    // Convolved objects
    if( saveMask & SF_SAVE_COBJ ){
        od.cobj.resize(imgCount, otfSize, otfSize );
        od.cobj.zero( );
        Array<float> view(od.cobj,0,0,0,otfSize-1,0,otfSize-1 );
        for( shared_ptr<SubImage>& si: images ) {
            view.assign( si->convolveImage( od.img) );
            view.shift(0,1 );
        }
    } else {
        od.cobj.clear( );
    }

    // Residuals
    if( saveMask & SF_SAVE_RESIDUAL ){
        od.res.resize( imgCount, otfSize, otfSize );
        od.res.zero( );
        Array<float> view( od.res, 0, 0, 0, otfSize-1, 0, otfSize-1 );
        if( od.cobj.sameSizes(od.res ) ){
            Array<float> cview(od.cobj,0,0,0,otfSize-1,0,otfSize-1 );
            for( shared_ptr<SubImage>& si: images ) {
                view.assign( si->convolvedResidual(cview) );
                view.shift( 0, 1 );
                cview.shift( 0, 1 );
            }
        } else {
            for( shared_ptr<SubImage>& si: images ) {
                view.assign( si->residual(od.img) );
                view.shift( 0, 1 );
            }
        }
    } else {
        od.res.clear( );
    }
    
    if( otfSize > patchSize ) {
        size_t offset = (otfSize - patchSize)/2;
        od.img.setLimits( offset, offset+patchSize-1, offset, offset+patchSize-1 );
        od.img.trim();
        if( nPSF ) {
            od.psf.setLimits( 0, nPSF-1, offset, offset+patchSize-1, offset, offset+patchSize-1 );
            od.psf.trim(false);
        }
        if( saveMask & SF_SAVE_COBJ ) {
            od.cobj.setLimits( 0, imgCount-1, offset, offset+patchSize-1, offset, offset+patchSize-1 );
            od.cobj.trim(false);
        }
        if( saveMask & SF_SAVE_RESIDUAL ) {
            od.res.setLimits( 0, imgCount-1, offset, offset+patchSize-1, offset, offset+patchSize-1 );
            od.res.trim(false);
        }
    }

    // Note: re-adding the plane has to be done *after* calculating the convolved objects and residuals
    if( fittedPlane.sameSize(od.img) ) {
        LOG_DEBUG << "Object " << to_string(ID) << ": re-adding fitted plane to result." << ende;
        od.img += fittedPlane;
    } else if( !fittedPlane.empty() ) {
        LOG_WARN << "Size mismatch when re-adding fitted plane." << ende;
    }
    
    // Mode coefficients
    if( saveMask & SF_SAVE_ALPHA ){
        od.alpha.resize( imgCount, nModes );
        float* aPtr = od.alpha.get();
        for( shared_ptr<SubImage>& si: images ) {
            std::copy( si->wfAlpha, si->wfAlpha+nModes, aPtr );
            aPtr += nModes;
        }
    } else {
        od.alpha.clear( );
    }

    // Diversity
    int nCh = channels.size( );
    if( nCh && (saveMask & SF_SAVE_DIVERSITY) ) {
        od.div.resize( nCh, pupilPixels, pupilPixels );
        Array<float> view(od.div,0,0,0,pupilPixels-1,0,pupilPixels-1);
        for( auto& ch: channels ){
            view = ch->phi_fixed;
            view.shift(0,1);
        }
    } else {
        od.div.clear( );
    }

}


void Object::getInit(ObjectData& od, double* alpha ){

    lock_guard<mutex> lock( mtx );
    // Mode coefficients
    size_t nAlpha = od.alpha.nElements( );
    if( nAlpha > 0 ){
        std::copy( od.alpha.get(), od.alpha.get()+nAlpha, alpha );
    }


}


void Object::initPQ( void ){
    std::fill_n( P.get(), otfSize2, complex_t(0) );
    std::fill_n( Q.get(), otfSize2, reg_gamma );
}


void Object::addRegGamma( double rg ){
    lock_guard<mutex> lock( mtx );
    reg_gamma += 0.10*rg/nObjectImages;
}


void Object::addToFT( const complex_t* ft ){
    lock_guard<mutex> lock( mtx );
    double* ftSumPtr = ftSum.get( );
    transform(ftSumPtr, ftSumPtr+otfSize2, ft, ftSumPtr,
              [](const double& a, const complex_t& b ){ return a+norm(b ); }
             );

}


void Object::addDiffToFT( const complex_t* newFT, const complex_t* oldFT ){
    
    lock_guard<mutex> lock( mtx );
    double* ftSumPtr = ftSum.get( );
    for(size_t ind=0; ind<otfSize2; ++ind ){
        ftSumPtr[ind] +=( norm(newFT[ind])-norm(oldFT[ind]) );
    }
    
}


void Object::addDiffToPQ(const redux::image::FourierTransform& ft, const Array<complex_t>& otf, const Array<complex_t>& oldotf ){

    lock_guard<mutex> lock( mtx );
    double *qPtr = Q.get( );
    complex_t *pPtr = P.get( );
    const complex_t *ftPtr = ft.get( );
    const complex_t *otfPtr = otf.get( );
    const complex_t *ootfPtr = oldotf.get( );

    for( auto& ind: pupil->otfSupport ){
        qPtr[ind] += norm(otfPtr[ind] )- norm(ootfPtr[ind] );
        pPtr[ind] += conj(ftPtr[ind] )*( otfPtr[ind] - ootfPtr[ind] );
    }
}


void Object::addToPQ(const complex_t* pp, const double* qq ){

    lock_guard<mutex> lock( mtx );
    double *qPtr = Q.get( );
    complex_t *pPtr = P.get( );
    for( auto& ind: pupil->otfSupport ){
        qPtr[ind] += qq[ind];
        pPtr[ind] += pp[ind];
    }
    
}


void Object::addAllPQ(void ){
    for( auto& ch : channels ){
        for( auto& im : ch->subImages ){
            lock_guard<mutex> lock( mtx );
            im->addPQ( P.get(),Q.get() );
        }
    }
}


void Object::calcHelpers(void ){
    
    lock_guard<mutex> lock( mtx );
    std::fill_n( PQ.get(), otfSize2, complex_t(0) );
    std::fill_n( PS.get(), otfSize2, 0.0 );
    std::fill_n( QS.get(), otfSize2, 0.0 );
    
    const double *qPtr = Q.get( );
    const complex_t *pPtr = P.get( );
    complex_t *pqPtr = PQ.get( );
    double *psPtr = PS.get( );
    double *qsPtr = QS.get( );

    for( auto& ind: pupil->otfSupport ){
        pqPtr[ind] = pPtr[ind] * qPtr[ind];
        psPtr[ind] = norm( pPtr[ind] );
        qsPtr[ind] = qPtr[ind] * qPtr[ind];
    }
    
}


void Object::fitAvgPlane( redux::util::Array<float>& plane, const vector<uint32_t>& wf ){
    
    set<uint32_t> wfSet( wf.begin(), wf.end() );
    vector<shared_ptr<SubImage>> images;
    for( const shared_ptr<Channel>& ch : channels ) {
        if( !ch ) continue;
        for( size_t i=0; i<ch->waveFrontList.size(); ++i) {
            if( wfSet.count( ch->waveFrontList[i] ) && ch->subImages[i] ) {
                images.push_back( ch->subImages[i] );
            }
        }
    }
    
    if( (myJob.runFlags & RF_FIT_PLANE) && !images.empty() ){
        size_t count(0);
        plane.zero();
        for( const shared_ptr<SubImage>& im: images ){
            if( !im )continue;
            if( !im->sameSize(plane) ){
                LOG_ERR << "Size mismatch when fitting average plane for object #" << ID
                        << printArray(im->dimensions(),"  imdims") << printArray(plane.dimensions(),"  pdims") << ende;
                plane.clear( );
                return;
            }
            plane += *im;
            count++;
        }
        if( count == 0 ){
            LOG_ERR << "No images present when fitting average plane for object #" << ID << ende;
            plane.clear( );
            return;
        }
        
        plane /= count;
        vector<double> coeffs(3,0.0);
        plane = fitPlane( plane, true, coeffs.data( ) );          // fit plane to the average image, and subtract average
        LOG_DEBUG << "Fitting average plane:  p = " << coeffs[0] << "x + " << coeffs[1] << "y + " << coeffs[2] << ende;
    } else {
        plane.clear( );
    }
    
}


void Object::calcMetric( void ){

    if( weight == 0.0 )return;
    
    const double* ftsPtr = ftSum.get( );
    const complex_t* pPtr = P.get( );
    const double* qPtr = Q.get( );
    
    lock_guard<mutex> lock( mtx );
    currentMetric = 0;
    //size_t N = 4*pupilPixels*pupilPixels;
    //for( size_t ind=0; ind<N; ++ind ){
    for( auto & ind : pupil->otfSupport ){
        currentMetric +=( ftsPtr[ind] - norm( pPtr[ind] )/ qPtr[ind] );
    }

    currentMetric /=( patchSize*patchSize );

}


bool Object::checkImageScale( double& F, double& A, double& P ) {

    double rad2asec = 180.0 * 3600.0 / M_PI;
    size_t count = F > 0 ? 1 : 0;
    count += A > 0 ? 1 : 0;
    count += P > 0 ? 1 : 0;
    if( count > 2 ){
        LOG_WARN << "Too many parameters specified: replacing telescope focal length (=" << F
                    << ") with computed value (=" <<( P * rad2asec / A )<< ")" << ende;
        F = P * rad2asec / A;
        return true;
    } else if( count < 2 ){
        LOG_ERR << "At least two of the parameters \"TELESCOPE_F\", \"ARCSECPERPIX\" and \"PIXELSIZE\" has to be provided." << ende;
    } else {    // count == 2
        if( F <= 0 ){
            F = P * rad2asec / A;
        } else if( A <= 0 ){
            A = P * rad2asec / F;
        } else if( P <= 0 ){
            P = F * A / rad2asec;
        }
        return true;
    }
    return false;
}


bool Object::checkCfg( void ){

    if( ( saveMask & SF_SAVE_PSF )&&( saveMask & SF_SAVE_PSF_AVG) ){
        LOG_WARN << "Both GET_PSF and GET_PSF_AVG mode specified." << ende;
    }
    if( channels.empty() ){
        LOG_FATAL << "Each object must have at least 1 channel specified." << ende;
        return false;
    }
    const float minWavelength = 1E-7;
    const float maxWavelength = 3E-6;
    if( wavelength < minWavelength || wavelength > maxWavelength ) {
        LOG_FATAL << "The wavelength value (=" << wavelength << ") of Object #" << ID << " is out of bounds. [" << minWavelength << "," << maxWavelength << "]" << ende;
        return false;
    }
    if( !checkImageScale( telescopeF, arcSecsPerPixel, pixelSize ) ){
        return false;
    }

    for( auto& ch : channels ){
        if( !ch->checkCfg() ) return false;
    }

    if( outputFileName.empty() ){   // TODO: clean this up
        string tpl = channels[0]->imageTemplate;
        size_t p = tpl.find_first_of( '%' );
        if( p != string::npos ){
            string tmpString = boost::str( boost::format( tpl )% 1 );
            p = 0;
            size_t i = std::min( tmpString.length(), tpl.length() );
            while( p < i && tmpString[p] == tpl[p] )p++;
            size_t ii = tmpString.length( )- 1;
            i = tpl.length( )- 1;
            while( ii && i && tmpString[ii] == tpl[i] ){
                ii--;
                i--;
            }
            tmpString.replace( p, ii - p + 1, "%d..%d" );
            if( count( tmpString.begin(), tmpString.end(), '%' )== 2 ){
                outputFileName = boost::str( boost::format( tmpString )% *channels[0]->fileNumbers.begin( )% *channels[0]->fileNumbers.rbegin() );
            } else {
                LOG_FATAL << boost::format( "failed to generate output filename from \"%s\" ( ->\"%s\")." )% tpl % tmpString << ende;
                return false;
            }
        } else LOG_FATAL << boost::format( "first filename template \"%s\" does not contain valid format specifier." )% tpl << ende;
    }
    
    /*if( !modeFile.empty( )&& !modeNumbers.empty( ) ){
        LOG_WARN << "A modefile was specified together with mode-numbers. BE AWARE that these numbers will be interpreted as indices in the supplied file and NOT Zernike/KL indices!!" << ende;
        return false;
    }*/
    

    return true;
}


bool Object::checkData( bool verbose ) {

    bfs::path outDir( myJob.info.outputDir );
    bfs::path tmpOF(outputFileName+".ext" );
    
    if( isRelative(tmpOF) && !outDir.empty() ){
        tmpOF = outDir / tmpOF;
    }
    
    for( int i = 1; i & FT_MASK; i <<= 1 ){
        if( i & myJob.outputFileType ) {  // this filetype is specified.
            tmpOF.replace_extension(FileTypeExtensions.at( ( FileType )i) );
            if( bfs::exists(tmpOF ) && !(myJob.runFlags & RF_FORCE_WRITE ) ){
                LOG_FATAL << boost::format( "output file %s already exists! Use -f (or cfg-keyword OVERWRITE) to replace file." ) % tmpOF << ende;
                return false;
            } else {
                if( verbose ) LOG_DETAIL << "Object #" << ID << ", filename: " << tmpOF << ende;
            }
        }
    }
    
    set<uint32_t> tmpWf;
    for( shared_ptr<Channel>& ch: channels ) {
        if( !ch->checkData() ) return false;
        tmpWf.insert( ch->waveFrontList.begin(), ch->waveFrontList.end() );
    }
    waveFrontList.assign( tmpWf.begin(), tmpWf.end() );
    string wfStr = redux::util::uIntsToString( waveFrontList );
    LOG << "Object " << ID << " waveFronts: " << wfStr << ende;
    
    bfs::path tmpPath = tmpOF.parent_path( );
    
    try {
        struct Writable {}; 
        auto isWritable = Cache::get().getMap<string,Writable>();
        auto it = isWritable.second.find( tmpPath.string() );
        if( it == isWritable.second.end() ) {
            if( !tmpPath.empty() && !bfs::exists(tmpPath) ) {
                if( !bfs::create_directories(tmpPath) ) {
                    LOG_FATAL << boost::format( "failed to create directory for output: %s" ) % tmpPath << ende;
                    return false;
                } else LOG_TRACE << boost::format( "create output directory %s" ) % tmpPath << ende;
            }
            bfs::path slask = tmpPath / bfs::path("_test_writability_");
            bfs::create_directory( slask );
            bfs::remove_all( slask );
            isWritable.second.emplace( tmpPath.string(), Writable() );
        }
    } catch( bfs::filesystem_error& e ){
        LOG_FATAL << boost::format( "output directory %s not writable: %s" ) % tmpPath % e.what( )<< ende;
        return false;
    }

    if( !pupilFile.empty() ){
        if( !bfs::is_regular_file(pupilFile) ){
            bfs::path fn = bfs::path(imageDataDir )/ bfs::path(pupilFile );
            if( ! bfs::is_regular_file(fn) ){
                //logAndThrow("Pupil-file " + pupilFile + " not found!" );
                LOG_FATAL << boost::format( "Pupil-file %s not found!" )% pupilFile << ende;
                return false;
            } else pupilFile = fn.string( );
        }
    }

    if( !modeFile.empty() ){
        if( !bfs::is_regular_file(modeFile) ){
            bfs::path fn = bfs::path(imageDataDir )/ bfs::path(modeFile );
            if( !bfs::is_regular_file(fn) ){
                //logAndThrow("Mode-file " + modeFile + " not found!" );
                LOG_FATAL << boost::format( "Mode-file %s not found!" )% modeFile << ende;
                return false;
            } else modeFile = fn.string( );
        }
    }


    return true;
}


void Object::initObject( void ){
    
    Pupil::calculatePupilSize( frequencyCutoff, pupilRadiusInPixels, pupilPixels, wavelength, patchSize, myJob.telescopeD, arcSecsPerPixel );
    
    if( !modes ) modes.reset( new ModeSet() );
    if( !pupil ) pupil.reset( new Pupil() );

    myJob.patchSize = patchSize;     // TODO: fulhack until per-channel sizes is implemented
    myJob.pupilPixels = pupilPixels;
    if( bfs::is_regular_file(pupilFile) ){
        PupilInfo info( pupilFile, pupilRadiusInPixels, pupilPixels );
        shared_ptr<Pupil> ret = myJob.globalData->get(info);
        lock_guard<mutex> lock(ret->mtx );
        if( ret->empty() ){    // this set was inserted, so it is not loaded yet.
            if( ret->load( pupilFile, pupilPixels, pupilRadiusInPixels ) ){
                LOG << "Loaded Pupil-file (" << info << ")" << ende;
                ret->info = info;
                pupil = ret;
            } else LOG_ERR << "Failed to load Pupil-file " << info << ende;
        } else {
            if( ret->nPixels && ret->nPixels == pupilPixels ){    // matching pupil
                //LOG_DEBUG << "Using pre-calculated pupil: ( " << pupilPixels << "x" << pupilPixels << "  radius=" << pupilRadiusInPixels << ")" << ende;
                pupil = ret;
            } else {
                LOG_ERR << "The Cache returned a non-matching Pupil. This might happen if a loaded Pupil was rescaled( which is not implemented yet)." << ende;
            }
        }
    }
    
    if( pupil->empty() ) {
        double co_radius = 0.0;
        if( myJob.telescopeCO > 0.0 ) {
            co_radius = pupilRadiusInPixels * myJob.telescopeCO/myJob.telescopeD;
        }
        PupilInfo info( pupilPixels, pupilRadiusInPixels, co_radius );
        shared_ptr<Pupil> ret = myJob.globalData->get(info );
        lock_guard<mutex> lock(ret->mtx );
        if( ret->empty() ){    // this set was inserted, so it is not generated yet.
            ret->generate( pupilPixels, pupilRadiusInPixels, co_radius );
            if( ret->nDimensions() != 2 || ret->dimSize(0) != pupilPixels || ret->dimSize(1) != pupilPixels ){    // mismatch
                LOG_ERR << "Generated Pupil does not match. This should NOT happen!!" << ende;
            } else {
                LOG_DETAIL << "Generated pupil: " << info << ende;
                pupil = ret; 
            }
        } else {
            if( ret->nPixels && ret->nPixels == pupilPixels ){    // matching pupil
                //LOG_DEBUG << "Using pre-calculated Pupil: ( " << pupilPixels << "x" << pupilPixels << "  radius=" << pupilRadiusInPixels << ")" << ende;
                pupil = ret;
            } else {
                LOG_ERR << "The Cache returned a non-matching Pupil (" << info << "). PupilPixels = " << pupilPixels
                        << "\n\t\tThis should NOT happen!!" << ende;
            }
        }
    }

    ModeInfo info( modeFile, pupilPixels );
    if( modeFile.empty() ) {
        info = ModeInfo( myJob.klMinMode, myJob.klMaxMode, myJob.modeList, pupilPixels, pupilRadiusInPixels, rotationAngle, myJob.klCutoff );
        if( myJob.modeBasis == ZERNIKE ){
            info.firstMode = info.lastMode = 0;
        }
    }
    
    modes = myJob.globalData->get( info );          // get normalized modes

    bool doPrint(false);    // just to reduce spamming the log, print only when modes are generated (i.e. on the master)
    lock_guard<mutex> lock( modes->mtx );
    if( modes->empty() ) {
        if( bfs::is_regular_file(modeFile) ){
            LOG_DEBUG << "initObject(" << to_string(ID) << "):   Empty ModeSet: loading \"" << modeFile << "\"" << ende;
            modes->load( modeFile, pupilPixels );
        } else if( info.firstMode == info.lastMode ) {
            LOG_DEBUG << "initObject(" << to_string(ID) << "):   Empty ModeSet: generating Zernike set." << ende;
            modes->generate( pupilPixels, pupilRadiusInPixels, rotationAngle, myJob.modeList, Zernike::NORMALIZE );
        } else {
            LOG_DEBUG << "initObject(" << to_string(ID) << "):   Empty ModeSet: generating KL set." << ende;
            modes->generate( pupilPixels, pupilRadiusInPixels, rotationAngle, myJob.klMinMode, myJob.klMaxMode, myJob.modeList, myJob.klCutoff, Zernike::NORMALIZE );
        }
        if( pupil ) modes->getNorms( *pupil );
        else modes->getNorms();
        LOG_DEBUG << "initObject(" << to_string(ID) << "): " << printArray( modes->norms,"measured norms" ) << ende;
        if( !bfs::is_regular_file(modeFile) || modeFileNormalize ) modes->normalize();
        doPrint = true;
    }
    
    modes->measureJacobian( *pupil, (0.5*frequencyCutoff)/util::pix2cf( arcSecsPerPixel, myJob.telescopeD ) );
    
    if( doPrint ) {
        LOG_DETAIL << "Alpha-to-Shift: " << modes->alphaToShift( PointF(1,1) )
                   << "    Shift-to-Alpha: " << modes->shiftToAlpha( PointF(1,1) ) << ende;
    }
    
    for( shared_ptr<Channel>& ch: channels ){
        ch->initChannel();
    }
    
}


void Object::reInitialize( boost::asio::io_service& service, ProgressWatch& pw, bool doReset ) {
    
    progWatch.clear();
    if( imgShifted ) {
        fitAvgPlane();
        if( doReset ) initPatch();
        for( const shared_ptr<Channel>& c: channels ) {
            for( const shared_ptr<SubImage>& im: c->getSubImages() ) {
                service.post( [this,&pw,im,doReset](){
                    //im->initialize(true);   // always re-calculate noise statistics
                    im->initialize( *this, doReset );
                    ++pw;
                } );
            }
        }
    } else {
        pw.increase( nImages() );
    }

}


void Object::loadData( boost::asio::io_service& service, Array<PatchData::Ptr>& patches ){
    
    startT = bpx::pos_infin;
    endT = bpx::neg_infin;
    
    progWatch.set( nImages( ) );

    loadInit( service, patches );
    progWatch.setTicker(nullptr );
//     progWatch.setTicker([&](){
//         LOG_WARN << "Object" << to_string(ID )<< ")::loadData( ) otick: " << progWatch.dump( )<< ende;
//     } );

    progWatch.setHandler([this](){        // this will be triggered after all images in this object are loaded/pre-processed
        vector<double> channelNorms;
        normalizeTo = std::numeric_limits<double>::lowest( );
        for( auto& ch : channels ){
            if( myJob.normType == NORM_OBJ_MAX_MEAN ) channelNorms.push_back( ch->getMaxMean() );
            else if( myJob.normType == NORM_OBJ_MAX_MEDIAN ) channelNorms.push_back( ch->getMaxMedian() );
            else if( myJob.normType == NORM_OBJ_MEDIAN_MEDIAN ) channelNorms.push_back( ch->getMedianMedian() );
            if(startT.is_special() )startT = ch->startT;
            else startT = std::min(startT,ch->startT );
            if(endT.is_special() )endT = ch->endT;
            else endT = std::max(endT,ch->endT );
        }
        std::sort( channelNorms.begin(), channelNorms.end() );
        if( myJob.normType == NORM_OBJ_MAX_MEAN || myJob.normType == NORM_OBJ_MAX_MEDIAN ) {
            normalizeTo = *(channelNorms.rbegin());
        } else if( myJob.normType == NORM_OBJ_MEDIAN_MEDIAN ) {
            size_t nM = channelNorms.size();
            size_t nHalf = (nM>>1) - !(nM&1);
            auto mid = channelNorms.begin() + nHalf;
            std::nth_element( channelNorms.begin(), mid, channelNorms.end());
            normalizeTo = *mid;
            if( !(nM&1) ) {
                std::nth_element( mid, mid+1, channelNorms.end());
                normalizeTo += *(mid+1);
                normalizeTo /= 2;
            }
        }
        LOG_DEBUG << "Object " << ID << ": The images will be normalized to " << normalizeTo << ende;
    } );
    
    for( shared_ptr<Channel>& ch: channels ){
        ch->loadData( service, patches );
    }
 
}


void Object::loadInit( boost::asio::io_service& service, Array<PatchData::Ptr>& patches ){
    
    //if( initFile == "" )return;
    
    bfs::path tmpFile( initFile );
    if( isRelative(tmpFile ) ){
        tmpFile = bfs::path(myJob.info.outputDir )/ tmpFile;
    }

    if( !bfs::exists(tmpFile ) ){
        tmpFile = bfs::path(tmpFile.string( )+ ".momfbd" );
    }

    if( !bfs::exists(tmpFile ) )return;
    
    ifstream file;
    std::shared_ptr<FileMomfbd> info( new FileMomfbd( ) );

    try {
        file.open( tmpFile.string( ) );
        info->read( file );
        if( !(info->dataMask&MOMFBD_ALPHA ) ){
            LOG_ERR << "Initilazation file: " << tmpFile <<  " does not contain alphas." << ende;
            return;
        }
        
        size_t nPatchesX = info->nPatchesX;
        size_t nPatchesY = info->nPatchesY;
        size_t patchSize = info->getPatchSize( MOMFBD_ALPHA );
        size_t offset = info->getPatchSize( 0 );

        if(( patches.dimSize(0 )!= nPatchesY )||( patches.dimSize(1 )!= nPatchesX ) ){
            LOG_ERR << "Initialization file: " << tmpFile <<  " has the wrong number of patches." << ende;
            LOG_ERR << "nPatchesX = " << nPatchesX << " - " << patches.dimSize(1)<< ende;
            LOG_ERR << "nPatchesY = " << nPatchesY << " - " << patches.dimSize(0)<< ende;
            return;
        }
        size_t nModes = info->patches(0,0).nm;
        if( nModes != myJob.nModes ){
            LOG_ERR << "Initilazation file: " << tmpFile <<  " has the wrong number of modes." << ende;
            LOG_ERR << "init->nModes = " << nModes << " != " << myJob.nModes << ende;
            return;
        }
        size_t nImgs = ( patchSize-offset)/(nModes*sizeof(float) );
        if( nImgs != nImages( ) ){
            LOG_ERR << "Initilazation file: " << tmpFile <<  " has the wrong number of alphas." << ende;
            LOG_ERR << "nImgs = " << nImgs << " != " << nImages( )<< ende;
            return;
        }

        LOG << "Loading initialization from file: " << tmpFile << ende;
        shared_ptr<char> tmpData = rdx_get_shared<char>(patchSize);
        for( size_t y=0; y<nPatchesY; ++y ){
            for( size_t x=0; x<nPatchesX; ++x ){
                info->patches(y,x).load( file, tmpData.get(), info->swapNeeded, info->version, MOMFBD_ALPHA );
                auto oData = patches(y,x)->getObjectData(ID);
                if( !oData ) throw runtime_error("patches(y,x)->getObjectData(ID) returned a null pointer !");
                oData->alpha.resize( nImgs, nModes );
                oData->alpha.copyFrom<float>( reinterpret_cast<float*>(tmpData.get()+offset ) );
            }
        }
        
    } catch( exception& e ){
        LOG_ERR << "Failed to read initilazation file: " << tmpFile <<  "\n\t Reason: " << e.what( )<< ende;
        return;
    }
/*
    if( verbosity > 1 ){
        cout << "File Version:       \"" << info->versionString << "\"" << endl;
    }

    uint8_t loadMask = 0;

    if( kw.img )    loadMask |= MOMFBD_IMG;
    if( kw.psf )    loadMask |= MOMFBD_PSF;
    if( kw.obj )    loadMask |= MOMFBD_OBJ;
    if( kw.res )    loadMask |= MOMFBD_RES;
    if( kw.alpha )  loadMask |= MOMFBD_ALPHA;
    if( kw.div )    loadMask |= MOMFBD_DIV;
    if( kw.modes )  loadMask |= MOMFBD_MODES;
    if( kw.names )  loadMask |= MOMFBD_NAMES;
    if( kw.all )    loadMask = info->dataMask;

    IDL_KW_FREE;

    if( !loadMask && !checkData ){
        loadMask = info->dataMask;
    } else {
        loadMask &= info->dataMask;
    }
*/
    
    
    //LOG << "Loading initialization from file: " << tmpFile << ende;
    //++progWatch;
}


size_t Object::getResultSize( void ) {
    size_t ret = patchSize*patchSize;
    if( saveMask & (SF_SAVE_PSF|SF_SAVE_PSF_AVG) ) {
        uint16_t nPSF = ( saveMask&SF_SAVE_PSF_AVG)? 1 : nObjectImages;
        ret += nPSF*patchSize*patchSize;
    }
    if( saveMask & SF_SAVE_COBJ ){
        ret += nObjectImages*patchSize*patchSize;
    }
    if( saveMask & SF_SAVE_RESIDUAL ){
        ret += nObjectImages*patchSize*patchSize;
    }
    if( saveMask & SF_SAVE_ALPHA ){
        ret += nObjectImages*myJob.nModes;
    }
    vector<shared_ptr<Channel>>& objChannels = channels;
    if( objChannels.empty() ) {
        objChannels = myJob.getChannels(ID);
    }
    int nCh = objChannels.size();
    if( nCh && (saveMask & SF_SAVE_DIVERSITY) ){
        ret += nCh*pupilPixels*pupilPixels;
    }
    return ret;
}


void Object::maybeInitializeStorage( void ) {
    
    lock_guard<mutex> lock(mtx);
    if( results.nElements() ) return;        // already allocated/loaded
    
    size_t nPatchesX = myJob.subImagePosX.size();
    size_t nPatchesY = myJob.subImagePosY.size();
#ifdef RDX_DO_TRANSPOSE
    std::swap( nPatchesX, nPatchesY );
#endif
    
    size_t resultSize = getResultSize();

    if( !(myJob.runFlags&RF_NOSWAP) ) {     // should we mmap?
        if( bfs::exists( bfs::path(cacheFile) ) ) {
            LOG_DEBUG << "opening: Object: " << cacheFile << ende;
            results.openMmap( cacheFile, nPatchesY, nPatchesX, resultSize );
        } else {
            LOG_DEBUG << "creating: Object: " << cacheFile << ende;
            results.createMmap( cacheFile, nPatchesY, nPatchesX, resultSize );
        }
    } else {
        results.resize(  nPatchesY, nPatchesX, resultSize );
    }
    
}


void Object::getStorage( PatchData& pData, shared_ptr<ObjectData> oData ) {

    maybeInitializeStorage();
    
    if( !oData ) return;
    
    float* resPtr = results.ptr( pData.index.y, pData.index.x, 0 );
    
    oData->img.wrap( resPtr, patchSize, patchSize );
    resPtr += patchSize*patchSize;
    if( saveMask & (SF_SAVE_PSF|SF_SAVE_PSF_AVG) ) {
        uint16_t nPSF = ( saveMask&SF_SAVE_PSF_AVG)? 1 : nObjectImages;
        oData->psf.wrap( resPtr, nPSF, patchSize, patchSize );
        resPtr += nPSF*patchSize*patchSize;
    }
    if( saveMask & SF_SAVE_COBJ ){
        oData->cobj.wrap( resPtr, nObjectImages, patchSize, patchSize );
        resPtr += nObjectImages*patchSize*patchSize;
    }
    if( saveMask & SF_SAVE_RESIDUAL ){
        oData->res.wrap( resPtr, nObjectImages, patchSize, patchSize );
        resPtr += nObjectImages*patchSize*patchSize;
    }
    if( saveMask & SF_SAVE_ALPHA ){
        oData->alpha.wrap( resPtr, nObjectImages, myJob.nModes );
        resPtr += nObjectImages*myJob.nModes;
    }
    vector<shared_ptr<Channel>>& objChannels = channels;
    if( objChannels.empty() ) {
        objChannels = myJob.getChannels(ID);
    }
    int nCh = objChannels.size();
    if( nCh && (saveMask & SF_SAVE_DIVERSITY) ){
        oData->div.wrap( resPtr, nCh, pupilPixels, pupilPixels );
        //resPtr += nCh*pupilPixels*pupilPixels; // resPtr is not used again
    }

}


void Object::doMozaic( float**& img, size_t& imgRows, size_t& imgCols, const redux::util::Array<PatchData::Ptr>& patches ) {

    size_t nPatches = patches.nElements( );
    vector<shared_ptr<float*>> patchPtrs;
    vector<float**> patchData;
    vector<int32_t> xpos,ypos;
    uint16_t maxPosX(0);
    uint16_t minPosX( std::numeric_limits<uint16_t >::max() );
    uint16_t maxPosY = maxPosX;
    uint16_t minPosY = minPosX;

    for( unsigned int y = 0; y < patches.dimSize(0 ); ++y ){
        for( unsigned int x = 0; x < patches.dimSize(1 ); ++x ){
            const Point16& first = patches(y,x)->roi.first;
            if( first.x > maxPosX ) maxPosX = first.x;
            if( first.x < minPosX ) minPosX = first.x;
            if( first.y > maxPosY ) maxPosY = first.y;
            if( first.y < minPosY ) minPosY = first.y;
            auto oData = patches(y,x)->getObjectData(ID);
            if( !oData ) throw runtime_error("patches(y,x)->getObject(ID) returned a null pointer !");
            auto pPtr = oData->img.reshape(patchSize,patchSize );
            patchPtrs.push_back( pPtr );
            patchData.push_back( pPtr.get( ) );
            xpos.push_back( first.x );
            ypos.push_back( first.y );
        }
    }

    imgCols = maxPosX+patchSize+1;
    imgRows = maxPosY+patchSize+1;
    
    delArray( img );
    img = newArray<float>( imgRows, imgCols );
    int margin = patchSize/8;
    int blend = (patchSize-2*margin)/3;

    
    bool doTranspose(false);
#ifdef RDX_DO_TRANSPOSE
    doTranspose = true;
#endif
    
    mozaic( img, imgRows, imgCols, const_cast<const float***>(patchData.data()), nPatches, patchSize, patchSize, ypos.data(), xpos.data(), blend, margin, doTranspose );

    if( !(myJob.runFlags&RF_NO_CLIP) ) {
        img_trim( img, imgRows, imgCols, 1E-15 );
    }


}


void Object::writeAna( const redux::util::Array<PatchData::Ptr>& patches ) {

    bfs::path fn = bfs::path(outputFileName + ".f0" );
    if( isRelative(fn ) ){
        fn = bfs::path(myJob.info.outputDir )/ fn;
    }
    LOG << "Writing output to file: " << fn << ende;
    
    try {

        size_t imgCols;
        size_t imgRows;
        float** tmpImg = nullptr;
        doMozaic( tmpImg, imgRows, imgCols, patches );
        
        shared_ptr<redux::file::Ana> hdr( new redux::file::Ana() );
        boost::posix_time::ptime now = boost::posix_time::microsec_clock::universal_time();

        string timeString = "N/A";
        if( !startT.is_special( ) && !endT.is_special() ){
            bpx::time_duration obs_interval = ( endT - startT );
            boost::posix_time::ptime avgtime = (startT+obs_interval/2);
            if( !avgtime.is_special() ) {
                timeString = bpx::to_simple_string( avgtime.time_of_day() );
            }
        } else if( !endT.is_special() ){
            timeString = bpx::to_simple_string( endT.time_of_day() );
        } else if( !startT.is_special() ){

            timeString = bpx::to_simple_string( startT.time_of_day() );
        }
        hdr->m_ExtendedHeader = "TIME_OBS=" + timeString + " DATE_OBS=" + myJob.observationDate
                            + " NX=" + to_string(imgCols) + " NY=" + to_string(imgRows) + " DATE=" + to_iso_extended_string( now );

        
        if( myJob.outputDataType == DT_F32T ){
            Array<float> wrap(*tmpImg, imgRows, imgCols );
            Ana::write( fn.string(), wrap, hdr );
        } else {
            Array<int16_t> wrap(imgRows, imgCols );
            wrap.copyFrom<float>( *tmpImg );
            Ana::write( fn.string(), wrap, hdr, 5 );
        }

        delArray( tmpImg );
        
        if( saveMask & SF_SAVE_ALPHA ){
            bfs::path fn = bfs::path(outputFileName + ".alpha.f0" );
            if( isRelative(fn ) ){
                fn = bfs::path(myJob.info.outputDir )/ fn;
            }
            LOG << "Saving alpha-coefficients to: " << fn << ende;
            Array<float> alpha(patches.dimSize(0), patches.dimSize(1), nObjectImages, myJob.nModes );
            for( auto& patch: patches ){
                Array<float> subalpha(alpha, patch->index.y, patch->index.y, patch->index.x, patch->index.x, 0, nObjectImages-1, 0, myJob.nModes-1 );
                auto oData = patch->getObjectData(ID);
                if( !oData ) throw runtime_error("patches(y,x)->getObject() returned a null pointer !");
                oData->alpha.copy( subalpha );
            }
            Ana::write(fn.string(), alpha );
        }
        
    } catch( const exception& e ) {
        LOG_ERR << "Error when writing output to file: " << fn << " what: " << e.what() << ende;
        myJob.setFailed();
    }

    ++progWatch;
    
}


void Object::writeFits( const redux::util::Array<PatchData::Ptr>& patches ) {
    
    bfs::path fn = bfs::path( outputFileName + ".fits" );
    if( isRelative( fn ) ) {
        fn = bfs::path( myJob.info.outputDir ) / fn;
    }
    LOG << "Writing output to file: " << fn << ende;
    
    try {

        size_t imgCols;
        size_t imgRows;
        float** tmpImg = nullptr;
        doMozaic( tmpImg, imgRows, imgCols, patches );
  
        if( myJob.outputDataType == DT_F32T ){
            Array<float> wrap(*tmpImg, imgRows, imgCols );
            Fits::write( fn.string(), wrap );
        } else {
            Array<int16_t> wrap(imgRows, imgCols );
            wrap.copyFrom<float>( *tmpImg );
            Fits::write( fn.string(), wrap );
        }

        delArray( tmpImg );
        
        if( saveMask & SF_SAVE_ALPHA ){
            bfs::path fn = bfs::path(outputFileName + ".alpha.fits" );
            if( isRelative(fn ) ){
                fn = bfs::path(myJob.info.outputDir ) / fn;
            }
            LOG << "Saving alpha-coefficients to: " << fn << ende;
            Array<float> alpha(patches.dimSize(0), patches.dimSize(1), nObjectImages, myJob.nModes );
            for( auto& patch: patches ){
                Array<float> subalpha(alpha, patch->index.y, patch->index.y, patch->index.x, patch->index.x, 0, nObjectImages-1, 0, myJob.nModes-1 );
                auto oData = patch->getObjectData(ID);
                if( !oData ) throw runtime_error("patches(y,x)->getObject() returned a null pointer !");
                oData->alpha.copy( subalpha );
            }
            Fits::write( fn.string(), alpha );
        }
        
    } catch( const exception& e ) {
        LOG_ERR << "Error when writing output to file: " << fn << " what: " << e.what() << ende;
        myJob.setFailed();
    }

    ++progWatch;
}



void Object::writeMomfbd( const redux::util::Array<PatchData::Ptr>& patchesData ) {

    bfs::path fn = bfs::path(outputFileName + ".momfbd" );      // TODO: fix storage properly

    if( isRelative(fn ) ){
        fn = bfs::path(myJob.info.outputDir )/ fn;
    }
    LOG << "Writing output to file: " << fn << ende;
    
    uint16_t oID = ID;
    if( traceID >= 0 ) {    // for trace-objects, update settings from the reference object (was modified in pre-processing)
        oID = traceID;
        shared_ptr<Object> myCfg = myJob.getObject(oID);
        if( myCfg ) {
            channels = myCfg->getChannels();
            modes = myCfg->modes;
            pupil = myCfg->pupil;
        }
    }
    
    try {
        
        std::shared_ptr<FileMomfbd> info( new FileMomfbd() );

        // Extract date/time from the git commit.
        int day, month, year, hour;
        char buffer [15];
        sscanf (reduxCommitTime, "%4d-%2d-%2d %2d", &year, &month, &day, &hour);
        snprintf( buffer, 15, "%4d%02d%02d.%01d", year, month, day, hour );
        info->versionString = buffer;
        info->version = atof( info->versionString.c_str( ) );

        info->dateString = myJob.observationDate;
        if( startT.is_special() && endT.is_special() ){
            info->timeString = "N/A";
        } else if( startT.is_special() ){
            info->timeString = bpx::to_simple_string( endT.time_of_day() );
        } else if( endT.is_special() ){
            info->timeString = bpx::to_simple_string( startT.time_of_day() );
        } else {
            bpx::time_duration obs_interval =( endT - startT );
            info->timeString = bpx::to_simple_string((startT+obs_interval/2).time_of_day() );
        }

        int32_t nChannels = info->nChannels = channels.size( );
        info->clipStartX = sharedArray<int16_t>( nChannels );
        info->clipEndX = sharedArray<int16_t>( nChannels );
        info->clipStartY = sharedArray<int16_t>( nChannels );
        info->clipEndY = sharedArray<int16_t>( nChannels );
        info->fileNames.clear( );

        for( int i = 0; i < nChannels; ++i ){
            channels[i]->getFileNames( info->fileNames, waveFrontList );
            if( channels[i]->alignClip.empty() ){
                Point16 sz = channels[i]->getImageSize( );
                info->clipStartX.get()[i] = info->clipStartY.get()[i] = 1;
                info->clipEndX.get()[i] = sz.x;
                info->clipEndY.get()[i] = sz.y;
            } else {
                info->clipStartX.get()[i] = channels[i]->alignClip[0]+1;
                info->clipEndX.get()[i] = channels[i]->alignClip[1]+1;
                info->clipStartY.get()[i] = channels[i]->alignClip[2]+1;
                info->clipEndY.get()[i] = channels[i]->alignClip[3]+1;
            }
        }

        info->nPH = pupilPixels;
        uint8_t writeMask = MOMFBD_IMG;                                                 // always output image
        int64_t imgSize = patchSize*patchSize*sizeof(float );
   
        if( info->fileNames.size( ) ) writeMask |= MOMFBD_NAMES;
        if( saveMask & (SF_SAVE_PSF|SF_SAVE_PSF_AVG) ) writeMask |= MOMFBD_PSF;
        if( saveMask & SF_SAVE_MODES && (info->nPH > 0) ) writeMask |= MOMFBD_MODES;
        if( saveMask & SF_SAVE_COBJ ) writeMask |= MOMFBD_OBJ;
        if( saveMask & SF_SAVE_RESIDUAL ) writeMask |= MOMFBD_RES;
        if( saveMask & SF_SAVE_ALPHA ) writeMask |= MOMFBD_ALPHA;
        if( nChannels && (saveMask & SF_SAVE_DIVERSITY) ) writeMask |= MOMFBD_DIV;
       
        Array<float> tmpModes;
        if( writeMask&MOMFBD_MODES ){     // copy modes from local cache
            //double pupilRadiusInPixels = pupilPixels / 2.0;
            //if( objChannels.size() )pupilRadiusInPixels = objChannels[0]->pupilRadiusInPixels;
            if( !modes || !pupil ) {
                writeMask &= ~MOMFBD_MODES;
            } else {
                tmpModes.resize( myJob.nModes+1, info->nPH, info->nPH );            // +1 to also fit pupil in the array
                tmpModes.zero( );
                Array<double> mode_wrap(reinterpret_cast<Array<double>&>(*modes), 0, myJob.nModes-1, 0, info->nPH-1, 0, info->nPH-1 );
                Array<float> tmp_slice(tmpModes, 0, 0, 0, info->nPH - 1, 0, info->nPH - 1 );     // subarray
                tmp_slice.assign(reinterpret_cast<const Array<double>&>(*pupil) );
                info->pix2cf = modes->shiftToAlpha( PointD(1,1) ).avg();
                info->cf2pix = modes->alphaToShift( PointD(1,1) ).avg();
                info->phOffset = 0;
                tmp_slice.wrap(tmpModes, 1, myJob.nModes, 0, info->nPH - 1, 0, info->nPH - 1 );
                tmp_slice.assign(mode_wrap);
                tmp_slice *= 1.0/wavelength;
                if( myJob.nModes ){
                    //ModeInfo id( myJob.klMinMode, myJob.klMaxMode, 0, pupilPixels, pupilRadiusInPixels, rotationAngle, myJob.klCutoff );
                    info->nModes = myJob.nModes;
                    info->modesOffset = pupilPixels * pupilPixels * sizeof( float );
                }
            }
        }

        Point16 nP;                         // number of patches in output file.
        nP.y = patchesData.dimSize(0);      // nPatchesY is the slowest dimension in the .momfbd file.
        nP.x = patchesData.dimSize(1);
        
        info->nPatchesY = nP.y;             // nPatchesY is the slowest dimension in the .momfbd file.
        info->nPatchesX = nP.x;
#ifdef RDX_DO_TRANSPOSE
        std::swap( info->nPatchesY, info->nPatchesX );      // old MvN code saved the output transposed
#endif
        info->patches.resize( info->nPatchesY, info->nPatchesX );
        info->nPoints = patchSize;
        
        info->region[0] = info->region[2] = numeric_limits<int32_t>::max();
        info->region[1] = info->region[3] = numeric_limits<int32_t>::min();

        size_t modeSize = tmpModes.nElements()*sizeof( float );
        size_t blockSize = modeSize;
        shared_ptr<ObjectData> oData( make_shared<ObjectData>() );
        for( uint16_t y(0); y < nP.y; ++y ) {                    // Loops are in the standard order
            for( uint16_t x(0); x < nP.x; ++x ) {
                FileMomfbd::PatchInfo& pi = info->patches
#ifdef RDX_DO_TRANSPOSE
                ( x, y );
#else
                ( y, x );
#endif
                PatchData::Ptr thisPatchData = patchesData( y, x );
                if( thisPatchData ){
                    getStorage( *thisPatchData, oData );
                    shared_ptr<ObjectData> refData = thisPatchData->getObjectData( oID );
                    if( refData ) {
                        oData->channels = refData->channels;
                    }
                    pi.region[0] = thisPatchData->roi.first.x+1;         // store as 1-based indices
                    pi.region[1] = thisPatchData->roi.last.x+1;
                    pi.region[2] = thisPatchData->roi.first.y+1;
                    pi.region[3] = thisPatchData->roi.last.y+1;
                    pi.nChannels = nChannels;
      
                    info->region[0] = std::min( info->region[0], pi.region[0] );
                    info->region[1] = std::max( info->region[1], pi.region[1] );
                    info->region[2] = std::min( info->region[2], pi.region[2] );
                    info->region[3] = std::max( info->region[3], pi.region[3] );

                    pi.nim = sharedArray<int32_t>( nChannels );
                    pi.dx = sharedArray<int32_t>( nChannels );
                    pi.dy = sharedArray<int32_t>( nChannels );
                    for( int i=0; i < nChannels; ++i ){
                        pi.nim.get()[i] = channels[i]->nImages( waveFrontList );
                        pi.dx.get()[i] = oData->channels[i]->channelOffset.x;
                        pi.dy.get()[i] = oData->channels[i]->channelOffset.y;
                    }
                    blockSize += imgSize;
                    if( writeMask&MOMFBD_PSF ){
                        if(oData->psf.nDimensions()>1 ){
                            pi.npsf = oData->psf.dimSize(0);
                            blockSize += pi.npsf*imgSize;
                        }
                    }
                    if( writeMask&MOMFBD_OBJ ){
                        if(oData->cobj.nDimensions()>1 ){
                            pi.nobj = oData->cobj.dimSize(0);
                            blockSize += pi.nobj*imgSize;
                        }
                    }
                    if( writeMask&MOMFBD_RES ){
                        if(oData->res.nDimensions()>1 ){
                            pi.nres = oData->res.dimSize(0);
                            blockSize += pi.nres*imgSize;
                        }
                    }
                    if( writeMask&MOMFBD_ALPHA ){
                        if(oData->alpha.nDimensions()==2 ){
                            pi.nalpha = oData->alpha.dimSize(0);
                            pi.nm = oData->alpha.dimSize(1 );
                            blockSize += pi.nalpha*pi.nm*sizeof(float );
                        }
                    }
                    if( writeMask&MOMFBD_DIV ){
                        if(oData->div.nDimensions()>1 ){
                            pi.ndiv = oData->div.dimSize(0);
                            pi.nphx = info->nPH;
                            pi.nphy = info->nPH;
                            blockSize += pi.ndiv*pi.nphx*pi.nphy*sizeof(float );
                        }
                    }
                }
            }   // x-loop
        }   // y-loop

        auto tmp = sharedArray<char>( blockSize );
        memcpy(tmp.get(), tmpModes.get(), modeSize );
        char* tmpPtr = tmp.get( );
        int64_t offset = modeSize;
        results.copyTo<float>(tmpPtr+offset);
        for( int y(0); y < info->nPatchesY; ++y ) {         // N.B. Loops are in "file order, which was transposed in the old code.
            for( int x(0); x < info->nPatchesX; ++x ) {
                FileMomfbd::PatchInfo& pi = info->patches( y, x );
                PatchData::Ptr thisPatchData =
#ifdef RDX_DO_TRANSPOSE
                patchesData( x, y );
#else
                patchesData( y, x );
#endif
                if( thisPatchData ){
                    pi.imgPos = offset;
                    offset += imgSize;
                    pi.psfPos = offset;
                    offset += pi.npsf*imgSize;
                    pi.objPos = offset;
                    offset += pi.nobj*imgSize;
                    pi.resPos = offset;
                    offset += pi.nres*imgSize;
                    pi.alphaPos = offset;
                    offset += pi.nalpha*pi.nm*sizeof(float);
                    pi.diversityPos = offset;
                    offset += pi.ndiv*pi.nphx*pi.nphy*sizeof(float);
                }
            }
        }
        
        info->write( fn.string(), reinterpret_cast<char*>( tmp.get()), writeMask );

    } catch( const exception& e ) {
        LOG_ERR << "Error when writing output to file: " << fn << " what: " << e.what() << ende;
        myJob.setFailed();
    }
    
    ++progWatch;

}

void Object::writeResults( boost::asio::io_service& service, const redux::util::Array<PatchData::Ptr>& patches ){
    
    progWatch.set( 1 );
    progWatch.setHandler([this](){
        ++myJob.progWatch;
    } );

    if( myJob.outputFileType & FT_ANA ) {
        progWatch.increaseTarget( 1 );
        service.post( std::bind( &Object::writeAna, this, std::ref( patches) ) );
    }
    if( myJob.outputFileType & FT_FITS ) {
        progWatch.increaseTarget( 1 );
        service.post( std::bind( &Object::writeFits, this, std::ref( patches) ) );
    }
    if( myJob.outputFileType & FT_MOMFBD ){
        progWatch.increaseTarget(1 );
        service.post( std::bind( &Object::writeMomfbd, this, std::ref(patches) ) );
    }
    ++progWatch;
    
}

void Object::writeResults( redux::util::Array<PatchData::Ptr>& patches ) {
    
    if( myJob.outputFileType & FT_ANA ) {
        writeAna(patches);
    }
    if( myJob.outputFileType & FT_FITS ) {
        writeFits(patches);
    }
    if( myJob.outputFileType & FT_MOMFBD ) {
        writeMomfbd(patches);
    }

    ++myJob.progWatch;
}


size_t Object::estimateOutputSizeANA( void) {
    size_t sz( 0);
    
//     uint16_t maxPosX( 0);
//     uint16_t minPosX( std::numeric_limits<uint16_t >::max() );
//     uint16_t maxPosY = maxPosX;
//     uint16_t minPosY = minPosX;
// 
//     for( unsigned int y = 0; y < patches.dimSize( 0 ); ++y ) {
//         for( unsigned int x = 0; x < patches.dimSize( 1 ); ++x ) {
//             const Point16& first = patches( y, x) ->roi.first;
//             if( first.x > maxPosX ) maxPosX = first.x;
//             if( first.x < minPosX ) minPosX = first.x;
//             if( first.y > maxPosY ) maxPosY = first.y;
//             if( first.y < minPosY ) minPosY = first.y;
//             auto pPtr = patches( y, x) ->objects[ID]->img.reshape( patchSize, patchSize );
//             patchPtrs.push_back( pPtr );
//             patchData.push_back( pPtr.get() );
//             xpos.push_back( first.x );
//             ypos.push_back( first.y );
//         }
//     }
//    
//     size_t imgCols = maxPosX+patchSize;
//     size_t imgRows = maxPosY+patchSize;
    return sz;
}

size_t Object::estimateOutputSizeFITS( void) {
    return 0;
}
size_t Object::estimateOutputSizeMOMFBD( void) {
    return 0;
}

size_t Object::estimateOutputSize( void) {
    size_t sz( 0);
    if( myJob.outputFileType & FT_ANA ) {
        sz += estimateOutputSizeANA();
    }
    if( myJob.outputFileType & FT_FITS ) {
        sz += estimateOutputSizeFITS();
    }
    if( myJob.outputFileType & FT_MOMFBD ) {
        sz += estimateOutputSizeMOMFBD();
    }
    return sz;
}



Point16 Object::getImageSize( bool force ){
    if( force || (imgSize == 0) ) {
        imgSize = 0;
        for( auto& ch : channels ){
            Point16 tmp = ch->getImageSize( force );
            if( imgSize == 0 ){
                imgSize = tmp;
            } else if( tmp != imgSize ){
                throw std::logic_error( "The images have different sizes for the different channels, please verify the ALIGN_CLIP values." );
            }
        }
    }
    return imgSize;
}


void Object::dump( std::string tag ){
    
    tag += "_"+to_string(ID );

    Ana::write( tag + "_ftsum.f0", Array<double>( ftSum.get(), otfSize, otfSize ) );
    Ana::write( tag + "_q.f0", Array<double>( Q.get(), otfSize, otfSize ) );
    Ana::write( tag + "_p.f0", Array<complex_t>(P.get(), otfSize, otfSize) );
    Ana::write( tag + "_fittedplane.f0", fittedPlane );
    if( pupil ) Ana::write( tag + "_pupil.f0", *pupil );
    if( modes ) Ana::write( tag + "_modes.f0", *modes );
    for( auto& ch : channels ){
        ch->dump(tag );
    }

}
