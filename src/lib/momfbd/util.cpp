#include "redux/momfbd/util.hpp"

#include "redux/file/fileio.hpp"
#include "redux/logger.hpp"
#include "redux/constants.hpp"

#include <boost/format.hpp>
#include <boost/filesystem.hpp>

using namespace redux::momfbd::util;
using namespace redux::momfbd;
using namespace redux::file;
using namespace redux;
using namespace std;

namespace bfs = boost::filesystem;

#define lg Logger::mlg
namespace {
    const string thisChannel = "util";
}


int64_t redux::momfbd::util::extractTime( const char* header ) {
    
    if( header ) {
        if( const char *q=strstr( header,"Ts=" ) ) { // try new style
            if( const char *r=strstr( header,"Te=" ) ) { //second entry
                int yy,mm,dd,hh,mn,ss,us;
                if( sscanf( q+3,"%04d.%02d.%02d %02d:%02d:%02d.%06d",&yy,&mm,&dd,&hh,&mn,&ss,&us )>=7 ) {
                    int yye,mme,dde,hhe,mne,sse,use;
                    if( sscanf( r+3,"%04d.%02d.%02d %02d:%02d:%02d.%06d",&yye,&mme,&dde,&hhe,&mne,&sse,&use )>=7 ) {
                       // io.msg( MFBD_XNFO,"new style header contained time stamps [%04d.%02d.%02d %02d:%02d:%02d.%06d,%04d.%02d.%02d %02d:%02d:%02d.%06d]\n",yy,mm,dd,hh,mn,ss,us,yye,mme,dde,hhe,mne,sse,use );
                        LOG_DEBUG << boost::format( "new style header contained time stamps [%04d.%02d.%02d %02d:%02d:%02d.%06d,%04d.%02d.%02d %02d:%02d:%02d.%06d]" )
                        % yy % mm % dd % hh % mn % ss % us % yye % mme % dde % hhe % mne % sse % use;
                        return ( ( int64_t )( 3600000*hh+60000*mn+1000*ss+( us/1000 ) )+( int64_t )( 3600000*hhe+60000*mne+1000*sse+( use/1000 ) ) )/( int64_t )2; // full succes
                    } else {
                        //io.msg( MFBD_WARN,"new style header contained only one usable time stamp (Ts)\n" );
                        LOG_WARN << "new style header contained only one usable time stamp (Ts)";
                        return ( int64_t )( 3600000*hh+60000*mn+1000*ss+( us/1000 ) ); // half succes
                    }
                } else {
                    int yye,mme,dde,hhe,mne,sse,use;
                    if( sscanf( r+3,"%04d.%02d.%02d %02d:%02d:%02d.%06d",&yye,&mme,&dde,&hhe,&mne,&sse,&use )>=7 ) {
                        //io.msg( MFBD_WARN,"new style header contained only one usable time stamp (Te)\n" );
                        LOG_WARN << "new style header contained only one usable time stamp (Te)";
                        return ( int64_t )( 3600000*hhe+60000*mne+1000*sse+( use/1000 ) ); // half succes
                    }
                    LOG_ERR << "new style header contained no usable time stamp!";
                    //io.msg( MFBD_ERROR,"new style header contained no usable time stamp!\n" );
                    return -1;
                }
            } else { // second entry is missing: what to do?
                LOG_WARN << "new style header contained only starting time stamp";
                //io.msg( MFBD_WARN,"new style header contained only starting time stamp\n" );
                int yy,mm,dd,hh,mn,ss,us;
                if( sscanf( q+3,"%04d.%02d.%02d %02d:%02d:%02d.%06d",&yy,&mm,&dd,&hh,&mn,&ss,&us )>=7 )
                    return ( int64_t )( 3600000*hh+60000*mn+1000*ss+us/1000 ); // succes
                LOG_ERR << "new style header contained one time stamp (Ts) which was not usable!";
                //io.msg( MFBD_ERROR,"new style header contained one time stamp (Ts) which was not usable!\n" );
                return -1;
            }
        } else {
            if( const char *r=strstr( header,"Te=" ) ) { //second entry
                //io.msg( MFBD_WARN,"new style header contained only ending time stamp\n" );
                LOG_WARN << "new style header contained only ending time stamp";
                int yy,mm,dd,hh,mn,ss,us;
                if( sscanf( r+3,"%04d.%02d.%02d %02d:%02d:%02d.%06d",&yy,&mm,&dd,&hh,&mn,&ss,&us )>=7 )
                    return ( int64_t )( 3600000*hh+60000*mn+1000*ss+us/1000 ); // succes
                LOG_ERR << "new style header contained one time stamp (Te) which was not usable!";
                //io.msg( MFBD_ERROR,"new style header contained one time stamp (Te) which was not usable!\n" );
                return -1;
            }
        } // nothing: try old style

        const char *p=strchr( header,':' );
        if( !p ) return -1;
        int h,m,s,ms;
        if( sscanf( p-2,"%d:%d:%d.%d",&h,&m,&s,&ms )!=4 )
            return -1; // read error!
        return ( int64_t )( 3600000*h+60000*m+1000*s+ms ); // succes
    } // no header

    return -1;

}

void redux::momfbd::util::loadPupil( const std::string& filename, redux::util::Array<double>& pupil, uint32_t pupilSize ) {
    
    if( bfs::exists( filename ) ) {

        readFile( filename, pupil );
       
        if( pupil.nDimensions() != 2 ) {
            LOG_ERR << "Pupil file \"" << filename << "\" not 2 dimensional.";
            pupil.resize();
            return;
        }
        
        uint32_t tmpInt1 = pupil.dimSize(0);
        uint32_t tmpInt2 = pupil.dimSize(1);
        if( tmpInt1 != tmpInt2 ) {
            LOG_ERR << "Pupil file \"" << filename << "\" not square (" << tmpInt1 << "x" << tmpInt2 << ")";
            pupil.resize();
            return;
        }
        
       if( pupilSize == 0 || tmpInt1 == pupilSize ) {
           LOG_DETAIL << "Pupil \"" << filename << "\" loaded.";
           return;
       } else if( tmpInt1 > pupilSize ) {
            LOG_WARN << "Pupil image \"" << filename << "\" needs to be rebinned " << tmpInt1 << " -> " << pupilSize << ".";
            // pupil=rebin(myJob.pupil,1,myJob.pupilSize,1,myJob.pupilSize,r_c,1,nph,1,nph);
        } else {
            LOG_WARN << "Pupil image \"" << filename << "\" needs to be interpolated " << tmpInt1 << " -> " << pupilSize << ".";
            // pupil=resample(myJob.pupil,1,myJob.pupilSize,1,myJob.pupilSize,r_c,1,nph,1,nph);
        }
        
    }

    pupil.resize();

}


double redux::momfbd::util::cf2pix(double arcsecperpix, double telescope_d) {
    return 1.0 / pix2cf(arcsecperpix,telescope_d);
}


double redux::momfbd::util::pix2cf(double arcsecperpix, double telescope_d) {
    static const double a = (redux::PI*redux::PI)/648000.0;
    return a * (arcsecperpix*telescope_d);
}


double redux::momfbd::util::def2cf(double pd_defocus,double telescope_r) { // defocus distance in meters
    static const double a = -redux::PI/(8.0*sqrt(3.0));
    return a * pd_defocus*sqr(telescope_r)/(8.0*sqrt(3.0));
}


double redux::momfbd::util::cf2def(double alpha,double telescope_r) {
    return 1.0 / def2cf(1.0/alpha,telescope_r);
}
