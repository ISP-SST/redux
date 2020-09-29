
#include "idlutil.hpp"

#include "redux/util/astro_comp.hpp"
#include "redux/util/solar_system_geo.hpp"

#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"
#include "redux/util/point.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

#include <sys/time.h>

using namespace redux::util;
using namespace redux;
using namespace std;


// Coordinate system
//
// Pixels: dx, dy,   unit = pixel
// Raw:    dx, dy,   unit = solar radius
// SOHO:   dx, dy,   unit = arcsec
// Helio:  lon, lat, unit = radians
//
// Positive direction is west and north
// Negative direction is east and south

using namespace std;

namespace {
    
    double julday(struct tm tm) {
        
        static const long IGREG = 15 + 31L * (10 + 12L * 152);
        long jul, jm;

        tm.tm_year += 1900;

        if( tm.tm_year < 0 ) {
            tm.tm_year++;
        }
        
        if(tm.tm_mon > 3) {
            jm = tm.tm_mon + 2;
        } else {
            jm = tm.tm_mon + 14;
            tm.tm_year--;
        }

        jul = (long)(floor(365.25 * tm.tm_year) + floor(30.6001 * jm) + tm.tm_mday + 1 + 1720995);

        if( tm.tm_mday + 1 + 31L * (tm.tm_mon + 1 + 13L * tm.tm_year) >= IGREG ) {
            long ja = (long)(0.01* tm.tm_year);
            jul += 2 - ja + (long)(0.25*ja);
        }

        return jul + tm.tm_hour / 24.0 + tm.tm_min / 1440.0 + tm.tm_sec / 86400.0;
        
    }

    double julday( time_t t = time(0) ) {
        struct tm tm;
        gmtime_r(&t, &tm);
        return julday(tm);
    }


    void sunpb( double *p0, double *b0, double* r0 = nullptr, double jd = std::numeric_limits<float>::infinity() ) {
        
        if( !isfinite(jd) ) {
            jd = julday();
        }
        
        // correction for Ephemeris date
        double dt = -15 + pow(jd - 2382148, 2) / 41048480;
        jd += dt / 86400;

        // Inclination of solar equator
        static const double IS = 0.12653637;

        // lon. of asc. node of solar equator
        double K = 1.2857260 + 0.024361886 * (jd - 2396758) / 36525;

        // Jul. centuries since 2000
        double T = (jd - 2451545) / 36525;
        double T2 = T * T;

        // mean lon. of sun
        double L0 = 4.8950630 + 628.3319667 * T + 5.2918383e-6 * T2;

        // Eccentricity of earth orbit
        double e = 0.016708617 - 0.000042037 * T - 0.0000001236 * T2;

        // mean anomaly of sun
        double M = 6.24006 + 628.30195532 * T - 2.7209683e-6 * T2 - 8.37758e-9 * T * T2;

        // suns equation of center
        double C = (0.03341607  - 8.40725e-5 * T - 2.4434e-7 * T2) * sin(M)
          + (0.00034894  - 1.76278e-6 * T) * sin(2 * M)
          +  5.0614548e-6 * sin(3 * M);

        // sun's true longitude
        double TL = L0 + C;

        // True anomaly of sun
        double nu = M + C;

        // Radius vector of earth
        double R = (1.000001018 *(1 - e * e)) / (1 + e * cos(nu));

        // mean obliquity of ecliptic
        double eps = 0.40909280  - 2.2696551e-4 * T - 2.8604e-9 * T2 + 8.7897e-9 * T * T2;

        // Correction for aberration
        double ab = -20.4898 / 3600.0 / R;

        double om = 2.1824 - 33.757041 * T;
        double lam = TL - 9.9309234e-5 - 8.34267e-5 * sin(om) + ab;
        double lamp = lam + 15 / 3600.0;

        double x = atan(-cos(lamp) * tan(eps));
        double y = -cos(lam - K) * tan(IS);

        *p0 = x + y;
        *b0 = -asin(sin(lam - K) * tan(IS));

        if( r0 ) {
            *r0 = 959.63 / R;
        }
        
    }
    
    // Taken from SST turret code
    double solarLocalTilt( const Observatory& observatory, const struct tm *pTime = 0, int msec = 0, double* derotation = nullptr ) {
        
        struct tm tm;

        if( !pTime ) {
            struct timeval tv;
            gettimeofday(&tv, 0);
            gmtime_r(&tv.tv_sec, &tm);
            pTime = &tm;
            msec = tv.tv_usec / 1000;
        }

        AstroTime astrotime;
        double nut_long;
        double true_obliquity;
        double mean_obliquity;
        double nut_obl;
        double parallax_object;
        double hour_angle;
        double apparent_geo_long;
        double radius_vec;
        double aberration;
        double refrac;
        double apparent_geo_lat;
        double parallactic_angle;
        double azimuth;
        double elevation;
        double ra;
        double decl;

        //  Date, Time and Julian Day at 0 hours UT (jd0)
        astrotime.year = pTime->tm_year + 1900;
        astrotime.month = pTime->tm_mon + 1;
        astrotime.day = pTime->tm_mday;
        astrotime.hour = pTime->tm_hour;
        astrotime.min = pTime->tm_min;
        astrotime.sec = pTime->tm_sec;
        astrotime.day_frc = astrotime.hour / 24.0 + astrotime.min / 1440.0 + astrotime.sec / 86400.0 + msec / 86400000.0;
        astrotime.jd0 = redux::util::julday( astrotime.year, astrotime.month, astrotime.day, 0.0 );

        //  Difference between Dynamical Time and Universal Time in
        //  fractions of a day (td_corr):
        astrotime.td_corr = time_dyn( astrotime.year, astrotime.month ) / 86400.0;

        //  Nutation in longitude and latitude. Mean and true obliquity:
        nutation_obliquity(astrotime.jd0 + 0.5, &nut_long, &nut_obl, &mean_obliquity, &true_obliquity);

        //  Greenwich Apparent Sidereal Time at 0 hours (gast0):
        astrotime.gast0 = greenwich_app_sid_time_0(astrotime.jd0, nut_long, true_obliquity);

        //  Get solar aberration and solar semi diameter:
        //  REFERENCE: Astronomical Algorithms, by Jean Meeus, p. 155 (solar
        //             aberration) and chapter 53 (solar semi diameter).
        aberration = 20.4898 / 3600.0 * DtoR;

        redux::util::sun_geo(astrotime.jd0 + 0.5, nut_long,  aberration, true_obliquity, &apparent_geo_long, &apparent_geo_lat, &radius_vec, &ra, &decl);

        aberration = aberration / radius_vec;

        //  Get right ascension (ra_td0) and declination (decl_td0) at
        //  0 hours Dynamical Time for the Sun:
        redux::util::sun_geo(astrotime.jd0, nut_long, aberration, true_obliquity, &apparent_geo_long, &apparent_geo_lat, &radius_vec, &ra, &decl);

        parallax_object = asin((sin(8.794 / 3600.0) * DtoR) / radius_vec);

        astrotime.ut = astrotime.day_frc * 2.0 * M_PI;

        local_position( observatory, parallax_object, astrotime.gast0, astrotime.ut, &hour_angle, &azimuth, &elevation, &ra, &decl, &astrotime.last, &refrac);
        //fprintf(stderr, "az = %8.3lf  el = %8.3lf\n", r2d(azimuth), r2d(elevation));
        if( derotation ) {
            *derotation = -(azimuth - elevation);
        }
        //fprintf(stderr, "%le %le %le %le\n", parallax_object, astrotime.gast0, astrotime.ut, r2d(elevation));

        parallactic_angle = atan2(sin(hour_angle ), tan(observatory.latitude) * cos(decl) - sin(decl) * cos(hour_angle));

    #if 0
        double axis;
        double stony_lat_cen;
        //  Get solar axis position angle (axis) and the
        //  heliographic latitude (stony_lat_cen) and longitude
        //  (stony_long_cen) for the center of solar disk:
        //  Equivalent to p0.
        phys_sun(astrotime.jd0 + 0.5, apparent_geo_long, nut_long, true_obliquity, &axis, &stony_lat_cen);

        return -axis + parallactic_angle;
    #else
        return parallactic_angle;
    #endif
        
    }
    
    
    // Convert position in solar radii to SOHO coordinates
    PointF pos2soho( const PointF& pos, double jd = julday() ) {
        
        double p0, b0, r0;
        sunpb( &p0, &b0, &r0, jd );
        return pos * (float)r0;
        
    }
        
        
    PointF pos2helio( const PointF& pos, double jd = julday() ) {
        
        double p0, b0, r0;

        sunpb( &p0, &b0, &r0, jd );

        // Calculate heliocentric coordinates
        double rho = hypot(pos.x, pos.y);
        double z = sqrt(1 - rho * rho);

        // Rotate in Y-Z plane
        double z_ = z * cos(b0) + pos.y * sin(b0);
        double y_ = pos.y * cos(b0) - z * sin(b0);

        // Convert to longitude/lattitude
        PointF helio;
        helio.y = asin(y_);
        helio.x = asin( pos.x / cos(helio.y));

        if( z_ < 0 ) {
            helio.x = -M_PI - helio.x;
            if(helio.x < -M_PI)
                helio.x += 2 * M_PI;
        }

        return helio;
    }


}


/*
class config config(string(getenv("HOME")) + "/.dotguiderrc");

static double fovw = 89.6; // field of view width in arcsec
static double fovh = 71.68; // height
static double fovoverlap = 0.2; // desired fraction of overlap for mosaics

static bool fliph;

static double longitude;
static double latitude;
static double air_pressure;
static double air_temp;
static bool altaz;

static double r0; // theoretical solar radius in arcsec
static double p0; // solar orientation in radians
static double b0;
static double l0;
static double tilt; // angle of alt/az axes
static double derotation; // rotation of image due to alt/az mount

static glm::vec2 offset;
static double radius;
static double scale;
static double rotation;

enum orientation {
	HELIOCENTRIC = 0,
	GEOCENTRIC,
	EQUATORIAL,
	ALTAZ,
	IMAGE,
	MAXORIENTATIONS
} orientation;

const char *orientation_name[MAXORIENTATIONS] = {
	"heliocentric",
	"geocentric",
	"equatorial",
	"altaz",
	"image",
};

static double d2r(double d) {
	return d / 180.0 * M_PI;
}

static double r2d(double d) {
	return 180.0 * d / M_PI;
}

static float d2rf(float d) {
	return d / 180.0 * M_PI;
}

static float r2df(float d) {
	return 180.0 * d / M_PI;
}



// Convert position in solar radii to SOHO coordinates
static glm::vec2 pos2soho(glm::vec2 pos, double jd = julday()) {
	double p0, b0, r0;
	glm::vec2 soho;

	sunpb(&p0, &b0, &r0, jd);

	return pos * (float)r0;
}

// From SOHO to position in solar radii
static glm::vec2 soho2pos(glm::vec2 soho, double jd = julday()) {
	double p0, b0, r0;
	glm::vec2 pos;

	sunpb(&p0, &b0, &r0, jd);

	return soho / (float)r0;
}

static glm::vec2 pos2helio(glm::vec2 pos, double jd = julday()) {
	double p0, b0, r0;
	glm::vec2 helio;

	sunpb(&p0, &b0, &r0, jd);

	// Calculate heliocentric coordinates
	double rho = hypot(pos.x, pos.y);
	double z = sqrt(1 - rho * rho);

	// Rotate in Y-Z plane
	double z_ = z * cos(b0) + pos.y * sin(b0);
	double y_ = pos.y * cos(b0) - z * sin(b0);

	// Convert to longitude/lattitude
	helio.y = asin(y_);
	helio.x = asin(pos.x / cos(helio.y));

	if(z_ < 0) {
		helio.x = -M_PI - helio.x;
		if(helio.x < -M_PI)
			helio.x += 2 * M_PI;
	}

	return helio;
}

static glm::vec2 soho2helio(glm::vec2 soho, double jd = julday()) {
	double p0, b0, r0;
	glm::vec2 helio;

	sunpb(&p0, &b0, &r0, jd);

	// Convert from arcsec to fraction of solar radius
	soho /= r0;

	// Calculate heliocentric coordinates
	double rho = hypot(soho.x, soho.y);
	double z = sqrt(1 - rho * rho);

	// Rotate in Y-Z plane
	double z_ = z * cos(b0) + soho.y * sin(b0);
	double y_ = soho.y * cos(b0) - z * sin(b0);

	// Convert to longitude/lattitude
	helio.y = asin(y_);
	helio.x = asin(soho.x / cos(helio.y));

	if(z_ < 0) {
		helio.x = -M_PI - helio.x;
		if(helio.x < -M_PI)
			helio.x += 2 * M_PI;
	}

	return helio;
}

static glm::vec2 helio2soho(glm::vec2 helio, double jd = julday()) {
	double p0, b0, r0;
	glm::vec2 soho;

	sunpb(&p0, &b0, &r0, jd);
	// Convert to heliocentric coordinates
	soho.x = cos(helio.y) * sin(helio.x);
	soho.y = sin(helio.y) * cos(b0) + cos(helio.y) * cos(helio.x) * sin(b0);

	soho *= r0;

	return soho;
}

static glm::vec2 helio2pos(glm::vec2 helio, double jd = julday()) {
	return soho2pos(helio2soho(helio, jd), jd);
}

// Apply rotation compensation to heliographic coordinates
static glm::vec2 rotcomp(glm::vec2 helio, double start, double end = julday()) {
	// Taken from http://en.wikipedia.org/wiki/Solar_rotation
	static const double A = 14.713;
	static const double B = -2.396;
	static const double C = -1.787;
	static const double synodic = 24.47 / 26.24;
	double sinphi = sin(helio.y);
	double sin2phi = sinphi * sinphi;
	double sin4phi = sin2phi * sin2phi;
	double omega = A + B * sinphi + C * sin4phi;
	helio.x += d2r((end - start) * omega * synodic);
	return helio;
}

void show_help(const char *argv0) {
	fprintf(stderr, "Usage: %s [options] [input]\n", argv0);
	fprintf(stderr, "\n");
	fprintf(stderr, " -a, --average <N>         Average every N frames.\n");
	fprintf(stderr, " -x, --dx <pixels>         Set the x offset of the solar image.\n");
	fprintf(stderr, " -y, --dy <pixels>         Set the y offset of the solar image.\n");
	fprintf(stderr, " -s, --scale <arcs/pixel>  Set the scale of the solar image.\n");
	fprintf(stderr, " -R, --rotation <degrees>  Set rotation of PIG camera.\n");
	fprintf(stderr, " -h, --human-readable      Output time in human readable format.\n");
	fprintf(stderr, " -H, --helio               Output heliographic coordinates.\n");
	fprintf(stderr, " --help                    Display this help and exit.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Report bugs to Guus Sliepen <Guus.Sliepen@astro.su.se>.\n");
}
*/


/*
crispex cmd:
convertlog --dx 31.92 --dy 14.81 --rotation 84.87 --scale 4.935 -a 16

~mats/bin/convertlog --dx 31.92 --dy 14.81 --rotation 84.87 --scale 4.935 -a 16 rmslog_guidercams.1 > rmslog_guidercams_old

.dotguiderrc:
offset.x = 23.79
offset.y = 29.90
radius = 195.79
scale = 4.32
rotation = 84.60
*/

typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_INT average;
    float x;
    float y;
    float radius;
    float rotation;
    float scale;
    IDL_INT helio;
    IDL_INT help;
    IDL_INT human;
    IDL_INT verbose;
} CVTLOG_KW;

// NOTE:  The keywords MUST be listed in alphabetical order !!
static IDL_KW_PAR cvtlog_pars[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "AVERAGE",    IDL_TYP_INT, 1, IDL_KW_ZERO,   0, (char*) IDL_KW_OFFSETOF2(CVTLOG_KW, average) },
    { (char*) "DX",       IDL_TYP_FLOAT, 1,           0,   0, (char*) IDL_KW_OFFSETOF2(CVTLOG_KW, x) },
    { (char*) "DY",       IDL_TYP_FLOAT, 1,           0,   0, (char*) IDL_KW_OFFSETOF2(CVTLOG_KW, y) },
    { (char*) "HELIO",      IDL_TYP_INT, 1, IDL_KW_ZERO,   0, (char*) IDL_KW_OFFSETOF2(CVTLOG_KW, helio) },
    { (char*) "HELP",       IDL_TYP_INT, 1, IDL_KW_ZERO,   0, (char*) IDL_KW_OFFSETOF2(CVTLOG_KW, help) },
    { (char*) "HUMAN",      IDL_TYP_INT, 1, IDL_KW_ZERO,   0, (char*) IDL_KW_OFFSETOF2(CVTLOG_KW, human) },
    { (char*) "RADIUS",   IDL_TYP_FLOAT, 1,           0,   0, (char*) IDL_KW_OFFSETOF2(CVTLOG_KW, radius) },
    { (char*) "ROTATION", IDL_TYP_FLOAT, 1,           0,   0, (char*) IDL_KW_OFFSETOF2(CVTLOG_KW, rotation) },
    { (char*) "SCALE",    IDL_TYP_FLOAT, 1,           0,   0, (char*) IDL_KW_OFFSETOF2(CVTLOG_KW, scale) },
    { (char*) "VERBOSE",    IDL_TYP_INT, 1, IDL_KW_ZERO,   0, (char*) IDL_KW_OFFSETOF2(CVTLOG_KW, verbose) },
    { (char*) "X",        IDL_TYP_FLOAT, 1,           0,   0, (char*) IDL_KW_OFFSETOF2(CVTLOG_KW, x) },
    { (char*) "Y",        IDL_TYP_FLOAT, 1,           0,   0, (char*) IDL_KW_OFFSETOF2(CVTLOG_KW, y) },
    { NULL }
};
/*
    static struct option const long_options[] = {
        {"help", no_argument, NULL, 1},
        {"dx", required_argument, NULL, 'x'},
        {"dy", required_argument, NULL, 'y'},
        {"scale", required_argument, NULL, 's'},
        {"rotation", required_argument, NULL, 'R'},
        {"soho", no_argument, NULL, 's'},
        {"human-readable", no_argument, NULL, 'h'},
        {"helio", no_argument, NULL, 'H'},
        {"average", required_argument, NULL, 'a'},
        {NULL, 0, NULL, 0}
    };
*/
string cvtlog_info( int lvl ) {
    
    string ret = "RDX_CONVERTLOG";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"     ");          // newline if lvl>1
        ret += "   Syntax:   output = rdx_convertlog( input, /KEYWORDS )\n";
        ret += "             or       rdx_convertlog, logfile, outputfile, /KEYWORDS\n";
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      AVERAGE             Calculate averages over specified number of entries. (integer)\n"
                    "      DX or X             Specify x-offset solar image. (float)\n"
                    "      DY or Y             Specify y-offset solar image. (float)\n"
                    "      HELIO               Convert to Heliocentric coordinates. (bool)\n"
                    "      HUMAN               Convert timestamps to human readable form. (bool)\n"
                    "      RADIUS              Rotation. (float)\n"
                    "      ROTATION            Rotation. (float, degrees)\n"
                    "      SCALE               Scale. (float)\n"
                    "      HELP                Display this info. (bool)\n"
                    "      VERBOSE             More detailed output. -1 suppresses all output (integer).\n";
        }
    } else ret += "\n";

    return ret;
    
}


IDL_VPTR convertlog( int argc, IDL_VPTR* argv, char* argk ) {
    
    CVTLOG_KW kw;
    kw.x = std::numeric_limits<float>::infinity();
    kw.y = std::numeric_limits<float>::infinity();
    kw.rotation = std::numeric_limits<float>::infinity();
    kw.scale = std::numeric_limits<float>::infinity();
    int nPlainArgs = IDL_KWProcessByOffset( argc, argv, argk, cvtlog_pars, (IDL_VPTR*)0, 255, &kw );
    
    if( kw.help ) {
        cout << cvtlog_info(2) << endl;
        return IDL_StrToSTRING( (char*)"" );
    }

    if( nPlainArgs < 1 ) {
        return IDL_StrToSTRING( (char*)"" );
    }

    if( !isfinite(kw.x) || !isfinite(kw.y) ) {
        if( kw.verbose > -1 ) cout << "WARNING: no x/y offset given!" << endl;
        kw.x = kw.y = 0;
    }

    if( !isfinite(kw.rotation) ) {
        if( kw.verbose > -1 ) cout << "WARNING: no rotation given!" << endl;
        kw.rotation = 0;
    } else kw.rotation *= DtoR;     // convert to radians

    if( !isfinite(kw.radius) && !isfinite(kw.scale) ) {
        if( kw.verbose > -1 ) cout << "WARNING: no image scale given!" << endl;
        kw.radius = 0;
        kw.scale = 0;
    }

    IDL_VPTR input = argv[0];
    string tmpS;

    IDL_ENSURE_STRING( input );
    IDL_ENSURE_SIMPLE( input );

    if ( !(input->flags & IDL_V_ARR) ) {
        tmpS = string(IDL_VarGetString( input ));
    } else {
        IDL_ARRAY* strarr = input->value.arr;
        IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>(strarr->data);
        for( int i=0; i<strarr->n_elts; ++i ) {
            tmpS += string( strptr[i].s ) + "\n";
        }
    }
    
    stringstream iss( tmpS );
    stringstream oss;
    
    double j0 = -1;
    double sum_t = 0, sum_cx = 0, sum_cy = 0, sum_cr = 0;
    int sum = 0;

    if( kw.human ) {
        oss << "#timestamp         ";
    } else {
        oss << "#unixtime        ";
    }
    
    if( kw.helio ) {
        oss << " helio_lon helio_lat";
    } else {
        oss << " soho_x     soho_y";
    }
    oss << endl;

    const double past = 946681200;
    const double future = time(NULL) + 157762800;
    
    // SST specific parameters
    Observatory sst;
    sst.latitude = 28.7597 * (DtoR);
    sst.longitude = 17.8756 * (DtoR);
    sst.air_pressure = 700;
    sst.air_temp = 20;
    sst.refraction_coeff = 0.00115;

    
    string line;
    char buf[1024];
    while( iss.good() ) {
        
        getline( iss, line );
        if( line.empty() || line[0] == '#' ) continue;
        
        double t, cx, cy, cr;
        if( sscanf( line.data(), "%lf %*f %*f %lf %lf %lf", &t, &cx, &cy, &cr) != 4 ) {
            if( kw.verbose > 1 ) cout << "Error parsing line: \"" << line << "\"" << endl;
            continue;
        }

        if( !t || !cr || !isfinite(t) || !isfinite(cx) || !isfinite(cy) || !isfinite(cr) ) {
            continue;
        }
        
        if( t < past || t > future ) {
            continue;
        }
        
        if( kw.average ) {
            sum_t += t;
            sum_cx += cx;
            sum_cy += cy;
            sum_cr += cr;
            sum++;
            if( sum < kw.average ) {
                continue;
            }
            t = sum_t / sum;
            cx = sum_cx / sum;
            cy = sum_cy / sum;
            cr = sum_cr / sum;
            sum_t = sum_cx = sum_cy = sum_cr = sum = 0;
        }

        time_t sec = lrint(t);
        static struct tm tm;
        if( !gmtime_r(&sec, &tm) ) {
            return IDL_StrToSTRING( (char*)"" );
        }
        
        double jd = julday(sec);
        if( j0 < 0 ) {
            j0 = jd;
        }
        
        double p0, b0, r0;
        double derotation;
        sunpb( &p0, &b0, &r0, jd );

        double tilt = solarLocalTilt( sst, &tm, 0, &derotation );
        if( isfinite( kw.scale ) ) {
            kw.radius = r0 / kw.scale;
        }

        double angle = kw.rotation + derotation - p0 + tilt;
        PointF cen;
        cen.x = (cx - kw.x) * cos(angle) - (cy - kw.y) * sin(angle);
        cen.y = (cx - kw.x) * sin(angle) + (cy - kw.y) * cos(angle);
        cen.x /= kw.radius;
        cen.y /= kw.radius;

        PointF soho = pos2soho(cen, jd);
        PointF helio = pos2helio(cen, jd);
        
        if( kw.human ) {
            if( !strftime( buf, sizeof buf, "%F %T ", &tm ) ) {
                return IDL_StrToSTRING( (char*)"" );
            }
            oss << buf;
        } else {
            snprintf( buf, sizeof buf, "%lf ", t);
            oss << buf;
        }
        
        if( kw.helio ) {
            helio *= RtoD;
            snprintf( buf, sizeof buf, "%lf %lf\n", helio.x, helio.y );
        } else {
            snprintf( buf, sizeof buf, "%lf %lf\n", soho.x, soho.y );
        }
        
        oss << buf;

    }
    
    return IDL_StrToSTRING( (char*)oss.str().c_str() );

}


void convertlog_proc( int argc, IDL_VPTR* argv, char* argk ) {

    CVTLOG_KW kw;
    kw.x = std::numeric_limits<float>::infinity();
    kw.y = std::numeric_limits<float>::infinity();
    kw.rotation = std::numeric_limits<float>::infinity();
    kw.scale = std::numeric_limits<float>::infinity();
    int nPlainArgs = IDL_KWProcessByOffset( argc, argv, argk, cvtlog_pars, (IDL_VPTR*)0, 255, &kw );
    
    if( kw.help ) {
        cout << cvtlog_info(2) << endl;
        return;
    }

    if( nPlainArgs < 1 ) {
        cout << "No input specified." << endl;
        return;
    }

    IDL_ENSURE_STRING( argv[0] );
    IDL_ENSURE_SIMPLE( argv[0] );
    IDL_ENSURE_STRING( argv[1] );
    IDL_ENSURE_SIMPLE( argv[1] );
    
    string infile( IDL_VarGetString( argv[0] ) );
    string outfile( IDL_VarGetString( argv[1] ) );
    
    std::string content;
    try {
        std::ifstream in( infile, std::ios::in | std::ios::binary );
        if( in ) {
            in.seekg( 0, std::ios::end );
            content.resize( in.tellg() );
            in.seekg( 0, std::ios::beg );
            in.read( &content[0], content.size() );
            in.close();
        } else {
            throw std::ios_base::failure( "Failed to open file: " + infile );
        }
    } catch( std::exception& e ) {
        cout << e.what() << endl;
        return;
    }
    
    IDL_VPTR tmpIn = IDL_StrToSTRING( (char*)content.c_str() );
    std::swap( argv[1], tmpIn );
    IDL_VPTR tmpOut = convertlog( argc-1, argv+1, argk+1 ); // skip first argv and call convertlog for processing
    std::swap( argv[1], tmpIn );
    
    std::string ocontent( IDL_VarGetString( tmpOut )  );
    
    IDL_Freetmp(tmpIn);
    IDL_Freetmp(tmpOut);
    
    if( outfile.empty() ) {
        outfile = infile + "_converted";
    }
    
    try {
        std::ofstream out( outfile, std::ofstream::trunc );
        if( out ) {
            out << ocontent;
        } else {
            throw std::ios_base::failure( "Failed to write to file: " + outfile );
        }
    } catch( std::exception& e ) {
        cout << e.what() << endl;

    }
    
    
}


namespace {
    static int dummy RDX_UNUSED =
    IdlContainer::registerRoutine( {{(IDL_SYSRTN_GENERIC)convertlog_proc}, (char*)"RDX_CONVERTLOG", 0, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 0, cvtlog_info ) +
    IdlContainer::registerRoutine( {{(IDL_SYSRTN_GENERIC)convertlog}, (char*)"RDX_CONVERTLOG", 0, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, 0 );
}
