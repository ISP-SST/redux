#ifndef REDUX_UTIL_ASTRO_COMP_HPP
#define REDUX_UTIL_ASTRO_COMP_HPP

#include "redux/util/datautil.hpp"

#include <math.h>

/**
 * astronomical computations   
 */


/* the planets */
#define MERCURY 1
#define VENUS 2
#define EARTH 3
#define MARS 4
#define JUPITER 5
#define SATURN 6
#define URANUS 7
#define NEPTUNE 8

#include <time.h>

namespace redux {
    
    namespace util {
        
        static double DtoR RDX_UNUSED = ( M_PI / 180.0 );
        static double RtoD RDX_UNUSED = ( 180.0 / M_PI );

        /**
         * structure for storing observatory particulars
         */
        typedef struct Observatory {

            char location[50];
            double latitude;          // LP: 28.7597
            double longitude;         // LP: 17.8756
            double air_pressure;      // LP: 700.0
            double air_temp;          // LP: 20.0
            double refraction_coeff;  // LP: 0.00115

            /* distance from rotation axis and height above equatorial plane (in 
               Earth radii) */
            double parallax_sin;      // LP: 0.47844447
            double parallax_cos;      // LP: 0.87766283

        } Observatory;

        /**
         * structure for storing RA/decl coordinates
         */
        typedef struct RaDecl {
          double ra;                /* right ascension */
          double decl;              /* declination */
        } RaDecl;


        /**
         * structure for storing broken down RA/decl coords
         */
        typedef struct HoursAndDegrees {
          int ra_h;                 /* RA hours (i hr = 15 degrees) */
          int ra_m;                 /* RA minutes */
          double ra_s;              /* RA seconds */
          int decl_d;               /* Decl degrees */
          int decl_m;               /* Decl minutes */
          int decl_s;               /* Decl seconds */
        } HoursAndDegrees;


        /**
         * structure for equatorial coordinates
         */
        typedef struct EqCoords {
          RaDecl original;          /* orginal coords */
          RaDecl corrected;         /* corrected coords */
          double prop_ra;           /* proper motion of star in ra */
          double prop_decl;         /* proper motion of star in decl */
          double vmagn;             /* visual magnitude */
          double epoch;             /* initial Epoch for star position */
          char constellation[4];    /* name of constellation */
          char object_name[100];    /* name of object */
          int star_number;          /* star number */
          double screen_x;          /* x screen coord */
          double screen_y;          /* y screen coord */
        } EqCoords;


        /**
         * structure for astronomical time
         */
        typedef struct AstroTime {
          int year, month, day;  /* Year, Month, Day */
          int hour ,min, sec;    /* Hours, Minutes, Seconds of time */
          double day_frc;        /* Fractions of day */
          double td_corr;        /* Dynamical Time Correction */
          double jd;             /* Julian Day, including fractions */
          double jd0;            /* Julian Day (UT 0 hours) */
          double jde;            /* Julian Ephemeris Day, incl. fractions */
          double ut;             /* Universal Time (UT) */
          double gast0;          /* Greenwich Apparent Sidereal Time 0h */
          double last;           /* Local Apparent Sidereal Time */
          time_t astrocalc_time; /* ??? */
          double hour_angle;
          double local_sidereal;
          double greenwich_sidereal;
        } AstroTime;


        void get_current_time( AstroTime *current, struct tm *offset, double obs_longitude );
                                                             
        double time_dyn( int year, int month=0 );

        double greenwich_app_sid_time_0 (double jd0, double nut_long, double true_obliquity);

        void nutation_obliquity (double jd, double *nut_long, double *nut_obl, double *mean_obliquity, double *true_obliquity);

        void phys_sun (double jd, double sun_apparent_geo_long,
                              double nut_long, double true_obliquity,
                              double *solar_axis, double *stony_lat_cen);



        double julday (int year, int month, int day, double day_frc);


        double refraction( const Observatory& observatory, double elevation );

        void local_position( const Observatory& observatory,
                                    double parallax_object, double gast0, double ut,
                                    double *ha, double *azimuth, double *elevation,
                                    double *ra, double *decl, double *last,
                                    double *refrac );

        void star_geo (double jd, double nut_long, double nut_obl,
                              double mean_obliquity, double sun_apparent_geo_long,
                              double solar_aberration, EqCoords *eqCoords );

    }   // util
    
}   // redux


#endif
