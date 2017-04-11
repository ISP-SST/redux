#include "redux/util/astro_comp.hpp"

#include "redux/util/solar_system_geo.hpp"

/**
 * astronomical computations
 */

#include <map>


using namespace redux::util;


namespace {

    /**
     * adds an offset in years, months, days etc.
     */
    time_t add_time_offset( time_t time, struct tm *offset ) {

        struct tm *tm = localtime( &time );
        time_t new_time;

        tm->tm_year += offset->tm_year;
        tm->tm_mon += offset->tm_mon;
        tm->tm_mday += offset->tm_mday;
        tm->tm_hour += offset->tm_hour;
        tm->tm_min += offset->tm_min;

        if( ( new_time = mktime( tm ) ) == -1 ) {
            return time;
        }

        return new_time;
    }



    /**
       calculates greenwich sidereal time

       NOTES:

       Calculates the sidereal time at Greenwich.

    */
    double get_greenwich_sidereal_time( AstroTime *astrotime ) {

        double UT;
        double T;
        double GST;

        /* get universal time in decimal hours */
        UT =
            (double)astrotime->hour +
            (double)astrotime->min / 60.0 +
            (double)astrotime->sec / 3600.0;

        T = ( astrotime->jd0 - 2451545.0 ) / 36525.0;

        /* calculate greenwich sidereal time (GST) */
        GST =
            6.697374558 + ( 2400.051336*T ) +
            ( 0.000025862*T*T )+( UT*1.0027379093 );

        /* adjust to be in range 0 - 24 */
        if( GST > 24.0 )
            while( GST > 24.0 )
                GST = GST - 24.0;

        if( GST < 0.0 )
            while( GST < 0.0 )
                GST = GST + 24.0;

        return GST;
    }


    /**
       calculates local sidereal time

       NOTES:

       Local sidereal time is the right ascension of a star on the observers
       meridian. One sidereal day corresponds to the time taken for the Earth
       to rotate once with respect to the stars and lasts approximately
       23 h 56 min.
    */
    double get_local_sidereal_time( double greenwich_sidereal,
                                    double obs_longitude ) {
        return greenwich_sidereal + ( ( -obs_longitude*RtoD ) / 15 );
    }


    /**
       Corrects solar system object ra/decl for parallax.

       NOTES:

       the observatory paralax constant (held in the Observatory structure) are:

       rho*sin( phi' )
       rho*cos( phi' )

       the values of which can be obtained (at the time of writing) from:

       http://cfa-www.harvard.edu/iau/lists/ObsCodes.html

       for each observatory.

       REFERENCES:

       Astronomical Algorithms, by Jean Meeus, chapter 39

    */
    __attribute__((unused))
    void correct_for_parallax( const Observatory& observatory,
                               double *ra, double *decl,
                               double parallax_object,
                               double *ha ) {

        double ra_corr;

        ra_corr = atan( -observatory.parallax_cos *
                        sin( parallax_object ) * sin( *ha ) /
                        ( cos( *decl ) - observatory.parallax_cos *
                          sin( parallax_object ) *
                          cos( *ha ) ) );

        *ra = *ra + ra_corr;
        *ra = fmod( *ra, 2.0*M_PI );

        if( *ra < 0.0 ) *ra = *ra + 2.0*M_PI;

        *decl = atan( ( sin( *decl ) - observatory.parallax_sin *
                        sin( parallax_object ) ) *
                      cos( ra_corr ) /
                      ( cos( *decl ) - observatory.parallax_cos *
                        sin( parallax_object ) *
                        cos( *ha ) ) );

    }


    /**
       correct ra and decl for nutation

       nutation is the short-period oscillation superimposed on the
       precession of the rotational axis of the earth.

       it amounts to a maximum of 15 arcsecs during a period of about
       18.6 years and is caused by changes in the moon's orbit over this period.
    */
    void correct_for_nutation( double nut_long,
                               double nut_obl,
                               double mean_obliquity,
                               double *ra, double *decl ) {
        *ra = *ra + nut_long *
              ( cos( mean_obliquity ) +
                sin( mean_obliquity ) *
                sin( *ra ) * tan( *decl ) ) -
              nut_obl * ( cos( *ra ) * tan( *decl ) );

        *decl = *decl + nut_long *
                ( sin( mean_obliquity ) * cos( *ra ) ) +
                nut_obl * sin(*ra);
    }


    /**
       correct for annual aberration

       NOTES:

       annual aberration is the aberration caused by the earth's
       orbital motion around the sun
    */
    void correct_for_annual_aberration( double t,
                                        double nut_long,
                                        double mean_obliquity,
                                        double sun_apparent_geo_long,
                                        double solar_aberration,
                                        double *ra,
                                        double *decl ) {

        double k, e, p;
        double sun_true_geometric_long;

        k = 20.49552 / 3600.0 * DtoR;
        e = 0.016708617 - 0.000042037 * t - 0.0000001236 * t * t;
        p = ( 102.93735 + 0.71953 * t + 0.00046 * t * t ) * DtoR;
        sun_true_geometric_long = sun_apparent_geo_long -
                                  nut_long + solar_aberration;

        // correct right ascension
        *ra = *ra - k * (cos( *ra ) * cos( sun_true_geometric_long ) *
                         cos( mean_obliquity ) +
                         sin( *ra ) * sin( sun_true_geometric_long ) ) /
              cos( *decl ) +
              e * k * ( cos( *ra ) * cos( p ) * cos( mean_obliquity ) +
                        sin( *ra ) * sin( p ) ) / cos( *decl );

        *ra = fmod( *ra, 2.0 * M_PI );

        if( *ra < 0.0 ) *ra = *ra + 2.0 * M_PI;

        // correct declination
        *decl = *decl - k * ( cos( sun_true_geometric_long ) *
                              cos( mean_obliquity ) *
                              ( tan( mean_obliquity ) *
                                cos( *decl ) -
                                sin( *ra ) * sin( *decl ) ) +
                              cos( *ra ) * sin( *decl ) *
                              sin( sun_true_geometric_long ) ) +
                e * k * ( cos( p ) * cos( mean_obliquity ) *
                          ( tan( mean_obliquity ) * cos( *decl ) - sin( *ra ) *
                            sin( *decl ) ) +
                          cos( *ra ) * sin( *decl ) * sin( p ) );
    }


    // List of delta-T with half-year intervals (data is for january & july)
    const std::map<double,double> dynTime = { { 1972.0, 42.23 }, { 1972.5, 42.80 }, { 1973.0, 43.37 }, { 1973.5, 43.93 },
        { 1974.0, 44.49 }, { 1974.5, 44.99 }, { 1975.0, 45.48 }, { 1975.5, 45.97 }, { 1976.0, 46.46 }, { 1976.5, 46.99 },
        { 1977.0, 47.52 }, { 1977.5, 48.03 }, { 1978.0, 48.53 }, { 1978.5, 49.06 }, { 1979.0, 49.59 }, { 1979.5, 50.07 },
        { 1980.0, 50.54 }, { 1980.5, 50.96 }, { 1981.0, 51.38 }, { 1981.5, 51.78 }, { 1982.0, 52.17 }, { 1982.5, 52.57 },
        { 1983.0, 52.96 }, { 1983.5, 53.38 }, { 1984.0, 53.79 }, { 1984.5, 54.07 }, { 1985.0, 54.34 }, { 1985.5, 54.61 },
        { 1986.0, 54.87 }, { 1986.5, 55.10 }, { 1987.0, 55.32 }, { 1987.5, 55.57 }, { 1988.0, 55.82 }, { 1988.5, 56.06 },
        { 1989.0, 56.30 }, { 1989.5, 56.58 }, { 1990.0, 56.86 }, { 1990.5, 57.22 }, { 1991.0, 57.57 }, { 1991.5, 57.94 },
        { 1992.0, 58.31 }, { 1992.5, 58.72 }, { 1993.0, 59.12 }, { 1993.5, 59.55 }, { 1994.0, 59.98 }, { 1994.5, 60.38 },
        { 1995.0, 60.78 }, { 1995.5, 61.20 }, { 1996.0, 61.63 }, { 1996.5, 61.96 }, { 1997.0, 62.29 }, { 1997.5, 62.63 },
        { 1998.0, 62.97 }, { 1998.5, 63.22 }, { 1999.0, 63.47 }, { 1999.5, 63.66 }, { 2000.0, 63.82 }, { 2000.5, 63.98 },
        { 2001.0, 64.09 }, { 2001.5, 64.20 }, { 2002.0, 64.30 }, { 2002.5, 64.41 }, { 2003.0, 64.47 }, { 2003.5, 64.55 },
        { 2004.0, 64.57 }, { 2004.5, 64.65 }, { 2005.0, 64.68 }, { 2005.5, 64.80 }, { 2006.0, 64.85 }, { 2006.5, 64.99 },
        { 2007.0, 65.15 }, { 2007.5, 65.34 }, { 2008.0, 65.45 }, { 2008.5, 65.63 }, { 2009.0, 65.78 }, { 2009.5, 65.95 },
        { 2010.0, 66.07 }, { 2010.5, 66.24 }, { 2011.0, 66.32 }, { 2011.5, 66.47 }, { 2012.0, 66.60 }, { 2012.5, 66.77 },
        { 2013.0, 66.91 }, { 2013.5, 67.13 }, { 2014.0, 67.28 }, { 2014.5, 67.49 }, { 2015.0, 67.64 }, { 2015.5, 67.86 },
        { 2016.0, 68.10 }, { 2016.5, 68.40 },
        { 2017.0, 68.50 },    // predicted
        { 2017.5, 68.85 },    // predicted
        { 2018.0, 69.20 },    // predicted
        { 2018.5, 69.50 },    // predicted
        { 2019.0, 69.80 },    // predicted
        { 2020.0, 70 },       // predicted
        { 2021.0, 71 },       // predicted
        { 2022.0, 71 },       // predicted
        { 2023.0, 72 },       // predicted
        { 2024.0, 72 },       // predicted
        { 2025.0, 73 },       // predicted
        { 2026.0, 73 }        // predicted

    };

}


/**
   gets the current universal time, and applies an offset
   if one is supplied
*/
void redux::util::get_current_time( AstroTime *current, struct tm *offset,
                                    double obs_longitude ) {

    time_t timebuf[128];
    struct tm *ltime;

    /* get current time */
    timebuf[0] = time(&timebuf[0]);

    if( offset != NULL )
        timebuf[0] = add_time_offset( timebuf[0], offset );

    current->astrocalc_time = timebuf[0];

    ltime = localtime( &timebuf[0] );

    current->year = 1900 + ltime->tm_year;
    current->month = 1 + ltime->tm_mon;
    current->day = ltime->tm_mday;
    current->hour = ltime->tm_hour;
    current->min = ltime->tm_min;
    current->sec = ltime->tm_sec;

    current->day_frc =
        current->hour/24.0 +
        current->min/1440.0 +
        current->sec/86400.0;

    current->jd0 = julday( current->year, current->month, current->day, 0.0 );
    current->td_corr = time_dyn( current->year ) / 86400.0;
    current->jde = current->jd0 + current->td_corr + current->day_frc;
    current->ut = current->day_frc * 2.0 * M_PI;

    /* get sidereal time */
    current->greenwich_sidereal = get_greenwich_sidereal_time( current );
    current->local_sidereal =
        get_local_sidereal_time( current->greenwich_sidereal, obs_longitude );
}


/**
   calculates the julian day from a given date at 0 hours universal
   time. function is valid between 1900 Mar. 1 and 2100 Feb. 28

   NOTES:

   julian day is the number of days since 4713 B.C. (-4712) at
   greenwich mean noon (12h UT).

   REFERENCE:

   Astronomical Algorithms, by Jean Meeus, chapter 7
   (Practical Ephemeris Calculations, by Oliver Montenbruck,
   chapter 2. Telescope Pointing and Tracking, page 3.)
*/
double redux::util::julday( int year, int month, int day, double day_frc ) {

    double jd;

    if( month <= 2 ) {

        jd = (int)(365.25 * (year + 4715.0)) +
             (int)(30.6001 * (month + 13.0)) + day - 1537.5 + day_frc;
    }
    else {
        jd = (int)(365.25 * (year + 4716.0)) +
             (int)(30.6001 * (month + 1.0)) + day - 1537.5 + day_frc;
    }

    return jd;
}


/**
   returns the difference between universal time and
   dynamical time (in seconds).

   NOTES:

   universal time lags behind dynamical time. in 1902
   they were roughly the same, but by 1996 universal time was
   lagging by about 62 seonds.

   the difference is due to variations in earth's orbit
   and cannot be predicted precisely because the rotation
   rate is too unpredictable.

   This correction is used only for objects in the solar
   system.

   An error in TD_CORR of 10 seconds of time
   gives a position error of about 0.25 arc seconds in
   Right Ascension for the position of the Sun.

   REFERENCE:

   Astronomical Almanac 1995 or later, p. K9.
*/
double redux::util::time_dyn( int year, int month ) {

    if( month > 9 ) {
        year++;
        month = 0;
    }
    month = (month+3)/6;

    double now = year + 0.5*month;

    dynTime.find( now );
    auto it = dynTime.find( now );
    if( it != dynTime.end() ) {
        return it->second;
    } else {
        return 68.5;    // return value for 2017.01 as default
    }

}


/**
   calculates nutation of the earth's rotational axis

   NOTES:

   nutation is a (relatively) short-period oscillation superimposed
   on the precession of the earth. it is caused by changes in the moon's
   orbit. it has a period of 18.6 years and maximum of about 15 arc
   seconds.

   this function finds:

   - nutation in longitude (nut_long)
   - nutation in obliquity (nut_obl)
   - mean obliquity (mean_obliquity)
   - true obliquity (true_obliquity)

   these values are needed for solar/planet/star positions.

   REFERENCE: Astronomical Algorithms, by Jean Meeus, chapter 21.

*/
void redux::util::nutation_obliquity( double jd, double *nut_long, double *nut_obl,
                                      double *mean_obliquity, double *true_obliquity ) {

    double t, d, m, mp, f, omega;

    t = (jd - 2451545.0) / 36525.0;

    d = (297.85036 + 445267.111480 * t - 0.0019142 *t*t +
         t*t*t / 189474.0) * DtoR;
    m = (357.52772 + 35999.050340 * t - 0.0001603 * t*t -
         t*t*t / 300000.0) * DtoR;
    mp = (134.96298 + 477198.867398 * t + 0.0086972 * t*t +
          t*t*t / 56250.0) * DtoR;
    f = (93.27191 + 483202.017538 * t - 0.0036825 * t*t +
         t*t*t / 327270.0) * DtoR;
    omega = (125.04452 - 1934.136261 * t + 0.0020708 * t*t +
             t*t*t / 450000.0) * DtoR;

    *nut_long = (-17.1996 * sin(omega)
                 - 1.3187 * sin(-2.0*d + 2.0*f + 2.0*omega)
                 - 0.2274 * sin(2.0*f + 2.0*omega)
                 + 0.2062 * sin(2.0*omega)
                 + 0.1426 * sin(m)
                 +  0.0712 * sin(mp)
                 - 0.0517 * sin(-2.0*d + m + 2.0*f + 2.0*omega)
                 - 0.0386 * sin(2.0*f + omega)
                 - 0.0301 * sin(mp + 2.0*f + 2.0*omega))
                / 3600.0 * DtoR;

    *nut_obl = (9.2025 * cos(omega)
                + 0.5736 * cos(-2.0*d + 2.0*f + 2.0*omega)
                + 0.0977 * cos(2.0*f + 2.0*omega)
                - 0.0895 * cos(2.0*omega))
               / 3600.0 * DtoR;

    *mean_obliquity = (84381.448 - 46.8150 * t - 0.00059 * t*t +
                       0.001813 * t*t*t) / 3600.0 * DtoR;

    *true_obliquity = *mean_obliquity + *nut_obl;
}


/**
   computes Greenwich Apparent Sidereal Time at 0 hours

   REFERENCE:

   Astronomical Algorithms, by Jean Meeus, chapter 21.
*/
double redux::util::greenwich_app_sid_time_0( double jd0, double nut_long,
        double true_obliquity ) {
    double tm, gmst0, gast0;

    tm = ( jd0 - 2451545.0 ) / 36525.0;

    gmst0 = ( 100.46061837 + 36000.770053608 * tm ) * DtoR;
    gast0 = gmst0 + nut_long * cos( true_obliquity );
    gast0 = fmod( gast0, 2.0*M_PI );
    if( gast0 < 0.0 ) gast0 = gast0 + 2.0*M_PI;

    return gast0;
}


/**
   gets the following:

   - solar axis position angle (solar_axis)
   - heliographic latitude (stony_lat_cen)
   - longitude for the center of solar disk (stony_long_cen)

   REFERENCE:

   Astronomical Algorithms, by Jean Meeus, chapter 28.
   (Astronomical Almanac 1995, p. C3. Textbook on
   Spherical Astronomy, by Smart, p. 169 - 174).

*/
void redux::util::phys_sun( double jd, double sun_apparent_geo_long,
                            double nut_long, double true_obliquity,
                            double *solar_axis, double *stony_lat_cen ) {

    double incl, selan, longitude;
    
    incl = 7.25 * DtoR;
    selan = ( 73.6667 + 1.3958333 * ( jd + 0.5 - 2396758.0 ) /
              36525.0 ) * DtoR;

    /* get geocentric longitude of the Sun taking into account
       aberration but not nutation */
    longitude = sun_apparent_geo_long - nut_long;

    *solar_axis = atan( - cos( sun_apparent_geo_long ) *
                        tan( true_obliquity ) ) +
                  atan( - cos( longitude - selan ) * tan( incl ) );

    *stony_lat_cen = asin( sin( longitude - selan ) * sin( incl ) );

    //  Heliographic longitude for the center of the solar disk (not used
    //  for the moment. Needed for computing Carrington coordinates):
    //    double theta = ( ( jd - 2398220.0 ) * 360.0/25.38 ) * DtoR;
    //    eta = atan2(- sin(longitude - selan) * cos(incl) ,
    //               - cos(longitude - selan));
    //    stony_long_cen = eta - theta;
    //    stony_long_cen = fmod(stony_long_cen, 2.0 * M_PI);
    //    if (stony_long_cen < 0.0) stony_long_cen =
    //        stony_long_cen + 2.0 * M_PI;
}



/**
   computes Az/El for an object taking into account the current
   Universal Time, observatory location, parallax and atmospheric
   refraction

   REFERENCES:

   Astronomical Algorithms, by Jean Meeus,
   chapter 12, p. 89 (to get azimuth and elevation),
   chapter 39 (to correct for parallax in ra and decl,
   in case of solar system objects).
   (Practical Ephemeris Calculations, by Oliver
   Montenbruck, p. 12 - 13. Telescope Pointing and
   Tracking, p. 5.)

*/
void redux::util::local_position( const Observatory& observatory,
                                  double parallax_object, double gast0, double ut,
                                  double *ha, double *azimuth, double *elevation,
                                  double *ra, double *decl, double *last,
                                  double *refrac ) {

    /* compute hour angle for object */
    *last = gast0 - observatory.longitude + 1.0027379093 * ut;
    *last = fmod( *last, 2.0*M_PI );
    if( *last < 0.0 ) *last = *last + 2.0*M_PI;
    *ha = *last - *ra;

    /* for solar system objects, add parallax correction */
    //Does really strange stuff?
    //correct_for_parallax( observatory, ra, decl, parallax_object, ha );

    *ha = *last - *ra;
    *ha = fmod( *ha, 2.0*M_PI );

    if( *ha < 0.0 ) *ha = *ha + 2.0*M_PI;

    /* compute azimuth */
    *azimuth = atan2( sin( *ha ), cos( *ha ) * sin( observatory.latitude ) -
                      tan( *decl ) * cos( observatory.latitude ) ) + M_PI;

    /* computer elevation */
    *elevation = asin( sin( observatory.latitude ) * sin( *decl ) +
                       cos( observatory.latitude ) * cos( *decl ) *
                       cos( *ha ) );

    /* correct for the effect of atmospheric refraction
       (for elevation < 3 degrees, keep atmospheric refraction constant) */
    if( *elevation < ( 3.0*DtoR ) )
        *refrac = refraction( observatory, 3.0*DtoR );
    else
        *refrac = refraction( observatory, *elevation );

    *elevation = *elevation + *refrac;
}


/**
   computes atmospheric refraction (in radians) for a specified
   zenith distance, taking into account the atmospheric pressure
   and temperature. (For La Palma, height: 2350 m)

   METHOD:

   Using parameter determined at the Carlsberg
   Automatic Transit Circle, La Palma (0.00115)

   REFERENCE:

   Telescope Pointing and Tracking, p. 6.

*/
double redux::util::refraction( const Observatory& observatory, double elevation ) {

    double true_z, refr_par, r1;

    true_z = M_PI/2.0 - elevation;

    /* Formulae not valid for zenith distance > about 1.5 radians */
    if( true_z > 1.5 ) true_z = 1.5;

    refr_par = 60.4 * observatory.air_pressure /
               ( 760.0 * ( 1.0 + observatory.air_temp / 273.0 ) ) / 3600.0 * DtoR;

    r1 = refr_par * ( tan( true_z ) - observatory.refraction_coeff *
                      tan( true_z ) *
                      tan( true_z ) * tan( true_z ) );

    return( refr_par * ( tan( true_z - r1 ) - observatory.refraction_coeff *
                         tan( true_z - r1 ) * tan( true_z - r1 ) *
                         tan( true_z - r1 ) ) );

}



/**
   corrects the equatorial coordinates of a star, taking into account:

   - proper motion (if any) of the star
   - precession of earth's rotational axis
   - nutation of earth's rotational axis
   - annual aberration

   original coords are held in eqCoords->original
   corrected coords are held in eqCoords->corrected
*/
void redux::util::star_geo( double jd, double nut_long, double nut_obl,
                            double mean_obliquity, double sun_apparent_geo_long,
                            double solar_aberration, EqCoords *eqCoords ) {


    double year_diff, zeta, zi, theta;
    double a, b, c, d, t;
    double jd_epoch;
    double corrected_prop_ra;
    double corrected_prop_decl;


    /* calculate julian day of this epoch */
    jd_epoch = 365.25 * ( eqCoords->epoch + 4715.0 ) - 1108.5;

    /* N.B. 2451545.0 was hardcoded where shown below */

    /* difference in years between initial and final epoch of a star */
    year_diff = (jd - /* 2451545.0 */ jd_epoch ) / 365.25;

    /* reduction for proper motion of a star (NB proper motions are
       in milli-arsecs/year. we convert to arcsecs/year)
    */
    corrected_prop_ra =
        eqCoords->prop_ra / 3600000.0 * 15.0 * DtoR * year_diff;

    corrected_prop_decl =
        eqCoords->prop_decl/3600000.0 * DtoR * year_diff;

    eqCoords->corrected.ra =
        eqCoords->original.ra + corrected_prop_ra;

    eqCoords->corrected.decl =
        eqCoords->original.decl + corrected_prop_decl;

    /* correct for precession of the earth's rotational axis */
    t = (jd - /*2451545.0*/jd_epoch) / 36525.0;
    zeta = (2306.2181 * t + 0.30188 * t * t +
            0.017998 * t * t * t) / 3600.0 * DtoR;
    zi = (2306.2181 * t + 1.09468 * t * t +
          0.018203 * t * t * t ) / 3600.0 * DtoR;
    theta = (2004.3109 * t - 0.42665 * t * t -
             0.041833 * t * t * t) / 3600.0 * DtoR;

    a = cos(eqCoords->corrected.decl) * sin(eqCoords->corrected.ra + zeta);
    b = cos(theta) * cos(eqCoords->corrected.decl) *
        cos(eqCoords->corrected.ra + zeta) -
        sin(theta) * sin(eqCoords->corrected.decl);
    c = sin(theta) * cos(eqCoords->corrected.decl) *
        cos(eqCoords->corrected.ra + zeta) +
        cos(theta) * sin(eqCoords->corrected.decl);
    d = atan2(a, b);
    eqCoords->corrected.ra = d + zi;
    eqCoords->corrected.decl = asin(c);

    /* correct for nutation of earth's rotational axis */
    correct_for_nutation( nut_long,
                          nut_obl,
                          mean_obliquity,
                          &eqCoords->corrected.ra, &eqCoords->corrected.decl );

    /* correct for annual aberration */
    correct_for_annual_aberration( t,
                                   nut_long,
                                   mean_obliquity,
                                   sun_apparent_geo_long,
                                   solar_aberration,
                                   &eqCoords->corrected.ra,
                                   &eqCoords->corrected.decl );
}

