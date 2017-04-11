#ifndef REDUX_UTIL_SOLAR_SYSTEM_GEO_HPP
#define REDUX_UTIL_SOLAR_SYSTEM_GEO_HPP


namespace redux {
    
    namespace util {
                
        //void comet_pos (double jde, double *ra, double *decl);


        /**
         * PURPOSE: To compute the following parameters for the Sun
         *          at a specific dynamical time (jde) instant:
         *
         *          a: The Apparent Geocentric Longitude
         *             (sun_apparent_geo_long) and Latitude
         *             (sun_apparent_geo_lat) in radians.
         *          b: The Apparent Geocentric Right Ascension (ra)
         *             and Declination (decl) in radians.
         *          c: Radius Vector in Astronomical Units (radius_vec).
         *
         * METHOD: Using the most important terms from the VSOP87
         *         planetary theory by Bretagnon & Francou.
         *
         * REFERENCE: Astronomical Algorithms, by Jean Meuus, chapter 24.
         *            (For a theoretical explanation, see Astronomy and 
         *            Astrophysics, 202, p. 309 - 315, 1988.
         */
        void sun_geo (double jde, double nut_long, double solar_aberration,
                             double true_obliquity,
                             double *sun_apparent_geo_long, double *sun_apparent_geo_lat,
                             double *radius_vec, double *ra, double *decl);


        /**
         * calculate the apparent Right Ascension and Declination
         * of the moon
         */
        void moon_geo (double jde, double nut_long, double true_obliquity,
                              double *dist_moon_earth, double *ra, double *decl);

        /**
         * computes the apparent geocentric position for the planets
         */
        void planets_geo (int planet, double jde, double nut_long,
                                 double solar_aberration, double true_obliquity,
                                 double *dist_planet_earth, double *ra, double *decl);


    }   // util
    
}   // redux

#endif
