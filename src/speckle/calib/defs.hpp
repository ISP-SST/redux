#ifndef REDUX_SPECKLE_CALIB_DEFS_HPP
#define REDUX_SPECKLE_CALIB_DEFS_HPP

#ifndef SPECKLE_EPS
  #define SPECKLE_EPS       ((double)(1.0e-6))    // a small number
#endif
#ifndef SPECKLE_EPS2
  #define SPECKLE_EPS2      ((double)(1.0e-12))   // a smaller number
#endif
#ifndef SPECKLE_IMAX
  #define SPECKLE_IMAX      231                   // truncation index for inf sum (Noll)
#endif
#ifndef SPECKLE_FREQSTEPS
  #define SPECKLE_FREQSTEPS 100                   // stepwidth in frequency
#endif
#ifndef SPECKLE_NOGO
  #define SPECKLE_NOGO      SPECKLE_FREQSTEPS     // NOGO tag
#endif


#endif  // REDUX_SPECKLE_CALIB_DEFS_HPP

