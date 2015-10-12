#ifndef REDUX_ALGO_GENETIC_PARAMETERRANGE_HPP
#define REDUX_ALGO_GENETIC_PARAMETERRANGE_HPP

#include "redux/util/boundvalue.hpp"

#include <algorithm>
#include <cstring>
#include <stdio.h>

namespace redux {

    namespace algo {

        namespace genetic {

            /*!  @ingroup genetic
             *  @{
             */

            /*!  @class     ParameterRange
             *   @brief     Defining the limits for the parameters the genetic algorithm is fitting.
             *   @author    Tomas Hillberg <hillberg@astro.su.se>
             *   @date      2011
             *   @todo      Cleanup & document.
             */

            /*! @} */

            template <class T = float>
            class ParameterRange : public redux::util::BoundValue<T> {

                T mask;
                uint8_t diff, precision;
                bool isDynamic;

            public:

                ParameterRange(T min=0, T max=1, T val=0) : redux::util::BoundValue<T>(val,min,max), diff(0), precision(0), isDynamic(false) {
                    calculateDiff();
                    calculateMask();
                }
                
                void setRange ( T min, T max ) {
                    this->setLimits(min,max);
                    calculateDiff();
                }
                
                void stepMin ( T val ) {
                    this->setMin(this->minVal_ + val);
                    calculateDiff();
                }

                void stepMax ( T val ) {
                    this->setMax(this->maxVal_ + val);
                    calculateDiff();
                }
                

                void setDynamic ( bool d ) {
                    isDynamic = d;
                };
                bool dynamic ( void ) const {
                    return isDynamic;
                };

                void setPeriodic ( bool val ) {
                    this->trim_ = this->selectTrimFunction ( redux::util::detail::WRAP );
                }

                void setPrecision ( unsigned char val ) {
                    if ( val <= ( ( sizeof ( T ) - 1 ) << 3 ) ) {
                        precision = val;
                    }
                };
                unsigned char getPrecision() {
                    return diff + precision;
                };
                void increasePrecision ( unsigned int val = 1 ) {
                    if ( precision + val < ( ( sizeof ( T ) - 1 ) << 3 ) ) {
                        precision += val;
                    }
                };
                void decreasePrecision ( unsigned int val = 1 ) {
                    if ( precision - val >= 0 ) {
                        precision -= val;
                    }
                };

                T getMask ( void ) const {
                    return mask;
                };
                T getMask ( uint8_t m ) {
                    calculateMask ( m );
                    return mask;
                };
                T getDigits ( uint8_t minBits = 0 ) {
                    T tmp;
                    unsigned int maxN = ( sizeof ( T ) << 3 );

                    if ( diff + precision > maxN ) {
                        precision = maxN - diff;
                    }

                    //cout << " maxN + diff + prec + m = " << maxN << "  " << (int)diff << "  " << (int)precision << "  " << (int)m << endl;
                    unsigned int n = std::max ( (uint8_t)( diff + precision ), minBits );
                    uint8_t* ptr = ( uint8_t* ) &tmp;
                    memset ( ptr, 0, sizeof ( T ) );
#if CPU_ARCHITECTURE == BIG_ENDIAN

                    while ( n >= 8 ) {
                        *ptr = 255;
                        ptr++;
                        n -= 8;
                    }

                    if ( n > 0 ) {
                        *ptr = ( ( uint8_t ) 255 << ( 8 - n ) );
                    }

#else
                    ptr += ( sizeof ( T ) - 1 );

                    while ( n >= 8 ) {
                        *ptr = 255;
                        ptr--;
                        n -= 8;
                    }

                    if ( n > 0 ) {
                        *ptr = ( ( uint8_t ) 255 << ( 8 - n ) );
                    }

#endif
                    /*    if ( __BYTE_ORDER == __LITTLE_ENDIAN ) {     // if little endian
                    ptr += (sizeof(T)-1);
                    // *ptr = 255;
                    //ptr--;
                    //cout << colorString(printBits(mask), RED) << endl;
                    while ( n >= 8 ) {
                      *ptr = 255;
                      ptr--;
                      n -= 8;
                    }

                    if ( n > 0 ) {
                      *ptr = ((uint8_t)255 << (8-n));
                    }
                    //cout << "precision = " << (int)precision << "   diff = " << (int)diff << endl;
                    //cout << colorString(printBits(mask), GREEN) << endl;
                        } else {                       // big endian
                    // *ptr = 255;
                    //ptr++;
                    //cout << colorString(printBits(mask),RED) << endl;
                    while ( n >= 8 ) {
                      *ptr = 255;
                      ptr++;
                      n -= 8;
                    }
                    if ( n > 0 ) {
                      *ptr = ((uint8_t)255 << (8-n));
                    }
                        }
                    */
                    return tmp;
                };

                void calculateMask ( uint8_t m = 255 ) {
                    //calculateDiff();
                    //diff = 21;
                    //precision = 11;
                    unsigned int maxN = ( ( sizeof ( T ) - 1 ) << 3 );

                    //if ( diff > maxN ) diff = maxN;
                    if ( diff + precision > maxN ) {
                        precision = maxN - diff;
                    }

                    //cout << " maxN + diff + prec = " << maxN << "  " << (int)diff << "  " << (int)precision << endl;
                    unsigned int n = ( unsigned int ) ( ( diff + precision ) >> 1 );
                    n += ( ( unsigned int ) m ) % ( maxN - n + 8 );
                    //if (qqq >= maxN-2) cout << "bla    " << qqq << endl;
                    //unsigned int n = ((diff + precision) >> 1);// + m % ((diff + precision) >> 1); //(maxN >> 1);  // testing to semi-randomize the shift...
                    uint8_t* ptr = ( uint8_t* ) &mask;
                    memset ( ptr, 0, sizeof ( T ) );
#if CPU_ARCHITECTURE == BIG_ENDIAN
                    *ptr = m;

                    while ( n >= 8 ) {
                        ( * ( ( unsigned short* ) ptr++ ) ) >>= 8;
                        n -= 8;
                    }

                    if ( n ) {
                        if ( ( ptr - ( uint8_t* ) &mask ) > ( sizeof ( T ) - 1 ) ) {
                            ( * ( ( unsigned short* ) ptr ) ) >>= n;
                        }
                        else if ( ( ptr - ( uint8_t* ) &mask ) == ( sizeof ( T ) - 1 ) ) {
                            *ptr >>= n;
                        }
                    }

#else
                    ptr += ( sizeof ( T ) - 1 );
                    *ptr = m;

                    while ( n >= 8 ) {
                        ( * ( ( unsigned short* )--ptr ) ) >>= 8;
                        n -= 8;
                    }

                    if ( n > 0 ) {
                        if ( ( ptr - ( uint8_t* ) &mask ) > 0 ) {
                            ( * ( ( unsigned short* )--ptr ) ) >>= n;
                        }
                        else if ( ( ptr - ( uint8_t* ) &mask ) == 0 ) {
                            *ptr >>= n;
                        }
                    }

#endif
                    /*    if ( __BYTE_ORDER == __LITTLE_ENDIAN ) {     // if little endian
                    ptr += (sizeof(T)-1);
                    *ptr = m;
                    //cout << colorString(printBits(mask), RED) << endl;
                    while ( n >= 8 ) {
                      (*((unsigned short*)--ptr)) >>= 8;
                      n -= 8;
                    }
                    if ( n > 0 ) {
                      if( (ptr - (uint8_t*)&mask) > 0 ) (*((unsigned short*)--ptr)) >>= n;
                      else if( (ptr - (uint8_t*)&mask) == 0 ) *ptr >>= n;
                    }
                    //cout << "precision = " << (int)precision << "   diff = " << (int)diff << endl;
                    //cout << colorString(printBits(mask), GREEN) << endl;
                        } else {                       // big endian
                    *ptr = m;
                    //cout << colorString(printBits(mask),RED) << endl;
                    while ( n >= 8 ) {
                      (*((unsigned short*)ptr++)) >>= 8;
                      n -= 8;
                    }
                    if ( n ) {
                      if( (ptr - (uint8_t*)&mask) > (sizeof(T)-1) ) (*((unsigned short*)ptr)) >>= n;
                      else if( (ptr - (uint8_t*)&mask) == (sizeof(T)-1) ) *ptr >>= n;
                    }
                        }
                        */
                }

                void calculateDiff ( void ) {
                    //unsigned int i = 0;
                    diff = 0;
//                     uint8_t* minPtr = ( ( uint8_t* ) &minValue ) + ( sizeof ( T ) - 1 );
//                     uint8_t* maxPtr = ( ( uint8_t* ) &maxValue ) + ( sizeof ( T ) - 1 );
// 
//                     while ( minPtr > ( uint8_t* ) &minValue && *minPtr == *maxPtr ) {
//                         diff += 8;
//                         minPtr--;
//                         maxPtr--;
//                         i++;
//                     }
// 
//                     if ( minPtr > ( uint8_t* ) &minValue ) {
//                         for ( unsigned int j = 0; j < 8; ++j ) {
//                             if ( ( *minPtr ^ *maxPtr ) & ( 1 << ( 7 - j ) ) ) {
//                                 break;
//                             }
//                             else {
//                                 diff++;
//                             }
//                         }
//                     }
                }

                /*void trimParam( T* val ) {
                  if ( !(*val <= maxValue) || !(*val >= minValue) ) {
                T bla = *val;
                T cnt = floor((*val-minValue)/(maxValue-minValue));

                if (ABS(cnt) < 10) {
                  if (isPeriodic) {
                    //cout << "periodic (" << minValue << "," << maxValue << ") and adjusting:  old=" << *val;
                    //float integral_part;
                    //float fractional_part = modf((*val-minValue)/(maxValue-minValue), &integral_part);
                    // *val = *val - (maxValue-minValue)*integral_part;
                    *val = *val - (maxValue-minValue)*cnt;
                    //cout << "  new=" << *val << endl;
                  } else *val = getRandom();

                  if ( *val > maxValue || *val < minValue ) {
                    cout << "failed fixing " << *val << " (" << bla << "," << cnt << ")";
                    cout << "MIN: " << minValue << "  MAX: " << maxValue << "   -> trimming ";

                    *val = getRandom(); //maxValue;

                    cout << "  new random value = " << *val << endl;
                    //else if (*val < minValue) *val = minValue;
                    //cout << *val << "   cnt = " << endl;
                    //cout << "  new2=" << *val << endl;
                  }

                } else {
                  *val = getRandom();
                  return;
                }
                }*/


                bool checkParam ( T* val ) {


                    return false;
                };


                ParameterRange& operator= ( const ParameterRange &rhs ) {
                    // Check for self-assignment!
                    if ( this == &rhs ) {
                        return *this;
                    }

                    redux::util::BoundValue<T>::operator= ( reinterpret_cast<const redux::util::BoundValue<T>&>(rhs) );

                    mask = rhs.mask;
                    diff = rhs.diff;
                    precision = rhs.precision;
                    isDynamic = rhs.isDynamic;

                    return *this;
                };


                template <class U>
                friend class Population;
                template <class U>
                friend class Individual;

            };        // ParameterRange

        }       // genetic
    }       // algo
}       // pallax

#endif      // REDUX_ALGO_GENETIC_PARAMETERRANGE_HPP
