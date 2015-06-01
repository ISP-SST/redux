#ifndef REDUX_ALGO_GENETIC_POPULATION_HPP
#define REDUX_ALGO_GENETIC_POPULATION_HPP

#include "redux/algo/genetic/parameterrange.hpp"
#include "redux/algo/genetic/individual.hpp"

#include "redux/util/datautil.hpp"
#include "redux/util/bitoperations.hpp"
// #include "plx/util/random.hpp"

#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <float.h>

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>

namespace redux {

    namespace algo {

        namespace genetic {

            /*!  @ingroup genetic
             *  @{
             */

            /*!  @class     Population
             *   @brief     Container class for a set of individuals in the genetic algorithm.
             *   @author    Tomas Hillberg <hillberg@astro.su.se    >
             *   @date      2011
             *   @todo      Cleanup & document.
             */


            /*! @} */

            template <class T>
            class Genetic;

            struct stat_struct {
                double mean, median, std, min, max;
            };

            template <class T = float>
            class Population {

                Genetic<T>* myGenetic;

                uint popSize;
                uint genomeSize;
                uint nInserted;

                //uint maxGenerations;

                stat_struct* stats;

                Individual<T>** individuals;
                ParameterRange<T>* ranges;

                float crossoverProbability, mutationProbability;
                float *rndArray, *rndPtr;
                int *rndIntArray, *rndIntPtr;
                uint rndSize;
                
                boost::uniform_real<double> uni_real;
                boost::uniform_int<int> uni_int;

                void randomizeGenome ( Individual<T>* );
                void calculateFitness ( int ind = -1 );

                int estimateRandomUsage ( void );

                bool compare ( Individual<T>*, Individual<T>* );

                uint selectParent ( float );
                void crossOver ( Individual<T>*, Individual<T>* );
                void mutate ( Individual<T>* );

                void replaceAll ( void );
                void replaceRandom ( void );
                void replaceWorst ( void );

                void refineCoords ( void );
                void getStats ( void );

            public:
                Population ( Genetic<T>* );
                ~Population();

                void setPopSize ( const uint& value );
                uint getPopSize() const {
                    return popSize;
                }

                void setCrossoverProbability ( float val ) {
                    crossoverProbability = val;
                };
                void setMutationProbability ( float val ) {
                    mutationProbability = val;
                };

                //void setMaxGenerations(const uint& value)   { maxGenerations = value; };
                //uint getMaxGenerations() const { return maxGenerations; }

                void populate ( void );
                void loadRanges ( void );
                void setDynamicRange ( uint n, bool val );
                void setDynamicRanges ( bool val );
                void evacuate ( void );
                void reset ( void );
                void keep ( uint );

                double evolve ( uint& steps ) {
                    return evolveTo ( DBL_MAX, steps );
                };
                double evolveTo ( double goal, uint& steps );
                void sort ( void );

                void modifyPopSize ( int n );

                //      int compareFitness (const void*, const void*);
                //      void setFitness( FitnessFunction f ) { FF = f; };

                void insert ( T* );

                std::string dump ( void ) const;

                template <class U>
                friend class Genetic;
                template <class U>
                friend class Individual;

            };


            template <class T>
            Population<T>::Population ( Genetic<T>* gen ) : individuals ( NULL ),
                                                            ranges ( NULL ),
                                                            rndArray ( NULL ),
                                                            rndIntArray ( NULL ),
                                                            uni_real(0,1),
                                                            uni_int(0,gen->popSize) {
                myGenetic = gen;
                nInserted = 0;
                popSize = gen->popSize;
                genomeSize = gen->genomeSize;
                crossoverProbability = gen->crossoverProbability;
                mutationProbability = gen->mutationProbability;
                stats = new stat_struct[ genomeSize ];
                loadRanges();
                populate();
                // Sort over "entire" array, i.e. the starting population is twice as large.
                redux::util::quickSort ( &individuals[0], popSize << 1, true );
            }


            template <class T>
            Population<T>::~Population() {
                evacuate();

                if ( ranges ) {
                    delete[] ranges;
                }

                if ( stats ) {
                    delete[] stats;
                }

                if ( rndArray ) {
                    delete[] rndArray;
                }

                if ( rndIntArray ) {
                    delete[] rndIntArray;
                }
            }

            template <class T>
            void Population<T>::setPopSize ( const uint& value ) {
                if ( value != popSize ) {
                    modifyPopSize ( value - popSize );
                    popSize = value;
                }
            }

            template <class T>
            void Population<T>::populate() {
                if ( individuals ) {
                    evacuate();
                }

                // enough random numbers for 300 generations initially (estimateRandomUsage() gives ~ 100% margin !!), then reload.
//                 rndSize = estimateRandomUsage() * 150;
// 
//                 if ( rndArray ) {
//                     delete[] rndArray;
//                 }
// 
//                 if ( rndIntArray ) {
//                     delete[] rndIntArray;
//                 }
// 
//                 rndArray = new float[ rndSize ];
//                 rndPtr = &rndArray[ rndSize - 1 ];
//                 rndIntArray = new int[ rndSize ];
//                 rndIntPtr = &rndIntArray[ rndSize - 1 ];
//                 myGenetic->rand.rand ( rndArray, rndSize );
//                 myGenetic->rand.rand ( rndIntArray, rndSize );
                individuals = new Individual<T>* [ popSize << 1 ];

                for ( uint i = 0; i < ( popSize << 1 ); ++i ) {
                    individuals[i] = new Individual<T> ( this, myGenetic );
                    randomizeGenome ( individuals[i] );
                    individuals[i]->getFitness();
                }
            }

            template <class T>
            void Population<T>::evacuate() {
                if ( ! individuals ) {
                    return;
                }

                for ( uint i = 0; i < popSize << 1; ++i ) {
                    delete individuals[i];
                }

                delete[] individuals;
                individuals = NULL;
            }

            template <class T>
            void Population<T>::reset ( void ) {
                loadRanges();

                for ( uint i = 0; i < popSize << 1; ++i ) {
                    randomizeGenome ( individuals[i] );
                    individuals[i]->getFitness();
                }

                redux::util::quickSort ( &individuals[0], popSize << 1, true );
            }

            template <class T>
            void Population<T>::loadRanges ( void ) {
                if ( ranges ) {
                    delete[] ranges;
                }

                ranges = new ParameterRange<T>[ genomeSize ];

                for ( uint i = 0; i < genomeSize; ++i ) {
                    ranges[i] = myGenetic->ranges[i];
                }
            }

            template <class T>
            void Population<T>::setDynamicRange ( uint n, bool val ) {
                ranges[n].setDynamic ( val );
            }

            template <class T>
            void Population<T>::setDynamicRanges ( bool val ) {
                for ( uint i = 0; i < genomeSize; ++i ) {
                    ranges[i].setDynamic ( val );
                }
            }

            template <class T>
            void Population<T>::keep ( uint n ) {
                uint m = std::min ( popSize, n );

                if ( ( rndPtr - rndArray ) < popSize * genomeSize << 1 ) {
                    myGenetic->rand.rand ( rndArray, rndSize );
                    rndPtr = &rndArray[rndSize - 1];
                }

                for ( uint i = 0; i < popSize << 1; ++i ) {
                    if ( i >= m ) {
                        randomizeGenome ( individuals[i] );
                    }

                    individuals[i]->getFitness();
                }

                loadRanges();
                redux::util::quickSort ( &individuals[0], popSize << 1, true );
            }

            template <class T>
            void Population<T>::modifyPopSize ( int n ) {
                if ( individuals ) {
                    Individual<T>** newIndividuals = new Individual<T>* [ ( popSize + n ) << 1 ];

                    for ( uint i = 0; i < ( popSize + n ) << 1; ++i ) {
                        if ( i < popSize << 1 ) {
                            newIndividuals[i] = individuals[i];
                        }
                        else {
                            newIndividuals[i] = new Individual<T> ( this, myGenetic );
                            randomizeGenome ( newIndividuals[i] );
                        }
                    }

                    if ( n < 0 ) {
                        for ( uint i = ( popSize + n ) << 1; i < popSize << 1; ++i ) {
                            delete individuals[i];
                        }
                    }

                    delete individuals;
                    individuals = newIndividuals;
                }

                popSize += n;
            }

            template <class T>
            void Population<T>::randomizeGenome ( Individual<T>* ind ) {
                for ( uint i = 0; i < genomeSize; ++i ) {
                    ind->genome[i] = ranges[i].getRandom ( *rndPtr-- );
                }

                ind->reCalculate = true;
            }

            template <class T>
            void Population<T>::calculateFitness ( int ind ) {
                if ( ind < 0 ) {
                    for ( uint i = 0; i < popSize; i++ ) {
                        individuals[i]->getFitness();
                    }

                    //        sort();
                }
                else {
                    individuals[ind]->getFitness();
                    //        sort();
                }
            }

            template <class T>
            uint Population<T>::selectParent ( float dice ) {
                double* low = myGenetic->breedProb;
                double* high = & ( myGenetic->breedProb[popSize - 1] );
                double* mid;

                if ( dice < *low ) {
                    return 0;
                }

                while ( high - low > 1 ) {
                    mid = low + ( ( high - low ) >> 1 );

                    if ( dice > *mid ) {
                        low = mid;
                    }
                    else {
                        high = mid;
                    }
                }

                return ( uint ) ( high - myGenetic->breedProb );
            }

            template <class T>
            int Population<T>::estimateRandomUsage() {
                double* low = myGenetic->breedProb;
                double* high = & ( myGenetic->breedProb[popSize - 1] );
                double usage = 0.0, tmp = 0.0;
                double inbreedProb = ( *low ) * ( *low );

                while ( low < high ) {
                    inbreedProb += ( * ( low + 1 ) - * ( low ) ) * ( * ( low + 1 ) - * ( low ) );
                    low++;
                }

                usage = genomeSize * ( 1 + ( 1 - myGenetic->crossoverProbability ) * ( 1 - myGenetic->mutationProbability ) );
                tmp = myGenetic->crossoverProbability;

                if ( myGenetic->crossOverType == Genetic<T>::X_2POINT ) {
                    tmp += 1;
                }

                usage += tmp;
                usage += ( 1 + inbreedProb );
                usage *= ( popSize << 1 ); // 100% margin
                return ( int ) usage;
            }

            template <class T>
            bool Population<T>::compare ( Individual<T>* ind1, Individual<T>* ind2 ) {
                //uint nBytes = (sizeof(T));
                //T mask;
                //u_int8_t *ptr,*ptr1,*ptr2;
                //bool ret = true;
                for ( uint i = 0; i < genomeSize; ++i ) {
                    if ( ind1->genome[i] != ind2->genome[i] ) {
                        return false;
                    }

                    /*memset( &mask, 0, sizeof(T) );
                    mask = ranges[i].getDigits( nBytes<<2 );
                    ptr = (u_int8_t*)&mask;
                    ptr1 = (u_int8_t*)&(ind1->genome[i]);
                    ptr2 = (u_int8_t*)&(ind2->genome[i]);
                    for (uint j=0; j<nBytes; ++j) {
                    //(*ptr2) = (*ptr1)>>1;
                    if( *ptr & ( (*ptr1) ^ (*ptr2) ) ) {
                    return false;
                    }
                    ptr++; ptr1++; ptr2++;
                    }*/
                }

                return true;
            }

            template <class T>
            void Population<T>::crossOver ( Individual<T>* ind1, Individual<T>* ind2 ) {
                uint nBits = genomeSize * ( sizeof ( T ) << 3 ); // N.B. left-shift bits by 3 is a faster way of multiplying by 8.
                uint firstPoint = ( *rndIntPtr-- ) % nBits;
                uint secondPoint = ( nBits - 1 );

                if ( myGenetic->crossOverType == Genetic<T>::X_1POINT ) {
                    // One point crossover is already defined, so nothing to do here...
                }
                else if ( myGenetic->crossOverType == Genetic<T>::X_2POINT || *rndPtr-- > 0.5 ) {
                    // if "random" it will use 50% 1-point & 50% 2-point crossover
                    // Select second point, swap them if it is located before the first point.
                    secondPoint = ( *rndIntPtr-- ) % nBits;

                    while ( secondPoint == firstPoint ) {
                        secondPoint = ( *rndIntPtr-- ) % nBits;
                    }

                    if ( secondPoint < firstPoint ) {
                        std::swap ( secondPoint, firstPoint );
                    }
                }

                // size is the the number of (whole or partial) bytes involved in the swap
                int size = ( secondPoint >> 3 ) - ( firstPoint >> 3 ) + 1;
                // Set the pointers to the first byte to be swapped
                u_int8_t *ptr1 = ( ( u_int8_t* ) ind1->genome ) + ( firstPoint >> 3 );
                u_int8_t *ptr2 = ( ( u_int8_t* ) ind2->genome ) + ( firstPoint >> 3 );
                u_int8_t mask;

                if ( firstPoint % 8 ) {                               // Check if there are any individual bits to swap in the first byte.
                    mask = ( ( u_int8_t ) 255 ) >> ( 8 - firstPoint % 8 ); // Set all bits to 1 and then shift in zeroes before "firstPoint"

                    if ( mask & ( *ptr1 ^ *ptr2 ) ) { // does this check save or waste time ??!!
                        redux::util::swapBits ( ptr1, ptr2, mask );              // Swap the remaining bits located before the chunk.
                    }

                    ptr1++;                                             // move pointers to next byte
                    ptr2++;
                    size--;                                             // one byte less to swap below
                }

                if ( ( secondPoint + 1 ) % 8 ) {                      // Check if there is any individual bits to swap in the last byte.
                    mask = ( ( u_int8_t ) 255 ) << ( 8 - secondPoint % 8 ); // Set all bits to 1 and then shift in zeroes after "secondPoint"
                    size--;                                             // one byte less to swap below

                    if ( mask & ( * ( ptr1 + size ) ^ * ( ptr2 + size ) ) ) { // does this check save or waste time ??!!
                        redux::util::swapBits ( ptr1 + size, ptr2 + size, mask );  // Swap the remaining bits located before the chunk.
                    }
                }

                if ( size > 0 ) {                   // if there is a continuous chunk consisting of at least 1 byte, swap it.
                    char* tmp = new char[size];
                    memcpy ( tmp, ptr1, size );
                    memcpy ( ptr1, ptr2, size );
                    memcpy ( ptr2, tmp, size );
                    delete[] tmp;
                }

                // signal that the genome has been changed and we need to re-calculate the fitness for these individuals.
                ind1->reCalculate = true;
                ind2->reCalculate = true;
            }

            template <class T>
            void Population<T>::mutate ( Individual<T>* ind ) {
                //    uint genomeBits = sizeof(T)<<3;
                T mask;
                uint rnd;
                u_int8_t count = 0;

                for ( uint i = 0; i < genomeSize; ++i ) {
                    if ( *rndPtr-- < mutationProbability ) {
                        if ( count == 0 ) {
                            rnd = *rndIntPtr--;
                            count = 4;
                        }

                        count--;
                        mask = ranges[i].getMask ( * ( ( ( u_int8_t* ) &rnd ) + count ) );
                        //uint mask = (((u_int8_t)255) << 24) & rnd;
                        //ranges[i].calculateMask();
                        //T mask2 = ranges[i].getMask();
                        //uint* ptr = (uint*)&(ind->genome[i]);
                        //u_int8_t offset = rnd%(genomeBits-7-ranges[i].precision);
                        //mask >>= (offset+ranges[i].precision);

                        if ( mask > 0 ) {
                            u_int8_t* ptr = ( u_int8_t* ) & ( ind->genome[i] );
                            u_int8_t* maskPtr = ( u_int8_t* ) ( &mask );

                            for ( uint j = 0; j < sizeof ( T ); ++j ) {
                                if ( * ( maskPtr + j ) ) {
                                    * ( ptr + j ) ^= * ( maskPtr + j );
                                }
                            }

                            ind->reCalculate = true;
                        }
                    }
                }
            }

            template <class T>
            void Population<T>::replaceAll() {
                Individual<T> **tmp = new Individual<T>*[ popSize ];           // Temporary storage while moving individuals
                memcpy ( tmp, & ( individuals[popSize] ), popSize * sizeof ( Individual<T>* ) );    // copy last half of array to tmp
                memmove ( & ( individuals[popSize] ), individuals, popSize * sizeof ( Individual<T>* ) ); // move the upper half to the lower half
                memcpy ( individuals, tmp, popSize * sizeof ( Individual<T>* ) );                   // copy tmp to first half of array
                delete[] tmp;
            }

            template <class T>
            void Population<T>::replaceRandom() {  // not yet implemented !!!
                return;
            }

            template <class T>
            void Population<T>::replaceWorst() {
                bool done = false;
                Individual<T>* tmp;                    // Temporary storage while moving individual
                Individual<T>** oldInd;                // Iterator used to locate and point to the insertion-point of new individual.
                Individual<T>** newInd = & ( individuals[popSize] ); // Iterator to step through all individuals to be inserted.

                while ( !done ) {
                    oldInd = & ( individuals[popSize - 1] );
                    bool skip = false;

                    // since they are sorted, we are done if this "newInd" is worse than the last one in the population,
                    // or if we have reached the end of the array.
                    if ( newInd >= & ( individuals[ popSize << 1 ] ) || **newInd < **oldInd ) {
                        done = true;
                        break;
                    }

                    oldInd++;  // just so we don't miss checking against the last individual in the population

                    // find the highest ranking individual with lower fitness that the new one (i.e. where to insert the new one)
                    while ( oldInd > individuals && **newInd >= ** ( oldInd - 1 ) ) {
                        oldInd--;

                        // if some individuals have the same fitness as the new one, step through them
                        if ( ! ( **newInd != **oldInd ) ) { //&& compare( *newInd, *(oldInd-1) ) ) {    // ...and check if the two individuals are identical, i.e. the same genome.
                            skip = true;                              // ...if they are, skip it and move on to next one,
                            //(*newInd)->reCalculate = true;
                            ( *newInd )->randomizeGenome(); // just to avoid identical individuals to be "sorted" into the population
                            ( *newInd )->getFitness();
                            break;
                        }
                    }

                    if ( !skip && ( oldInd < & ( individuals[popSize] ) ) ) { // are we skipping this individual ??
                        tmp = *newInd;                              // ...save the new individual in a temporary place,
                        memmove ( oldInd + 1, oldInd, ( newInd - oldInd ) *sizeof ( Individual<T>* ) ); // push the whole block "down" 1 position
                        *oldInd = tmp;                              // and then insert the new one in the right place.
                    }

                    newInd++;                                     // proceed to the next individual to be inserted.
                }
            }


            template <class T>
            void Population<T>::refineCoords() {
                for ( uint i = 0; i < genomeSize; ++i ) {
                    if ( 20 * stats[i].std < ranges[i].getRange() ) {
                        if ( ranges[i].dynamic() ) {
                            ranges[i].stepMin ( ( stats[i].mean - 3 * stats[i].std - ranges[i].getMin() ) / 5 );
                            ranges[i].stepMax ( ( ranges[i].getMax() - stats[i].mean - 3 * stats[i].std ) / 5 );
                            //ranges[i].calculateMask();
                            //ranges[i].increasePrecision(); //sizeof(T)/2+1);
                        }
                        else {
                            ranges[i].increasePrecision();
                        }

                        ranges[i].increasePrecision( );
                        //cout << colorString(std::to_string(i) + string("  ") + std::to_string((int)ranges[i].precision) + string("\n"),MAGENTA);
                        /*T mymask = ranges[i].getMask();
                          char* ptr1 = ((char*)&mymask);
                          cout << i << " ";
                          for (uint i=0; i<sizeof(T); ++i)
                          cout << colorString(bitString(ptr1[i]) + std::to_string(" "),MAGENTA);
                          cout << alignRight(mymask,20,-10) << endl;*/
                    }
                }
            }


            template <class T>
            void Population<T>::getStats() {

                for ( uint i = 0; i < genomeSize; ++i ) {
                    // calculate mean
                    double sum = 0;
                    stats[i].min = stats[i].max = individuals[0]->genome[i];

                    for ( uint j = 0; j < popSize; j++ ) {
                        sum += individuals[j]->genome[i];

                        if ( individuals[j]->genome[i] > stats[i].max ) {
                            stats[i].max = individuals[j]->genome[i];
                        }
                        else if ( individuals[j]->genome[i] < stats[i].min ) {
                            stats[i].min = individuals[j]->genome[i];
                        }
                    }

                    stats[i].mean = ( sum / popSize );
                    // calculate standard deviation
                    sum = 0;

                    // NOTE: this method fails because of floating-point precision limits.
                    //for (uint j=0; j<popSize; j++) sum += individuals[j]->genome[i]*individuals[j]->genome[i];
                    //sum -= popSize*stats[i].mean*stats[i].mean;
                    for ( uint j = 0; j < popSize; j++ ) {
                        sum += ( individuals[j]->genome[i] - stats[i].mean ) * ( individuals[j]->genome[i] - stats[i].mean );
                    }

                    stats[i].std = ( sqrt ( sum / popSize ) );
                }
            }

            /*! @fn float Population<T>::evolveTo( double goal, int& steps )
             *  @brief Run the optimization for \em at \em most the given number of iterations, but if a fitness value > goal is reached earlier, don't continue.
             *  @details This is the main routine for the genetic algorithm.
             *  @param goal  Fitness limit which will trigger conditional return.
             *  @param steps Number of iterations before returning.
             *  @returns The best fitness-value found.
             */
            template <class T>
            double Population<T>::evolveTo ( double goal, uint& nGen ) {
                if ( !nGen ) {
                    nGen = myGenetic->maxGenerations;
                }

                bool done ( false );
                uint nIter ( 0 );
                uint limit = estimateRandomUsage();
                uint remain;

                for ( uint ii = 0; ii < genomeSize; ++ii ) {
                    ranges[ii].setPrecision ( 0 );
                }

                /*    for ( uint j=0; j<popSize<<1; ++j ) {
                  individuals[j]->getFitness();
                  }*/

                while ( !done ) {
                    // check if we need to reload the random numbers...
                    /*remain = ( rndPtr - rndArray );

                    if ( remain < limit ) {
                        myGenetic->rand.rand ( rndPtr, rndSize - remain );
                        rndPtr = &rndArray[rndSize - 1];
                    }

                    remain = ( rndIntPtr - rndIntArray );

                    if ( remain < limit ) {
                        myGenetic->rand.rand ( rndIntPtr, rndSize - remain );
                        rndIntPtr = &rndIntArray[rndSize - 1];
                    }*/

                    // create "clones" of the parents in the second half of the array...
                    for ( uint child = ( popSize << 1 ) - 1; child >= popSize; --child ) {
                        *individuals[ child ] = *individuals[ selectParent ( *rndPtr-- ) ];
                    }

                    // Then do CrossOver to mix the genome of the parents,
                    for ( uint child = ( popSize << 1 ) - 1; child >= popSize; child -= 2 ) {
                        if ( *rndPtr-- < crossoverProbability ) {
                            while ( ! ( *individuals[ child ] != *individuals[ child - 1 ] ) ) { // make sure they are not identical.
                                *individuals[ child ] = *individuals[ selectParent ( *rndPtr-- ) ];
                            }

                            crossOver ( individuals[child], individuals[child - 1] );
                        }
                    }

                    // ...and finally Random Mutation of the children
                    for ( uint child = ( popSize << 1 ) - 1; child >= popSize; --child ) {
                        mutate ( individuals[ child ] );
                    }

                    uint count ( 0 );

                    // New individuals that are pure copies of the parents are randomized.
                    // ...then we calculate the fitness and the ones that are "better" than
                    // the worst old one are moved to the beginning of the array.
                    for ( uint j = popSize; j < popSize << 1; ++j ) {
                        // bit-swapping floats/doubles or signed integers can cause the variables to jump outside their defined intervals,
                        // so we need to check for that by calling checkGenome(). Randomize unmodified children.
                        if ( individuals[j]->reCalculate ) {
                            individuals[j]->checkGenome();
                        }
                        else {
                            randomizeGenome ( individuals[j] );
                        }

                        individuals[j]->getFitness();

                        if ( individuals[j]->fitness > individuals[ popSize - 1 ]->fitness ) {
                            std::swap ( individuals[j], individuals[ popSize + count ] );
                            count++;
                        }
                    }

                    // Then we do a proper sort of the ones that are good enough and insert them into the population
                    if ( myGenetic->replaceMethod == Genetic<T>::REPLACE_ALL ) {
                        redux::util::quickSort ( &individuals[popSize], popSize, true );
                        replaceAll();
                    }
                    else if ( myGenetic->replaceMethod == Genetic<T>::REPLACE_RANDOM && count > 0 ) {
                        redux::util::quickSort ( &individuals[popSize], count, true );
                        replaceAll();
                    }
                    else if ( myGenetic->replaceMethod == Genetic<T>::REPLACE_WORST && count > 0 ) {
                        redux::util::quickSort ( &individuals[popSize], count, true );
                        replaceWorst();
                    }

                    //float bla[ popSize ];
                    //for( int ii=0; ii<popSize; ++ii) {
                    //if (individuals[ii+1]->fitness > individuals[ii]->fitness)
                    //  std::cout << "fitness decrease !!!!" << std::endl;
                    //bla[ii] = individuals[ii]->fitness;
                    //}
                    //std::cout << dumpArray(bla,popSize,"fitness") << std::endl << "plot,fitness" << std::endl;
                    // getStats() calculates min,max, mean & standard deviation. Which we can use to refine
                    // coordinate limits and mutation mask.
                    nIter++;

                    if ( nIter % myGenetic->statFrequency == 0 ) {
                        getStats();
                        refineCoords();
                    }

                    if ( nIter >= nGen || nIter >= myGenetic->maxGenerations ||
                            ( individuals[0]->fitness > goal && nIter >= myGenetic->minGenerations ) ) {
                        done = true;
                    }
                } // end while

                nGen = nIter;
                return individuals[0]->fitness;
            }

            template <class T>
            void Population<T>::sort ( void ) {
                redux::util::quickSort ( &individuals[0], popSize, true );
                redux::util::quickSort ( &individuals[popSize], popSize, true );
                replaceWorst();
                nInserted = 0;
            }

            template <class T>
            void Population<T>::insert ( T* genome ) {
                if ( nInserted > popSize ) {
                    return;
                }

                nInserted++;
                individuals[ ( popSize << 1 ) - nInserted]->setGenome ( genome );
                //individuals[(popSize<<1)-nInserted]->trimParams();
                individuals[ ( popSize << 1 ) - nInserted]->getFitness();
                sort();
            }

            template <class T>
            std::string Population<T>::dump ( void ) const {
                std::string tmp = "\nPopulation: (" + redux::util::hexString( this ) + ")\n";
                tmp += "     Population Size:              " + std::to_string(popSize) + "\n";
                tmp += "     CrossOver Probability:        " + std::to_string(crossoverProbability) + "\n";
                tmp += "     Mutation Rate:                " + std::to_string(mutationProbability) + "\n";
                tmp += "     Dynamic Mutation Rate ?       ";

                if ( false ) {
                    tmp += "YES  ";
                }
                else {
                    tmp += "NO   ";
                }

                tmp += "\n     Parameters:          Min         Max       Dynamic?   Precision          Mean         Standard deviation\n";

                for ( uint i = 0; i < genomeSize; ++i ) {
                    tmp += "    Parameter    " + std::to_string ( i );
                    tmp += std::to_string( ranges[i].getMin() ) + std::to_string(ranges[i].getMax());

                    if ( ranges[i].dynamic() ) {
                        tmp += "    YES     ";
                    }
                    else {
                        tmp += "     NO     ";
                    }

                    tmp += std::to_string((int)ranges[i].getPrecision()) + std::to_string(stats[i].mean);
                    tmp += std::to_string(stats[i].std / ranges[i].getRange()) + "\n";
                }

                tmp += "\n";
                tmp += "     " +  individuals[0]->dump();
                tmp += "     " +  individuals[1]->dump();
                tmp += "     " +  individuals[10]->dump();
                /*for( uint i=0; i<2; ++i ) { //popSize; ++i) {
                  tmp += string( 5, ' ' ) +  individuals[i]->dump();
                }

                for( uint i=popSize-1; i<popSize; ++i ) { //popSize; ++i) {
                  tmp += string( 5 ) + individuals[i]->dump();
                }*/
                return tmp;
            }
        }       // genetic
    }       // algo
}       // pallax

#endif  // REDUX_ALGO_GENETIC_POPULATION_HPP
