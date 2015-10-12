#ifndef REDUX_ALGO_GENETIC_INDIVIDUAL_HPP
#define REDUX_ALGO_GENETIC_INDIVIDUAL_HPP

#include <string>

#include "redux/util/stringutil.hpp"

namespace redux {

    namespace algo {

        namespace genetic {

            /*!  @ingroup genetic
             *  @{
             */

            /*!  @class     Individual
             *   @brief     Container class for an "individual" in the genetic algorithm.
             *   @author    Tomas Hillberg <hillberg@astro.su.se    >
             *   @date      2011
             *   @todo      Cleanup & document.
             */


            /*! @} */


            template <class T>
            class Genetic;

            template <class T>
            class Population;

            template <class T = float>
            class Individual {

                Population<T>* myPopulation;
                Genetic<T>* myGenetic;

                Individual<T> *fittestParent;

                T* genome;

                bool reCalculate;

                void setGenome ( T* );
                void checkGenome ( void );
                void randomizeGenome ( void );
                void getFitness ( void );

            public:
                Individual ( Population<T>*, Genetic<T>* );
                Individual ( const Individual& );
                ~Individual();

                double fitness;

                bool operator< ( const Individual& ) const;
                bool operator> ( const Individual& ) const;
                bool operator<= ( const Individual& ) const;
                bool operator>= ( const Individual& ) const;
                bool operator== ( const Individual& ) const;
                bool operator!= ( const Individual& ) const;
                Individual<T>& operator= ( const Individual<T>& rhs );

                //      inline friend bool operator<=(const Individual& lhs, const Individual* rhs)

                template <class U>
                friend class Population;

                template <class U>
                friend class Genetic;

                std::string dump ( void ) const;

            };

            template <class T>
            Individual<T>::Individual ( Population<T>* pop, Genetic<T>* gen ) {
                myPopulation  = pop;
                myGenetic     = gen;
                fittestParent = NULL;
                fitness = 0.0;
                genome = new T [ pop->genomeSize ];
                reCalculate = true;
                //randomizeGenome();
                //getFitness();
            }

            template <class T>
            Individual<T>::Individual ( const Individual& rhs ) {
                myPopulation  = rhs.myPopulation;
                myGenetic     = rhs.myGenetic;
                fittestParent = rhs.fittestParent;
                fitness = rhs.fitness;
                genome = new T [ myPopulation->genomeSize ];
                memcpy ( genome, rhs.genome, sizeof ( T ) *myPopulation->genomeSize );
                reCalculate = rhs.reCalculate;
                //randomizeGenome();
                //getFitness();
            }

            template <class T>
            Individual<T>::~Individual() {
                if ( genome ) {
                    delete[] genome;
                }
            }

            template <class T>
            void Individual<T>::setGenome ( T* gen ) {
                for ( unsigned int i = 0; i < myPopulation->genomeSize; ++i ) {
                    genome[i] = gen[i];
                    //if ( ! myPopulation->ranges[i].checkParam( &(genome[i]) ) )
                    //genome[i] = myPopulation->ranges[i].getRandom(*(myPopulation->rndPtr--));
                    //myPopulation->ranges[i].trimParam( &(genome[i]) );
                }
            }

            template <class T>
            void Individual<T>::checkGenome ( void ) {
                for ( unsigned int i = 0; i < myPopulation->genomeSize; ++i ) {
                    if ( ! myPopulation->ranges[i].checkParam ( & ( genome[i] ) ) ) {
                        genome[i] = myPopulation->ranges[i].getRandom ( * ( myPopulation->rndPtr-- ) );
                    }

                    //myPopulation->ranges[i].trimParam( &(genome[i]) );
                }
            }

            template <class T>
            void Individual<T>::randomizeGenome() {
               // boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
                for ( unsigned int i = 0; i < myPopulation->genomeSize; ++i ) {
               //     genome[i] = myPopulation->ranges[i].getRandom ( myGenetic->rand() );
                }
            }

            template <class T>
            void Individual<T>::getFitness() {
                fitness = myGenetic->FF ( myPopulation->genomeSize, genome, myGenetic->fitnessArg );
                //fitness = rand()/RAND_MAX;
                reCalculate = false;
            }


            template <class T>
            bool Individual<T>::operator< ( const Individual& rhs ) const {
                return ( fitness < rhs.fitness );
            }

            template <class T>
            bool Individual<T>::operator> ( const Individual& rhs ) const {
                return ( fitness > rhs.fitness );
            }

            template <class T>
            bool Individual<T>::operator<= ( const Individual& rhs ) const {
                //if (!(rhs.fitness == rhs.fitness)) return false;
                //else if (!(this->fitness == this->fitness)) return true;
                return ( fitness <= rhs.fitness );
            }

            template <class T>
            bool Individual<T>::operator>= ( const Individual& rhs ) const {
                //if (!(rhs.fitness == rhs.fitness)) return true;
                //else if (!(this->fitness == this->fitness)) return false;
                return ( fitness >= rhs.fitness );
            }

            template <class T>
            bool Individual<T>::operator== ( const Individual& rhs ) const {
                return ( fitness == rhs.fitness );
                /*   if ( this->fitness != rhs.fitness ) return false;
                return true;
                for ( unsigned int i=0; i<myPopulation->genomeSize; ++i )
                  if ( genome[i] != rhs.genome[i] ) return false;

                return true;*/
            }

            template <class T>
            bool Individual<T>::operator!= ( const Individual& rhs ) const {
                if ( this->fitness != rhs.fitness ) {
                    return true;
                }

                //return false;
                for ( unsigned int i = 0; i < myPopulation->genomeSize; ++i )
                    if ( genome[i] != rhs.genome[i] ) {
                        return true;
                    }

                return false;
            }

            template <class T>
            Individual<T>& Individual<T>::operator= ( const Individual<T> &rhs ) {
                // Check for self-assignment!
                if ( this == &rhs ) {
                    return *this;
                }

                //int chunkSize = sizeof(genome[0])*myPopulation->genomeSize;
                memcpy ( genome, rhs.genome, sizeof ( T ) *myPopulation->genomeSize );
                fitness = rhs.fitness;
                reCalculate = rhs.reCalculate;
                myPopulation = rhs.myPopulation;
                myGenetic = rhs.myGenetic;
                return *this;
            }

            template <class T>
            std::string Individual<T>::dump ( void ) const {
                std::string tmp = "Individual: ";
                tmp += "   Fitness: " + std::to_string(fitness) + "\n";
                tmp += "            Genome: \n";

                if ( myPopulation ) {
                    for ( unsigned int i = 0; i < myPopulation->genomeSize; ++i ) {
                        tmp += std::to_string(genome[i]) + "    " + redux::util::bitString ( genome[i] ) + "\n";
                    }
                }

                return tmp;
            }
        }       // genetic
    }       // algo
}       // pallax

#endif          // REDUX_ALGO_GENETIC_INDIVIDUAL_HPP
