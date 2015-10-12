#ifndef REDUX_ALGO_GENETIC_HPP
#define REDUX_ALGO_GENETIC_HPP

#include <boost/random/mersenne_twister.hpp>

#include "redux/algo/genetic/population.hpp"
#include "redux/algo/genetic/parameterrange.hpp"

namespace redux {

    namespace algo {

        /*! @defgroup genetic Genetic Algorithm
         *  @-page GeneticAlgorithm
         *  @{
         */

        /*!  @class     Genetic
         *   @brief     Class providing a genetic fitting algorithm.
         *   @details
         *   @code
         *   Usage: Genetic [-?] [-h] [-debug] InputFile
         *   @endcode
         *   <table>
         *   <tr> <td>   -? or -h       </td> <td> Show a little help text explaining the command line options.                                                                                                                                                                                                                     </td> </tr>
         *   <tr> <td>   -v or -verbose </td> <td> Switch on the debug mode. With -verbose=0xFF you will get all messages. </td></tr>
         *   <tr> <td>   -bla           </td> <td> Name of input file.                                                                                                                                                                                                                                                     </td> </tr>
         *   </table>
         *   @author    Tomas Hillberg <hillberg@astro.su.se>
         *   @date      2011
         *   @todo  Properly implement support for multiple populations.
         *   @todo      Cleanup & document.
         */


        /*! @} */

        class GeneticBase {

        protected:
            
            typedef double ( *BreedingPDF ) ( uint, unsigned int );

            enum crossover_t { X_RAND = 0, X_1POINT, X_2POINT };
            enum replace_t { REPLACE_ALL = 0, REPLACE_RANDOM, REPLACE_WORST };

            crossover_t crossOverType;
            replace_t replaceMethod;

            unsigned int maxGenerations, minGenerations, nGenerations;
            unsigned int popSize;
            unsigned int nParameters;
            unsigned int statFrequency;
            unsigned int totalEffort;

            float crossoverProbability, mutationProbability;

            BreedingPDF bPDF;
            std::unique_ptr<float[]> breedingProbability;

            bool allowPopFork, allowCrossBreed;

            //plx::Random rand;
            //boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
            typedef boost::random::mt19937 rng_type;
            rng_type rng;


            void calcProb ( void );

        public:
            GeneticBase ( unsigned int nParams, unsigned int popSize );
            
            void setPopSize ( const size_t& n );

            /*! @fn void setBreedingPDF( BreedingPDF f )
             *  @brief Specify the breeding probability distribution to be used.
             *  @details The default bPDF is a simple linear weighting according to fitness.
             *  @param f Breeding PDF
             */
            void setBreedingPDF ( BreedingPDF f ) { bPDF = f; calcProb(); };
            
            void setStatFrequency ( const uint& value ) { statFrequency = value; };
            unsigned int getStatFrequency() const { return statFrequency; }

            void setCrossoverProbability ( float val ) { crossoverProbability = val; }
            void setMutationProbability ( float val ) { mutationProbability = val; }


        };


        template <class T = float>
        class Genetic : public GeneticBase {

            struct Individual {

                T* genome;
                bool reCalculate;
                /*
                                    void setGenome ( T* );
                                    void checkGenome ( void );
                                    void randomizeGenome ( void );
                                    void getFitness ( void );

                                    Individual ( Population<T>*, Genetic<T>* );
                                    Individual ( const Individual& );
                                    ~Individual();*/

                double fitness;

            };

            typedef double ( *FitnessFunction ) ( uint, T*, void* );
            FitnessFunction FF;

            void* fitnessArg;

            std::vector<Individual> population;
            std::vector<genetic::ParameterRange<T>> ranges;


        public:

            Genetic ( unsigned int nVars = 2, unsigned int nInd = 100 );

            void calcFitness ( void );
            double evolve ( uint& steps );
            double evolveTo ( double goal, uint& steps );
            T* getFittest ( void );
            T* getIndividual ( unsigned int );
            void setParameterRange ( unsigned int n, T low, T high, bool per = false );
            void setPeriodicParameter ( unsigned int n, bool val = true );
            void setDynamicRange ( unsigned int n, bool val );
            void setDynamicRanges ( bool val );

            void clear ( void );
            void reset ( void );
            /*! @fn resetCounter( void )
             *  @brief Reset the iteration counter.
             */
            void resetCounter ( void ) {
                totalEffort = 0;
            };
            void sort ( void );
            void keep ( unsigned int );

            /*! @fn void setFitnessFunction( FitnessFunction f, void* arg=NULL )
             *  @brief Specify the fitness function to be used.
             *  @param f Fitness function.
             *  @param arg Additional argument passed to the fitness function upon evaluation.
             */
            void setFitnessFunction ( FitnessFunction f, void* arg = NULL ) {
                FF = f;
                fitnessArg = arg;
            }; //calcFitness(); };
            void setFitnessArg ( void* ar ) {
                fitnessArg = ar;
            }
            void* getFitnessArg() {
                return fitnessArg;
            }

            std::string dump ( void ) const;


        };

        /*! @fn double defaultMerit( unsigned int n, T* data, void* arg )
         *  @brief Default merit-function. A simple example, fitting each of the parameters to the value 0.5.
         */
        template <class T>
        double defaultMerit ( unsigned int n, T* data, void* arg ) {
            double chisq = 0;

            for ( unsigned int i = 0; i < n; ++i ) {
                chisq += ( data[i] - 0.5 ) * ( data[i] - 0.5 );
            }

            return ( double ) 1.0 / chisq;
        }

        /*! @fn defaultBreedingPDF(unsigned int n, unsigned int rank)
         *  @brief Default breeding pdf.
         */
        inline double defaultBreedingPDF ( unsigned int n, unsigned int rank ) {
            double probability ( 0 );

            for ( unsigned int i = n; i >= rank; --i ) {
                probability += ( n - i + 1 );
            }

            return probability;
        }

        template <class T>
        Genetic<T>::Genetic ( unsigned int nParameters, unsigned int popSize ) : GeneticBase ( nParameters,popSize ) {
            //MTRandom::uint32 seed[ MTRandom::N ];
            //for (unsigned int i=0; i<MTRandom::N; ++i) seed[i] = 127*i;
            //rand.seed( (long)this );
            //rand.seed ( ::rand() );
            ranges.resize ( nParameters );
//                 setFitnessFunction ( defaultMerit );
//                 setBreedingPDF ( defaultBreedingPDF );
//                 init();
//                 calcProb();
//                 calcFitness();
        }


        /*! @fn void Genetic<T>::clear( void )
         *  @brief Delete the population(s).
         */
        template <class T>
        void Genetic<T>::clear ( void ) {
//                 if ( !populations ) {
//                     return;
//                 }
//
//                 for ( unsigned int i = 0; i < nPopulations; ++i ) {
//                     delete populations[i];
//                 }
//
//                 delete[] populations;
//                 populations = NULL;
        }

        /*! @fn void Genetic<T>::reset( void )
         *  @brief Reset the population(s).
         */
        template <class T>
        void Genetic<T>::reset ( void ) {
//                 for ( unsigned int i = 0; i < nPopulations; ++i ) {
//                     populations[i]->reset();
//                 }

            totalEffort = 0;
        }


        /*! @fn void Genetic<T>::keep( unsigned int n )
         *  @brief Keep the best n individuals for next run, randomize the rest.
         */
        template <class T>
        void Genetic<T>::keep ( unsigned int n ) {
//                 for ( unsigned int i = 0; i < nPopulations; ++i ) {
//                     populations[i]->keep ( n );
//                 }
        }

        /*! @fn void Genetic<T>::sort( void )
         *  @brief Sort the population(s).
         */
        template <class T>
        void Genetic<T>::sort ( void ) {
//                 for ( unsigned int i = 0; i < nPopulations; ++i ) {
//                     populations[i]->sort();
//                 }
        }


        /*! @fn void Genetic<T>::calcFitness( void )
         *  @brief Calculate the fitness for the population(s).
         */
        template <class T>
        void Genetic<T>::calcFitness ( void ) {
//                 for ( unsigned int i = 0; i < nPopulations; ++i ) {
//                     populations[i]->calculateFitness();
//                 }
        }

        /*! @fn double Genetic<T>::evolve( int& steps )
         *  @brief Run the genetic algorithm for the given number of iterations.
         *  @details This is the main routine for the genetic algorithm. It will call "evolve()" for each of the populations.
         *  @param steps Number of iterations before exiting
         *  @returns The best fitness-value found.
         */
        template <class T>
        double Genetic<T>::evolve ( uint& steps ) {
            if ( steps == 0 || steps > maxGenerations ) {
                steps = maxGenerations;
            }

            //unsigned int thisEffort = steps;
            steps = 0;
            double fn = 0.0;//, tmp = 0.0;

//                 for ( unsigned int j = 0; j < nPopulations; ++j ) {
//                     tmp = populations[j]->evolve ( thisEffort );
//
//                     if ( tmp > fn ) {
//                         fn = tmp;
//                     }
//
//                     steps += thisEffort;
//                 }

            totalEffort += steps;
            return fn;
        }

        /*! @fn double Genetic<T>::evolveTo( double goal, int& steps )
         *  @brief Run the genetic algorithm for \em at \em most the given number of iterations, but if a fitness value of goal is reached earlier, don't continue.
         *  @details This is the main routine for the genetic algorithm. It will call "evolveTo()" for each of the populations.
         *  @param goal  Fitness limit which will trigger conditional return.
         *  @param steps Number of iterations before returning.
         *  @returns The best fitness-value found.
         */
        template <class T>
        double Genetic<T>::evolveTo ( double goal, uint& steps ) {
            if ( steps == 0 || steps > maxGenerations ) {
                steps = maxGenerations;
            }

            unsigned int thisEffort = steps;
            steps = 0;
            double fn = 0.0; //, tmp = 0.0;

//                 for ( unsigned int j = 0; j < nPopulations; ++j ) {
//                     tmp = populations[j]->evolveTo ( goal, thisEffort );
//
//                     if ( tmp > fn ) {
//                         fn = tmp;
//                     }
//
//                     steps += thisEffort;
//                 }

            totalEffort += steps;
            return fn;
        }

        /*! @fn T* Genetic<T>::getFittest( void )
         *  @brief Get the parameters for the best fit.
         *  @returns A pointer to the "genome" (i.e. the parameters) of the fittest individual.
         */
        template <class T>
        T* Genetic<T>::getFittest ( void ) {
            return nullptr; //populations[0]->individuals[0]->genome;
        }

        /*! @fn T* Genetic<T>::getIndividual( int i  )
         *  @brief Get the parameters for the best fit.
         *  @returns A pointer to the "genome" (i.e. the parameters) of the fittest individual.
         */
        template <class T>
        T* Genetic<T>::getIndividual ( unsigned int i ) {
//                 if ( i < popSize ) {
//                     return populations[0]->individuals[i]->genome;
//                 }

            return NULL;
        }



//             /*-! @fn void Genetic<T>::setGenomeSize(const uint& n)
//              *  @brief Set the genome size, i.e. how many parameters we want to fit.
//              */
//             template <class T>
//             void Genetic<T>::setGenomeSize ( const uint& n ) {
// //                 if ( n == genomeSize || n == 0 ) {
// //                     return;
// //                 }
// //
// //                 ParameterRange<T>* newRanges = new ParameterRange<T>[n];
// //
// //                 for ( unsigned int i = 0; i < n; ++i ) {
// //                     if ( i < genomeSize ) {
// //                         newRanges[i] = ranges[i];
// //                     }
// //                 }
// //
// //                 delete[] ranges;
// //                 ranges = newRanges;
// //                 genomeSize = n;
//                 init();
//             }

        /*! @fn void Genetic<T>::setParameterRange( unsigned int n, T low, T high, bool per )
         *  @brief Set the limits for parameter number n (n >= 0)
         *  @param n Index of the parameter to modify (starting at 0 for the first parameter).
         *  @param low Smallest allowed value for this parameter.
         *  @param high Largest allowed value for this parameter.
         *  @param per Is this parameted periodic ? (i.e. should a value slightly outside the high limit be re-mapped to value slightly
         *  inside the low limit?)
         */
        template <class T>
        void Genetic<T>::setParameterRange ( unsigned int n, T low, T high, bool per ) {
            ranges[n].setRange ( low, high );
            ranges[n].setPeriodic ( per );
            //ranges[n].calculatePrecision();
            ranges[n].calculateMask();
            /*
                            for ( unsigned int i = 0; i < nPopulations; ++i ) {
                                populations[i]->loadRanges();
                            }*/
        }

        /*! @fn void Genetic<T>::setPeriodicParameter( unsigned int n, bool val )
         *  @brief Specify if parameter n is periodic.
         *  @param n Index of the parameter to modify (starting at 0 for the first parameter).
         *  @param val true/false
         */
        template <class T>
        void Genetic<T>::setPeriodicParameter ( unsigned int n, bool val ) {
            ranges[n].setPeriodic ( val );

//                 for ( unsigned int i = 0; i < nPopulations; ++i ) {
//                     populations[i]->ranges[n].setPeriodic ( val );
//                 }
        }

        /*! @fn void Genetic<T>::setDynamicRange( unsigned int n, bool val )
         *  @brief Specify if parameter n is dynamic.
         *  @details Should the allowed range be modified according to the standard deviation of current values of the population ? \n
         *  The allowed ranges will be modified at a frequency specified by the statFrequency parameter.
         *  @param n Index of the parameter to modify (starting at 0 for the first parameter).
         *  @param val true/false
         */
        template <class T>
        void Genetic<T>::setDynamicRange ( unsigned int n, bool val ) {
            ranges[n].setDynamic ( val );

//                 for ( unsigned int i = 0; i < nPopulations; ++i ) {
//                     populations[i]->setDynamicRange ( n, val );
//                 }
        }

        /*! @fn void Genetic<T>::setDynamicRanges( bool val )
         *  @brief Convenience function for (un)setting the dynamic feature for \em all parameters.
         *  @param val true/false
         */
        template <class T>
        void Genetic<T>::setDynamicRanges ( bool val ) {
            for ( unsigned int i = 0; i < nParameters; ++i ) {
                ranges[i].setDynamic ( val );
            }

//                 for ( unsigned int i = 0; i < nPopulations; ++i ) {
//                     populations[i]->setDynamicRanges ( val );
//                 }
        }


        /*! @fn string Genetic<T>::dump( void ) const
         *  @brief Generate a string with stats/debugging info of the genetic object.
         */
        template <class T>
        std::string Genetic<T>::dump ( void ) const {
            std::string tmp = "Genetic:    (" + redux::util::hexString ( this ) + ")\n";
            tmp += "     Population size:        " + std::to_string ( popSize ) + "\n";
            tmp += "     Number of Parameters:         " + std::to_string ( nParameters ) + "\n";
            tmp += "     Number of Generations:        " + std::to_string ( totalEffort ) + "\n";
            tmp += "     Replacement Method:           ";
            if ( replaceMethod == REPLACE_ALL ) {
                tmp += "Full Generational Replacement ";
            } else if ( replaceMethod == REPLACE_RANDOM ) {
                tmp += "Replace Random                ";
            } else if ( replaceMethod == REPLACE_WORST ) {
                tmp += "Replace Worst                 ";
            } else tmp += "?????????????                 ";
            tmp +=  "\n" ;

//                 for ( unsigned int i = 0; i < nPopulations; ++i ) {
//                     tmp += populations[i]->dump();
//                 }

            return tmp;
        }


    }       // algo

}       // redux

#endif  //  REDUX_ALGO_GENETIC_HPP
