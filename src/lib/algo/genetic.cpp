#include "redux/algo/genetic.hpp"

using namespace redux::algo::genetic;
using namespace redux::algo;
using namespace redux::util;

GeneticBase::GeneticBase ( uint nParameters, uint popSize ) : crossOverType ( X_RAND ),
    replaceMethod ( REPLACE_WORST ),
    maxGenerations ( 500 ),
    minGenerations ( 5 ),
    nGenerations ( 100 ),
    popSize ( popSize ),
    nParameters ( nParameters ),
    statFrequency ( 10 ),
    totalEffort ( 0 ),
    crossoverProbability ( 0.85 ),
    mutationProbability ( 0.1 ),
    breedingProbability ( nullptr ),
    rng ( ::rand() ) {

//    setFitnessFunction ( defaultMerit );
//    setBreedingPDF ( defaultBreedingPDF );
//    init();
//    calcProb();
//    calcFitness();
}

void GeneticBase::calcProb ( void ) {

    float* breedPtr = new float[popSize];
    breedingProbability.reset ( breedPtr );

    breedPtr[0] = bPDF ( popSize, 1 );

    // cumulative probability.
    for ( uint i = 1; i < popSize; ++i ) {
        breedPtr[i] = breedPtr[i - 1] + bPDF ( popSize, i + 1 );
    }

    // normalize probability to 1.
    for ( uint i = 0; i < popSize; ++i ) {
        breedPtr[i] /= breedPtr[ popSize - 1 ];
    }
}


void GeneticBase::setPopSize ( const size_t& n ) {
    if ( n != popSize ) {
        popSize = n;
//        init();
    }
}

/*
template <class T>
Genetic<T>::~Genetic() {
    clear();
    delete[] ranges;

    if ( breedProb ) {
        delete[] breedProb;
    }
}

template <class T>
void Genetic<T>::init ( void ) {
    if ( populations ) {
        clear();
    }

    populations = new Population<T>*[ nPopulations ];

    for ( uint i = 0; i < nPopulations; ++i ) {
        populations[i] = new Population<T> ( this );
    }

    totalEffort = 0;
}

template <class T>
void Genetic<T>::clear ( void ) {
    if ( !populations ) {
        return;
    }

    for ( uint i = 0; i < nPopulations; ++i ) {
        delete populations[i];
    }

    delete[] populations;
    populations = NULL;
}

template <class T>
void Genetic<T>::reset ( void ) {
    for ( uint i = 0; i < nPopulations; ++i ) {
        populations[i]->reset();
    }

    totalEffort = 0;
}


template <class T>
void Genetic<T>::keep ( uint n ) {
    for ( uint i = 0; i < nPopulations; ++i ) {
        populations[i]->keep ( n );
    }
}

template <class T>
void Genetic<T>::sort ( void ) {
    for ( uint i = 0; i < nPopulations; ++i ) {
        populations[i]->sort();
    }
}

template <class T>
void Genetic<T>::calcProb ( void ) {
    if ( breedProb ) {
        delete[] breedProb;
    }

    breedProb = new double[ popSize ];
    breedProb[0] = bPDF ( popSize, 1 );

    // cumulative probability.
    for ( uint i = 1; i < popSize; ++i ) {
        breedProb[i] = breedProb[i - 1] + bPDF ( popSize, i + 1 );
    }

    // normalize probability to 1.
    for ( uint i = 0; i < popSize; ++i ) {
        breedProb[i] /= breedProb[ popSize - 1 ];
    }
}

template <class T>
void Genetic<T>::calcFitness ( void ) {
    for ( uint i = 0; i < nPopulations; ++i ) {
        populations[i]->calculateFitness();
    }
}

template <class T>
double Genetic<T>::evolve ( uint& steps ) {
    if ( steps == 0 || steps > maxGenerations ) {
        steps = maxGenerations;
    }

    uint thisEffort = steps;
    steps = 0;
    double fn = 0.0, tmp = 0.0;

    for ( uint j = 0; j < nPopulations; ++j ) {
        tmp = populations[j]->evolve ( thisEffort );

        if ( tmp > fn ) {
            fn = tmp;
        }

        steps += thisEffort;
    }

    totalEffort += steps;
    return fn;
}

template <class T>
double Genetic<T>::evolveTo ( double goal, uint& steps ) {
    if ( steps == 0 || steps > maxGenerations ) {
        steps = maxGenerations;
    }

    uint thisEffort = steps;
    steps = 0;
    double fn = 0.0, tmp = 0.0;

    for ( uint j = 0; j < nPopulations; ++j ) {
        tmp = populations[j]->evolveTo ( goal, thisEffort );

        if ( tmp > fn ) {
            fn = tmp;
        }

        steps += thisEffort;
    }

    totalEffort += steps;
    return fn;
}

template <class T>
T* Genetic<T>::getFittest ( void ) {
    return populations[0]->individuals[0]->genome;
}

template <class T>
T* Genetic<T>::getIndividual ( uint i ) {
    if ( i < popSize ) {
        return populations[0]->individuals[i]->genome;
    }

    return NULL;
}


template <class T>
void Genetic<T>::setCrossoverProbability ( float val ) {
    crossoverProbability = val;

    for ( uint i = 0; i < nPopulations; ++i ) {
        populations[i]->setCrossoverProbability ( val );
    }
}

template <class T>
void Genetic<T>::setMutationProbability ( float val ) {
    mutationProbability = val;

    for ( uint i = 0; i < nPopulations; ++i ) {
        populations[i]->setMutationProbability ( val );
    }
}

template <class T>
void Genetic<T>::setPopSize ( const int& value ) {
    if ( value != popSize && value != 0 ) {
        clear();
        popSize = value;
        init();
    }
}

template <class T>
void Genetic<T>::setGenomeSize ( const uint& n ) {
    if ( n == genomeSize || n == 0 ) {
        return;
    }

    ParameterRange<T>* newRanges = new ParameterRange<T>[n];

    for ( uint i = 0; i < n; ++i ) {
        if ( i < genomeSize ) {
            newRanges[i] = ranges[i];
        }
    }

    delete[] ranges;
    ranges = newRanges;
    genomeSize = n;
    init();
}

template <class T>
void Genetic<T>::setParameterRange ( uint n, T low, T high, bool per ) {
    ranges[n].setRange ( low, high );
    ranges[n].setPeriodic ( per );
    //ranges[n].calculatePrecision();
    ranges[n].calculateMask();

    for ( uint i = 0; i < nPopulations; ++i ) {
        populations[i]->loadRanges();
    }
}

template <class T>
void Genetic<T>::setPeriodicParameter ( uint n, bool val ) {
    ranges[n].setPeriodic ( val );

    for ( uint i = 0; i < nPopulations; ++i ) {
        populations[i]->ranges[n].setPeriodic ( val );
    }
}

template <class T>
void Genetic<T>::setDynamicRange ( uint n, bool val ) {
    ranges[n].setDynamic ( val );

    for ( uint i = 0; i < nPopulations; ++i ) {
        populations[i]->setDynamicRange ( n, val );
    }
}

template <class T>
void Genetic<T>::setDynamicRanges ( bool val ) {
    for ( uint i = 0; i < genomeSize; ++i ) {
        ranges[i].setDynamic ( val );
    }

    for ( uint i = 0; i < nPopulations; ++i ) {
        populations[i]->setDynamicRanges ( val );
    }
}

template <class T>
void Genetic<T>::insert ( T* genome ) {
    for ( uint i = 0; i < nPopulations; ++i ) {
        populations[i]->insert ( genome );
    }
}

template <class T>
std::string Genetic<T>::dump ( void ) const {
    std::string tmp = "Genetic:    (" + redux::util::hexString ( this ) + ")\n";
    tmp += "     Number of Populations:        " + std::to_string ( nPopulations ) + "\n";
    tmp += "     Number of Parameters:         " + std::to_string ( genomeSize ) + "\n";
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

    for ( uint i = 0; i < nPopulations; ++i ) {
        tmp += populations[i]->dump();
    }

    return tmp;
}
*/
