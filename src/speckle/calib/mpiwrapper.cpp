#include "mpiwrapper.hpp"

#include "defs.hpp"
#include "wam.hpp"

#include <functional>   // std::bind
#include <iostream>
#include <iomanip>
#include <complex>
#include <thread>

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

using namespace redux::speckle;
using namespace std;

namespace {

    double xl1[2] = { 0, 0 };
    double xu1[2] = { ( M_PI ) * 2.0, 1 };
    double xl2[4] = { 0, 0, 0, 0 };
    double xu2[4] = { ( M_PI ) * 2.0, 1, ( M_PI ) * 2.0, 1 };

    const gsl_rng_type *rngType;

    struct ThreadArg {
        ThreadArg( size_t nI, size_t nC ) : calls( nC ), iterations( nI ) {
            s1 = gsl_monte_vegas_alloc( 2 );
            s2 = gsl_monte_vegas_alloc( 4 );
            
            r = gsl_rng_alloc( rngType );
            gsl_rng_set( r, time( 0 ) );
            G = { &wam, 2, &p };
            H = { &wam2, 4, &p };
        }
        ~ThreadArg( void ) {
            gsl_monte_vegas_free( s1 );
            gsl_monte_vegas_free( s2 );
            gsl_rng_free(r);
        }
        
        bool init(void) {
            s2->stage = s1->stage = 0;
            s2->iterations = s1->iterations = 1;
            double dummy1,dummy2;
            gsl_monte_vegas_integrate (&G, xl1, xu1, 2,  8*calls, r, s1, &dummy1, &dummy2);
            gsl_monte_vegas_integrate (&H, xl2, xu2, 4, 16*calls, r, s2, &dummy1, &dummy2);
            s2->stage = s1->stage = 1;
            s2->iterations = s1->iterations = iterations;
            return false;
        }
        
        void run( const function<bool(double*&)>& next ) { 

            double* data;
            bool doInit(true);
            while( next(data) ) {
                p.q_five_thirds = data[0];                // current frequency    (actually: pow( q, 5.0/3.0 ) )
                p.alpha_five_thirds = data[1];            // current alpha        (actually: pow( 1.0 / alpha, 5.0/3.0 ) )
                p.delta = {data[2], data[3]};
                
                if(doInit) doInit = init();

                double ltf_result, stf_result, ltf_sigma, stf_sigma;

                // long exposure
                uint16_t loop = 0;
                do {
                    gsl_monte_vegas_integrate (&G, xl1, xu1, 2, 2 * calls, r, s1, &ltf_result, &ltf_sigma);
                    loop++;
                } while ( (fabs(gsl_monte_vegas_chisq(s1) - 1.0) > 0.5) && (loop <= iterations));

                // speckle transfer function
                loop = 0;

                do {
                    gsl_monte_vegas_integrate (&H, xl2, xu2, 4, 4 * calls, r, s2, &stf_result, &stf_sigma);
                    loop++;
                } while ( (fabs(gsl_monte_vegas_chisq(s2) - 1.0) > 0.5) && (loop <= iterations));

                memset(data,0,4*sizeof(double));
                if( finite(ltf_result) && finite(ltf_sigma) ) {
                    data[0] = ltf_result;
                    data[1] = ltf_sigma;
                }
                if( finite(stf_result) && finite(stf_sigma) ) {
                    data[2] = stf_result;
                    data[3] = stf_sigma;
                }
            }
        }
        size_t calls, iterations;
        struct parameters p;
        gsl_monte_vegas_state *s1, *s2;
        gsl_rng *r;
        gsl_monte_function G, H;
    };

    void arrayDeleter( double* p ) {
        delete[] p;
    }
}

MpiWrapper::MpiWrapper( int argc, char *argv[] ) : running( true ), zd( ZernikeData::get() ) {

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
    MPI_Comm_size( MPI_COMM_WORLD, &worldSize );

    if( argc != 6 ) return;

    int threadCount = 4; //std::max(std::thread::hardware_concurrency(),1);      // some implementations (e.g. gcc 4.4) seems to return 0, default to 1 thread in that case.

    nThreads.push_back( threadCount );

    nCalls = atol( argv[4] );
    nIterations = 10;
    thread initThread;
    
    if( myRank == 0 ) {     // master
        
        nThreads.resize( worldSize, 1 );
        nTotalThreads = nThreads[0];

        efficiencies = ZernikeData::load( argv[5] );
        
        // Start initializing radial/angular lookup-tables in the background while establishing the MPI group.
        initThread = std::thread( (void(ZernikeData::*)(const vector<float>&,uint32_t,float))&ZernikeData::init,
                                         &zd,efficiencies, SPECKLE_IMAX, SPECKLE_COV_CUTOFF);

        for( int i = worldSize - 1; i; --i ) {
            int slaveThreads;
            MPI_Recv( &slaveThreads, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
            nThreads[status.MPI_SOURCE] = slaveThreads;
            nTotalThreads += slaveThreads;
            //cout << "node #" << status.MPI_SOURCE << " reported " << slaveThreads << " threads.  Total = " << nTotalThreads << endl;
        }

        double as = atof( argv[1] );
        double ae = atof( argv[2] );
        double ast = atof( argv[3] );

        for( double j = as; j <= ae + SPECKLE_EPS; j += ast ) {
            alpha.push_back( j );
        }

        double freqstep = ( double )( 1.0 / SPECKLE_FREQSTEPS );
        for( int i = 0; i < SPECKLE_FREQSTEPS; ++i ) {
            freq.push_back( i * freqstep );
        }

        makeCube();

        int nParts = freq.size() * alpha.size();
        int nPartsPerThread = nParts / nTotalThreads;

        blockSizes.resize( worldSize, 0 );
        blockOffsets.resize( worldSize, 0 );

        for( int i = 0; i < worldSize; ++i ) {
            blockSizes[i] = nThreads[i] * nPartsPerThread;
            nParts -= blockSizes[i];
            blockSizes[i] *= 4;

            if( i ) blockOffsets[i] = blockOffsets[i - 1] + blockSizes[i - 1];
        }

        if( nParts > 0 ) *( blockSizes.rbegin() ) += 4 * nParts; // just add leftover parts to the last slave for now.

        int sz = efficiencies.size();
        MPI_Bcast( &sz, 1, MPI_INT, 0, MPI_COMM_WORLD );

        if( sz ) {
            MPI_Bcast( efficiencies.data(), sz, MPI_FLOAT, 0, MPI_COMM_WORLD );
        }

        MPI_Bcast( blockSizes.data(), worldSize, MPI_INT, 0, MPI_COMM_WORLD );

    } else {        // slave
        int tmp;
        MPI_Send( nThreads.data(), 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
        MPI_Bcast( &tmp, 1, MPI_INT, 0, MPI_COMM_WORLD );
        if( tmp ) {
            efficiencies.resize( tmp );
            MPI_Bcast( efficiencies.data(), tmp, MPI_FLOAT, 0, MPI_COMM_WORLD );
        }
        // Start initializing radial/angular lookup-tables in the background while establishing the MPI group.
        initThread = std::thread( (void(ZernikeData::*)(const vector<float>&,uint32_t,float))&ZernikeData::init,
                                         &zd,efficiencies, SPECKLE_IMAX, SPECKLE_COV_CUTOFF);        // Start initializing radial/angular lookup-tables.
        blockSizes.resize( worldSize );
        blockOffsets.resize( worldSize, 0 );
        MPI_Bcast( blockSizes.data(), worldSize, MPI_INT, 0, MPI_COMM_WORLD );

    }

    initThread.join();

}


MpiWrapper::~MpiWrapper( void ) {

    MPI_Finalize();

}


bool MpiWrapper::next(const double* base, double *& ptr) {
    unique_lock<mutex> lock(mtx);
    if( currentOffset < 4 ) return false;
    currentOffset -= 4;
    ptr = const_cast<double*>(base+currentOffset);
    return true;
}


void MpiWrapper::run( void ) {
    
    using namespace std::placeholders;

    currentOffset = 0;
    int myThreads = nThreads.empty() ? 1 : nThreads[0];

    if( myRank < (int)blockSizes.size() ) {
        currentOffset = blockSizes[myRank];
    }

    if( currentOffset == 0 ) return;

    if( currentOffset % 4 ) {
        cout << myRank << ": Weird blocksize (" << currentOffset << "), bailing out." << endl;
        return;
    }

    unique_ptr<double[]> buf( new double[currentOffset] );
    double* ptr = buf.get();
    MPI_Scatterv( data.get(), blockSizes.data(), blockOffsets.data(), MPI_DOUBLE, ptr, currentOffset, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    function<bool(double*&)> next_fn = std::bind(&MpiWrapper::next,this,ptr,_1);

    gsl_rng_env_setup();
    rngType = gsl_rng_default;

    vector<shared_ptr<ThreadArg>> threadArgs;
    vector<thread> threads;

    for( int i = myThreads - 1; i >= 0; --i ) {
        shared_ptr<ThreadArg> ta( new ThreadArg( nIterations, nCalls ) );
        threadArgs.push_back( ta );
        threads.push_back( std::thread( &ThreadArg::run, ta.get(), next_fn ) );
    }

    for( auto it = threads.begin(); it != threads.end(); it++ ) {
        it->join();
    }

    MPI_Gatherv( ptr, blockSizes[myRank], MPI_DOUBLE, data.get(), blockSizes.data(), blockOffsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD );

}

void MpiWrapper::makeCube( void ) {


    uint32_t nParts = freq.size() * alpha.size();
    uint32_t dataSize = 4 * nParts;

    data.reset( new double[dataSize], arrayDeleter );           // data[alpha][q][param]

    // input:  [0]=pow( freq, 5.0/3.0 ) ), [1]=pow( 1.0 / alpha, 5.0/3.0 )
    //         [2]=delta_x                 [3]= delta_y

    memset( data.get(), 0, dataSize * sizeof( double ) );
    double* ptr = data.get();

    for( uint i = 0; i < alpha.size(); ++i ) {
        double a_53 = pow( 1.0 / alpha[i], 5.0 / 3.0 );
        for( uint j = 0; j < freq.size(); ++j ) {
            complex_t delta = polar( freq[j], M_PI_4 );
            ptr[0] = pow( freq[j], 5.0 / 3.0 );
            ptr[1] = a_53;
            ptr[2] = delta.real();
            ptr[3] = delta.imag();
            ptr += 4;
        }
    }

}


void MpiWrapper::writeResults(std::ostream& strm) const {
    
    if( myRank != 0 || !data ) return;
    
    streamsize oldPrec = strm.precision(8);
    
    strm << "# (partially) corrected Zernikes: " << efficiencies.size() << "   efficiencies=[";
    for(uint i=0; i<efficiencies.size(); ++i) { strm << (i?", ":" ") << freq[i]; }
    strm << "]\n# spatial frequency\n";
    for(uint i=0; i<freq.size(); ++i) strm << freq[i] << endl;
    strm << "# seeing alpha\n";
    for(uint i=0; i<alpha.size(); ++i) strm << alpha[i] << endl;
    strm << "# functions\n";
    double* ptr = data.get();
    double sr, sr_sigma;

    for( uint i = 0; i < alpha.size(); ++i ) {
        strm << "letf           s_letf         stf            s_stf          sr             s_sr\n";
        for( uint j = 0; j < freq.size(); ++j ) {
            if(ptr[0] > 0 && ptr[2] && finite(ptr[1]) && finite(ptr[3]) ) {
                sr = ptr[0]*ptr[0]/ptr[2];
                sr_sigma = sr*sqrt(ptr[1] * ptr[1] / (ptr[0] * ptr[0]) + ptr[3] * ptr[3] / (ptr[2] * ptr[2]));
            } else sr = sr_sigma = 0;
            strm << setw(15) << left << ptr[0];             // LTF
            strm << setw(15) << left << ptr[1];             // LTF sigma
            strm << setw(15) << left << ptr[2];             // STF
            strm << setw(15) << left << ptr[3];             // STF sigma
            strm << setw(15) << left << sr;                 // Spectral ratio
            strm << setw(15) << left << sr_sigma << endl;     // Spectral ratio sigma
            ptr += 4;
        }
    }
    strm.precision(oldPrec);

}


