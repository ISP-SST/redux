#include "defs.hpp"
#include "wam.hpp"
#include "zernike.hpp"

#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

using namespace redux::speckle;
using namespace std;

int main (int argc, char *argv[]) {
    
    // Declare important MPI variables here
    int myRank;                    // Rank of processor
    int worldSize;                  // Number of available processors
    MPI_Status *status;          // status variable for communication

    // parameters
    double params[4];
    vector<float> eff;
    
    // Initialize the (singleton) container.
    ZernikeData& zd = ZernikeData::get();
    
    // result variables
    double *result;

    // Initialize MPI
    MPI_Init (&argc, &argv);
    // Get rank
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    // Get the total number of processors
    MPI_Comm_size (MPI_COMM_WORLD, &worldSize);
    
    // Get memory for status variable
    status = (MPI_Status *) malloc (sizeof (MPI_Status));

    if (myRank == 0) {
        // helpers
        double step = (double) (1.0 / (SPECKLE_FREQSTEPS));
        int l;
        double tmp[7];
        
        // deal with the command line arguments
        if (argc != 6) {
            printf ("\nIncorrect number of arguments.\n");
            printf ("\nArguments:\n");
            printf ("  alp_start  alp_end  alp_step  #calls  efile\n");
            printf ("where:\n");
            printf ("\nalp_start  alp_end  alp_step:   ");
            printf ("start, end and stepwidth of alpha\n");
            printf ("#calls:                         ");
            printf ("number of calls (increases accuracy and computational time!)\n");
            printf ("efile:                          ");
            printf ("Filename with efficency coefficents\n\n");
            printf ("Recommended: Output redirection with ... 2> \"filename\"\n");
            printf ("\nThanks for the attention :)\n");
            printf ("F. Woeger, KIS\n\n");
            MPI_Finalize();
            exit (0);
        }

        // alpha settings
        double as = atof (argv[1]);
        double ae = atof (argv[2]);
        double ast = atof (argv[3]);

        // get memory
        result = (double *) malloc (7 * SPECKLE_FREQSTEPS * sizeof (double));

        // initialize using provided efficency file
        eff = zd.init(argv[5]);
        int sz = eff.size();
        string effString = "eff=[";
        for(int i=0; i<sz; ++i) {
            if(i) effString += " ,";
            effString += to_string((long double)eff[i]);
        }
        effString += "]";
        MPI_Bcast (&sz, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if(sz) {
            MPI_Bcast (eff.data(), sz, MPI_FLOAT, 0, MPI_COMM_WORLD);
        }
        fprintf (stderr, "# (partially) corrected Zernikes: %d   %s\n", sz-1, effString.c_str() );
        fprintf (stderr, "# spatial frequency\n");

        for (double i = 0.0; i < 1.0; i += step) fprintf (stderr, "%.7lf\n", i);

        fprintf (stderr, "# seeing alpha\n");

        for (double j = as; j <= ae + SPECKLE_EPS; j += ast) fprintf (stderr, "%.7lf\n", j);

        fprintf (stderr, "# functions\n");

        int nActive(1);     // start at 1 to simplify loops below.
        
        for (double j = as; j <= ae + SPECKLE_EPS; j += ast) {      // loop over alpha
            int frequencyIndex(0);
            int nCompleted(0);
            fprintf (stdout, "Starting alpha=%f  [%f, %f]\n", j, as, ae);
            for (double i = 0.0; i < 1.0; i += step, frequencyIndex++) {      // loop over frequencies

                params[0] = pow( i, 5.0/3.0 );                // current frequency   (actually: pow( q, 5.0/3.0 ) )
                params[1] = pow( 1.0 / j, 5.0/3.0 );          // current alpha       (actually: pow( 1.0 / alpha, 5.0/3.0 ) )
                complex_t delta = polar(i,M_PI_4);
                params[2] = delta.real();
                params[3] = delta.imag();
                if ( nActive < worldSize ) { // first send a job to each slave (skipping 0)
                    MPI_Send (params, 4, MPI_DOUBLE, nActive++, frequencyIndex, MPI_COMM_WORLD);
                } else {    // all slaves got a job, wait for result and then send next job. (NB: nActive remains the same)
                    MPI_Recv (tmp, 7, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status);
                    memcpy(result+7*status->MPI_TAG,tmp,7*sizeof(double));
                    nCompleted++;
                    MPI_Send (params, 4, MPI_DOUBLE, status->MPI_SOURCE, frequencyIndex, MPI_COMM_WORLD);
                }
                
                fprintf (stdout, "%d percent completed\r", (int) (nCompleted * step * 100));
                fflush (stdout);
            }

            // collect remaining jobs
            while (nActive > 1) {
                MPI_Recv (tmp, 7, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status);
                memcpy(result+7*status->MPI_TAG,tmp,7*sizeof(double));
                nCompleted++;
                fprintf (stdout, "%d percent completed\r", (int) (nCompleted * step * 100));
                fflush (stdout);
                nActive--;
            }

            // printout
            fprintf (stderr, "letf      s_letf    stf       s_stf     sr        s_sr\n");

            for (l = 0; l < SPECKLE_FREQSTEPS; l++) {
                fprintf (stderr, "%.7lf  %.7lf  %.7lf  %.7lf  %.7lf  %.7lf\n",
                         result[l * 7 + 1], result[l * 7 + 2], result[l * 7 + 3],
                         result[l * 7 + 4], result[l * 7 + 5], result[l * 7 + 6]);
            }
        }

        fprintf (stderr, "# end");

        // quit slaves
        for (int r = 1; r < worldSize; r++) {
            MPI_Send (params, 4, MPI_DOUBLE, r, SPECKLE_NOGO, MPI_COMM_WORLD);
        }

    } else {

        // setup function
        struct parameters p;
        gsl_monte_function G = { &wam, 2, &p };
        gsl_monte_function H = { &wam2, 4, &p };

        double xl1[2] = { 0, 0 };
        double xu1[2] = { (M_PI) * 2.0, 1 };
        double xl2[4] = { 0, 0, 0, 0 };
        double xu2[4] = { (M_PI) * 2.0, 1, (M_PI) * 2.0, 1 };

        // setup integrator
        size_t calls;
        double r1, r2, e1, e2;
        int loop, maxloop = 10;
        gsl_monte_vegas_state *s1 = gsl_monte_vegas_alloc (2);
        gsl_monte_vegas_state *s2 = gsl_monte_vegas_alloc (4);
        // Random Number Generator
        const gsl_rng_type *T;
        gsl_rng *r;
        gsl_rng_env_setup ();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);

        if (argc != 6) {
            MPI_Finalize();
            exit (0);
        }

        calls = (size_t) atol (argv[4]);

        // get memory
        result = (double *) malloc (7 * sizeof (double));

        // initialise efficency factors
        int sz;
        MPI_Bcast (&sz, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if( sz ) {

            eff.resize(sz);
            MPI_Bcast (eff.data(), sz, MPI_FLOAT, 0, MPI_COMM_WORLD);
            zd.init(eff);

            // begin computation
            while (1) {
                // receive the parameters from master
                MPI_Recv (params, 4, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status);

                // break, if tag says NOGO
                if (status->MPI_TAG == SPECKLE_NOGO) break;

                // set the parameters for the function
                result[0] = params[0];                      // current frequency
                p.q_five_thirds = params[0];                // current frequency    (actually: pow( q, 5.0/3.0 ) )
                p.alpha_five_thirds = params[1];            // current alpha        (actually: pow( 1.0 / alpha, 5.0/3.0 ) )
                p.delta = {params[2],params[3]};

                gsl_monte_vegas_init(s1);                   // reset what was learned in previous run
                gsl_monte_vegas_init(s2);                   // reset what was learned in previous run

                // do some work
                // long exposure
                gsl_monte_vegas_integrate (&G, xl1, xu1, 2, 4 * calls, r, s1, &r1, &e1);
                loop = 0;

                do {
                    gsl_monte_vegas_integrate (&G, xl1, xu1, 2, 2 * calls, r, s1, &r1, &e1);
                    loop++;
                } while ( (fabs (s1->chisq - 1.0) > 0.5) && (loop <= maxloop));

                if (!finite (r1) || !finite (e1)) {
                    result[1] = 0.0;
                    result[2] = 0.0;
                } else {
                    result[1] = r1;
                    result[2] = e1;
                }

                // speckle transfer function
                gsl_monte_vegas_integrate (&H, xl2, xu2, 4, 8 * calls, r, s2, &r2, &e2);
                loop = 0;

                do {
                    gsl_monte_vegas_integrate (&H, xl2, xu2, 4, 4 * calls, r, s2, &r2, &e2);
                    loop++;
                } while ( (fabs (s2->chisq - 1.0) > 0.5) && (loop <= maxloop));

                if (!finite (r2) || !finite (e2)) {
                    result[3] = 0.0;
                    result[4] = 0.0;
                } else {
                    result[3] = r2;
                    result[4] = e2;
                }

                // spectral ratio
                if (!finite (r1) || !finite (e1) || !finite (r2) || !finite (e2)) {
                    result[5] = 0.0;
                    result[6] = 0.0;
                } else {
                    result[5] = (r1 * r1) / r2;
                    result[6] = result[5] * sqrt (e1 * e1 / (r1 * r1) + e2 * e2 / (r2 * r2));

                    if (!finite (result[6])) result[6] = 0.0;
                }

                // send results back to master
                MPI_Send (result, 7, MPI_DOUBLE, 0, status->MPI_TAG, MPI_COMM_WORLD);
            }
        }

        // free integrator workspace
        gsl_monte_vegas_free (s1);
        gsl_monte_vegas_free (s2);
    }

    MPI_Finalize();
    return 0;
}
