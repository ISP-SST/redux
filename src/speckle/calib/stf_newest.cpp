#include "mpiwrapper.hpp"

using namespace redux::speckle;
using namespace std;

void showUsage(void) {
    
    cout << "\nIncorrect number of arguments.\n"
         << "\nArguments:\n"
         << "  alp_start  alp_end  alp_step  #calls  efile\n"
         << "where:\n"
         << "\nalp_start  alp_end  alp_step:   "
         << "start, end and stepwidth of alpha\n"
         << "#calls:                         "
         << "number of calls (increases accuracy and computational time!)\n"
         << "efile:                          "
         << "Filename with efficency coefficents\n\n"
         << "Recommended: Output redirection with ... 2> \"filename\"\n"
         << "\nThanks for the attention :)\n"
         << "F. Woeger, KIS\n" << endl;
         
}


int main (int argc, char *argv[]) {
    
    MpiWrapper myMPI(argc,argv);

    if (argc != 6) {
        showUsage();
        return 0;
    }

    myMPI.run();
    myMPI.writeResults(cerr);
    
}
