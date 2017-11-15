#include "redux/util/stopwatch.hpp"

#include <sstream>

using namespace redux::util;
using namespace std;

StopWatch::StopWatch() {
    
    start();

}


void StopWatch::analyze(void) {

    deltaT = subtract(endT, beginT);
    deltaR = subtract(endR, beginR);

    totRes = add( deltaR.ru_utime, deltaR.ru_stime );

    prePrint();

}


void StopWatch::prePrint(void) {

    trimTime( deltaT );
    trimTime( totRes );
    trimTime( deltaR.ru_utime );
    trimTime( deltaR.ru_stime );

    setSuffix( deltaT, suffixT, lapsedT );
    setSuffix( totRes, suffixTot, lapsedTot );
    setSuffix( deltaR.ru_utime, suffixRU, lapsedRU );
    setSuffix( deltaR.ru_stime, suffixRS, lapsedRS );

}


void StopWatch::trimTime( timeval& tv ) {
    
    if ( tv.tv_usec < 0 ) {
        tv.tv_sec--;
        tv.tv_usec += 1000000;
    }

}


void StopWatch::setSuffix( timeval& tv, std::string& str, double& res ) {
    
    if ( tv.tv_sec >= 0 ) {
        if ( tv.tv_sec == 0 ) {
            res = (double)tv.tv_usec/1000;
            str = "ms";
        } else {
            res = (double)((int)(tv.tv_usec/1000))/1000 + tv.tv_sec;
            str = "s";
        }
    } else {
        res = 0.0;
        str = "s";
    }
    
}


timeval StopWatch::add(const timeval& tv1, const timeval& tv2 ) const {

    timeval tmp;

    tmp.tv_sec = tv1.tv_sec + tv2.tv_sec;
    tmp.tv_usec = tv1.tv_usec + tv2.tv_usec;

    return tmp;
}


timeval StopWatch::subtract(const timeval& tv1, const timeval& tv2 ) const {

    timeval tmp;

    tmp.tv_sec = tv1.tv_sec - tv2.tv_sec;
    tmp.tv_usec = tv1.tv_usec - tv2.tv_usec;

    return tmp;
}


rusage StopWatch::add(const rusage& ru1, const rusage& ru2 ) const {

    rusage tmp;

    tmp.ru_utime = add( ru1.ru_utime, ru2.ru_utime );  // Total amount of user time used
    tmp.ru_stime = add( ru1.ru_stime, ru2.ru_stime );  // Total amount of system time used

    /*  tmp.ru_maxrss   = ru1.ru_maxrss   + ru2.ru_maxrss;      // Maximum resident set size (in kilobytes).
        tmp.ru_ixrss    = ru1.ru_ixrss    + ru2.ru_ixrss;       // Amount of sharing of text segment memory (kilobyte+seconds).
        tmp.ru_idrss    = ru1.ru_idrss    + ru2.ru_idrss;       // Amount of data segment memory used (kilobyte+seconds).
        tmp.ru_isrss    = ru1.ru_isrss    + ru2.ru_isrss;       // Amount of stack memory used (kilobyte+seconds).
        tmp.ru_minflt   = ru1.ru_minflt   + ru2.ru_minflt;      // Number of soft page faults.
        tmp.ru_majflt   = ru1.ru_majflt   + ru2.ru_majflt;      // Number of hard page faults.
        tmp.ru_nswap    = ru1.ru_nswap    + ru2.ru_nswap;       // Number of times a process was swapped out of physical memory.
        tmp.ru_inblock  = ru1.ru_inblock  + ru2.ru_inblock;     // Number of input operations via the file system.
        tmp.ru_oublock  = ru1.ru_oublock  + ru2.ru_oublock;     // Number of output operations via the file system.
        tmp.ru_msgsnd   = ru1.ru_msgsnd   + ru2.ru_msgsnd;      // Number of IPC messages sent..
        tmp.ru_msgrcv   = ru1.ru_msgrcv   + ru2.ru_msgrcv;      // Number of IPC messages received.
        tmp.ru_nsignals = ru1.ru_nsignals + ru2.ru_nsignals;    // Number of signals delivered.
        tmp.ru_inblock  = ru1.ru_inblock  + ru2.ru_inblock;     // Number of input operations via the file system..
        tmp.ru_nvcsw    = ru1.ru_nvcsw    + ru2.ru_nvcsw;       // Number of voluntary context switches.
        tmp.ru_nivcsw   = ru1.ru_nivcsw   + ru2.ru_nivcsw;      // Number of involuntary context switches.
    */

    return tmp;
}


rusage StopWatch::subtract(const rusage& ru1, const rusage& ru2 ) const {

    rusage tmp;
    tmp.ru_utime = subtract( ru1.ru_utime, ru2.ru_utime );  // Total amount of user time used
    tmp.ru_stime = subtract( ru1.ru_stime, ru2.ru_stime );  // Total amount of system time used
    return tmp;
    
}


void StopWatch::start() {
    
    gettimeofday(&beginT,0);
    getrusage( RUSAGE_SELF, &beginR );
    
}


void StopWatch::stop() {
    
    gettimeofday(&endT,0);
    getrusage( RUSAGE_SELF, &endR );
    analyze();
    
}


float StopWatch::getLoad( void ) {
    
    stop();
    return (lapsedTot/lapsedT);
   
}


float StopWatch::getSeconds( void ) {
    
    stop();
    if( suffixT == "ms" ) return lapsedT/1000.0;
    return lapsedT;
   
}


float StopWatch::getCPUSeconds( void ) {
    
    stop();
    if( suffixTot == "ms" ) return lapsedTot/1000.0;
    return lapsedTot;
   
}


string StopWatch::print(void) {
    
    stop();
    float usage = (100*lapsedTot/lapsedT);
    if( suffixT == "ms" && suffixTot == "s") usage *= 1000;
    else if( suffixT == "s" && suffixTot == "ms") usage /= 1000;
    std::stringstream ss;
    ss.setf( std::ios::fixed, std:: ios::floatfield );
    ss.precision(3);
    ss << "Elapsed:  " << lapsedT << suffixT << " wall, ";
    ss << lapsedRU << suffixRU << " user + ";
    ss << lapsedRS << suffixRS << " system = ";
    ss << lapsedTot << suffixTot << " CPU (";
    ss.precision(2);
    ss << usage << "%)";
   
    return ss.str();
    
}
