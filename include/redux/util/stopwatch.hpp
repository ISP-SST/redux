#ifndef REDUX_UTIL_STOPWATCH_HPP
#define REDUX_UTIL_STOPWATCH_HPP

#include <sys/time.h>
#include <sys/resource.h>

#include <string>

namespace redux {
        
    namespace util {

      class StopWatch {

          timeval beginT,endT,deltaT,totRes;
          rusage beginR,endR,deltaR;

          double lapsedT,lapsedRS,lapsedRU,lapsedTot;
          std::string suffixT,suffixRS,suffixRU,suffixTot;

          void analyze(void);
          void prePrint(void);
          void trimTime( timeval& tv ); 
          void setSuffix( timeval& tv, std::string& str, double& res );

          timeval add(const timeval& tv1, const timeval& tv2 ) const;
          timeval subtract(const timeval& tv1, const timeval& tv2 ) const;
          rusage add(const rusage& ru1, const rusage& ru2 ) const;
          rusage subtract(const rusage& ru1, const rusage& ru2 ) const;

        public:

          void start(void);
          void stop(void);
          
          float getSeconds( void );

          std::string print( void );

        }; // end class



    } // util



} // redux



#endif // REDUX_UTIL_STOPWATCH_HPP
