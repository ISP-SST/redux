#ifndef REDUX_LOGGING_LOGTOFILE_HPP
#define REDUX_LOGGING_LOGTOFILE_HPP

#include "redux/logging/logtostream.hpp"

namespace redux {

    namespace logging {

        
        class LogToFile : public LogToStream {
            
        public:
            LogToFile( const std::string &filename, uint8_t m=LOG_MASK_ANY, bool replace=false, unsigned int flushPeriod=1);
            ~LogToFile() { fout.close(); }
            
        private:
            
            std::ofstream fout;

            friend class Log;
            friend class Logger;
            
        };

    } // end namespace logging

} // end namespace redux



#endif // REDUX_LOGGING_LOGTOFILE_HPP
