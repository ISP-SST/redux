#ifndef REDUX_LOGGING_LOGTOSTREAM_HPP
#define REDUX_LOGGING_LOGTOSTREAM_HPP

#include "redux/logging/logoutput.hpp"

#include <fstream>
#include <mutex>
#include <vector>


namespace redux {

    namespace logging {

        
        class LogToStream : public LogOutput {


        public:
            LogToStream( std::ostream &os, uint8_t m=LOG_MASK_ANY, unsigned int flushPeriod=1);
            ~LogToStream();

            void flushBuffer( void );

        private:
            void writeFormatted( const LogItem& );

            std::ostream& out;
            bool color,localtime;

        };

    } // end namespace logging

} // end namespace redux



#endif // REDUX_LOGGING_LOGTOSTREAM_HPP
