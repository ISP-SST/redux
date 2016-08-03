#ifndef REDUX_LOGGING_LOGTONETWORK_HPP
#define REDUX_LOGGING_LOGTONETWORK_HPP

#include <redux/logging/logoutput.hpp>
#include <redux/network/tcpconnection.hpp>

#include <fstream>
#include <mutex>
#include <vector>


namespace redux {

    namespace logging {

        
        class LogToNetwork : public LogOutput {

        public:
            LogToNetwork( const network::TcpConnection::Ptr&, uint32_t id, uint8_t m=LOG_MASK_ANY, unsigned int flushPeriod=5);
            ~LogToNetwork();

            void flushBuffer( void );

        private:

            network::TcpConnection::Ptr conn;

        };

    } // end namespace logging

} // end namespace redux



#endif // REDUX_LOGGING_LOGTONETWORK_HPP
