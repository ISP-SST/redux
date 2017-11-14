#ifndef REDUX_LOGGING_LOGTONETWORK_HPP
#define REDUX_LOGGING_LOGTONETWORK_HPP

#include <redux/logging/logoutput.hpp>
#include <redux/network/host.hpp>
#include <redux/network/tcpconnection.hpp>

#include <fstream>
#include <mutex>
#include <vector>


namespace redux {

    namespace logging {

        
        class LogToNetwork : public LogOutput {

        public:
            LogToNetwork( boost::asio::io_service&, const network::Host::Ptr&, uint32_t id, uint8_t m=LOG_MASK_ANY, unsigned int flushPeriod=5);
            ~LogToNetwork();

            void connect(void);
            void flushBuffer( void );

        private:
            boost::asio::io_service& service;
            network::Host::Ptr host;
            network::TcpConnection::Ptr conn;
            uint32_t id;
            std::vector<LogItem> sendBuffer;
        };

    } // end namespace logging

} // end namespace redux



#endif // REDUX_LOGGING_LOGTONETWORK_HPP
