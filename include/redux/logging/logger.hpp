#ifndef REDUX_LOGGING_LOGGER_HPP
#define REDUX_LOGGING_LOGGER_HPP

#include <redux/logging/logitem.hpp>
#include <redux/logging/logoutput.hpp>
#include <redux/network/host.hpp>
#include <redux/network/tcpconnection.hpp>
#include <redux/util/datautil.hpp>


#include <mutex>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

namespace bpo = boost::program_options;


namespace redux {

    namespace logging {
        
        class Logger : public LogOutput {
        public:

            Logger(void);
            explicit Logger( bpo::variables_map& );
            ~Logger();
            
            void append( LogItem& );
            void flushBuffer( void ) override;
            void flushAll( void );

            void addLogger( Logger& );      // forward output to another Logger instance.
            void addStream( std::ostream&, uint8_t m=0, unsigned int flushPeriod=1 );
            void addFile( const std::string &filename, uint8_t m=0, bool replace=false, unsigned int flushPeriod=1 );
            void addNetwork( boost::asio::io_service&, const network::Host::Ptr, uint32_t id, uint8_t m=0, unsigned int flushPeriod=5 );
            void removeOutput( const std::string& );
            void removeAllOutputs( void );
            void addConnection( network::TcpConnection::Ptr conn, network::Host::Ptr host );
            void removeConnection( network::TcpConnection::Ptr conn );
            void netReceive( network::TcpConnection::Ptr conn );
            void setContext( const std::string& c ) { context = c; };
            
            LogItem& getItem( LogMask m=LOG_MASK_NORMAL ) {
                threadItem.setLogger( this );
                threadItem.entry.setMask(m);
                threadItem.context = context;
                return threadItem;
            }

            LogEntry& getEntry( LogMask m=LOG_MASK_NORMAL ) {
                threadItem.setLogger( this );
                threadItem.entry.setMask(m);
                return threadItem.entry;
            }
                
            static inline void setDefaultLevel( int lvl ) { defaultLevelMask = LOG_UPTO(lvl); }
            static int getDefaultLevel(void);
            static inline void setDefaultMask( uint8_t m ) { defaultLevelMask = m; }
            static inline uint8_t getDefaultMask(void) { return defaultLevelMask; }
            
            static std::pair<std::string, std::string> customParser( const std::string& s );
            static std::string environmentMap( const std::string& );
            static bpo::options_description getOptions( const std::string& application_name );

        private:
            std::string context;
            
            typedef std::map<std::string, LogOutputPtr> OutputMap;
            OutputMap outputs;
            std::mutex outputMutex;
            
            std::map<network::TcpConnection::Ptr, network::Host::Ptr,
                     redux::util::PtrCompare<network::TcpConnection>> connections;
            
            static uint8_t defaultLevelMask;
            static thread_local LogItem threadItem;
        };

    }

}

#define LLOG_TRACE(mylog) mylog.getItem(redux::logging::LOG_MASK_TRACE)
#define LLOG_DEBUG(mylog) mylog.getItem(redux::logging::LOG_MASK_DEBUG)
#define LLOG_DETAIL(mylog) mylog.getItem(redux::logging::LOG_MASK_DETAIL)
#define LLOG_NOTICE(mylog) mylog.getItem(redux::logging::LOG_MASK_NOTICE)
#define LLOG(mylog) mylog.getItem()
#define LLOG_WARN(mylog) mylog.getItem(redux::logging::LOG_MASK_WARNING)
#define LLOG_ERR(mylog) mylog.getItem(redux::logging::LOG_MASK_ERROR)
#define LLOG_FATAL(mylog) mylog.getItem(redux::logging::LOG_MASK_FATAL)

#define LOG_TRACE LLOG_TRACE(logger)
#define LOG_DEBUG LLOG_DEBUG(logger)
#define LOG_DETAIL LLOG_DETAIL(logger)
#define LOG_NOTICE LLOG_NOTICE(logger)
#define LOG  LLOG(logger)
#define LOG_WARN  LLOG_WARN(logger)
#define LOG_ERR  LLOG_ERR(logger)
#define LOG_FATAL  LLOG_FATAL(logger)




#endif //   REDUX_LOGGING_LOGGER_HPP
