#ifndef REDUX_LOGGING_LOGENTRY_HPP
#define REDUX_LOGGING_LOGENTRY_HPP

#include <redux/util/stringutil.hpp>

#include <string>
#include <ostream>

#include <boost/io/ios_state.hpp>
#include <boost/thread/xtime.hpp>

namespace redux {
    
    namespace logging {

        class LogEntry;

        enum severity_level {
            sev_none,
            sev_trace,
            sev_debug,
            sev_detail,
            sev_normal,
            sev_warning,
            sev_error,
            sev_critical
        };

        enum LogLevel {
            LOG_LEVEL_NONE,
            LOG_LEVEL_FATAL,
            LOG_LEVEL_ERROR,
            LOG_LEVEL_WARNING,
            LOG_LEVEL_NOTICE,
            LOG_LEVEL_NORMAL,
            LOG_LEVEL_DETAIL,
            LOG_LEVEL_DEBUG,
            LOG_LEVEL_TRACE
        };

        
        enum LogMask {
            LOG_MASK_NONE    = 0,
            LOG_MASK_FATAL   = 1,
            LOG_MASK_ERROR   = 2,
            LOG_MASK_WARNING = 4,
            LOG_MASK_NOTICE  = 8,
            LOG_MASK_NORMAL  = 16,
            LOG_MASK_DETAIL  = 32,
            LOG_MASK_DEBUG   = 64,
            LOG_MASK_TRACE   = 128,
            LOG_MASK_ANY     = 255
        };

        
        class LogEntry {
        public:
            
            LogEntry(void) : mask(LOG_MASK_ERROR), settings(buffer) {}
            LogEntry( LogMask m, const std::string &msg ) : mask(m), message(msg), settings(buffer) {}
            LogEntry( LogMask m, const std::string &msg, const boost::xtime &t ) : mask(m), message(msg),
                    entryTime(t), settings(buffer) {}

            LogEntry( const LogEntry& rhs) : mask(rhs.mask), message(rhs.message),
                    entryTime(rhs.entryTime), settings(buffer) {
            }

            LogEntry &operator=( const LogEntry& rhs ) {
                mask = rhs.mask;
                message = rhs.message;
                entryTime = rhs.entryTime;
                return *this;
            }
            
            inline void now(void){
                entryTime = boost::posix_time::microsec_clock::local_time();
                //entryTime = boost::posix_time::microsec_clock::universal_time();
            }
            inline void setMask( LogMask m ) { mask = m; }
            inline uint8_t getMask(void) const { return mask; }
            inline const char *getMessage(void) const { return message.c_str(); }
            inline const boost::posix_time::ptime &getTime(void) const { return entryTime; }

            template <typename T>
            inline LogEntry &operator<<(const T &v) {
                buffer << v;
                return *this;
            }

            // forward iostream manipulators (i.e. endl)
            inline LogEntry &operator<<(std::ostream &(*f)(std::ostream &)) {
                buffer << f;
                return *this;
            }

            // forward ios manipulators
            inline LogEntry &operator<<(std::ios &(*f)(std::ios &)) {
                buffer << f;
                return *this;
            }

            // forward ios_base manipulators
            inline LogEntry &operator<<(std::ios_base &(*f)(std::ios_base &)) {
                buffer << f;
                return *this;
            }

            // User defined manipulators for the LogEntry class
            inline LogEntry &operator<<(LogEntry &(*f)(LogEntry &)) {
                return f(*this);
            }
            
            inline void finalize(void) {
                now();
                message = buffer.str();
                buffer.str( std::string() );
                settings.restore();
            }

            inline void prepend(const std::string& m) { message = m+message; }

            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);


        private:
            uint8_t mask;
            std::string message;
            boost::posix_time::ptime entryTime;
            std::stringstream buffer;
            boost::io::ios_all_saver  settings;

        };


    }   // log

}   // redux

std::ostream &operator<<(std::ostream &os, const redux::logging::LogEntry &le);

#endif //  REDUX_LOGGING_LOGENTRY_HPP
