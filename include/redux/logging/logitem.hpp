#ifndef REDUX_LOGGING_LOGITEM_HPP
#define REDUX_LOGGING_LOGITEM_HPP

#include "logentry.hpp"


namespace redux {

    namespace logging {
        
        class Logger;

        struct LogItem {
            
            inline void setLogger( Logger* lg ) { logger = lg; }
            
            inline LogItem& operator<<(LogItem &(*f)(LogItem &)) {
                return f(*this);
            }

            inline LogItem& operator<<(const LogMask& m) {
                entry.setMask(m);
                return *this;
            }

            template <typename T>
            inline LogItem& operator<<(const T &v) {
                entry << v;
                return *this;
            }

            // forward iostream manipulators (i.e. endl)
            inline LogItem& operator<<(std::ostream &(*f)(std::ostream &)) {
                entry << f;
                return *this;
            }

            // forward ios manipulators
            inline LogItem &operator<<(std::ios &(*f)(std::ios &)) {
                entry << f;
                return *this;
            }

            // forward ios_base manipulators
            inline LogItem &operator<<(std::ios_base &(*f)(std::ios_base &)) {
                entry << f;
                return *this;
            }

            // apply LogEntry manipulators
            inline LogItem& operator<<(LogEntry &(*f)(LogEntry &)) {
                f(entry);
                return *this;
            }
            
            void endEntry(void);
            
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);

            Logger* logger;
            LogEntry entry;
            std::string context;
            
        };
        
        // mark the end of the log entry and trigger publishing.
        inline LogItem& ende( LogItem &i ) {
            i.endEntry();
            return i;
        }
        typedef std::shared_ptr<LogItem> LogItemPtr;


    }   // log

}   // redux

#endif // REDUX_LOGGING_LOGITEM_HPP
