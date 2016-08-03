#ifndef REDUX_LOGGING_LOGOUTPUT_HPP
#define REDUX_LOGGING_LOGOUTPUT_HPP

#include "redux/logging/logitem.hpp"

#include <atomic>
#include <deque>
#include <mutex>
#include <vector>


namespace redux {

    namespace logging {


#define LOG_PRINT_HEADER        0x0100
#define LOG_PRINT_TIME          0x0200
#define LOG_PRINT_PRIO          0x0400
#define LOG_PRINT_COLORS        0x0800
#define LOG_PRINT_ORIGIN        0x1000
#define LOG_PRINT_FUNC          0x2000

#define FILE_MODE_READONLY      (uint8_t)0
#define FILE_MODE_INCREMENT     (uint8_t)1
#define FILE_MODE_APPEND        (uint8_t)2
#define FILE_MODE_REPLACE       (uint8_t)3


#ifdef LOG_UPTO
#undef LOG_UPTO
#endif
#define LOG_UPTO(lvl)   ((1<<lvl)-1)  // all levels up to lvl


        class LogOutput {
            
            unsigned int flushPeriod;
            void maybeFlush( void );
            
        protected:

            std::mutex flushMutex;

            uint8_t mask;
            std::mutex queueMutex;
            std::deque<LogItemPtr> itemQueue;
            std::atomic<unsigned int> itemCount;
            
            std::string name_;

            virtual void flushBuffer( void ) {};
            
            LogOutput( uint8_t m=LOG_MASK_ANY, unsigned int flushPeriod=1 );

        public:
            virtual ~LogOutput();

            
            inline void setName( const std::string& n ) { name_ = n; }
            inline std::string name(void) const { return name_; }
            
            inline void setLevel( uint8_t l ) { mask = LOG_UPTO(l); }
            inline void setMask( uint8_t m ) { mask = m; }
            inline uint8_t getMask(void) const { return mask; }

            void addItem( LogItemPtr );
            void addItems( const std::vector<LogItemPtr>& );

            friend class Logger;

        };
        typedef std::shared_ptr<LogOutput> LogOutputPtr;




    } // end namespace logging

} // end namespace redux



#endif // REDUX_LOGGING_LOGOUTPUT_HPP
