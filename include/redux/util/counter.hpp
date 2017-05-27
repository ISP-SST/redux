#ifndef REDUX_UTIL_COUNTER_HPP
#define REDUX_UTIL_COUNTER_HPP

#include <condition_variable>
#include <mutex>
#include <string>

#define TRIGGER_OUTSIDE_RANGE   0x000100
#define TRIGGER_INSIDE_RANGE    0x000200
#define TRIGGER_ON_RANGE        0x000400
#define TRIGGER_RELEASES_ALL    0x000800
#define INCREMENT_ON_WAKEUP     0x001000
#define DECREMENT_ON_WAKEUP     0x002000
#define RESET_ON_TRIGGER        0x004000
#define WRAP_ON_TRIGGER         0x008000
#define LOCK_DURING_CALLBACK    0x010000

namespace redux {

    namespace util {

        class Counter {

            //uint32_t m_Flags;
            typedef void ( *CB_t )( int64_t, int64_t, int64_t );
            CB_t m_CallBack;
            //string m_Name;

            /*! @var Conditional m_Cond
            *  @brief pthread conditional. Used for signaling.  */
            std::condition_variable m_Cond, m_Tick;
            /*! @var Mutex m_Mutex
                @brief pthread mutex. Providing thread-safety for local variables. */
            std::mutex m_Mutex;

            /*********** Local variables, protected by m_Mutex ***********/
            int64_t m_Value, m_StartValue, m_Range[2];
            uint32_t m_Settings;
            /**************************************************************/

            bool conditionMet( void );

        public:
            bool check( void );

            //typedef void (*callback_t)( Socket* );
            /********** Constructors/Destructor **********/
            Counter( int64_t high = 100, int64_t low = 0, int64_t init = 0, std::string nm = "" );
            ~Counter();
            /*********************************************/

            /*********** run() is the entry-point for the threadpool, overloading PoolTask::run() ***********/
            void run( void );
            /************************************************************************************************/

            /********* Methods to access local variables via jobMutex *********/
            void set( int64_t high, int64_t low, int64_t init );
            void set( uint32_t f );
            void unset( uint32_t f ) { m_Settings &= ~f; }
            uint32_t settings(void) const { return m_Settings; }

            void add( int64_t );
            void subtract( int64_t );
            bool compare( int64_t );
            //void lock( void );
            //void unlock( void );
            int64_t high( void );
            int64_t low( void );
            int64_t init( void );
            int64_t range( void );
            void setRange( int64_t high, int64_t low );
            int64_t value( void );
            void setValue( int64_t );
            //std::string name( void );
            //void setName( std::string );
            /******************************************************************/

            void setCallBack( CB_t );

            void trigger( void );
            void wait( void );
            void waitForChange( void );

            virtual size_t size( void );
            virtual char* pack( char* ptr );
            virtual char* unpack( char* ptr, bool swap );

            Counter& operator++( int );
            Counter& operator--( int );

        };

       /*! @class Counter counter.hpp "redux/counter.hpp"
        *  @brief Counter-class to be used for triggering, waiting for certain values and similar.
        *  @details   Providing a "wait()" function and a callback-facility to allow triggering when a specified target-value is reached.
        *  @author    Tomas Hillberg (hillberg@astro.su.se)
        *  @date      2011
        */


    }  // end namespace util

}  // end namespace redux


#endif // REDUX_UTIL_COUNTER_HPP
