#include "redux/util/counter.hpp"

#include "redux/util/endian.hpp"

// #include "redux/util/datatools.hpp"
// #include "redux/log.hpp"
// #include "redux/threadpool.hpp"

#include <string.h>

using namespace redux::util;
using namespace std;

/*! @fn Counter::Counter( int64_t high, int64_t low, int64_t init )
  * @brief Constructor
  * @param high Upper limit of range.
  * @param low Lower limit of range
  * @param init Initial value
  */
Counter::Counter( int64_t high, int64_t low, int64_t init, string nm ) : m_CallBack( NULL ), m_Value( init ), m_StartValue( init ) { //, m_Name(nm) {

////    setFlag( TRIGGER_OUTSIDE_RANGE | TRIGGER_RELEASES_ALL | RESET_ON_TRIGGER );

    m_Range[0] = low;
    m_Range[1] = high;

//    Msg( nm + string("  arg(")+toString(high)+string(",")+toString(low)+string(",")+toString(init)+string(") (")+hexString(this)+string(")"), MASK_DEBUG, m_TaskLogMask );
//m_Name = nm;
    //setTaskType(TT_OTHER);

}

Counter::~Counter() {

    //myCout("Counter::~Counter(1) "+ name()+string(" ")+hexString(this));
    //Msg( name()+string("  (")+hexString(this)+string(")"), MASK_ERR, m_TaskLogMask );
    //setExitAction(EA_IGNORE);
    //unsetConfig( TP_REINSERT_ON_COMPLETION|TP_DELETE_ON_COMPLETION );
    //Msg( string("2(")+hexString(this)+string(")"), MASK_DEBUG, m_TaskLogMask );
    unique_lock<mutex> lock(m_Mutex);

    m_CallBack = NULL;
    m_Cond.notify_all();
    
    //exit();
    //Msg( string("3(")+hexString(this)+string(")"), MASK_DEBUG, m_TaskLogMask );
    //m_Cond.broadcast();
    //Msg( name()+string("  4(")+hexString(this)+string(")"), MASK_ERR, m_TaskLogMask );
    //PoolTask::waitForState( TP_STATE_DONE );
    //Msg( name()+string("  Done. (")+hexString(this)+string(")"), MASK_ERR, m_TaskLogMask );
}


/*! @fn void Counter::run( void )
  * @brief Entry-point for thread.
  * @details Wait on a conditional until the condition is met, then do callback
  */
void Counter::run( void ) {
/*
    //m_Mutex.lock();

    int64_t v, r[2];
    uint32_t f;
    bool lockDuring;
    CB_t tmp;
    {
        unique_lock<mutex> lock(m_Mutex);
        tmp = m_CallBack;
    }
    //Msg( dumpArray( m_Range, 2, "range" )+string(" value=")+toString(m_Value)+string("   (")+hexString(this)+string(")"), MASK_ERR, m_TaskLogMask );
    while ( tmp && !checkState( TP_STATE_EXIT ) ) {
        m_Mutex.lock();
        Msg( name() + dumpArray( m_Range, 2, "  range" ) + string( " value=" ) + toString( m_Value ) + string( "   (" ) + hexString( this ) + string( ")" ), MASK_DEBUG, getLogMask() );

        if ( !conditionMet() )
            m_Cond.wait( m_Mutex );

        v = m_Value;
        f = flags();
        lockDuring = f & LOCK_DURING_CALLBACK;
        r[0] = m_Range[0];
        r[1] = m_Range[1];

        if ( f & RESET_ON_TRIGGER )
            m_Value = m_StartValue;
        else
            if ( f & WRAP_ON_TRIGGER ) {
                if ( ( m_Value - m_Range[1] ) > ( m_Value - m_Range[0] ) )
                    m_Value = m_Range[1];
                else
                    m_Value = m_Range[0];
            }

        if ( f & INCREMENT_ON_WAKEUP ) {
            m_Value++;
            //check();
        }
        else
            if ( f & DECREMENT_ON_WAKEUP ) {
                m_Value--;
                //check();
            }

        tmp = m_CallBack;

        if ( tmp ) {
            Msg( name() + string( "  Callback=(" ) + hexString( tmp ) + string( ")  (" ) + hexString( this ) + string( ")" ), MASK_DEBUG2, getLogMask() );
        }
        else {
            Msg( name() + string( "  No callback defined when triggered !!!  (" ) + hexString( this ) + string( ")" ), MASK_WARNING, getLogMask() );
        }

        if ( tmp &&  lockDuring )
            tmp( v, r[1], r[0] );

        m_Mutex.unlock();

        if ( tmp && !lockDuring )
            tmp( v, r[1], r[0] );

        //m_Mutex.lock();
        //Msg( string("Triggered:  ")+dumpArray( r, 2, "range" )+string(" value=")+toString(v)+string("   (")+hexString(this)+string(")"), MASK_DETAIL, m_TaskLogMask );
        //Msg(  string("Done:  ")+dumpArray( m_Range, 2, "range" )+string(" value=")+toString(m_Value)+string("   (")+hexString(this)+string(")"), MASK_ERR, m_TaskLogMask );
        //m_Mutex.unlock();

    }

    //Msg( name()+ string("   EXITED:  ")+dumpArray( r, 2, "range" )+string(" value=")+toString(v)+string("   (")+hexString(this)+string(")"), MASK_EMERG, m_TaskLogMask );
    //m_Mutex.unlock();
*/
}

/*! @fn bool Counter::conditionMet( void )
  * @brief Check if the condition is met.
  */
bool Counter::conditionMet( void ) {

    uint32_t f( 0 ); //flags() );
    return ( ( ( f & TRIGGER_OUTSIDE_RANGE ) && ( ( m_Value   < m_Range[0] ) || ( m_Value   > m_Range[1] ) ) ) ||
             ( ( f & TRIGGER_INSIDE_RANGE )  && ( ( m_Value   > m_Range[0] ) && ( m_Value   < m_Range[1] ) ) ) ||
             ( ( f & TRIGGER_ON_RANGE )      && ( ( m_Value  == m_Range[0] ) || ( m_Value  == m_Range[1] ) ) ) );

}

/*! @fn bool Counter::check( void )
  * @brief Check if the condition is met, trigger if true;
  */
bool Counter::check( void ) {

    //Msg( dumpArray( m_Range, 2, "range" )+string(" value=")+toString(m_Value)+string("   (")+hexString(this)+string(")"), MASK_ERR, m_TaskLogMask );
    uint32_t f( 0 ); //flags() );
    //Msg( name()+string("   flags=")+bitString(f), MASK_ERR, m_TaskLogMask );
    bool ret(false);
    {
        unique_lock<mutex> lock(m_Mutex);
        ret = ( ( f & TRIGGER_OUTSIDE_RANGE ) && ( ( m_Value   < m_Range[0] ) || ( m_Value   > m_Range[1] ) ) ) ||
              ( ( f & TRIGGER_INSIDE_RANGE )  && ( ( m_Value   > m_Range[0] ) && ( m_Value   < m_Range[1] ) ) ) ||
              ( ( f & TRIGGER_ON_RANGE )      && ( ( m_Value  == m_Range[0] ) || ( m_Value  == m_Range[1] ) ) );
    }

    if ( ret ) {
        if ( f & TRIGGER_RELEASES_ALL )
            m_Cond.notify_all();
        else
            m_Cond.notify_one();
    }

    //Msg( name()+string("   flags2=")+bitString(f), MASK_ERR, m_TaskLogMask );
    return ret;
}


void Counter::set( int64_t high, int64_t low, int64_t init ) {
    {
        unique_lock<mutex> lock(m_Mutex);
        //if ( (low != m_Range[0]) || (high != m_Range[1]) ) m_Tick.broadcast();
        m_Range[0] = low;
        m_Range[1] = high;
        m_Value = m_StartValue = init;
    }
    check();
}


void Counter::set( uint32_t f ) {
    //m_Mutex.lock();
    if ( (f&WRAP_ON_TRIGGER) && (f&RESET_ON_TRIGGER) ) {
////        Msg( string( "Incompatible flags: WRAP_ON_TRIGGER & RESET_ON_TRIGGER   (" ) + hexString( this ) + string( ")" ), MASK_ERR, getLogMask() );
        return;
    }
    
    if ( (f&INCREMENT_ON_WAKEUP) && (f&DECREMENT_ON_WAKEUP) ) {
////            Msg( string( "Incompatible flags: INCREMENT_ON_WAKEUP & DECREMENT_ON_WAKEUP   (" ) + hexString( this ) + string( ")" ), MASK_ERR, getLogMask() );
        return;
    }
        
    m_Settings |= f;
////            unsetFlag( 0xFFFFFFFF );
////            setFlag( t );


    //m_Mutex.unlock();
}


void Counter::add( int64_t v ) {
    {
        unique_lock<mutex> lock(m_Mutex);

        if ( v ) {
            //m_Tick.broadcast();
            m_Value += v;
        }

    }
    check();
}


void Counter::subtract( int64_t v ) {
    {
        unique_lock<mutex> lock(m_Mutex);

        if ( v ) {
            //m_Tick.broadcast();
            m_Value -= v;
        }

    }
    check();
}

/*! @fn bool Counter::compare( int64_t v )
  * @brief Compare v to current value.
  */
bool Counter::compare( int64_t v ) {
    unique_lock<mutex> lock(m_Mutex);
    return (m_Value == v);
}

/*! @fn void Counter::lock( void )
  * @brief Compare v to current value.
  */
//void Counter::lock( void ) {
//    m_Mutex.lock();
//}

/*! @fn void Counter::unlock( void )
  * @brief Compare v to current value.
  */
//void Counter::unlock( void ) {
//    m_Mutex.unlock();
//}


/*! @fn int64_t Counter::high( void )
  * @brief Get value for upper limit of counter.
  */
int64_t Counter::high( void ) {
    unique_lock<mutex> lock(m_Mutex);
    return m_Range[1];
}

/*! @fn int64_t Counter::low( void )
  * @brief Get value for lower limit of counter.
  */
int64_t Counter::low( void ) {
    unique_lock<mutex> lock(m_Mutex);
    return m_Range[0];
}

/*! @fn int64_t Counter::init( void )
  * @brief Get initial-value for the counter.
  */
int64_t Counter::init( void ) {
    unique_lock<mutex> lock(m_Mutex);
    return m_StartValue;
}

/*! @fn int64_t Counter::range( void )
  * @brief Get target-value.
  */
int64_t Counter::range( void ) {
    unique_lock<mutex> lock(m_Mutex);
    return ( m_Range[1] - m_Range[0] );
}

/*! @fn void Counter::setRange( int64_t high, int64_t low )
  * @brief Specify range of counter.
  */
void Counter::setRange( int64_t high, int64_t low ) {
    {
        unique_lock<mutex> lock(m_Mutex);
        m_Range[0] = low;
        m_Range[1] = high;
    }
    check();
}

/*! @fn int64_t Counter::value( void )
  * @brief Get current value.
  */
int64_t Counter::value( void ) {
    unique_lock<mutex> lock(m_Mutex);
    return m_Value;
}

/*! @fn void Counter::setValue( int64_t v )
  * @brief Set current value to v.
  */
void Counter::setValue( int64_t v ) {
    {
        unique_lock<mutex> lock(m_Mutex);
        //Msg( name()+string("  oldValue=")+toString(m_Value)+string("  newValue=")+toString(v), MASK_DEBUG2, m_TaskLogMask );
        if ( m_Value != v ) {
            //m_Tick.broadcast();
            m_Value = v;
        }

    }
    check();
}


/*! @fn string Counter::name( void )
  * @brief Get name.
  */
/*string Counter::name( void ) {
    m_Mutex.lock();
    string ret(m_Name);
    m_Mutex.unlock();
    return ret;
}*/

/*! @fn void Counter::setName( string n )
  * @brief Set name to n.
  */
/*void Counter::setName( string n ) {
    m_Mutex.lock();
    m_Name = n;
    m_Mutex.unlock();
}*/


void Counter::setCallBack( CB_t cb ) {
    unique_lock<mutex> lock(m_Mutex);
    m_CallBack = cb;
}


/*! @fn void Counter::trigger( void )
  * @brief Release waiting threads..
  */
void Counter::trigger( void ) {
    uint32_t f( 0 ); //flags() );

    //cout << (string("Counter::trigger()  f=") +bitString(f)) << endl;
    if ( f & TRIGGER_RELEASES_ALL )
        m_Cond.notify_all();
    else
        m_Cond.notify_one();

    //m_Tick.broadcast();
}

/*! @fn void Counter::wait( void )
  * @brief Wait on a conditional until counter gets triggered.
  */
void Counter::wait( void ) {

    // myCout("Counter::wait(1) "+ name()+string(" ")+hexString(this));
    unique_lock<mutex> lock(m_Mutex);

    if ( !conditionMet() ) {
        //cout << (string("Counter::wait()  about to wait...")) << endl;
        //Msg( string(" Waiting... ")+toString(m_Value)+string("  ")+dumpArray(m_Range,2,"range"), MASK_DEBUG, getLogMask() );
        m_Cond.wait( lock );
    }

    // myCout("Counter::wait(2)");
    //Msg( string(" Triggered: value=")+toString(m_Value)+string("  ")+dumpArray(m_Range,2,"range"), MASK_DEBUG, getLogMask() );

    uint32_t f( 0 ); //flags() );

    //cout << (string("Counter::wait()  after wait:  f=") +bitString(f)) << endl;
    //Msg( name()+string("  check=true  flags=")+bitString(f), MASK_EMERG, m_TaskLogMask );
    if ( f & INCREMENT_ON_WAKEUP ) {
        m_Value++;
        //check();
    }
    else
        if ( f & DECREMENT_ON_WAKEUP ) {
            m_Value--;
            //check();
        }

    // myCout("Counter::wait(E)");
    //cout << (string("Counter::wait()  returning.")) << endl;
}

/*! @fn void Counter::waitForChange( void )
  * @brief Wait on a conditional until counter gets triggered.
  */
void Counter::waitForChange( void ) {

    unique_lock<mutex> lock(m_Mutex);
    //Msg( name()+string("  Before tick: value=")+toString(m_Value)+string("  ")+dumpArray(m_Range,2,"range"), MASK_EMERG, m_TaskLogMask );
    m_Tick.wait( lock );
    //Msg( name()+string("  After tick:  value=")+toString(m_Value)+string("  ")+dumpArray(m_Range,2,"range"), MASK_EMERG, m_TaskLogMask );

}


size_t Counter::size( void ) {

    return 4 * sizeof( int64_t ); //int64_t m_Value, m_StartValue, m_Range[2];

}

char* Counter::pack( char* ptr ) {

    unique_lock<mutex> lock(m_Mutex);
    *( int64_t* )ptr = m_Value;
    ptr += sizeof( int64_t );
    *( int64_t* )ptr = m_StartValue;
    ptr += sizeof( int64_t );
    size_t sz = 2 * sizeof( int64_t );
    memcpy( ptr, m_Range, sz );
    ptr += sz;

    return ptr;

}

char* Counter::unpack( char* ptr, bool swap ) {

    unique_lock<mutex> lock(m_Mutex);
    m_Value = *( int64_t* )ptr;
    ptr += sizeof( int64_t );
    m_StartValue = *( int64_t* )ptr;
    ptr += sizeof( int64_t );
    size_t sz = 2 * sizeof( int64_t );
    memcpy( m_Range, ptr, sz );
    ptr += sz;

    if ( swap ) {
        swapEndian( m_Value );
        swapEndian( m_StartValue );
        swapEndian( m_Range, 2 );
    }

    return ptr;

}


/*! @fn Counter& Counter::operator++(int)
  * @brief Increment current value by 1.
  */
Counter& Counter::operator++( int ) {
    {
        unique_lock<mutex> lock(m_Mutex);
        m_Value++;
        //Msg( name()+string(":  value=")+toString(m_Value), MASK_EMERG, getLogMask() );
        m_Tick.notify_all();
        //Msg( string("value=")+toString(m_Value), MASK_EMERG, getLogMask() );
    }
    check();
    return *this;
}

/*! @fn Counter& Counter::operator--(int)
  * @brief Decrement current value by 1.
  */
Counter& Counter::operator--( int ) {
     {
        unique_lock<mutex> lock(m_Mutex);
        m_Value--;
        //Msg( name()+string(":  value=")+toString(m_Value), MASK_EMERG, getLogMask() );
        m_Tick.notify_all();
        //Msg( string("value=")+toString(m_Value), MASK_EMERG, getLogMask() );
     }
    check();
    return *this;
}

