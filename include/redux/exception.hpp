#ifndef REDUX_EXCEPTION_HPP
#define REDUX_EXCEPTION_HPP

#include <exception>
#include <sstream>
#include <string>

namespace redux {

    
    /*!  @file      exception.hpp
     *   @details   Any specialized exceptions should be derived from these ones rather than
     *              the std exceptions directly, so that global error handling/reporting can
     *              be done on the library level.
     *   @name      Exceptions
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2013
     */

    /*!  @class     Exception
     *   @brief     Base class for all exceptions.
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2013
     */
    class Exception : public std::exception {
    public:
        Exception( void ) : message( "Exception" ) {}
        Exception( const std::string &message ) : message( message ) {}
        virtual ~Exception( void ) throw() {}

        virtual const char *what( void ) const throw() {
            return message.c_str();
        }

    private:
        const std::string message; // Only used when constructed with a string.
    };


    /*!  @class     RecoverableException
     *   @brief     Exception that may have a temporary cause. The failed operation might work upon a retry.
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2013
     */
    class RecoverableException : public Exception {
    public:
        RecoverableException( void ) : Exception( "RecoverableException" ) {}
        RecoverableException( const std::string &message ) : Exception( message ) {}
        virtual ~RecoverableException( void ) throw() {}
    };


    /*!  @class     UnrecoverableException
     *   @brief     Exception that will rise again if the failed operation is retried.
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2013
     */
    class UnrecoverableException : public Exception {
    public:
        UnrecoverableException( void ) : Exception( "UnrecoverableException" ) {}
        UnrecoverableException( const std::string &message ) : Exception( message ) {}
        virtual ~UnrecoverableException( void ) throw() {}
    };


    /*!  @class     ProgramError
     *   @brief     Program error.
     *   @details   Exception that indicates that the program has reached a
     *              state it should not be able to reach. This generally means
     *              there is a bug.
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2013
     */
    class ProgramError : public Exception {
    public:
        ProgramError( void ) : Exception( "ProgramError" ) {}
        ProgramError( const std::string &message ) : Exception( message ) {}
        virtual ~ProgramError( void ) throw() {}
    };


    /*!  @class     BadArgument
     *   @brief     Exception that indicates that something has been invoked with
     *              bad arguments.
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2013
     */
    class BadArgument : public UnrecoverableException {
    public:
        BadArgument( void ) : UnrecoverableException( "BadArgument" ) {}
        BadArgument( const std::string &message ) : UnrecoverableException( message ) {}

        virtual ~BadArgument( void ) throw() {}
    };


    /*!  @class     IndexOutOfBounds
     *   @brief     Exception that indicates that an index argument was out of bounds.
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2013
     */
    class IndexOutOfBounds : public BadArgument {
    public:
        IndexOutOfBounds( size_t index, size_t maxAllowed )
            : index( index ), maxAllowed( maxAllowed ) {
            std::ostringstream ss;
            ss << "Index out of bounds (" << index << " > " << maxAllowed
               << ")";
            message = ss.str();
        }

        virtual ~IndexOutOfBounds( void ) throw() {}

        size_t getIndex( void ) const throw() {
            return index;
        }
        size_t getMaxAllowed( void ) const throw() {
            return maxAllowed;
        }

        virtual const char *what( void ) const throw() {
            return message.c_str();
        }

    private:
        size_t index;
        size_t maxAllowed;
        std::string message;
    };


    /*!  @class     NotImplemented
     *   @brief     Exception that indicates that the operation is not (yet) implemented.
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2013
     */
    class NotImplemented : public UnrecoverableException {
    public:
        NotImplemented( void ) : UnrecoverableException( "Operation is not implemented" ) {}
        NotImplemented( const std::string &message ) : UnrecoverableException( message ) {}
        virtual ~NotImplemented( void ) throw() {}
    };


    /*!  @class     NotSupported
     *   @brief     Exception that indicates that the operation is not supported on the class in question.
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2013
     */
    class NotSupported : public ProgramError {
    public:
        NotSupported( void ) : ProgramError( "Operation is not supported" ) {}
        NotSupported( const std::string &message ) : ProgramError( message ) {}
        virtual ~NotSupported( void ) throw() {}
    };



}

#endif // REDUX_EXCEPTION_HPP
