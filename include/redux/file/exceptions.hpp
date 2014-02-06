#ifndef REDUX_FILE_EXCEPTIONS_HPP
#define REDUX_FILE_EXCEPTIONS_HPP

#include "redux/exception.hpp"

namespace redux {

    namespace file {

        class DataIOException : public redux::Exception {
        public:
            DataIOException ( const std::string &msg ) : redux::Exception ( msg ) {}
            DataIOException ( const char* message = "DataIOException" ) : redux::Exception ( message ) {}
            virtual ~DataIOException ( void ) throw () {}
        };
    }

}

#endif // REDUX_FILE_EXCEPTIONS_HPP
