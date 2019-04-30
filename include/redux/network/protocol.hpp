#ifndef REDUX_NETWORK_PROTOCOL_HPP
#define REDUX_NETWORK_PROTOCOL_HPP

#include <cstdint>
#include <string>

namespace redux {

    namespace network {

        enum Command : uint8_t { CMD_OK = 0,       // NB: specific-type enums are supported by gcc/clang, but may not behave well with other compilers
                                 CMD_NOTICE,
                                 CMD_CONNECT,      //   check this !!! /THI
                                 CMD_ADD_JOB,
                                 CMD_MOD_JOB,
                                 CMD_DEL_JOB,
                                 CMD_GET_WORK,
                                 CMD_GET_JOBLIST,
                                 CMD_PUT_PARTS,
                                 CMD_STAT,
                                 CMD_JSTAT,
                                 CMD_PSTAT,
                                 CMD_SLV_CFG,
                                 CMD_SLV_IO,
                                 CMD_SLV_RES,
                                 CMD_SLV_REJ,
                                 CMD_DEL_SLV,
                                 CMD_AUTH,
                                 CMD_CFG,
                                 CMD_DISCONNECT,
                                 CMD_LOG_CONNECT,
                                 CMD_PUT_LOG,
                                 CMD_RESET,
                                 CMD_FORCE_RESET,
                                 CMD_WAKE,
                                 CMD_DIE,
                                 CMD_EXIT,
                                 CMD_INTERACTIVE,
                                 CMD_LISTEN,
                                 CMD_ERR = 255
                               };
                               
        uint8_t cmdFromString( const std::string& str );
        std::string cmdToString( uint8_t cmd );
        
    }   // network
}   // redux

#endif // REDUX_NETWORK_PROTOCOL_HPP
