#ifndef REDUX_NETWORK_PROTOCOL_HPP
#define REDUX_NETWORK_PROTOCOL_HPP

#include <cstdint>

namespace redux {

    namespace network {

        enum Command : uint8_t { CMD_OK = 0,       // NB: specific-type enums are supported by gcc/clang, but may not behave well with other compilers
                                 CMD_CONNECT,      //   check this !!! /THI
                                 CMD_NEW_JOB,
                                 CMD_MOD_JOB,
                                 CMD_DEL_JOB,
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
                                 CMD_ERR = 255
                               };
    }   // network
}   // redux

#endif // REDUX_NETWORK_PROTOCOL_HPP

/*

#define CMD_OK                      (uint8_t)0
#define CMD_NONE                    CMD_OK
#define CMD_ERR                     (uint8_t)0xFF

#define CMD_DIE                     (uint8_t)0xF

#define CMD_CONNECT                 (uint8_t)0x1
#define CMD_START                   (uint8_t)0x2
#define CMD_STOP                    (uint8_t)0x3
#define CMD_RESTART                 (uint8_t)0x4
#define CMD_RENICE                  (uint8_t)0x5
#define CMD_STATUS                  (uint8_t)0x6
#define CMD_SET_LOG                 (uint8_t)0x7
#define CMD_DUMP_LOG                (uint8_t)0x8
#define CMD_DISCONNECT              (uint8_t)0x9
//#define CMD_DIE                     0x1E

// "Slave" commands, 4-bit mask, i.e. space for 15 non-zero commands.
#define CMD_ADD_THREAD              (uint8_t)0x21
#define CMD_DEL_THREAD              (uint8_t)0x22
#define CMD_RUN                     (uint8_t)0x23
#define CMD_PAUSE                   (uint8_t)0x24
#define CMD_DYNAMIC                 (uint8_t)0x25

// "FileServer" commands, 4-bit mask, i.e. space for 15 non-zero commands.
//#define CMD_START_LOG               (uint8_t)0x41     // reserved for log. cf. Log.h
//#define CMD_STOP_LOG                (uint8_t)0x42     // reserved for log. cf. Log.h
#define CMD_OPEN_FILE               (uint8_t)0x43
#define CMD_CLOSE_FILE              (uint8_t)0x44
#define CMD_GET_BLOCK               (uint8_t)0x45
#define CMD_PUT_BLOCK               (uint8_t)0x46
#define CMD_GET_FILE                (uint8_t)0x47
#define CMD_PUT_FILE                (uint8_t)0x48
#define CMD_CREATE_FILE             (uint8_t)0x49
#define CMD_DELETE_FILE             (uint8_t)0x4A
#define CMD_VERIFY                  (uint8_t)0x4B

// "Master" commands, 4-bit mask, i.e. space for 15 non-zero commands.
#define CMD_CFG                     (uint8_t)0x80
#define CMD_GET_JOB                 (uint8_t)0x81
#define CMD_GET_JOBS                (uint8_t)0x82
#define CMD_ADD_JOB                 (uint8_t)0x83
#define CMD_DEL_JOB                 (uint8_t)0x84
#define CMD_MOD_JOB                 (uint8_t)0x85
#define CMD_START_JOB               (uint8_t)0x86
#define CMD_STOP_JOB                (uint8_t)0x87
#define CMD_GET_SYS                 (uint8_t)0x88
#define CMD_ADD_HOST                (uint8_t)0x89
#define CMD_DEL_HOST                (uint8_t)0x8A
#define CMD_SIGNOFF                 (uint8_t)0x8B
#define CMD_LAUNCH                  (uint8_t)0x8C
#define CMD_KILL                    (uint8_t)0x8D
#define CMD_NEXT_PART               (uint8_t)0x8E
*/