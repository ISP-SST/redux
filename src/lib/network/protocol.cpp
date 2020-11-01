#include "redux/network/protocol.hpp"

#include <redux/util/stringutil.hpp>

#include <boost/lexical_cast.hpp>

using namespace redux::util;
using namespace std;


uint8_t redux::network::cmdFromString( const string& str ) {
    
    try {
        uint8_t cmd = boost::lexical_cast<uint8_t>( str );
        return cmd;
    } catch( const boost::bad_lexical_cast& e ) {
            // ignore fail, try parsing as string
    }
    
    if( contains(str, "OK", true ) ) return CMD_OK;
    if( contains(str, "NOTICE", true ) ) return CMD_NOTICE;
    if( contains(str, "CONNECT", true ) ) return CMD_CONNECT;
    if( contains(str, "ADD_JOB", true ) ) return CMD_ADD_JOB;
    if( contains(str, "MOD_JOB", true ) ) return CMD_MOD_JOB;
    if( contains(str, "DEL_JOB", true ) ) return CMD_DEL_JOB;
    if( contains(str, "GET_WORK", true ) ) return CMD_GET_WORK;
    if( contains(str, "GET_JOBLIST", true ) ) return CMD_GET_JOBLIST;
    if( contains(str, "PUT_PARTS", true ) ) return CMD_PUT_PARTS;
    if( contains(str, "JSTAT", true ) ) return CMD_JSTAT;
    if( contains(str, "PSTAT", true ) ) return CMD_PSTAT;
    if( contains(str, "STAT", true ) ) return CMD_STAT;
    if( contains(str, "SLV_CFG", true ) ) return CMD_SLV_CFG;
    if( contains(str, "SLV_IO", true ) ) return CMD_SLV_IO;
    if( contains(str, "SLV_RES", true ) ) return CMD_SLV_RES;
    if( contains(str, "SLV_REJ", true ) ) return CMD_SLV_REJ;
    if( contains(str, "DEL_SLV", true ) ) return CMD_DEL_SLV;
    if( contains(str, "AUTH", true ) ) return CMD_AUTH;
    if( contains(str, "CFG", true ) ) return CMD_CFG;
    if( contains(str, "DISCONNECT", true ) ) return CMD_DISCONNECT;
    if( contains(str, "LOG_CONNECT", true ) ) return CMD_LOG_CONNECT;
    if( contains(str, "PUT_LOG", true ) ) return CMD_PUT_LOG;
    if( contains(str, "RESET", true ) ) return CMD_RESET;
    if( contains(str, "DIE", true ) ) return CMD_DIE;
    
    return CMD_ERR;
    
     
}


string redux::network::cmdToString( uint8_t cmd ) {
    
    switch( cmd ) {
        case CMD_OK: return "CMD_OK";
        case CMD_NOTICE: return "CMD_NOTICE";
        case CMD_CONNECT: return "CMD_CONNECT";
        case CMD_ADD_JOB: return "CMD_ADD_JOB";
        case CMD_MOD_JOB: return "CMD_MOD_JOB";
        case CMD_DEL_JOB: return "CMD_DEL_JOB";
        case CMD_GET_WORK: return "CMD_GET_WORK";
        case CMD_GET_JOBLIST: return "CMD_GET_JOBLIST";
        case CMD_PUT_PARTS: return "CMD_PUT_PARTS";
        case CMD_STAT: return "CMD_STAT";
        case CMD_JSTAT: return "CMD_JSTAT";
        case CMD_PSTAT: return "CMD_PSTAT";
        case CMD_SLV_CFG: return "CMD_SLV_CFG";
        case CMD_SLV_IO: return "CMD_SLV_IO";
        case CMD_SLV_RES: return "CMD_SLV_RES";
        case CMD_SLV_REJ: return "CMD_SLV_REJ";
        case CMD_DEL_SLV: return "CMD_DEL_SLV";
        case CMD_AUTH: return "CMD_AUTH";
        case CMD_CFG: return "CMD_CFG";
        case CMD_DISCONNECT: return "CMD_DISCONNECT";
        case CMD_LOG_CONNECT: return "CMD_LOG_CONNECT";
        case CMD_PUT_LOG: return "CMD_PUT_LOG";
        case CMD_RESET: return "CMD_RESET";
        case CMD_DIE: return "CMD_DIE";
        case CMD_ERR: return "CMD_ERR";
        default: return "unrecognized";
    }
    
}
