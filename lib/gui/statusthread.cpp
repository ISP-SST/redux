#include "redux/gui/statusthread.hpp"

#include "redux/network/protocol.hpp"

/*
#include <REDUX/Common.h>
#include <REDUX/DaemonCmd.h>*/

using namespace redux::gui;
using namespace redux;
using namespace std;

StatusThread::StatusThread( Job::JobSet& j ) : isRunning( false ), isConnected( false ), jobs(j) {

}

void StatusThread::run() {

    isRunning = true;
    //uint count = 0;
    //uint8_t cmd;
    while( isRunning ) {
        usleep( 2000 );
//         if( ! isConnected ) {
//             for( int i = 0; i < myNet->hosts.size(); ++i ) {
//                 try {
//                     //cout << "Trying to connect to:  " << myNet->hosts[i].hostName << ":" << myNet->hosts[i].listeningPort << endl;
//                     conn = new Socket( myNet->hosts[i]->hostName, myNet->hosts[i]->listeningPort );
//                     usleep( 5000 );
//                     *conn << CMD_CONNECT;
//                     *conn >> cmd;
//                     if( cmd != CMD_OK ) {
//                         HostInfo mstInfo, myInfo;
//                         myInfo.hostType = TYPE_MASTER;
//                         struct sockaddr_in addr;
//                         socklen_t len = sizeof( addr );
//                         int ret = getsockname( conn->getFD(), ( struct sockaddr* )&addr, &len );
//                         myInfo.IP = ntohl( addr.sin_addr.s_addr );
//                         //*mySocket << myInfo;
//                         //*mySocket >> mstInfo;
//                         myInfo.send( *mySocket );
//                         mstInfo.recieve( *mySocket );
//                         *conn >> cmd;
//                     }
//                     while( isRunning ) {
//                         cmd = count++ % 256;
//                         cout << "cmd: " << ( int )cmd << " " << cmdString( cmd ) << endl;
//                         if( cmd == CMD_DISCONNECT || cmd == CMD_DIE ) {
//                         }
//                         else if( cmd == CMD_START || cmd == CMD_STOP || cmd == CMD_RESTART ) {
//                             *conn << cmd;
//                             *conn << TYPE_ALL;
//                         }
//                         else if( cmd == CMD_ADD_THREAD || cmd == CMD_DEL_THREAD ) {
//                             *conn << cmd;
//                             *conn << TYPE_SLAVE;
//                             *conn << ( byte )6;
//                         }
//                         else *conn << cmd;
//                         usleep( 50000 );
//                     }
//                     *conn << CMD_DISCONNECT;
//                     delete conn;
//                 }
// //            catch(ROYAC::ROYAC_Exception &ex) {
// //               //cout << ex << endl;
// //            }
//                 catch( ... ) {
//                     //cout << "unknown exception thrown." << endl;
//                 }
//
//             }
//         }
    }

}
