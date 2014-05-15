#ifndef REDUX_GUI_STATUSTHREAD_HPP
#define REDUX_GUI_STATUSTHREAD_HPP

#include "redux/job.hpp"
#include "redux/network/tcpconnection.hpp"

#include <vector>

#include <QThread>

namespace redux {

    namespace gui {

        class StatusThread : public QThread {

            Q_OBJECT
            
            bool isRunning, isConnected;
            Job::JobSet& jobs;
            //SysConfig *myNet;
            network::TcpConnection *conn;

        public:
            StatusThread( Job::JobSet& );
            ~StatusThread() {};
            void run();
            void stop() { isRunning = false; };

        };

    }   // gui

}

#endif  // REDUX_GUI_STATUSTHREAD_HPP
