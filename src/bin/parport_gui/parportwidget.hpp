#ifndef REDUX_GUI_PARPORTWIDGET_HPP
#define REDUX_GUI_PARPORTWIDGET_HPP

#include <redux/util/parport.hpp>
#include "parportframe.hpp"

#include <chrono>

//#include <QWidget>
//#include <QString>
#include <QtGui>


namespace gui {
    

    class ParPortWidget : public QWidget, public redux::util::ParPort {

        Q_OBJECT

    public:
        ParPortWidget ( QWidget *parent = 0 );
        ~ParPortWidget();
        int y;


    public slots:
        void grabRelease();
        void writeState ( ParPort::regtype& );
        void stateChanged ( uint32_t );
        void pinChanged ( int );
        void dirChanged ( int );
        void irqChanged ( int );
        void poll ( void );
        void receive(void);
        void stop(void) { receiving = false; }
        void editSettings(void);
        void printMsg( const QString msg, QColor c = Qt::black ) { emit print(msg,c); };

    signals:
        void print( const QString msg, QColor c );

    private:

        QGridLayout *mainLayout;
        ParPortFrame *parportFrame;
        QLabel *addressLabel;
        QLineEdit *portAddress;
        QPushButton *grabButton,*settingsButton,*runButton,*receiveButton,*stopButton;
        QCheckBox *dirBox,*irqBox;
        QAction *settingsAct;
        QTimer *pollTimer;

        void createLayout ( void );
        void createActions ( void );
        void readSettings(void);
        void writeSettings(void);

        redux::util::ParPort::regtype cachedRegisters;
        std::chrono::duration<int,std::milli> receiveSleep;
        
        bool showVoltages;
        bool receiving;

    };

}   // gui

#endif // REDUX_GUI_PARPORTWIDGET_HPP
