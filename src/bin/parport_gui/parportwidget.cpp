#include "parportwidget.hpp"


#include <redux/util/stringutil.hpp>
//#include <redux/util/timer.hpp>

#include <iostream>
#include <thread>
//#include <cstdio>
//#include <sys/io.h>
//#include <string.h>


using namespace redux::util;
using namespace gui;
using namespace std;


ParPortWidget::ParPortWidget ( QWidget *parent ) : QWidget ( parent ), pollTimer ( new QTimer ( this ) ), cachedRegisters ( registers ),
         receiveSleep(chrono::milliseconds(100)), showVoltages ( false ) {

    if ( parent ) {
        connect ( this, SIGNAL ( print ( QString, QColor ) ), parent, SLOT ( printMsg ( QString, QColor ) ) );
    }

    createLayout();
    createActions();

    readSettings();


    connect ( pollTimer, SIGNAL ( timeout() ), this, SLOT ( poll() ) );
    pollTimer->start();
}


ParPortWidget::~ParPortWidget() {
    writeSettings();
}


void ParPortWidget::grabRelease ( void ) {

    if ( status == PAR_OK ) {
        printMsg ( "Released address (" + QString::number ( baseaddress, 16 ) + ")", Qt::blue );
        pollTimer->stop();
        setBaseAddress ( 0 );
        grabButton->setText ( tr ( "Grab" ) );
        portAddress->setEnabled ( true );
        runButton->setEnabled ( false );
        receiveButton->setEnabled ( false );
        stopButton->setEnabled ( false );
        dirBox->setEnabled ( false );
        irqBox->setEnabled ( false );
        status = PAR_NONE;
        return;
    }

    uint64_t addr = portAddress->text().toULong ( 0, 16 );
    setBaseAddress ( addr );

    if ( status == PAR_OK ) {
        printMsg ( "Acquired address (" + QString::number ( baseaddress, 16 ) + ")", Qt::blue );
        grabButton->setText ( tr ( "Release" ) );
        runButton->setEnabled ( true );
        receiveButton->setEnabled ( true );
        stopButton->setEnabled ( true );
        portAddress->setEnabled ( false );
        dirBox->setEnabled ( true );
        irqBox->setEnabled ( true );
        pollTimer->start();
    } else {
        if ( status == PAR_ERR ) {
            printMsg ( "Could not acquire address (" + QString::number ( addr, 16 ) + ") " + QString ( strerror ( errno ) ), Qt::red );
        }
    }

}


void ParPortWidget::writeState ( ParPort::regtype& state ) {

    for ( int i = 0; i < 3; ++i ) write_register ( i, state.b[i] );

}


void ParPortWidget::stateChanged ( uint32_t diff ) {

    uint32_t state = cachedRegisters.i;
    if ( cachedRegisters.b[3] ) {
        state ^= inv_mask.i;
    }

    uint32_t mask = 1;
    for ( int i = 0; i < 24; ++i ) {
        if ( mask & PAR_MASK_PIN ) {
            parportFrame->setPin ( mask, mask & state );
        }
        mask <<= 1;
    }

}


void ParPortWidget::pinChanged ( int pin ) {
//    printMsg("pin " + QString::number(pin) + " toggled, mask = " + QString::number(ParPort::PIN[pin], 16));
    ParPort::regtype newState = cachedRegisters;
    newState.i ^= ParPort::PIN[pin];
    for ( int i = 0; i < 3; ++i ) {
        if ( newState.b[i] != cachedRegisters.b[i] ) {
            write_register ( i, newState.b[i] );
        }
    }
}


void ParPortWidget::dirChanged ( int val ) {
    ParPort::regtype newState = cachedRegisters;
    if ( val ) {
        //set_control(1 << 5);
        newState.i |= PAR_INPUT_MODE;
    } else {
        //clear_control(1 << 5);
        newState.i &= ~PAR_INPUT_MODE;
    }
    if ( newState.i != cachedRegisters.i ) {
        writeState ( newState );
    }
}


void ParPortWidget::irqChanged ( int val ) {
    ParPort::regtype newState = cachedRegisters;
    if ( val ) {
        //set_control(1 << 4);
        newState.i |= PAR_IRQ_MODE;
    } else {
        //clear_control(1 << 4);
        newState.i &= ~PAR_IRQ_MODE;
    }
    if ( newState.i != cachedRegisters.i ) {
        writeState ( newState );
    }
}


void ParPortWidget::poll ( void ) {

    if ( status != PAR_OK ) return;

    read_registers();

    registers.b[3] = showVoltages;      // use unused byte for flagging
    //registers.b[0]++;
    //registers.b[1]++;
    //registers.b[2]++;

    uint32_t diff = registers.i ^ cachedRegisters.i;
    if ( diff ) {
#ifdef _DEBUG
        printMsg ( "Poll: (status=" + QString::number ( status ) + ")  newval=" + QString::number ( registers.i, 2 ), Qt::blue );
#endif
        cachedRegisters = registers;
        stateChanged ( diff );
    }
}

#define PP_READ_MODE  0x20 // (0x01<<5)
#define PP_NO_RESET   0x1 //(0x01<<0)
#define PP_CLOCK_LOW  0x2 //(0x01<<1)
//#define PP_CLOCK_LOW  (0x0F)
#define PP_NO_DATA    0x40 //(0x01<<6)


void ParPortWidget::receive ( void ) {

    if(receiving) return;
       
    receiving = true;
    std::thread ( [this]() {
        ParPort::regtype value;
        uint8_t i;
        printMsg("Starting CounterBox listener.", Qt::blue);
        write_register(2,PP_READ_MODE);                          // read mode, RESET, clock=high
        this_thread::sleep_for ( chrono::microseconds ( 1 ) );
        write_register(2,PP_READ_MODE|PP_NO_RESET);              // read mode, NO RESET, clock=high
        //cout << "." << flush;
        //this_thread::sleep_for ( chrono::seconds ( 1 ) );
        //cout << "|" << endl;
        do {
            //auto next = std::chrono::steady_clock::now() + settings.interval;

            value.i=0;
            i=0;
            while(receiving && !(read_register(1)&PP_NO_DATA)){               // repeat while more data is available
                write_register(2,PP_READ_MODE|PP_NO_RESET|PP_CLOCK_LOW);      // read mode, NO RESET, clock=low
                //this_thread::sleep_for ( chrono::microseconds ( 10 ) );
                value.b[i++] = read_register(0);                              // read data from the port
                //this_thread::sleep_for ( chrono::milliseconds ( 10 ) );
                write_register(2,PP_READ_MODE|PP_NO_RESET);                   // read mode, NO RESET, clock=high
                //this_thread::sleep_for ( chrono::milliseconds ( 10 ) );
            }
            if (i) printMsg("Received " + QString::number(value.i,16) + "  nBytes="+ QString::number(i));
  
            this_thread::sleep_for ( receiveSleep );
            //std::this_thread::sleep_until(next);
        } while ( receiving );
        printMsg("Stopping CounterBox listener.", Qt::blue);
    } ).detach();

}


void ParPortWidget::editSettings ( void ) {

    QDialog dialog ( this );
    dialog.setWindowTitle ( "Settings" );

    QFormLayout form ( &dialog );

    form.addRow ( new QLabel ( "Settings", &dialog ) );

    QLabel* showLabel = new QLabel ( "Show Voltages", &dialog );
    QCheckBox *showCB = new QCheckBox ( &dialog );
    showCB->setChecked ( showVoltages );
    form.addRow ( showLabel, showCB );

    QLabel* timeoutLabel = new QLabel ( "Polling timeout (ms)", &dialog );
    QLineEdit *timeoutLE = new QLineEdit ( &dialog );
    QValidator *val = new QIntValidator ( 1, 1000000, &dialog );
    timeoutLE->setValidator ( val );
    timeoutLE->setText ( QString::number ( pollTimer->interval() ) );

    form.addRow ( timeoutLabel, timeoutLE );

    QLabel* receiveLabel = new QLabel ( "Receive sleep time (ms)", &dialog );
    QLineEdit *receiveLE = new QLineEdit ( &dialog );
    receiveLE->setValidator ( val );
    receiveLE->setText ( QString::number ( receiveSleep.count() ) );

    form.addRow ( receiveLabel, receiveLE );

    QDialogButtonBox buttonBox ( QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Horizontal, &dialog );
    form.addRow ( &buttonBox );
    QObject::connect ( &buttonBox, SIGNAL ( accepted() ), &dialog, SLOT ( accept() ) );
    QObject::connect ( &buttonBox, SIGNAL ( rejected() ), &dialog, SLOT ( reject() ) );

    if ( dialog.exec() == QDialog::Accepted ) {
        if ( showVoltages != showCB->isChecked() ) {
            showVoltages = showCB->isChecked();
            //stateChanged(0xFFFFFFFF);
        }
        pollTimer->setInterval ( timeoutLE->text().toInt() );
        receiveSleep = chrono::milliseconds(receiveLE->text().toInt());
        pollTimer->start();
    }

}




void ParPortWidget::createLayout() {

    mainLayout = new QGridLayout ( this );

    parportFrame = new ParPortFrame ( this );

    addressLabel = new QLabel ( this );
    addressLabel->setText ( tr ( "Port Address (hex)" ) );
    portAddress = new QLineEdit ( this );

    grabButton = new QPushButton ( this );
    grabButton->setText ( tr ( "Grab" ) );

    settingsButton = new QPushButton ( this );
    settingsButton->setText ( tr ( "Settings" ) );

    runButton = new QPushButton ( this );
    runButton->setText ( tr ( "Run All" ) );
    runButton->setEnabled ( false );

    receiveButton = new QPushButton ( this );
    receiveButton->setText ( tr ( "Receive" ) );
    receiveButton->setEnabled ( false );

    stopButton = new QPushButton ( this );
    stopButton->setText ( tr ( "Stop" ) );
    stopButton->setEnabled ( false );

    dirBox = new QCheckBox ( "DATA IN", this );
    dirBox->setEnabled ( false );

    irqBox = new QCheckBox ( "Pin 10 <-> IRQ", this );
    irqBox->setEnabled ( false );

    mainLayout->addWidget ( addressLabel, 0, 0, 1, 1 );
    mainLayout->addWidget ( portAddress, 0, 1, 1, 1 );
    mainLayout->addWidget ( grabButton, 0, 2, 1, 1 );

    mainLayout->addWidget ( dirBox, 0, 3, 1, 1 );
    mainLayout->addWidget ( irqBox, 0, 4, 1, 1 );

    mainLayout->addWidget ( settingsButton, 0, 5, 1, 1 );
    mainLayout->addWidget ( runButton, 0, 6, 1, 1 );
    
    mainLayout->addWidget ( receiveButton, 0, 7, 1, 1 );
    mainLayout->addWidget ( stopButton, 0, 8, 1, 1 );

    mainLayout->addWidget ( parportFrame, 1, 0, 2, 9 );

    connect ( grabButton, SIGNAL ( clicked() ), this, SLOT ( grabRelease() ) );
    connect ( settingsButton, SIGNAL ( clicked() ), this, SLOT ( editSettings() ) );
    connect ( runButton, SIGNAL ( clicked() ), parportFrame, SLOT ( runAll() ) );
    connect ( receiveButton, SIGNAL ( clicked() ), this, SLOT ( receive() ) );
    connect ( stopButton, SIGNAL ( clicked() ), this, SLOT ( stop() ) );
    connect ( portAddress, SIGNAL ( returnPressed() ), this, SLOT ( grabRelease() ) );
    connect ( dirBox, SIGNAL ( stateChanged ( int ) ), this, SLOT ( dirChanged ( int ) ) );
    connect ( irqBox, SIGNAL ( stateChanged ( int ) ), this, SLOT ( irqChanged ( int ) ) );


}

void ParPortWidget::createActions ( void ) {

    settingsAct = new QAction ( tr ( "Settings" ), this );
    connect ( settingsAct, SIGNAL ( triggered() ), this, SLOT ( editSettings() ) );

}

void ParPortWidget::readSettings() {

#ifdef _DEBUG
    print ( "ParportWidget: Loading saved settings" );
#endif

    QSettings settings;
    settings.beginGroup ( "ParPort" );
    baseaddress = settings.value ( "address", QString::number ( this->baseaddress, 16 ) ).toString().toULong ( 0, 16 );
    showVoltages = settings.value ( "showVoltages", true ).toBool();
    pollTimer->setInterval ( settings.value ( "pollTime", 1000 ).toInt() );
    receiveSleep = chrono::milliseconds(settings.value ( "receiveSleepTime", receiveSleep.count() ).toInt());
    settings.endGroup();

}

void ParPortWidget::writeSettings() {

#ifdef _DEBUG
    print ( "ParportWidget: Saving settings" );
#endif

    QSettings settings;
    settings.beginGroup ( "ParPort" );
    settings.setValue ( "address", QString::number ( this->baseaddress, 16 ) );
    settings.setValue ( "showVoltages", showVoltages );
    settings.setValue ( "pollTime", pollTimer->interval() );
    settings.setValue ( "receiveSleepTime", receiveSleep.count() );
    settings.endGroup();

}
