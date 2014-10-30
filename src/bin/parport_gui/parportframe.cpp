#include "parportframe.hpp"

#include "parportwidget.hpp"

#include <redux/util/parport.hpp>

using namespace redux::util;
using namespace gui;


namespace {

    const QString PinNames[] = { "~STR", "D0", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "ACK", "T2"/*"~BUS"*/, "CL"/*"PAP"*/, "DA"/*"SEL"*/,
                                "~LF", "ERR", "RES", "~SP", "GND", "GND", "GND", "GND", "GND", "GND", "GND", "GND", "GND" };
}



ParPortFrame::ParPortFrame(ParPortWidget* ppW) : QFrame(ppW), ppWidget(ppW) {

   // connect(this, SIGNAL(print(QString, QColor)), parent, SLOT(printMsg(QString, QColor)));
    
    layout = new QGridLayout(this);

    for(int i = 0; i < 25; ++i) {
        uint32_t pinmask = ParPort::PIN[i + 1];
        ParPortPin* tmp = new ParPortPin(ppWidget, i + 1);
        tmp->setLabel(PinNames[i]);
        layout->addWidget(tmp, (i / 13), (i % 13) * 2 + (i / 13), 1, 1);
        if(pinmask) {
            pins[ParPort::PIN[i + 1]] = tmp;
        }
        else {
            tmp->setColor(Qt::gray);
            tmp->setEnabled(false);
        }
    }


    setFrameShape(QFrame::StyledPanel);
    //setFrameShape( QFrame::Box|QFrame::Sunken );
    QSizePolicy sp = sizePolicy();
    sp.setWidthForHeight(true);
    //sp.setHeightForWidth(true);
    setSizePolicy(sp);
    //setSizePolicy(QSizePolicy::Fixed);
    setMaximumHeight(1200);
    setMinimumHeight(50);
    resize(130, 130);
    //setLineWidth( 10 );
    setVisible(true);

}

void ParPortFrame::setPin(int i, bool v) {
    auto it = pins.find(i);
    if(it != pins.end()) {
        it->second->setHigh(v);
    }
}

void ParPortFrame::setPinEnabled(int i, bool v) {
    auto it = pins.find(i);
    if(it != pins.end()) {
        it->second->setEnabled(v);
    }
}

void ParPortFrame::runAll(void) {
    for(auto &it : pins) {
        it.second->run();
    }
}

