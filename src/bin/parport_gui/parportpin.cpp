#include "parportpin.hpp"

#include "parportwidget.hpp"

#include <thread>

using namespace gui;


ParPortPin::ParPortPin(ParPortWidget *ppW, int id, const QColor& color)
    : QWidget(ppW), ppWidget(ppW), pinID(id), m_color(color), state(false), running(false) {
    // QSizePolicy sp = sizePolicy();
    //  sp.setWidthForHeight(true);
    //sp.setHeightForWidth(true);
    //  setSizePolicy(sp);
    settingsAct = new QAction(tr("Settings"), this);
    startAct = new QAction(tr("Start"), this);
    stopAct = new QAction(tr("Stop"), this);
    connect(settingsAct, SIGNAL(triggered()), this, SLOT(editSettings()));
    connect(startAct, SIGNAL(triggered()), this, SLOT(run()));
    connect(stopAct, SIGNAL(triggered()), this, SLOT(stop()));
    connect(this, SIGNAL(toggled(int)), ppWidget, SLOT(pinChanged(int)),Qt::DirectConnection);
    connect(this, SIGNAL(print(QString, QColor)), ppWidget, SLOT(printMsg(QString, QColor)));
}


void ParPortPin::setHigh(bool st) {
    if(st == state) return;
    state = st;
    update();
}


void ParPortPin::editSettings(void) {

    QDialog dialog(this);
    dialog.setWindowTitle("PinSettings");

    QFormLayout form(&dialog);

    //form.addRow(new QLabel("Interval", &dialog));

    QLabel* countLabel = new QLabel("Count", &dialog);
    QLineEdit *countLE = new QLineEdit(&dialog);
    QValidator *val = new QIntValidator(0, 1000000, &dialog);
    countLE->setValidator(val);
    countLE->setText(QString::number(settings.count));
    form.addRow(countLabel, countLE);

    int unit(0);
    size_t count = settings.interval.count();
    while(count && count % 1000 == 0) {
        count /= 1000;
        unit++;
    }
    
    QLabel* intervalLabel = new QLabel("Interval", &dialog);
    form.addRow(intervalLabel);
    QLineEdit *intervalLE = new QLineEdit(&dialog);
    intervalLE->setValidator(val);
    intervalLE->setText(QString::number(count));
    QComboBox *intervalUnit = new QComboBox(&dialog);
    intervalUnit->addItem("microseconds");
    intervalUnit->addItem("milliseconds");
    intervalUnit->addItem("seconds");
    intervalUnit->setCurrentIndex(unit);
    form.addRow(intervalLE, intervalUnit);

    unit = 0;
    count = settings.duration.count();
    while(count && count % 1000 == 0) {
        count /= 1000;
        unit++;
    }
    
    QLabel* durationLabel = new QLabel("Duration", &dialog);
    form.addRow(durationLabel);
    QLineEdit *durationLE = new QLineEdit(&dialog);
    durationLE->setValidator(val);
    durationLE->setText(QString::number(count));
    QComboBox *durationUnit = new QComboBox(&dialog);
    durationUnit->addItem("microseconds");
    durationUnit->addItem("milliseconds");
    durationUnit->addItem("seconds");
    durationUnit->setCurrentIndex(unit);
    form.addRow(durationLE, durationUnit);

    QDialogButtonBox buttonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Horizontal, &dialog);
    form.addRow(&buttonBox);
    QObject::connect(&buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
    QObject::connect(&buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));

    if(dialog.exec() == QDialog::Accepted) {
        settings.count = countLE->text().toInt();
        switch(intervalUnit->currentIndex()) {
            case(0): settings.interval = std::chrono::microseconds(intervalLE->text().toInt()); break;
            case(1): settings.interval = std::chrono::milliseconds(intervalLE->text().toInt()); break;
            default: settings.interval = std::chrono::seconds(intervalLE->text().toInt()); break;
        }
        switch(durationUnit->currentIndex()) {
            case(0): settings.duration = std::chrono::microseconds(durationLE->text().toInt()); break;
            case(1): settings.duration = std::chrono::milliseconds(durationLE->text().toInt()); break;
            default: settings.duration = std::chrono::seconds(durationLE->text().toInt()); break;
        }
    }

}


void ParPortPin::run(void) {

    if(settings.count) {
        running = true;
        update();
        std::thread([this]() {
            size_t cnt = settings.count - 1;
            //printMsg("Running pin " + QString::number(pinID) + "  count="+ QString::number(settings.count));
            do {
                auto next = std::chrono::steady_clock::now() + settings.interval;
                //emit toggled(pinID);
                ppWidget->pinChanged(pinID);
                if(settings.duration.count()) {
                    std::this_thread::sleep_for(settings.duration);
                    ppWidget->pinChanged(pinID);
                }
                //toggle();
                // std::cout << "." << std::flush;
                //std::this_thread::sleep_for(settings.interval);
                std::this_thread::sleep_until(next);
            }
            while(running && cnt--);
            running = false;
            update();
        }).detach();
    }
}


void ParPortPin::mousePressEvent(QMouseEvent *event) {
    if(isEnabled()) {
        Qt::MouseButtons mouseButtons = event->buttons();
        if(mouseButtons == Qt::LeftButton) {
            ppWidget->pinChanged(pinID);    // so clicking always toggles the state.
            //emit toggled(pinID);    // so clicking always toggles the state.
        }
    }
}


void ParPortPin::contextMenuEvent(QContextMenuEvent* event) {

    QMenu menu(parentWidget());
    menu.addAction(settingsAct);
    menu.addAction(startAct);
    menu.addAction(stopAct);
    menu.exec(event->globalPos());

}


void ParPortPin::paintEvent(QPaintEvent *) {
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    int radius = std::min(width(), height()) / 2;
    QPoint center(width() / 2, height() - radius);
    if(running) {
        painter.setPen(Qt::yellow);
        painter.drawEllipse(center, radius , radius);
    }

    radius /= 1.05;
    if(state) painter.setBrush(Qt::green);
    else painter.setBrush(Qt::red);
    painter.drawEllipse(center, radius , radius);


    int fontSize = 20;
    QFont font("Times", fontSize);
    painter.setFont(font);
    QFontMetrics fm = painter.fontMetrics();
    while(fm.width(label) > width()) {
        font.setPointSize(--fontSize);
        painter.setFont(font);
        fm = painter.fontMetrics();
    }


    center = QPoint((width() - fm.width(label)) / 2, 0.65 * fm.height());

    painter.setPen(Qt::black);
    painter.drawText(center, label);

}
