#ifndef REDUX_GUI_PARPORTPIN_HPP
#define REDUX_GUI_PARPORTPIN_HPP

#include <redux/util/parport.hpp>

#include <chrono>
#include <map>

#include <QtGui>


namespace gui {

    class ParPortWidget;

    class ParPortPin : public QWidget {
        
        struct ToggleSettings {
            ToggleSettings() : count(0), interval(0), duration(0) {}
            size_t count;
            std::chrono::duration<size_t,std::micro> interval, duration;
        };
        
        Q_OBJECT

    public:
        ParPortPin(ParPortWidget *ppW, int id = 0, const QColor& color = Qt::cyan);

        bool isHigh() const { return state; }
        void setHigh(bool st);
        //void setMask(uint32_t m) { pinID = m; }
        void setColor(QColor c) { m_color = c; }
        void setLabel(QString l) { label = l; };

    public slots:
        //void on(void) { setHigh(true); }
        //void off(void) { setHigh(false); }
        void toggle(void) { emit toggled(pinID); } //setHigh(!m_state); }
        void editSettings(void);
        void run(void);
        void stop(void) { running = false; }

    signals:
        void toggled(int);
        void print( const QString msg, QColor c );

    protected:
        //int widthForHeight(int w) { return 20 * w; };
        void mousePressEvent(QMouseEvent *event);
        void contextMenuEvent(QContextMenuEvent* event);
        void paintEvent(QPaintEvent *);
        
    private:
        
        void printMsg( const QString msg, QColor c = Qt::black ) { emit print(msg,c); };
        
        ParPortWidget* ppWidget;
        QAction *settingsAct,*startAct,*stopAct;
        QString label;
        int pinID;
        QColor m_color;
        bool state;
        bool running;
        ToggleSettings settings;
        
    };

}   // gui

#endif // REDUX_GUI_PARPORTPIN_HPP
