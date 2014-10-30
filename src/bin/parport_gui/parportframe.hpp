#ifndef REDUX_GUI_PARPORTFRAME_HPP
#define REDUX_GUI_PARPORTFRAME_HPP

#include "parportpin.hpp"

#include <map>

#include <QtGui>


namespace gui {
    
    class ParPortWidget;

    class ParPortFrame : public QFrame {

        Q_OBJECT

    public:
        ParPortFrame ( ParPortWidget* );
        ~ParPortFrame() {};

        void setPin ( int i, bool v=true );
        void setPinEnabled ( int i, bool v=true );

    public slots:
        void runAll(void);
        
    protected:

        //int heightForWidth( int w ) { return 2*w; }
        int widthForHeight ( int w ) {
            return 2 * w;
        }
    private:

        ParPortWidget* ppWidget;
        QGridLayout *layout;
        std::map<uint32_t,ParPortPin*> pins;

    };

}   // gui

#endif // REDUX_GUI_PARPORTFRAME_HPP
