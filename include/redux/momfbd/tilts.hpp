#ifndef REDUX_MOMFBD_TILTS_HPP
#define REDUX_MOMFBD_TILTS_HPP

#include "redux/momfbd/subimage.hpp"


#include <memory>
#include <mutex>
#include <set>


namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */

        class Channel;
        class WorkSpace;
        
        /*! @brief Constrained tilts.
         *
         */
        struct Tilts : public std::enable_shared_from_this<Tilts> {
            
            Tilts(Channel& c, uint16_t m);
            ~Tilts();
            
            size_t nFreeParameters(void) const;
            size_t setPointers(double* a, double* g );

            void init(void);
            void addRelativeTilt(std::shared_ptr<Tilts>& t);

            void applyTiltToImages(boost::asio::io_service&);
            void addTiltToImages(boost::asio::io_service&, double*);
            void calcGradient(boost::asio::io_service&, uint8_t, const grad_t&);

            double getPartial(size_t i);
            void addPartials(void);
            size_t size(void) { return 1+relativeTilts.size(); }

            Channel& channel;
            std::vector<std::shared_ptr<Tilts>> relativeTilts;
            std::shared_ptr<Tilts> refTilt;
            std::shared_ptr<WorkSpace> ws;
            uint16_t mode;
            size_t nFreeAlpha;
            
            double *x_, *df_;                               // shortcuts to data.
            
            redux::util::Array<double> phi;                 // temporary arrays for gradient calculations
            redux::util::Array<complex_t> otf;
            double* partialGrad;
            bool anchorChannel;
            
        };
        

        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_TILTS_HPP
