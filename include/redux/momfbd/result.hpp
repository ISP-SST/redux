#ifndef REDUX_MOMFBD_RESULT_HPP
#define REDUX_MOMFBD_RESULT_HPP

#include "redux/momfbd/modes.hpp"

#include "redux/image/fouriertransform.hpp"
#include "redux/types.hpp"
#include "redux/work.hpp"
#include "redux/util/array.hpp"

#include <memory>

namespace redux {

    namespace momfbd {

        /*! @ingroup momfbd
         *  @{
         */
        
        class Channel;
        struct ChannelResult : public std::enable_shared_from_this<ChannelResult> {
            
            typedef std::shared_ptr<ChannelResult> Ptr;
            
            ChannelResult(const std::shared_ptr<Channel>& c) : cfg(c) {};
            //~ChannelResult(void);
            
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            
            std::shared_ptr<Channel> cfg;           //!< Channel configuration

        };

        
        class Object;
        struct ObjectResult : public std::enable_shared_from_this<ObjectResult> {
            
            typedef std::shared_ptr<ObjectResult> Ptr;
            
            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);
            
            void addToFT(const redux::image::FourierTransform&);
            void addToPQ(const redux::image::FourierTransform&, const redux::util::Array<complex_t>);
            
            std::mutex mtx;
            redux::image::FourierTransform ftSum;
            redux::util::Array<double> Q;
            redux::util::Array<complex_t> P;
            
            std::shared_ptr<Object> object;

            
        };

        
        class MomfbdJob;
        struct PatchResult : public Part {
            typedef std::shared_ptr<PatchResult> Ptr;
            const MomfbdJob& myJob;
            Point16 index;                      //! Patch-index in mozaic
            Point16 pos;                        //! Position of patch, coordinates in the "anchor channel"
            std::vector<ObjectResult::Ptr> objects;
            PatchResult(const MomfbdJob&);

            uint64_t size(void) const;
            uint64_t pack(char*) const;
            uint64_t unpack(const char*, bool);

        };

 
        /*! @} */


    }

}

#endif  // REDUX_MOMFBD_RESULT_HPP
