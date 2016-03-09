#ifndef REDUX_IMAGE_CACHEDFILE_HPP
#define REDUX_IMAGE_CACHEDFILE_HPP

#include "redux/file/fileio.hpp"
#include "redux/image/image.hpp"
#include "redux/util/cache.hpp"

#include <iostream>
#include <mutex>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
namespace bfs = boost::filesystem;

namespace redux {

    namespace image {
                
        struct CachedFile {
            
            CachedFile( const std::string& fn ) : filename(fn) {
                filename = bfs::canonical( filename ) ;
            }
            
            bfs::path filename;

            bool operator<(const CachedFile& rhs) const { return (filename < rhs.filename); }

            template <typename T>
            void loadImage( redux::image::Image<T>& img ) const { redux::file::readFile( filename.string(), img ); }

            template <typename T>
            static void load( redux::image::Image<T>& img, const std::string& fn ) {
                CachedFile cf( fn );
                redux::image::Image<T>& cimg = redux::util::Cache::get<CachedFile, redux::image::Image<T> >(cf);
                {
                    std::unique_lock<std::mutex> lock( cimg.imgMutex );
                    if( cimg.nElements() == 0 ) cf.loadImage( cimg );
                }
                img = cimg;
            }

            template <typename T>
            static void unload( const std::string& fn ) {
                CachedFile cf( fn );
                redux::util::Cache::erase<CachedFile, redux::image::Image<T>>(cf);
            }
            
            template <typename T>
            static redux::image::Image<T>& get( const std::string& fn ) {
                CachedFile cf( fn );
                redux::image::Image<T>& cimg = redux::util::Cache::get<CachedFile, redux::image::Image<T> >(cf);
                {
                    std::unique_lock<std::mutex> lock( cimg.imgMutex );
                    if( cimg.nElements() == 0 ) cf.loadImage( cimg );
                    else std::cout << std::endl << "CACHED: " << fn << "  = " << cf.filename << std::endl << std::endl;
                }
                return cimg;
            }

            
        };
               
    }   // image

}   // redux

#endif  // REDUX_IMAGE_CACHEDFILE_HPP
