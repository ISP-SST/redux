#ifndef REDUX_IMAGE_CACHEDFILE_HPP
#define REDUX_IMAGE_CACHEDFILE_HPP

#include "redux/file/fileio.hpp"
#include "redux/image/image.hpp"
#include "redux/util/cache.hpp"
#include "redux/util/stringutil.hpp"

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
                filename = redux::util::cleanPath( filename.string() ) ;
            }
            
            bfs::path filename;

            bool operator<(const CachedFile& rhs) const { return (filename < rhs.filename); }

            template <typename T>
            static void load( redux::image::Image<T>& img, const std::string& fn, bool metaOnly=false ) {
                CachedFile cf( fn );
                if( !bfs::exists( cf.filename ) || bfs::is_directory(cf.filename) ) throw std::runtime_error("Not a readable file: " + cf.filename.string() );
                auto m = redux::util::Cache::get().getMap< CachedFile, redux::image::Image<T> >();
                auto it = m.second.find(cf);
                if( it != m.second.end() && ((it->second.nElements() > 0) || (it->second.meta && metaOnly)) ) {  // exists and loaded
                    std::unique_lock<std::mutex> lock( it->second.imgMutex ); // sync-point: image might not be done loading.
                    img = it->second;
                } else {    // not already loaded, load the file.
                    //std::cout << ("Loading file: " + cf.filename.string()) << std::endl;
                    auto it2 = m.second.emplace(cf,img);
                    std::unique_lock<std::mutex> lock( it2.first->second.imgMutex ); // sync-point to prevent other threads from getting this item before it is finished loading
                    m.first.unlock();
                    redux::file::readFile( fn, it2.first->second, metaOnly );
                    img = it2.first->second;
                }
//                 return ret.first->second;
//                 redux::image::Image<T>& cimg = redux::util::Cache::get<CachedFile, redux::image::Image<T> >(cf);
//                 std::unique_lock<std::mutex> lock( cimg.imgMutex );
//                 if( !cimg.meta || (!metaOnly && cimg.nElements() == 0) ) {
//                     std::cout << "Loading file: " << cf.filename << std::endl;
//                     redux::file::readFile( fn, cimg, metaOnly );
//                 }
//                 img = cimg;
            }

            template <typename T>
            static void unload( const std::string& fn ) {
                CachedFile cf( fn );
                auto m = redux::util::Cache::get().getMap< CachedFile, redux::image::Image<T> >();
                auto it = m.second.find(cf);
                if( it != m.second.end() ) {
                    {
                        std::unique_lock<std::mutex> lock( it->second.imgMutex ); // sync-point: prevent deletion while loading
                    }
//                     long int uc = it->second.getData().use_count();
//                     //if( it->second.getData().unique() ) {
//                     //std::cout << ("UnLoading file: " + cf.filename.string()) << std::endl;
                    m.second.erase(it);
//                     //} else std::cout << ("UseCount: " + to_string(uc) + "  " + cf.filename.string()) << std::endl;
                }
//                 {
//                     redux::image::Image<T>& cimg = redux::util::Cache::get<CachedFile, redux::image::Image<T> >(cf);
//                     std::unique_lock<std::mutex> lock( cimg.imgMutex );
//                 }
//                 redux::util::Cache::erase<CachedFile, redux::image::Image<T>>(cf);
            }

            
        };
               
    }   // image

}   // redux

#endif  // REDUX_IMAGE_CACHEDFILE_HPP
