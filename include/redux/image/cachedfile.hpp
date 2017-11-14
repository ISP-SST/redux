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
                filename = redux::util::cleanPath( filename ) ;
            }
            
            std::string filename;

            bool operator<(const CachedFile& rhs) const { return (filename < rhs.filename); }

            template <typename T>
            static void load( redux::image::Image<T>& img, const std::string& fn, bool metaOnly=false ) {
                CachedFile cf( fn );
                bfs::path cfp(cf.filename);
                if( !bfs::exists(cfp) || bfs::is_directory(cfp) ) throw std::runtime_error("Not a readable file: " + cf.filename );
                auto m = redux::util::Cache::get().getMap< CachedFile, redux::image::Image<T> >();
                auto it = m.second.find(cf);
                if( it != m.second.end() && ((it->second.nElements() > 0) || (it->second.meta && metaOnly)) ) {  // exists and loaded
                    std::unique_lock<std::mutex> lock( it->second.imgMutex ); // sync-point to prevent other threads from getting this item before it is finished loading
                    img = it->second;
                    m.first.unlock();
                } else {    // not already loaded, load the file.
                    auto it2 = m.second.emplace(cf,img);
                    std::unique_lock<std::mutex> lock( it2.first->second.imgMutex ); // sync-point to prevent other threads from getting this item before it is finished loading
                    m.first.unlock();
                    redux::file::readFile( fn, it2.first->second, metaOnly );
                    img = it2.first->second;
                }
            }

            template <typename T>
            static void unload( const std::string& fn ) {
                CachedFile cf( fn );
                auto m = redux::util::Cache::get().getMap< CachedFile, redux::image::Image<T> >();
                auto it = m.second.find(cf);
                if( it != m.second.end() ) {
                    std::unique_lock<std::mutex> lock( it->second.imgMutex ); // sync-point: prevent deletion while loading
                    if( it->second.getData().use_count() < 2 ) {
                        lock.unlock();
                        m.second.erase(it);
                    }
                }
            }

            
        };
               
    }   // image

}   // redux

#endif  // REDUX_IMAGE_CACHEDFILE_HPP
