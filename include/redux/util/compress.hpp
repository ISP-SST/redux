#ifndef REDUX_UTIL_COMPRESS_HPP
#define REDUX_UTIL_COMPRESS_HPP

#include "redux/util/datautil.hpp"

#include <memory>
#include <zlib.h>


namespace redux {

    namespace util {


        /*!  @ingroup util
         *  @{
         */

        std::unique_ptr<Bytef[]> compress(const Bytef* inData, uint64_t uncompressedSz, uint64_t& compressedSz, int compressionLevel=Z_DEFAULT_COMPRESSION);
        std::unique_ptr<Bytef[]> decompress( const Bytef* inData, uint64_t compressedSz, uint64_t& uncompressedSz, bool swap_endian );
        
        template <class T, int LVL=Z_DEFAULT_COMPRESSION>
        class Compressed : public T {
            
        public:
            template <class... S>
            explicit Compressed(S&& ... s) : T(std::forward<S>(s)...),
                isCompressed(false), compressionLevel(LVL),
                compressedSize(0), uncompressedSize(0) { }
            
            static inline uint64_t metaSize(void) { return 2*sizeof(uLongf)+1; }
            uint64_t size(void) const { return T::size() + metaSize(); };
            uint64_t pack( char* ptr ) const {
                uint64_t count = metaSize();      // skip space for meta
                uncompressedSize = T::pack(ptr+count);
                compressedSize = compressBound(uncompressedSize);
                std::unique_ptr<Bytef[]> buf( new Bytef[compressedSize] );
                int ret = compress2( buf.get(), &compressedSize, reinterpret_cast<Bytef*>(ptr)+count, uncompressedSize, compressionLevel );
                if( ret == Z_OK && compressedSize < uncompressedSize ) {
                    isCompressed = true;
                } else {
                    isCompressed = false;
                    compressedSize = 0;
                }
                count = packMeta( ptr );
                if( isCompressed ) {
                    memcpy( ptr+count, buf.get(), compressedSize );
                    count += compressedSize;
                } else count += uncompressedSize;
                return count;
                
            };
            uint64_t unpack( const char* ptr, bool swap_endian ) {
                uint64_t count = unpackMeta( ptr, swap_endian );
                if ( isCompressed ) {
                    uLongf tmpSize = uncompressedSize;
                    std::unique_ptr<Bytef[]> buf( new Bytef[ uncompressedSize ] );
                    int ret = uncompress( buf.get(), &tmpSize, reinterpret_cast<const Bytef*>(ptr)+count, compressedSize );
                    if( ret != Z_OK ){
                        // TODO error handling
                    }
                    count += T::unpack( reinterpret_cast<const char*>(buf.get()), swap_endian );
                } else {
                    count += T::unpack( ptr+count, swap_endian );
                }
                return count;
                
            };
            
        private:
            
            uint64_t packMeta( char* ptr ) const {
                using redux::util::pack;
                uint64_t count = pack( ptr, isCompressed );
                count += pack( ptr+count, compressedSize );
                count += pack( ptr+count, uncompressedSize );
                return count;
            }
            uint64_t unpackMeta( const char* ptr, bool swap_endian ) {
                using redux::util::unpack;
                uint64_t count = unpack( ptr, isCompressed );
                count += unpack( ptr+count, compressedSize, swap_endian );
                count += unpack( ptr+count, uncompressedSize, swap_endian );
                return count;
            }
            
            mutable bool isCompressed;
            mutable int compressionLevel;
            mutable uLongf compressedSize;
            mutable uLongf uncompressedSize;
            
        };
        
        
        /*! @} */

    }

}

#endif // REDUX_UTIL_COMPRESS_HPP
