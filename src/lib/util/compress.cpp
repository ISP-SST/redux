#include "redux/util/compress.hpp"

#include "redux/util/datautil.hpp"
#include "redux/util/stringutil.hpp"

using namespace redux::util;
using namespace std;

unique_ptr<Bytef[]> redux::util::compress( const Bytef* inData, uint64_t uncompressedSz, uint64_t& compressedSz, int compressionLevel ) {
    
    uLongf tmpSz = compressBound(uncompressedSz);
    compressedSz = uncompressedSz;

    unique_ptr<Bytef[]> buf( new Bytef[ tmpSz + sizeof(uncompressedSz)] );
    int ret = compress2( buf.get()+sizeof(uncompressedSz), &tmpSz, inData, uncompressedSz, compressionLevel );
    switch(ret){
        case(Z_OK): compressedSz = tmpSz; break;
        case(Z_MEM_ERROR): cout << "compressing data: out of memory." << endl; break;
        case(Z_BUF_ERROR): cout << "compressing data: buffer not large enough" << endl; break;
        default: cout << "compressing data: unknown reason. ret=" << ret << endl; break;
    }
    pack( reinterpret_cast<char*>(buf.get()), uncompressedSz );

    return std::move(buf);
    
}


unique_ptr<Bytef[]> redux::util::decompress( const Bytef* inData, uint64_t compressedSz, uint64_t& uncompressedSz, bool swap_endian ) {
    
    unpack( reinterpret_cast<const char*>(inData), uncompressedSz, swap_endian );
    uLongf tmpSz = uncompressedSz+sizeof(uncompressedSz);
    unique_ptr<Bytef[]> buf( new Bytef[ tmpSz ] );
    int ret = uncompress( buf.get(), &tmpSz, inData+sizeof(uncompressedSz), compressedSz );
    
    switch(ret){
        case(Z_OK): break;
        case(Z_DATA_ERROR): cout << "decompressing data: corrupt buffer." << endl; break;
        case(Z_MEM_ERROR): cout << "decompressing data: out of memory." << endl; break;
        case(Z_BUF_ERROR): cout << "decompressing data: buffer not large enough." << endl; break;
        default: cout << "decompressing data: unknown reason. ret=" << ret << endl; break;
    }
    
    return std::move(buf);

}

