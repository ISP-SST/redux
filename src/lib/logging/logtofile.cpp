#include "redux/logging/logtofile.hpp"

using namespace redux::logging;


LogToFile::LogToFile( const std::string &filename, uint8_t m, bool replace, unsigned int flushPeriod )
  : LogToStream( fout, m, flushPeriod ),
    fout( filename.c_str(), (replace ? std::ios::out : (std::ios::out | std::ios::app)) )  {

   if( !fout.good() ) throw std::ios_base::failure("Failed to open file: " + filename );
       
}
