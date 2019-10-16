#ifndef REDUX_UTIL_RECURSEPATH_HPP
#define REDUX_UTIL_RECURSEPATH_HPP

#include <limits>

#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

namespace redux {

    namespace util {
        
       /*!  @ingroup util
        *  @{
        */


       /*! @enum EntryType
        *  @brief Identifiers for filesystem entries.s
        */
        enum EntryType { ET_UNDEF=0, ET_DIRECTORY=1, ET_FILE=2, ET_SYMLINK=4 };

        /*!  @class     RecursePath
         *   @brief     Wrapper class that calls a specified function for all files/folders in a tree.
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2016
         */
        class RecursePath  {
            
            typedef bool (*functype)( const bfs::path& );

        public:

            RecursePath( bfs::path& p, functype f, int nSubs=std::numeric_limits<int>::max() );
            RecursePath( const RecursePath& rhs, int sublevels );
            RecursePath( const RecursePath& rhs );

            void operator()( bfs::directory_entry& p ) const;
            void operator()( const bfs::path& p ) const;

        private:

            int nSubLevels;
            functype callBack;

        }; // end RecursePath


       /*! @} */

    } // end namespace util


} // end namespace redux



#endif // REDUX_UTIL_RECURSEPATH_HPP
