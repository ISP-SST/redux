#include "redux/util/stringutil.hpp"

#include "redux/translators.hpp"
#include "redux/util/datautil.hpp"

#include <cstdlib>
#include <mutex>
#include <pwd.h>
#include <unistd.h>

#ifdef __GNUG__
#   include <cxxabi.h>
#endif

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/version.hpp>

namespace bf = boost::filesystem;
namespace bpt = boost::property_tree;


using namespace redux::util;
using namespace std;


bool redux::util::onlyDigits(const string &s) {
    static const boost::regex e("[0-9]+");
    return boost::regex_match(s, e);
}

bool redux::util::onlyAlpha(const string &s) {
    static const boost::regex e("^[[:alpha:]]*$");
    return boost::regex_match(s, e);
}

bool redux::util::onlyAlnum(const string &s) {
    static const boost::regex e("^[[:alnum:]]*$");
    return boost::regex_match(s, e);
}

bool redux::util::onlyHex(const string &s) {
    static const boost::regex e("^[[:xdigit:]]*$");
    return boost::regex_match(s, e);
}

bool redux::util::isInteger(const string &s) {
    static const boost::regex e("^\\s*(-|0x|0X)?[[:xdigit:]]+\\s*$");
    return boost::regex_match(s, e);
}


bool redux::util::isHex(const string &s) {
    static const boost::regex e("^\\s*[[:xdigit:]]+\\s*$");
    return boost::regex_match(s, e);
}


bool redux::util::contains ( const string & haystack, const string & needle, bool ignoreCase ) {

    auto it = std::search (
                  haystack.begin(), haystack.end(),
                  needle.begin(),   needle.end(),
    [ignoreCase] ( char ch1, char ch2 ) {
        if ( ignoreCase ) return std::toupper ( ch1 ) == std::toupper ( ch2 );
        return ch1 == ch2;
    }
              );
    if ( it != haystack.end() ) return true;
    return false;

}


string redux::util::replace_n( std::string input, const std::string& loc, const std::string& replace, size_t n) {
    
    if(n == 0) return input;
    
    size_t count(0);
    for( size_t offset = input.find(loc); offset != std::string::npos; ) {
        count++;
        offset = input.find(loc, offset+loc.length());
    }

    count = std::min(n, count);
    
    while( count-- ) {
        boost::replace_first( input, loc, replace );
    }
    
    return input;
    
}


string redux::util::popword( std::string &line, const char *separator ) {
    
    size_t pos = line.find_first_not_of( separator );
    if( pos == string::npos ) line.clear();
    else if( pos ) line.erase(0, pos);
    
    pos = line.find_last_not_of( separator );
    if( pos != string::npos ) line.erase(pos+1);
    else line.clear();
    
    pos = line.find_first_of( separator );
    string result( line, 0, pos );
    line.erase( 0, pos );

    return result;
}


bool redux::util::nocaseLess(const string& lhs, const string& rhs) {

          return std::lexicographical_compare( lhs.begin (), lhs.end (), rhs.begin (), rhs.end (),
                                               [] ( char ch1, char ch2 ) { return std::toupper( ch1 ) < std::toupper( ch2 ); }
                                             );

}


bool redux::util::isRelative( const std::string &s ) {
    return (!s.empty() && s[0] != '/');
}


void redux::util::maybeCreateDir( const bfs::path& p ) {
    
    if( !p.empty() && !bfs::exists(p) ) {
        if( !bfs::create_directories(p) ) {
            cout << boost::format( "failed to create directory: %s" ) % p << endl;
        }
    }

}


vector<set<string>> redux::util::make_template( const vector<string>& list, string& out, string split_chars ) {
    
    vector<string> segments;
    vector< set<string> > seg_list;

    if( split_chars.empty() ) return seg_list;
    
    size_t n_segments(0);
    for( auto& str: list ) {
        boost::split(segments, str, boost::is_any_of(split_chars));
        if( !n_segments && (segments.size() >= 1) ) {    // TBD: how to treat size=1 (no split), return tpl="" or tpl="%1"
            n_segments = segments.size();
            seg_list.resize( n_segments );
        }
        if( segments.size() == n_segments ) {
            for( size_t j=0; j<n_segments; ++j ) {
                seg_list[j].insert(segments[j]);
            }
        } else if(n_segments) {
            //cout << "Multiple segment sizes!!  nS=" << n_segments << "  ss=" << segments.size() << endl;
        }
    }
    
    int arg_cnt(0);
    out = "";
    for( size_t j=0; j<n_segments; ++j ) {
        size_t nArgs = seg_list[j].size();
        if( nArgs > 1 ) {
            out += "%"+to_string(++arg_cnt);     // add a placeholder for this segment
        } else if( nArgs == 1 ) {
            out += *(seg_list[j].begin());       // a unique item, add it to template
        }
        if( j < n_segments-1 ) out += split_chars[0];       // separator    TBD: should this be a parameter
    }

    return std::move(seg_list);
    
}


string redux::util::alignCenter(const string& s, size_t n, unsigned char c) {

    if(s.length() > n) {
        return s.substr(0, n);
    }

    string tmp  = s + string((n - s.length()) >> 1, c);
    // if n is odd, the extra character goes on the left side
    return string(n - tmp.length(), c) + tmp;


}


string redux::util::alignLeft(const string& s, size_t n, unsigned char c) {

    if(s.length() > n) {
        return s.substr(0, n);
    }

    return s + string(n - s.length(), c);

}


string redux::util::alignRight(const string& s, size_t n, unsigned char c) {

    if(s.length() > n) {
        return s.substr(0, n);
    }

    return string(n - s.length(), c) + s;

}


string redux::util::getUname(__uid_t id) {
    
    static map<__uid_t,string> users;
    static mutex mtx;
    lock_guard<mutex> lock(mtx);
    
    if(!id) id = geteuid();
    
    auto it = users.find(id);
    if( it != users.end() ) return it->second;

    string tmp;
    struct passwd pwent;
    struct passwd *pwentp;
    char buf[1024];
    if( !getpwuid_r( id, &pwent, buf, 1024, &pwentp ) ) {
        tmp = pwent.pw_name;
    }
    else {
        tmp = std::to_string((int)id);
    }
    
    users.emplace( id, tmp );
    
    return tmp;
}


string expandTilde(string in) {
    if(in.empty() || in[0] != '~') return in;
    string tmp;
    size_t cut;
    struct passwd pwent;
    struct passwd *pwentp;
    char buf[1024];
    if(in.length() == 1 || in[1] == '/') {
        cut = 1;
        tmp = getenv("HOME");
        if(tmp.empty()) {
            if( !getpwuid_r( geteuid(), &pwent, buf, 1024, &pwentp ) ) {
                tmp = pwent.pw_dir;
            }
        }
    }
    else {
        cut = in.find_first_of('/');
        string user = in.substr(1, cut - 1);
        if( !getpwnam_r(user.c_str(), &pwent, buf, 1024, &pwentp) ) {
            tmp = pwent.pw_dir;
        }
    }
    if(tmp.empty()) return in;
    else return tmp + in.substr(cut);
}


string redux::util::cleanPath(string in, string base) {

    if(in.empty()) return in;
    bf::path fn, result, ain(in);
    if(!base.empty() && base[0] == '~') base = expandTilde(base);
    if(!base.empty() && base[0] != '/') result = bf::current_path() / bf::path(base);
    if(in[0] == '~') in = expandTilde(in);

    if(bf::is_regular_file(ain)) {
        fn = ain.filename();
        ain = ain.parent_path();
    }
    auto it = ain.begin();
    if(in[0] != '/' && !base.empty()) {
        if(!bf::is_directory(result)) return in;
    }
    else result = *it++;

    bool docanonical RDX_UNUSED = (result.string()[0] == '/');         // don't canonicalize relative paths
    for(; it != ain.end(); ++it) {
        if(*it == "..") result = result.parent_path();
        else if(*it != ".") {
#if BOOST_VERSION > 104800
            if(!exists(result / *it) && docanonical) {      // canonicalize the existing part (boost >= 1.48)
                result = canonical(result);
                docanonical = false;
            }
#else
            docanonical = false;
#endif
            result /= *it;
        }
    }

    return (result / fn).string();

}

void redux::util::printProgress( const string& text, float progress ) {
    static mutex mtx;
    unique_lock<mutex> lock(mtx);
    if( progress >= 0 ) {
        printf( "\r%s (%.1f%%)", text.c_str(), progress );
    } else printf( "\r%s", text.c_str() );
    fflush(stdout);
}


template <typename T>
vector<T> redux::util::stringToUInts(const string& str) {
    
    bpt::ptree tmpTree;                         // just to be able to use the VectorTranslator
    tmpTree.put( "tmp", str );
    return tmpTree.get<vector<T>>( "tmp", vector<T>() );

}
template vector<uint8_t> redux::util::stringToUInts(const string&);
template vector<uint16_t> redux::util::stringToUInts(const string&);
template vector<uint32_t> redux::util::stringToUInts(const string&);
template vector<uint64_t> redux::util::stringToUInts(const string&);


template <typename T>
std::string redux::util::uIntsToString(const std::vector<T>& ints) {
    
    bpt::ptree tmpTree;
    tmpTree.put("tmp", ints);
    return tmpTree.get<string>( "tmp", "" );
    
}
template string redux::util::uIntsToString(const vector<uint8_t>& );
template string redux::util::uIntsToString(const vector<uint16_t>& );
template string redux::util::uIntsToString(const vector<uint32_t>& );
template string redux::util::uIntsToString(const vector<uint64_t>& );


std::string redux::util::colorString( const std::string& in, StringColor col ) {

    static std::string a("\033[");
    static std::string b("m");
    static std::string c("\033[0m");
    return a + to_string(col) + b + in + c;

}


string redux::util::tvToString( const timeval& a, bool millis ) {


    char tmp[15];
    strftime( tmp, 14, "%H:%M:%S", gmtime( &a.tv_sec ) );
    string ret( tmp );

    if ( millis ) {
        sprintf( tmp, ".%.3u", ( uint )( a.tv_usec / 1000 ) );
        ret += tmp;
    }

    return ret;

}


string redux::util::tsToString( const timespec& a, bool millis ) {


    char tmp[15];
    strftime( tmp, 14, "%H:%M:%S", gmtime( &a.tv_sec ) );
    string ret( tmp );

    if ( millis ) {
        sprintf( tmp, ".%.6u", ( uint )( a.tv_nsec / 1000 ) );
        ret += tmp;
    }

    return ret;

}


#ifdef __GNUG__
    std::string redux::util::demangle_name( const string& name ) {
        int status(0);
        std::unique_ptr<char, void(*)(void*)> res {
            abi::__cxa_demangle( name.c_str(), 0, 0, &status ),
            std::free
        };
        return (status==0) ? res.get() : name ;
    }
#else
    std::string redux::util::demangle_name( const string& name ) { return name; }    // TODO: implement
#endif


string redux::util::demangle_symbol( const string& sym ) {
    // Example symbol: ./module(func_name+0x15c) [0x8048a6d]
    static const boost::regex sym_re("([[:alnum:]_]+)[+ ]+([0x[:xdigit:]]+)");      // match func_name & offset
    boost::regex re( "(\\d+)+" );
    boost::smatch match;
    string ret = sym;
    if( boost::regex_search( sym, match, sym_re ) ) {
        string sym_name = string(match[1]);
        ret = demangle_name( sym_name );
    }
    return ret;

}

