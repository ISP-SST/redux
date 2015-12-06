#include "redux/util/stringutil.hpp"

#include "redux/translators.hpp"

#include <cstdlib>
#include <pwd.h>
#include <unistd.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/regex.hpp>
#include <boost/property_tree/ptree.hpp>

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


bool redux::util::nocaseLess(const string& lhs, const string& rhs) {

          return std::lexicographical_compare( lhs.begin (), lhs.end (), rhs.begin (), rhs.end (),
                                               [] ( char ch1, char ch2 ) { return std::toupper( ch1 ) < std::toupper( ch2 ); }
                                             );

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
    if(!id) id = geteuid();
    struct passwd* pwd = getpwuid(id);

    string tmp;
    if(pwd) {
        tmp = pwd->pw_name;
    }
    else {
        tmp = std::to_string((int)id);
    }
    return tmp;
}


string expandTilde(string in) {
    if(in.empty() || in[0] != '~') return in;
    string tmp;
    size_t cut;
    if(in.length() == 1 || in[1] == '/') {
        cut = 1;
        tmp = getenv("HOME");
        if(tmp.empty()) {
            struct passwd* pwd = getpwuid(geteuid());
            if(pwd) tmp = pwd->pw_dir;
        }
    }
    else {
        cut = in.find_first_of('/');
        string user = in.substr(1, cut - 1);
        struct passwd* pwd = getpwnam(user.c_str());
        if(pwd) tmp = pwd->pw_dir;
    }
    if(tmp.empty()) return in;
    else return tmp + in.substr(cut);
}


string redux::util::cleanPath(string in, string base) {

    if(in.empty()) return in;
    if(!base.empty() && base[0] == '~') base = expandTilde(base);
    if(in[0] == '~') in = expandTilde(in);

    bf::path fn, result, ain(in);
    if(bf::is_regular_file(ain)) {
        fn = ain.filename();
        ain = ain.parent_path();
    }
    auto it = ain.begin();
    if(in[0] != '/' && !base.empty()) {
        result = bf::absolute(base);
        if(!bf::is_directory(result)) return in;
    }
    else result /= *it++;

    bool docanonical = (result.string()[0] == '/');         // don't canonicalize relative paths
    for(; it != ain.end(); ++it) {
        if(*it == "..") result = result.parent_path();
        else if(*it != ".") {
            if(!exists(result / *it) && docanonical) {      // canonicalize the existing part
                result = canonical(result);
                docanonical = false;
            }
            result /= *it;
        }
    }

    return (result / fn).string();

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
