#include "redux/util/stringutil.hpp"

#include <cstdlib>
#include <pwd.h>
#include <unistd.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/regex.hpp>

namespace bf = boost::filesystem;

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
    if(!bf::is_directory(ain)) {
        fn = ain.filename();
        ain = ain.parent_path();
    }
    bf::path::iterator it = ain.begin();
    if(!base.empty()) {
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

