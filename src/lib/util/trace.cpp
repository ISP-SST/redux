#include "redux/util/trace.hpp"

#ifdef __GNUG__
#   include <execinfo.h>
#   include <cxxabi.h>
#endif

#include "redux/util/cache.hpp"

#include <thread>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

//#define RDX_TRACE

#define BT_BUF_SIZE 100

using namespace redux;
using namespace redux::util;
using namespace std;


int redux::util::trace::BT::max_depth = 5;
std::map<size_t, Trace::trace_t> Trace::traces;
std::mutex Trace::mtx;

void redux::util::thread_trace( const char* file, int line ) {
    string mark = string(file) + ":" + to_string(line);
    std::thread::id thread_id = std::this_thread::get_id();
    auto m = Cache::get().getMap<std::thread::id,string>();        // m.first is a unique_lock for the map in m.second
    m.second[ thread_id ] = mark;

}


void redux::util::thread_untrace( void ) {
    
    std::thread::id thread_id = std::this_thread::get_id();
    auto m = Cache::get().getMap<std::thread::id,string>();        // m.first is a unique_lock for the map in m.second
    m.second.erase( thread_id );

}


string redux::util::thread_traces( bool all ) {

    auto m = Cache::get().getMap<std::thread::id,string>();        // m.first is a unique_lock for the map in m.second
    ostringstream ss;
    map<string,size_t> counts;
    for( auto& t: m.second ) {
        counts[t.second]++;
        ss <<  t.first << " - " << t.second << endl;
    }
    string ret;
    if( counts.size() ) {
        multimap<size_t, string, greater<size_t>> sorted_counts;
        std::transform( counts.begin(), counts.end(), inserter(sorted_counts, sorted_counts.begin()),
                        [](const pair<string,size_t>& in) { return pair<size_t,string>(in.second,in.first); });
        for( auto& t: sorted_counts ) {
            ret += to_string(t.first) + " - " + t.second + "\n";
        }
    }
    if( all ) ret += ss.str();
    return ret;
}


void getBT( redux::util::trace::BT& bt ) {

#ifdef __GNUG__
    bt.nSyms = backtrace( bt.mangled_syms, TRACE_BT_BUF_SIZE );
    if( bt.nSyms ) {
        char** syms = backtrace_symbols( bt.mangled_syms, bt.nSyms );
        for( int i(2); i<bt.nSyms; ++i ) bt.syms.push_back( syms[i] ); // Start at 2 in order to skip getBT() & Trace()
        bt.nSyms = bt.syms.size();
        free( syms );
    }
#else
// TODO: implement
#endif

}


trace::BT::BT(void) : nSyms(0) {

    getBT(*this);

}


string trace::BT::printBackTrace( size_t indent ) const {

    
    string ret; 
#ifdef RDX_TRACE

    if( !syms.empty() ) {
        int n = std::min( max_depth, nSyms );
        for( int i(0); i<n; ++i ) {
            ret += string(indent,' ') + to_string(i) + ": ";
            ret += demangle_symbol(syms[i]) + string("\n");
        }
    } else {    // no names, just print the pointers.
        for( int i=2; i<nSyms+2; ++i ) {
            ret += string( indent, ' ' ) + to_string(i) + ": ";
            ret += hexString( mangled_syms[i] ) + string("\n");
        }
    }

#endif
    return ret;
    
}


string Trace::getBackTraces(void) {
    
    string ret;
    lock_guard<mutex> lock(mtx);
    for( auto t: traces ) {
        ret += t.second.backtraces() + "\n";
    }
    return ret;
    
}


string Trace::getStats(void) {
    
    string ret = "Items       Size          Type\n";
    lock_guard<mutex> lock(mtx);
    size_t totalCount(0);
    size_t totalSize(0);
    for( auto t: traces ) {
        ret += t.second.stats() + "\n";
        totalSize += t.second.size();
        totalCount += t.second.count();
    }
    if( totalCount ) ret += "Total traced items: " + to_string( totalCount ) + "  size: " + to_string( totalSize );
    return ret;
    
}


Trace::trace_t& Trace::addTraceObject( size_t id, string_cb_t stats, string_cb_t bt, size_cb_t cnt, size_cb_t sz ) {
    trace_t tmp;
    tmp.stats = stats;
    tmp.backtraces = bt;
    tmp.count = cnt;
    tmp.size = sz;
    lock_guard<mutex> lock(mtx);
    auto it = traces.emplace(id,tmp);
    return it.first->second;
}


void Trace::removeTraceObject( size_t id ) {
    
    lock_guard<mutex> lock(mtx);
    auto it = traces.find(id);
    if( it != traces.end() ) {
        if( it->second.count() == 0 ) {
            traces.erase(id);
        }
    }
    
}


