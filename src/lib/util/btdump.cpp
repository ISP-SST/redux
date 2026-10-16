#include "redux/util/btdump.hpp"

#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <dirent.h>
#include <cstring>
#include <cstdlib>
#include <atomic>
#include <chrono>
#include <thread>
#include <sstream>
#include <vector>

using namespace std;

namespace {

    const int MAX_FRAMES = 64;
    const int MAX_SLOTS = 1024;

    struct Slot {
        std::atomic<pid_t> tid{0};
        std::atomic<int> nFrames{0};
        std::atomic<bool> ready{false};
        void* frames[MAX_FRAMES];
    };

    Slot slots[MAX_SLOTS];

    Slot* findOrClaimSlot( pid_t tid ) {
        Slot* free_slot = nullptr;
        for( auto& s: slots ) {
            pid_t cur = s.tid.load( std::memory_order_relaxed );
            if( cur == tid ) return &s;
            if( !free_slot && cur == 0 ) free_slot = &s;
        }
        if( free_slot ) {
            pid_t expected = 0;
            if( free_slot->tid.compare_exchange_strong( expected, tid ) ) return free_slot;
            // lost the race, fall back to a full re-scan for one that now matches.
            for( auto& s: slots ) {
                if( s.tid.load( std::memory_order_relaxed ) == tid ) return &s;
            }
        }
        return nullptr;    // table full
    }

    void sigusr1Handler( int ) {
        pid_t tid = (pid_t)syscall( SYS_gettid );
        Slot* s = findOrClaimSlot( tid );
        if( !s ) return;
        s->ready.store( false, std::memory_order_relaxed );
        int n = backtrace( s->frames, MAX_FRAMES );
        s->nFrames.store( n, std::memory_order_relaxed );
        s->ready.store( true, std::memory_order_release );
    }

    struct Installer {
        Installer( void ) {
            // Force glibc's lazy backtrace()-support init to happen here, outside of
            // signal context, so later calls from the handler don't risk doing it unsafely.
            void* warm[4];
            backtrace( warm, 4 );
            struct sigaction sa;
            memset( &sa, 0, sizeof(sa) );
            sa.sa_handler = sigusr1Handler;
            sa.sa_flags = SA_RESTART;
            sigaction( SIGUSR1, &sa, nullptr );
        }
    } installer;

}


string redux::util::dumpAllThreadBacktraces( void ) {

    vector<pid_t> tids;
    DIR* d = opendir( "/proc/self/task" );
    if( !d ) return "dumpAllThreadBacktraces: failed to open /proc/self/task\n";
    struct dirent* ent;
    while( (ent = readdir(d)) ) {
        if( ent->d_name[0] == '.' ) continue;
        tids.push_back( (pid_t)atoi( ent->d_name ) );
    }
    closedir(d);

    pid_t myTid = (pid_t)syscall( SYS_gettid );

    for( pid_t tid: tids ) {
        if( tid == myTid ) continue;
        syscall( SYS_tgkill, getpid(), tid, SIGUSR1 );
    }

    // Give signalled threads a moment to run the handler and fill their slot.
    std::this_thread::sleep_for( std::chrono::milliseconds(150) );

    ostringstream ss;
    for( pid_t tid: tids ) {
        void* frames[MAX_FRAMES];
        int n(0);
        if( tid == myTid ) {
            n = backtrace( frames, MAX_FRAMES );
        } else {
            Slot* s = nullptr;
            for( auto& sl: slots ) {
                if( sl.tid.load( std::memory_order_relaxed ) == tid ) { s = &sl; break; }
            }
            if( !s || !s->ready.load( std::memory_order_acquire ) ) {
                ss << "Thread " << tid << ": (no response)\n";
                continue;
            }
            n = s->nFrames.load( std::memory_order_relaxed );
            memcpy( frames, s->frames, n*sizeof(void*) );
        }
        ss << "Thread " << tid << ":\n";
        char** syms = backtrace_symbols( frames, n );
        if( syms ) {
            for( int i=0; i<n; ++i ) ss << "  " << syms[i] << "\n";
            free( syms );
        }
    }
    return ss.str();

}
