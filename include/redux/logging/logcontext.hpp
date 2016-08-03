#ifndef REDUX_LOGGING_LOGCONTEXT_HPP
#define REDUX_LOGGING_LOGCONTEXT_HPP

#include <string>
#include <stack>
#include <ostream>

namespace redux {

    namespace logging {

        class Logger;

        class LogContext {
        public:
            LogContext(void) : name(""), taskDepth(0) {}

            inline const char *getName(void) const { return name.c_str(); }
            inline int getTaskDepth(void) const { return taskDepth; }
            inline const char *getScope(void) const {
                if( scope.empty() ) {
                    return "global_scope";
                } else {
                    return scope.top().c_str();
                }
            }

            inline void setName(const std::string &s) { name = s; }
            inline void increaseTaskDepth(void) { ++taskDepth; }
            inline void decreaseTaskDepth(void) { --taskDepth; }
            inline void pushScope(const std::string &s) { scope.push(s); }
            inline void popScope(void) { scope.pop(); }

        private:

            std::string name;
            int taskDepth;
            std::stack<std::string> scope;
        };

        class LogContextSnapshot {
        public:
            LogContextSnapshot(void) : taskDepth(0) {}
            LogContextSnapshot(const LogContext &lc)
                : name(lc.getName()), taskDepth(lc.getTaskDepth()),
                  scope(lc.getScope()) {}

            uint64_t size(void) const { return 0; };
            uint64_t pack(char*) const { return 0; };
            uint64_t unpack(const char*, bool) { return 0; };
            
            std::string name;
            int taskDepth;
            std::string scope;
        };

    }

}

std::ostream &operator<<( std::ostream &os, const redux::logging::LogContext &lc );
std::ostream &operator<<( std::ostream &os, const redux::logging::LogContextSnapshot &lcs );

#endif // REDUX_LOGGING_LOGCONTEXT_HPP
