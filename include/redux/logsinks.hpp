#ifndef REDUX_LOGSINKS_HPP
#define REDUX_LOGSINKS_HPP

#include <string>
#include <vector>

#include <boost/log/sinks/async_frontend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sinks/bounded_fifo_queue.hpp>
#include <boost/log/sinks/drop_on_overflow.hpp>

namespace bsinks = boost::log::sinks;

#define SEV_MASK(sev) (1<<(sev))
#define SEV_OR_HIGHER(sev) ~((1<<(sev))-1)

namespace redux {

    enum severity_level {
        sev_none,
        sev_trace,
        sev_debug,
        sev_detail,
        sev_normal,
        sev_warning,
        sev_error,
        sev_critical
    };

    class LogSink {
    public:
        std::vector<std::string> channels;
        severity_level minSeverity;
        bool truncate;
        LogSink( const severity_level& sev = sev_normal );
        virtual ~LogSink( void ) {};
        void parseFilter( std::string );
    };

    class FileSink : public LogSink {
        typedef bsinks::asynchronous_sink < bsinks::text_ostream_backend,
                bsinks::bounded_fifo_queue< 100, bsinks::drop_on_overflow >
                > sink_type;

        void init( void );
        boost::shared_ptr< sink_type > sink;
    public:
        std::string fileName;
        FileSink( const std::string&, const severity_level& sev = sev_normal );
        ~FileSink( void );
        void parse( const std::string& s );
    };

    class StreamSink : public LogSink {
        typedef bsinks::asynchronous_sink < bsinks::text_ostream_backend,
                bsinks::bounded_fifo_queue< 100, bsinks::drop_on_overflow >
                > sink_type;

        void init( void );
        boost::shared_ptr< sink_type > sink;
        std::ostream& strm;
    public:
        StreamSink( std::ostream&, int = sev_normal );
        ~StreamSink( void );
    };


}

#endif  // REDUX_LOGSINKS_HPP
