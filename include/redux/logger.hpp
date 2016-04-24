#ifndef REDUX_LOGGER_HPP
#define REDUX_LOGGER_HPP

#include "redux/util/boundvalue.hpp"
#include "redux/logsinks.hpp"

#include <string>
#include <vector>

#include <boost/log/common.hpp>
#include <boost/log/sources/severity_channel_logger.hpp>
#include <boost/program_options.hpp>

namespace bsrc = boost::log::sources;
namespace bpo = boost::program_options;

namespace redux {


    class Logger {
    public:

        typedef bsrc::severity_channel_logger< severity_level, std::string > logger_t;
        typedef bsrc::severity_channel_logger_mt< severity_level, std::string > logger_mt;

        Logger( bpo::variables_map& );
        ~Logger( void );

        void addFileLog( const std::string& file );
        void addNullLog( void );
        
        static int getDefaultSeverity(void) { return (int)defaultSeverity; }
        
        static std::pair<std::string, std::string> customParser( const std::string& s );
        static bpo::options_description getOptions( const std::string& application_name );

        static logger_t lg;
        static logger_mt mlg;

    private:

        std::vector< boost::shared_ptr<LogSink> > logSinks;
        static redux::util::BoundValue<int> defaultSeverity;

    };

}


#define LOG_TRACE BOOST_LOG_CHANNEL_SEV(lg, logChannel, redux::sev_trace)
#define LOG_DEBUG BOOST_LOG_CHANNEL_SEV(lg, logChannel, redux::sev_debug)
#define LOG_DETAIL BOOST_LOG_CHANNEL_SEV(lg, logChannel, redux::sev_detail)
#define LOG BOOST_LOG_CHANNEL_SEV(lg, logChannel, redux::sev_normal)
#define LOG_WARN BOOST_LOG_CHANNEL_SEV(lg, logChannel, redux::sev_warning)
#define LOG_ERR BOOST_LOG_CHANNEL_SEV(lg, logChannel, redux::sev_error)
#define LOG_CRITICAL BOOST_LOG_CHANNEL_SEV(lg, logChannel, redux::sev_critical)

#define LOGC_TRACE(chan) BOOST_LOG_CHANNEL_SEV(lg, chan, redux::sev_trace)
#define LOGC_DEBUG(chan) BOOST_LOG_CHANNEL_SEV(lg, chan, redux::sev_debug)
#define LOGC_DETAIL(chan) BOOST_LOG_CHANNEL_SEV(lg, chan, redux::sev_detail)
#define LOGC(chan) BOOST_LOG_CHANNEL_SEV(lg, chan, redux::sev_normal)
#define LOGC_WARN(chan) BOOST_LOG_CHANNEL_SEV(lg, chan, redux::sev_warning)
#define LOGC_ERR(chan) BOOST_LOG_CHANNEL_SEV(lg, chan, redux::sev_error)
#define LOGC_CRITICAL(chan) BOOST_LOG_CHANNEL_SEV(lg, chan, redux::sev_critical)

#endif  // REDUX_LOGGER_HPP
