#include "redux/python/common.hpp"

#include <pybind11/stl.h>

#include "redux/job.hpp"
#include "redux/logging/logger.hpp"
#include "redux/momfbd/momfbdjob.hpp"

#include <boost/property_tree/info_parser.hpp>
#include <boost/program_options.hpp>
#include <sstream>



namespace redux {

    namespace python {

        namespace py = pybind11;
        namespace bpo = boost::program_options;
        using namespace redux::momfbd;
        using namespace std;


        static py::dict parse_cfg( const string& path, const string& name = "", const string& logfile = "" ) {

            require_file( path );
            string absPath = boost::filesystem::absolute( path ).string();

            string jobName = name.empty() ? boost::filesystem::path( absPath ).stem().string() : name;
            string logFile = logfile.empty() ? jobName + ".log" : logfile;
            string text = Job::filterOldCfg( absPath, jobName, logFile );

            boost::property_tree::ptree tree;
            std::istringstream iss( text );
            boost::property_tree::read_info( iss, tree );


            boost::filesystem::path cfgDir = boost::filesystem::path( absPath ).parent_path();
            boost::filesystem::path savedCwd = boost::filesystem::current_path();
            struct CwdGuard {
                boost::filesystem::path saved;
                ~CwdGuard() { boost::filesystem::current_path( saved ); }
            } cwdGuard{ savedCwd };
            boost::filesystem::current_path( cfgDir );

            bpo::variables_map vm;
            redux::logging::Logger logger;

            vector<Job::JobPtr> jobs = Job::parseTree( vm, tree, logger, false );

            py::list jobList;
            for( auto& job: jobs ) {
                auto mjob = dynamic_pointer_cast<MomfbdJob>( job );
                if( !mjob ) continue;

                py::dict jd;
                jd["name"] = mjob->info.name;
                jd["outputDir"] = mjob->info.outputDir;
                jd["logFile"] = mjob->info.logFile;

                py::list objList;
                for( auto& obj: mjob->getObjects() ) {
                    py::dict od;
                    od["outputFileName"] = obj->outputFileName;
                    od["wavelength"] = obj->wavelength;
                    od["telescopeF"] = obj->telescopeF;
                    od["arcSecsPerPixel"] = obj->arcSecsPerPixel;
                    od["pixelSize"] = obj->pixelSize;

                    py::list chanList;
                    for( auto& chan: obj->getChannels() ) {
                        py::dict cd;
                        cd["imageDataDir"] = chan->imageDataDir;
                        cd["imageTemplate"] = chan->imageTemplate;
                        cd["fileNumbers"] = chan->fileNumbers;
                        cd["waveFrontList"] = chan->waveFrontList;
                        cd["imageNumberOffset"] = chan->imageNumberOffset;
                        chanList.append( cd );
                    }
                    od["channels"] = chanList;
                    objList.append( od );
                }
                jd["objects"] = objList;
                jobList.append( jd );
            }

            py::dict result;
            result["jobs"] = jobList;
            return result;
        }

        void init_cfg( py::module_& m ) {
            m.def( "parse_cfg", &parse_cfg,
                   "Parse a redux cfg file the same way `rdx_sub -P` does. Returns "
                   "{jobs: [{name, outputDir, logFile, objects: [{outputFileName, wavelength, "
                   "telescopeF, arcSecsPerPixel, pixelSize, channels: [{imageDataDir, imageTemplate, "
                   "fileNumbers, waveFrontList, imageNumberOffset}, ...]}, ...]}, ...]}. "
                   "check=false means that the verification of on-disk data is skipped.",
                   py::arg( "path" ), py::arg( "name" ) = "", py::arg( "logfile" ) = "" );
        }

    }  // namespace python

}  // namespace redux
