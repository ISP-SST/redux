#include "idlutil.hpp"

#include <redux/file/fileana.hpp>
#include <redux/file/filefits.hpp>
#include <redux/util/stringutil.hpp>

#include <atomic>

#include <boost/asio.hpp>
#include <boost/thread/thread.hpp>
#include <boost/filesystem.hpp>

using namespace redux::file;
using namespace redux::util;
using namespace redux;

using namespace std;
namespace bfs = boost::filesystem;


namespace {

    typedef struct {
        IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
        IDL_INT help;
        UCHAR nthreads;
        IDL_VPTR times;
        IDL_INT verbose;
    } KW_LOADFILES;


    // NOTE:  The keywords MUST be listed in alphabetical order !!
    static IDL_KW_PAR kw_loadfiles_pars[] = {
        IDL_KW_FAST_SCAN,
        { (char*) "HELP",             IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(KW_LOADFILES,help) },
        { (char*) "NTHREADS",         IDL_TYP_BYTE,  1, 0,                      0, (char*) IDL_KW_OFFSETOF2(KW_LOADFILES,nthreads) },
        { (char*) "TIMES",            IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(KW_LOADFILES,times) },
        { (char*) "VERBOSE",          IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(KW_LOADFILES,verbose) },
        { NULL }
    };

}


string loadfiles_info( int lvl ) {
    string ret = "RDX_LOADFILES";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"      ");          // newline if lvl>1
        ret += "   Syntax:   out = rdx_loadfiles(file_list, /KEYWORDS)\n";
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      NTHREADS            Number of threads.\n"
                    "      TIMES               (output) Array with timestamps from headers.\n"
                    "      VERBOSE             Verbosity, default is 0 (only error output).\n";
        }
    } else ret += "\n";
    return ret;
}


IDL_VPTR loadfiles( int argc, IDL_VPTR* argv, char* argk ) {

    KW_LOADFILES kw;
    kw.nthreads = std::thread::hardware_concurrency();
    int nPlainArgs = IDL_KWProcessByOffset (argc, argv, argk, kw_loadfiles_pars, (IDL_VPTR*)0, 255, &kw);

    if( nPlainArgs < 1 ) {
        return IDL_GettmpInt(0);
    }

    IDL_VPTR filelist = argv[0];
    IDL_VPTR ret;
    IDL_ENSURE_STRING( filelist )
    IDL_ENSURE_SIMPLE( filelist );

    try {

        vector<string> filenames;
        if ( !(filelist->flags & IDL_V_ARR) ) {
            bfs::path fn( string(filelist->value.str.s) );
            if( bfs::is_regular_file(fn) ) {
                filenames.push_back( fn.string() );
            } else return IDL_GettmpInt(0);
        } else {
            IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>(filelist->value.arr->data);
            for( int i=0; i<filelist->value.arr->n_elts; ++i) {
                bfs::path fn( string(strptr[i].s) );
                if( bfs::is_regular_file(fn) ) {
                    filenames.push_back( fn.string() );
                }
            }
        }

        size_t nImages = filenames.size();
        if( nImages == 0 ) {
            return IDL_GettmpInt(0);
        }

        kw.nthreads = max<UCHAR>(1, min<UCHAR>(kw.nthreads, thread::hardware_concurrency()));

        if( kw.help ) {
            cout << loadfiles_info(2) << endl;
            return IDL_GettmpInt(0);
        }

        // Get size etc from first image.
        FileMeta::Ptr meta = getMeta( filenames[0] );
        if( !meta || (meta->nDims() != 2) ) {        // Only allow images for now
            cout << "rdx_loadfiles: Failed to get meta, or not 2D input." << endl;
            return IDL_GettmpInt(0);
        }

        string statusString;
        if( kw.verbose ) {
            statusString = "Loading " + to_string(nImages)
                           + " files using " +to_string((int)kw.nthreads) + string(" thread") + ((kw.nthreads>1)?"s.":".");
            cout << statusString << ((kw.verbose == 1)?"\n":"") << flush;
        }

        unique_ptr<double[]> times;
        double* timesPtr = nullptr;
        if( kw.times && nImages > 1 ) {
            times.reset( new double[nImages] );
            timesPtr = times.get();
        }

        IDL_MEMINT dims[] = { static_cast<IDL_MEMINT>(meta->dimSize(1)),
                              static_cast<IDL_MEMINT>(meta->dimSize(0)),
                              static_cast<IDL_MEMINT>(nImages)
                            };

        char* rawPtr;
        if( nImages == 1 ) {
            rawPtr = IDL_MakeTempArray( meta->getIDLType(), 2, dims, IDL_ARR_INI_NOP, &ret );
            readFile( filenames[0], rawPtr, meta );
            if( kw.times ) {
                IDL_ALLTYPES tmp;
                tmp.d = meta->getAverageTime().time_of_day().total_nanoseconds()*1E-9;
                IDL_StoreScalar( kw.times, IDL_TYP_DOUBLE, &tmp );
            }
            return ret;
        } else {
            if( kw.times ) {
                IDL_VPTR tmp;
                (void)IDL_MakeTempArray( IDL_TYP_DOUBLE, 1, dims+2, IDL_ARR_INI_ZERO, &tmp );
                IDL_VarCopy( tmp, kw.times );
                timesPtr = (double*)kw.times->value.arr->data;
            }
            rawPtr = IDL_MakeTempArray( meta->getIDLType(), 3, dims, IDL_ARR_INI_NOP, &ret );
        }

        size_t frameSize = meta->dataSize();

        atomic<size_t> nLoaded(0);
        postLoadCallback postLoad = [=,&nLoaded]( char* data, size_t i, FileMeta::Ptr& meta ) {
            if( timesPtr ) timesPtr[ i ] = meta->getAverageTime().time_of_day().total_nanoseconds()*1E-9;
            size_t nl = nLoaded++;
            if( kw.verbose > 1 ) printProgress( statusString, (nl*100.0/(nImages-1)));
        };

        loadFiles( filenames, rawPtr, frameSize, kw.nthreads, postLoad );

        if( kw.verbose > 1 ) {
            cout << endl;
        }

    } catch (const exception& e ) {
        cout << "rdx_loadfiles: unhandled exception: " << e.what() << endl;
        return IDL_GettmpInt(0);
    }

    return ret;

}


namespace {

typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_INT all;
    IDL_VPTR date_beg;
    IDL_VPTR header;
    IDL_INT help;
    IDL_VPTR status;
    IDL_INT raw;
    // IDL_STRING split_chars;
} KW_READDATA;

// NOTE:  The keywords MUST be listed in alphabetical order !!
static IDL_KW_PAR kw_readdata_pars[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "ALL",        IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(KW_READDATA,all) },
    { (char*) "DATE_BEG",   IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(KW_READDATA,date_beg) },
    { (char*) "HEADER",     IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(KW_READDATA,header) },
    { (char*) "HELP",       IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(KW_READDATA,help) },
    { (char*) "RAW",        IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(KW_READDATA,raw) },
    { (char*) "STATUS",     IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(KW_READDATA,status) },
    { NULL }
};

}


string readdata_info( int lvl ) {
    string ret = "RDX_READDATA";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"       ");          // newline if lvl>1
        ret += "   Syntax:   out = rdx_readdata(filename, /KEYWORDS)\n";
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      DATE_BEG            Extract tabulated frame timestamps.\n"
                    "      HELP                Display this info.\n"
                    "      HEADER              Return metadata.\n"
                    "      ALL                 Return all metadata.\n"
                    "      RAW                 Return metadata exactly as it is in the file, without manipulation.\n";
        }
    } else ret += "\n";
    return ret;
}


IDL_VPTR readdata( int argc, IDL_VPTR* argv, char* argk ) {

    IDL_VPTR filenames = argv[0];
    IDL_ENSURE_SIMPLE( filenames );

    if( filenames->type != IDL_TYP_STRING ) {
        return IDL_GettmpInt(0);
    }

    KW_READDATA kw;
    kw.header = nullptr;
    kw.status = nullptr;
    (void)IDL_KWProcessByOffset( argc, argv, argk, kw_readdata_pars, (IDL_VPTR*)0, 255, &kw );

    if( kw.help ) {
        cout << readdata_info(2) << endl;
        return IDL_GettmpInt(0);
    }

    //kw.nthreads = max<UCHAR>(1, min<UCHAR>(kw.nthreads, thread::hardware_concurrency()));
    if( kw.status ) {
        IDL_VarCopy( IDL_GettmpInt(-1), kw.status );
    }

    try {

        IDL_VPTR ret;
        shared_ptr<double> darkData, gainData, bsGainData;
        shared_ptr<uint8_t> maskData;
        shared_ptr<uint8_t*> mask2D;
        shared_ptr<fftw_complex> bsOtfData;

        vector<string> existingFiles;
        if ( !(filenames->flags & IDL_V_ARR) ) {
            bfs::path fn( string(filenames->value.str.s) );
            if( bfs::is_regular_file(fn) ) {
                existingFiles.push_back( fn.string() );
            } else return IDL_GettmpInt(0);
        } else {
            IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>(filenames->value.arr->data);
            for( int i=0; i<filenames->value.arr->n_elts; ++i ) {
                bfs::path fn( string(strptr[i].s) );
                if( bfs::is_regular_file(fn) ) {
                    existingFiles.push_back( fn.string() );
                }
            }
        }
        size_t nImages = existingFiles.size();

        if( !nImages ) {
            cout << "rdx_readdata: No input files." << endl;
            return IDL_GettmpInt(0);
        } else if ( nImages > 1 ) {
            cout << "rdx_readdata: Only a single file supported at the moment." << endl;
            return IDL_GettmpInt(0);
        }

        try {
            FileMeta::Ptr myMeta = getMeta( existingFiles[0] );
            if( kw.header ) {
                IDL_VPTR tmpHdr;
                myMeta->getAverageTime();
                myMeta->getEndTime();
                myMeta->getStartTime();
                bool raw = (kw.raw != 0);
                vector<string> hdrTexts = myMeta->getText( raw );
                int nTexts = hdrTexts.size();
                if( kw.all == 0 ) nTexts = 1;
                if( nTexts ) {
                    size_t charCount(0);
                    for( int i=0; i<nTexts; ++i ) charCount += hdrTexts[i].size();
                    if( charCount%80 ) {
                        throw logic_error("Header text-size is not a multiple of 80.");
                    }
                    IDL_MEMINT nKeys = charCount/80;
                    IDL_MEMINT dims[] = { nKeys };
                    IDL_MakeTempArray( IDL_TYP_STRING, 1, dims, IDL_ARR_INI_ZERO, &tmpHdr );
                    IDL_STRING* strPtr = reinterpret_cast<IDL_STRING*>(tmpHdr->value.arr->data);
                    size_t cnt(0);
                    for( int i=0; i<nTexts; ++i ) {
                        string hdr = hdrTexts[i];
                        while( !hdr.empty() ) {
                            string tmp = hdr.substr( 0, 80 );
                            hdr.erase( 0, 80 );
                            IDL_StrStore( &(strPtr[cnt++]), (char*)tmp.c_str() );
                        }
                    }
                } else tmpHdr = IDL_StrToSTRING((char*)"");
                IDL_VarCopy( tmpHdr, kw.header );
            }
            if( kw.date_beg ) {
                vector<boost::posix_time::ptime> date_beg = myMeta->getStartTimes();
                if( date_beg.empty() ) date_beg.push_back( myMeta->getStartTime() );
                IDL_VPTR tmp;
                IDL_MEMINT nDB = date_beg.size();
                IDL_MEMINT dims[] = { nDB };
                IDL_MakeTempArray( IDL_TYP_STRING, 1, dims, IDL_ARR_INI_ZERO, &tmp );
                IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>( tmp->value.arr->data );
                for( int j=0; j<nDB; ++j ) {
                    string tStr = bpx::to_iso_extended_string( date_beg[j] );
                    IDL_StrStore( strptr++, const_cast<char*>(tStr.c_str()) );
                }
                IDL_VarCopy( tmp, kw.date_beg );
            }
            size_t nDims = myMeta->nDims();
            if( !nDims ) {
                cerr << "rdx_readdata: file contains no data." << endl;
                return IDL_GettmpInt(0);
            }
            vector<IDL_MEMINT> dims;

            for( size_t i=0; i<nDims; ++i ) {
                dims.push_back( myMeta->dimSize(i) );
            }
            std::reverse( dims.begin(), dims.end() );
            char* data = (char*)IDL_MakeTempArray( myMeta->getIDLType(), nDims, dims.data(), IDL_ARR_INI_NOP, &ret ); //IDL_ARR_INI_ZERO
            readFile( existingFiles[0], data, myMeta );
        } catch (const exception& e ) {
            cout << "rdx_readdata: unhandled exception: " << e.what() << endl;
            return IDL_GettmpInt(0);
        }

        if( kw.status ) {
            IDL_VarCopy( IDL_GettmpInt(0), kw.status );
        }
        return ret;

    } catch (const exception& e ) {
        cout << "rdx_readdata: unhandled exception: " << e.what() << endl;
        return IDL_GettmpInt(0);
    }


}


namespace {

typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_INT all;
    IDL_INT help;
    IDL_VPTR status;
    IDL_INT raw;
    IDL_VPTR date_beg;
    IDL_VPTR framenumbers;
    // IDL_STRING split_chars;
} KW_READHEAD;

// NOTE:  The keywords MUST be listed in alphabetical order !!
static IDL_KW_PAR kw_readhead_pars[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "ALL",            IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(KW_READHEAD,all) },
    { (char*) "DATE_BEG",       IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(KW_READHEAD,date_beg) },
    { (char*) "FRAMENUMBERS",   IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(KW_READHEAD,framenumbers) },
    { (char*) "HELP",           IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(KW_READHEAD,help) },
    { (char*) "RAW",            IDL_TYP_INT,   1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(KW_READHEAD,raw) },
    { (char*) "STATUS",         IDL_TYP_UNDEF, 1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(KW_READHEAD,status) },
    { NULL }
};

}


string readhead_info( int lvl ) {
    string ret = "RDX_READHEAD";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"       ");          // newline if lvl>1
        ret += "   Syntax:   out = rdx_readhead(filename, /KEYWORDS)\n";
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      ALL                 Return all metadata.\n"
                    "      FRAMENUMBERS        (output) Get the framenumbers.\n"
                    "      DATE_BEG            (output) Get the start-times.\n"
                    "      RAW                 Return metadata exactly as it is in the file, without manipulation.\n"
                    "      STATUS              (output) Status flag.\n";
        }
    } else ret += "\n";
    return ret;
}


IDL_VPTR readhead( int argc, IDL_VPTR* argv, char* argk ) {

    IDL_VPTR filenames = argv[0];
    IDL_ENSURE_SIMPLE( filenames );

    if( filenames->type != IDL_TYP_STRING ) {
        cout << readhead_info(2) << endl;
        return IDL_GettmpInt(0);
    }

    KW_READHEAD kw;
    kw.status = nullptr;
    (void)IDL_KWProcessByOffset( argc, argv, argk, kw_readhead_pars, (IDL_VPTR*)0, 255, &kw );

    if( kw.help ) {
        cout << readhead_info(2) << endl;
        return IDL_GettmpInt(0);
    }

    //kw.nthreads = max<UCHAR>(1, min<UCHAR>(kw.nthreads, thread::hardware_concurrency()));
    if( kw.status ) {
        IDL_VarCopy( IDL_GettmpInt(-1), kw.status );
    }

    try {

        IDL_VPTR ret = IDL_GettmpInt(0);

        vector<string> existingFiles;
        if ( !(filenames->flags & IDL_V_ARR) ) {
            bfs::path fn( string(filenames->value.str.s) );
            if( bfs::is_regular_file(fn) ) {
                existingFiles.push_back( fn.string() );
            }
        } else {
            IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>(filenames->value.arr->data);
            for( int i=0; i<filenames->value.arr->n_elts; ++i ) {
                bfs::path fn( string(strptr[i].s) );
                if( bfs::is_regular_file(fn) ) {
                    existingFiles.push_back( fn.string() );
                }
            }
        }
        size_t nImages = existingFiles.size();
        //kw.nthreads = static_cast<UCHAR>( min<size_t>(kw.nthreads, nImages) );

        if( !nImages ) {
            cout << "rdx_readhead: No input files." << endl;
            return IDL_GettmpInt(0);
        } else if ( nImages > 1 ) {
            cout << "rdx_readhead: Only a single file supported at the moment." << endl;
            return IDL_GettmpInt(0);
        }

        try {
            FileMeta::Ptr myMeta = getMeta( existingFiles[0] );
            IDL_VPTR tmpHdr;
            myMeta->getAverageTime();
            myMeta->getEndTime();
            myMeta->getStartTime();
            bool raw = (kw.raw != 0);
            vector<string> hdrTexts = myMeta->getText( raw );
            int nTexts = hdrTexts.size();
            if( nTexts ) {
                if( kw.all == 0 ) nTexts = 1;
                if( raw ) {
                    string tmp;
                    for( auto& t: hdrTexts ) tmp += t;
                    tmpHdr = IDL_StrToSTRING( (char*)tmp.c_str() );
                } else {
                    size_t charCount(0);
                    for( int i=0; i<nTexts; ++i ) charCount += hdrTexts[i].size();
                    IDL_MEMINT nKeys = charCount/80;
                    if( charCount%80 ) {
                        nKeys++;
                        //throw logic_error("Header text-size is not a multiple of 80.");
                    }
                    IDL_MEMINT dims[] = { nKeys };
                    IDL_MakeTempArray( IDL_TYP_STRING, 1, dims, IDL_ARR_INI_ZERO, &tmpHdr );
                    IDL_STRING* strPtr = reinterpret_cast<IDL_STRING*>(tmpHdr->value.arr->data);
                    size_t cnt(0);
                    for( int i=0; i<nTexts; ++i ) {
                        string hdr = hdrTexts[i];
                        while( !hdr.empty() ) {
                            string tmp = hdr.substr( 0, 80 );
                            hdr.erase( 0, 80 );
                            IDL_StrStore( &(strPtr[cnt++]), (char*)tmp.c_str() );
                        }
                    }
                }
            } else tmpHdr = IDL_StrToSTRING((char*)"");
            IDL_VarCopy( tmpHdr, ret );
            
            if( kw.framenumbers ) {
                vector<size_t> frameNumbers = myMeta->getFrameNumbers();
                IDL_VPTR tmp;
                IDL_MEMINT nFN = frameNumbers.size();
                IDL_MEMINT dims[] = { nFN }; 
                int32_t* tmpData = (int32_t*)IDL_MakeTempArray( IDL_TYP_LONG, 1, dims, IDL_ARR_INI_NOP, &tmp );
                std::copy( frameNumbers.begin(), frameNumbers.end(), tmpData );
                IDL_VarCopy( tmp, kw.framenumbers );
            }
            
            if( kw.date_beg ) {
                vector<boost::posix_time::ptime> date_beg = myMeta->getStartTimes();
                if( date_beg.empty() ) date_beg.push_back( myMeta->getStartTime() );
                IDL_VPTR tmp;
                IDL_MEMINT nDB = date_beg.size();
                IDL_MEMINT dims[] = { nDB };
                IDL_MakeTempArray( IDL_TYP_STRING, 1, dims, IDL_ARR_INI_ZERO, &tmp );
                IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>( tmp->value.arr->data );
                for( int j=0; j<nDB; ++j ) {
                    string tStr = bpx::to_iso_extended_string( date_beg[j] );
                    IDL_StrStore( strptr++, const_cast<char*>(tStr.c_str()) );
                }
                IDL_VarCopy( tmp, kw.date_beg );
            }


        } catch (const exception& e ) {
            cout << "rdx_readhead: unhandled exception: " << e.what() << endl;
            return IDL_GettmpInt(0);
        }

        if( kw.status ) {
            IDL_VarCopy( IDL_GettmpInt(0), kw.status );
        }
        return ret;

    } catch (const exception& e ) {
        cout << "rdx_readhead: unhandled exception: " << e.what() << endl;
        return IDL_GettmpInt(0);
    }


}


namespace {
static int dummy RDX_UNUSED =
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)loadfiles, (char*)"RDX_LOADFILES", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, loadfiles_info ) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)readdata, (char*)"RDX_READDATA", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, readdata_info ) +
    IdlContainer::registerRoutine( {(IDL_SYSRTN_GENERIC)readhead, (char*)"RDX_READHEAD", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, readhead_info );
}
