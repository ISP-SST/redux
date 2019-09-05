#include "rdx.hpp"

#include "redux/file/fileio.hpp"
#include "redux/image/utils.hpp"
#include "redux/util/cache.hpp"
#include "redux/version.hpp"

#include <iomanip>
#include <fstream>
#include <typeinfo>
#include <zlib.h>
#include <omp.h>

using namespace redux;
using namespace redux::file;
using namespace redux::image;
using namespace redux::util;
using namespace std;


void redux::rdx_cache_gdl( EnvT* e ) {

    SizeT nParam RDX_UNUSED = e->NParam( 2 );
    
    DStringGDL* tags = e->GetParAs<DStringGDL>(0);
    SizeT nEl = tags->N_Elements();
    if ( nEl > 1 ) {    // TODO: implement support for multiple files
        e->Throw( "redux::rdx_cache_gdl: Only a single variable supported at the moment." );
        //return new DStringGDL("");
    }
    
    BaseGDL*& par = e->GetPar(1);
    if( par == nullptr ) {
        return;
    }

    auto map = Cache::get().getMap< string, std::shared_ptr<BaseGDL> >();
    for( SizeT i(0); i<nEl; ++i ) {
        BaseGDL* tmp = par->NewResult();
        tmp->AssignAt( par );
        map.second[ (*tags )[i] ].reset( tmp );
    }

}


void redux::rdx_cacheclear_gdl( EnvT* ) {
    Cache::cleanup();
}


void redux::rdx_cachedel_gdl( EnvT* e ) {
    
    SizeT nParam RDX_UNUSED = e->NParam( 1 );

    static int countIx = e->KeywordIx("COUNT");
    bool get_count = e->KeywordPresent( countIx );

    DStringGDL* tags = e->GetParAs<DStringGDL>(0);
    SizeT nEl = tags->N_Elements();

    int total_found(0);
    auto map = Cache::get().getMap< string, std::shared_ptr<BaseGDL> >();
    for( SizeT i(0); i<nEl; ++i ) {
        total_found += map.second.erase( (*tags )[i] );
    }
    
    if( get_count ) {
        DLongGDL* count = new DLongGDL(0);
        (*count)[0] = total_found;
        e->SetKW( countIx, count );
    }

}


void redux::rdx_cacheinfo_gdl( EnvT* e ) {
    SizeT nParam RDX_UNUSED = e->NParam( 0 );
    auto map = Cache::get().getMap< string, std::shared_ptr<BaseGDL> >();
    string retString = "rdx_cache contains " + to_string( map.second.size()) + " entries.\n";
    uint64_t sz(0);
    for( auto& p: map.second ) {
        sz += 1; //getVarSize(p.second.get());
    }
    retString += "  Total size: " + to_string(sz) + " bytes.";
    cout << retString << endl;
    //DStringGDL* ret = new DStringGDL(retString);
}


BaseGDL* redux::rdx_cacheget_gdl( EnvT* e ) {

    SizeT nParam RDX_UNUSED = e->NParam( 1 );
    
    static int countIx = e->KeywordIx("COUNT");
    bool get_count = e->KeywordPresent( countIx );

    DStringGDL* tags = e->GetParAs<DStringGDL>(0);
    SizeT nEl = tags->N_Elements();
    
    if ( nEl > 1 ) {    // TODO: implement support for multiple files
        e->Throw( "redux::rdx_cacheget_gdl: Only a single variable supported at the moment." );
        //return new DStringGDL("");
    }
    
    BaseGDL* ret = nullptr;
    int total_found(0);
    auto map = Cache::get().getMap< string, std::shared_ptr<BaseGDL> >();
    for( SizeT i(0); i<nEl; ++i ) {
        auto it = map.second.find( (*tags )[i] );
        int found = (it != map.second.end());
        if( found && it->second ) {
            ret = it->second.get()->NewResult();
            ret->AssignAt( it->second.get() );
            //ret = it->second.get()->Dup();
        } else if( !get_count ) {  // only warn if keyword found was not passed
            e->Throw( "redux::rdx_cacheget_gdl: "+(*tags )[i]+" not found." );
        }
        total_found += found;
    }
    
    if( get_count ) {
        DLongGDL* count = new DLongGDL(0);
        (*count)[0] = total_found;
        e->SetKW( countIx, count );
    }

    return ret;
    
}

/*string filetype_info( int lvl ) {
    
    string ret = "RDX_FILETYPE";
    if( lvl > 0 ) {
        ret += "          Syntax:   fmt = rdx_filetype()\n";
    } else ret += "\n";

    return ret;
    
}*/


BaseGDL* redux::rdx_filetype_gdl( EnvT* e ) {

    SizeT nParam RDX_UNUSED = e->NParam( 1 );
    DStringGDL* files = e->GetParAs<DStringGDL>(0);
    SizeT nEl = files->N_Elements();
    DStringGDL* ret = new DStringGDL( files->Dim(), BaseGDL::NOZERO );

    for( SizeT i(0); i<nEl; ++i ) {
        (*ret)[ i] = redux::file::filetype( (*files )[i] );
    }

    return ret;
    
}


namespace {
    const string fillpix_help[]={
        "RDX_FILLPIX:        out = rdx_fillpix(image(s)[,/KEYWORDS)]",
        "   Accepted Keywords:",
        "      HELP                Display this info.",
        "      MASK                Mask with non-zero values where filling is supposed to be performed.",
        "      NTHREADS            Number of threads.",
        "      THRESHOLD           Value below which a pixel qualifies for filling. (default: 0)",
        "      VERBOSE             Verbosity, default is 0 (only error output)."
    };
    constexpr int fillpix_help_sz = sizeof(fillpix_help) / sizeof(string);
}

BaseGDL* redux::rdx_fillpix_gdl( EnvT* e ) {

    static int helpIx = e->KeywordIx("HELP");
    static int maskIx = e->KeywordIx("MASK");
    static int ntIx = e->KeywordIx("NTHREADS");
    static int thrIx = e->KeywordIx("THRESHOLD");
    static int verbIx = e->KeywordIx("VERBOSE");
    
    if( e->KeywordSet( helpIx ) ) {
        e->Help( fillpix_help, fillpix_help_sz );
    }

    SizeT nParam RDX_UNUSED = e->NParam( 1 );

    DLong nThreads = std::thread::hardware_concurrency();
    e->AssureLongScalarKWIfPresent( ntIx, nThreads );

    DDouble threshold = 0;
    e->AssureDoubleScalarKWIfPresent( thrIx, threshold );

    DLong verbose = 0;
    e->AssureLongScalarKWIfPresent( verbIx, verbose );

    
    BaseGDL* data = e->GetPar(0);
    const dimension& dims = data->Dim();
    
    SizeT nDims = dims.Rank();
    if( nDims < 2 ) {
        e->Throw( "redux::rdx_fillpix_gdl: Input has to be one, or several 2D images." );
    }

    dimension imgDim = dims;
    SizeT nImages = 1;
    SizeT xSize = dims[0];
    SizeT ySize = dims[1];
    if( nDims == 3 ) {
        nImages = dims[2];
        imgDim.Remove(2);
    }
    size_t nPixels = xSize*ySize;
    
    BaseGDL* mask = nullptr;
    if( e->KeywordSet( maskIx ) ) {
        mask = e->GetKWAs<DByteGDL>( maskIx );
        if( mask && (mask->Dim() != imgDim) ){
            Warning("redux::rdx_fillpix_gdl: The size of the mask does not match the image(s). It will be ignored!");
            mask = nullptr;
        }
    }
    
    shared_ptr<uint8_t*> mask2D;
    if( mask ) {
        mask2D = reshapeArray( static_cast<uint8_t*>(mask->DataAddr()), ySize, xSize );
    }
    
    BaseGDL* ret = data->NewResult();
    ret->AssignAt( data );      // make a copy.
    void* retData = ret->DataAddr();
    namespace sp = std::placeholders;
    try {
        switch( data->Type() ) {
            case( GDL_BYTE ): {
                DByte* data = reinterpret_cast<DByte*>( retData );
#pragma omp parallel for
                for( SizeT i=0; i<nImages; ++i ) {
                    auto images2D = reshapeArray( data+i*nPixels, ySize, xSize );
                    function<DByte (SizeT, SizeT) > func = bind (horizontalInterpolation<DByte>, images2D.get(), ySize, xSize, sp::_1, sp::_2);
                    fillPixels( images2D.get(), ySize, xSize, func, std::bind2nd (std::less_equal<double>(), threshold), mask2D.get());
                }
                break;
            }
            case( GDL_INT ): {
                DInt* data = reinterpret_cast<DInt*>( retData );
#pragma omp parallel for
                for( SizeT i=0; i<nImages; ++i ) {
                    auto images2D = reshapeArray( data+i*nPixels, ySize, xSize );
                    function<DInt (SizeT, SizeT) > func = bind (horizontalInterpolation<DInt>, images2D.get(), ySize, xSize, sp::_1, sp::_2);
                    fillPixels( images2D.get(), ySize, xSize, func, std::bind2nd (std::less_equal<double>(), threshold), mask2D.get());
                }
                break;
            }
            case( GDL_LONG ): {
                DLong* data = reinterpret_cast<DLong*>( retData );
#pragma omp parallel for
                for( SizeT i=0; i<nImages; ++i ) {
                    auto images2D = reshapeArray( data+i*nPixels, ySize, xSize );
                    function<DLong (SizeT, SizeT) > func = bind (horizontalInterpolation<DLong>, images2D.get(), ySize, xSize, sp::_1, sp::_2);
                    fillPixels( images2D.get(), ySize, xSize, func, std::bind2nd (std::less_equal<double>(), threshold), mask2D.get());
                }
                break;
            }
            case( GDL_FLOAT ): {
                DFloat* data = reinterpret_cast<DFloat*>( retData );
#pragma omp parallel for
                for( SizeT i=0; i<nImages; ++i ) {
                    auto images2D = reshapeArray( data+i*nPixels, ySize, xSize );
                    function<DFloat (SizeT, SizeT) > func = bind (horizontalInterpolation<DFloat>, images2D.get(), ySize, xSize, sp::_1, sp::_2);
                    fillPixels( images2D.get(), ySize, xSize, func, std::bind2nd (std::less_equal<double>(), threshold), mask2D.get());
                }
                break;
            }
            case( GDL_DOUBLE ): {
                DDouble* data = reinterpret_cast<DDouble*>( retData );
#pragma omp parallel for
                for( SizeT i=0; i<nImages; ++i ) {
                    auto images2D = reshapeArray( data+i*nPixels, ySize, xSize );
                    function<DDouble (SizeT, SizeT) > func = bind (horizontalInterpolation<DDouble>, images2D.get(), ySize, xSize, sp::_1, sp::_2);
                    fillPixels( images2D.get(), ySize, xSize, func, std::bind2nd (std::less_equal<double>(), threshold), mask2D.get());
                }
                break;
            }
            default: throw std::runtime_error(" Type="+data->TypeStr()+" not implemented.");
        }

    } catch( const exception& ex ) {
        e->Throw(string("rdx_fillpix_gdl: unhandled exception: ") + ex.what());
    }
    
    return ret;

}


BaseGDL* redux::rdx_hasopencv_gdl( EnvT* ) {
    
#ifdef RDX_WITH_OPENCV
    return new DLongGDL( 1 );
#else
    return new DLongGDL( 0 );
#endif
    
}


namespace {
    const string readdata_help[]={
        "RDX_READDATA:       data = rdx_readdata(filename, /KEYWORDS)",
        "   Accepted Keywords:",
        "      ALL                 Return all metadata.",
        "      DATE_BEG            Extract tabulated frame timestamps.",
        "      FRAMENUMBERS        (output) Get the framenumbers.",
        "      HEADER              (output) Return metadata.",
        "      HELP                Display this info.",
        "      RAW                 Return metadata exactly as it is in the file, without manipulation.",
        "      STATUS              (output) Status flag."
    };
    constexpr int readdata_help_sz = sizeof(readdata_help) / sizeof(string);
}

BaseGDL* redux::rdx_readdata_gdl( EnvT* e ) {
    
    static int allIx = e->KeywordIx("ALL");
    static int dbIx = e->KeywordIx("DATE_BEG");
    static int fnIx = e->KeywordIx("FRAMENUMBERS");
    static int headIx = e->KeywordIx("HEADER");
    static int helpIx = e->KeywordIx("HELP");
    static int rawIx = e->KeywordIx("RAW");
    static int statIx = e->KeywordIx("STATUS");
    
    bool helpSet = e->KeywordSet( helpIx );
    if( helpSet ) {
        e->Help( readdata_help, readdata_help_sz );
    }

    SizeT nParam RDX_UNUSED = e->NParam( 1 );
    DStringGDL* files = e->GetParAs<DStringGDL>(0);
    SizeT nEl = files->N_Elements();

    if ( nEl > 1 ) {    // TODO: implement support for multiple files
        e->Throw( "redux::rdx_readdata_gdl: Only a single file supported at the moment." );
    }
    
    bool allSet = e->KeywordSet( allIx );
    bool dbSet = e->KeywordPresent( dbIx );
    bool fnSet = e->KeywordPresent( fnIx );
    bool headSet = e->KeywordPresent( headIx );
    bool rawSet = e->KeywordSet( rawIx );
    bool statusSet = e->KeywordPresent(statIx);
    
    DLongGDL* status = nullptr;
    if( statusSet ) {
        status = new DLongGDL(0);
        (*status)[0] = -1;
        e->SetKW( statIx, status );
    }

    BaseGDL* ret = nullptr;
    vector<string> texts;
    for( SizeT i(0); i<nEl; ++i ) {
        FileMeta::Ptr meta;
        if( headSet ) {
            if( !rawSet ) {
                texts = redux::file::getMetaTextAsCards( (*files)[i], meta, rawIx, allSet );
            } else {
                texts.push_back(redux::file::getMetaText( (*files)[i], meta, rawIx, allSet ));
            }
        } else {
            meta = getMeta( (*files)[i] );
        }
        if( fnSet && meta ) {
            vector<size_t> frameNumbers = meta->getFrameNumbers();
            SizeT nFN = frameNumbers.size();
            DLongGDL* frame_numbers = new DLongGDL( nFN, BaseGDL::NOZERO );
            for( SizeT i(0); i<nFN; ++i ) {
                (*frame_numbers)[i] = frameNumbers[i];
            }
            e->SetKW( fnIx, frame_numbers );
        }
        
        if( dbSet && meta ) {
            meta->getAverageTime();
            meta->getEndTime();
            vector<bpx::ptime> date_beg = meta->getStartTimes();
            if( date_beg.empty() ) date_beg.push_back( meta->getStartTime() );
            SizeT nDB = date_beg.size();
            DStringGDL* db = new DStringGDL( nDB, BaseGDL::NOZERO );
            for( SizeT i(0); i<nDB; ++i ) {
                //string tStr = bpx::to_iso_extended_string( date_beg[j] );
                (*db)[i] = bpx::to_iso_extended_string( date_beg[i] );
            }
            e->SetKW( dbIx, db );
        }

        size_t nDims = meta->nDims();
        if( !nDims ) {
            e->Throw( "redux::rdx_readdata_gdl: File contains no data: " + (*files)[i] );
            //Warning( "redux::rdx_readdata_gdl: File contains no data: " + (*files)[i] );
            //ret = new DLongGDL( 0 );
            //return ret;
        }
        vector<SizeT> fileDims;
        for( size_t i=0; i<nDims; ++i ) {
            fileDims.push_back( meta->dimSize(i) );
        }
        std::reverse( fileDims.begin(), fileDims.end() );   // IDL has reversed dimensions
        
        dimension gdlDims( fileDims.data(), nDims );

        int dataTYpe = meta->getIDLType();
        void* data = nullptr;
        switch( dataTYpe ) {
            case( GDL_BYTE ): ret = new DByteGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_INT ): ret = new DIntGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_LONG ): ret = new DLongGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_FLOAT ): ret = new DFloatGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_DOUBLE ): ret = new DDoubleGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_COMPLEX ): ret = new DComplexGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_COMPLEXDBL ): ret = new DComplexDblGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_UINT ): ret = new DUIntGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_ULONG ): ret = new DLongGDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_LONG64 ): ret = new DLong64GDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            case( GDL_ULONG64 ): ret = new DULong64GDL( gdlDims, BaseGDL::NOZERO ); data = ret->DataAddr(); break;
            default: e->Throw( "redux::rdx_readdata_gdl: Unrecognized dataType: " + (*files)[i] );
        }
        readFile( (*files)[i], reinterpret_cast<char*>(data), meta );
    }
    
    if( headSet ) {
        DStringGDL* hdr = new DStringGDL( texts.size(), BaseGDL::NOZERO );
        for( SizeT i(0); i<texts.size(); ++i ) {
            (*hdr)[i] = texts[i];
        }
        e->SetKW( headIx, hdr );
    }
    
    if( statusSet ) {
        (*status)[0] = 0;
    }

    return ret;    


}


namespace {
    const string readhead_help[]={
        "RDX_READHEAD:       hdr = rdx_readhead(filename, /KEYWORDS)",
        "   Accepted Keywords:",
        "      ALL                 Return all metadata.",
        "      DATE_BEG            Extract tabulated frame timestamps.",
        "      FRAMENUMBERS        (output) Get the framenumbers.",
        "      HELP                Display this info.",
        "      RAW                 Return metadata exactly as it is in the file, without manipulation.",
        "      STATUS              (output) Status flag."
    };
    constexpr int readhead_help_sz = sizeof(readhead_help) / sizeof(string);
}

BaseGDL* redux::rdx_readhead_gdl( EnvT* e ) {

    static int allIx = e->KeywordIx("ALL");
    static int dbIx = e->KeywordIx("DATE_BEG");
    static int fnIx = e->KeywordIx("FRAMENUMBERS");
    static int helpIx = e->KeywordIx("HELP");
    static int rawIx = e->KeywordIx("RAW");
    static int statIx = e->KeywordIx("STATUS");
    
    if( e->KeywordSet( helpIx ) ) {
        e->Help( readhead_help, readhead_help_sz );
    }
    SizeT nParam RDX_UNUSED = e->NParam( 1 );
    DStringGDL* files = e->GetParAs<DStringGDL>(0);
    SizeT nEl = files->N_Elements();

    bool allSet = e->KeywordSet( allIx );
    bool dbSet = e->KeywordPresent( dbIx );
    bool fnSet = e->KeywordPresent( fnIx );
    bool rawSet = e->KeywordSet( rawIx );
    bool statusSet = e->KeywordPresent(statIx);
    
    DLongGDL* status = nullptr;
    if( statusSet ) {
        status = new DLongGDL(0);
        (*status)[0] = -1;
        e->SetKW( statIx, status );
    }
    
    if ( nEl > 1 ) {    // TODO: implement support for multiple files
        e->Throw( "redux::readhead_gdl: Only a single file supported at the moment." );
        //return new DStringGDL("");
    }

    vector<string> texts;
    if( nEl == 1 ) {
        FileMeta::Ptr meta;
        if( !rawSet ) {
            texts = redux::file::getMetaTextAsCards( (*files)[0], meta, rawIx, allSet );
        } else {
            texts.push_back(redux::file::getMetaText( (*files)[0], meta, rawIx, allSet ));
        }

        if( fnSet && meta ) {
            vector<size_t> frameNumbers = meta->getFrameNumbers();
            SizeT nFN = frameNumbers.size();
            DLongGDL* frame_numbers = new DLongGDL( nFN, BaseGDL::NOZERO );
            for( SizeT i(0); i<nFN; ++i ) {
                (*frame_numbers)[i] = frameNumbers[i];
            }
            e->SetKW( fnIx, frame_numbers );
        }
        
        if( dbSet && meta ) {
            vector<bpx::ptime> date_beg = meta->getStartTimes();
            if( date_beg.empty() ) date_beg.push_back( meta->getStartTime() );
            SizeT nDB = date_beg.size();
            DStringGDL* db = new DStringGDL( nDB, BaseGDL::NOZERO );
            for( SizeT i(0); i<nDB; ++i ) {
                //string tStr = bpx::to_iso_extended_string( date_beg[j] );
                (*db)[i] = bpx::to_iso_extended_string( date_beg[i] );
            }
            e->SetKW( dbIx, db );
        }
    } else if( rawSet ) {   // for raw headers we can return results for multiple files.
        texts.resize(nEl);
#pragma omp parallel for
        for( SizeT i=0; i<nEl; ++i ) {
            FileMeta::Ptr meta;
            texts[i] = redux::file::getMetaText( (*files)[i], meta, rawIx, allSet );
        }
    }
    DStringGDL* ret = new DStringGDL( texts.size(), BaseGDL::NOZERO );
    for( SizeT i(0); i<texts.size(); ++i ) {
        (*ret)[i] = texts[i];
    }
    
    if( statusSet ) {
        (*status)[0] = 0;
    }

    return ret;
    

}


namespace {
    const string segment_help[]={
        "RDX_SEGMENT:        result = rdx_segment(first,last,segment-size[,min-overlap][,/KEYWORDS])",
        "   Accepted Keywords:",
        "      HELP                Display this info.",
        "      MOMFBD              Interpret arguments as pixels and do segmentation identical to the MOMFBD code."
    };
    constexpr int segment_help_sz = sizeof(segment_help) / sizeof(string);
}

BaseGDL* redux::rdx_segment_gdl( EnvT* e ) {

    static int helpIx = e->KeywordIx("HELP");
    static int momfbdIx = e->KeywordIx("MOMFBD");
    
    if( e->KeywordSet( helpIx ) ) {
        e->Help( segment_help, segment_help_sz );
    }
    
    SizeT nParam = e->NParam( 1 );
    DLong first, last, sz, minOverlap;
    e->AssureLongScalarPar( 0, first );
    e->AssureLongScalarPar( 1, last );
    e->AssureLongScalarPar( 2, sz );
    minOverlap = 0;
    if( nParam > 3 ) {
        e->AssureLongScalarPar( 3, minOverlap );
    }
    
    if( e->KeywordSet( momfbdIx ) /*&& ((last-first) > 2*sz)*/ ) {
        first += sz/2;
        last -= sz/2;
        minOverlap = 16 + sz/4;
    }
    
    /*if( 2*minOverlap >= sz ) {
        if( nParam > 3 ) {
            sz = 3*minOverlap;
        } else {
            minOverlap = sz/3;
        }

    }*/

    std::vector<DLong> resV = segment<DLong>( first, last, sz, minOverlap );
    
    DLongGDL* ret = nullptr;
    if( resV.empty() ) {
        ret = new DLongGDL(-1);
    } else {
        ret = new DLongGDL( dimension(resV.size()) );
    }
    for( SizeT i(0); i<resV.size(); ++i ) {
        (*ret)[i] = resV[i];
    }
    
    return ret;
    
}


namespace {
    const string ints2str_help[]={
        "RDX_INTS2STR:       list = rdx_ints2str(integer_array, /KEYWORDS)",
        "   Accepted Keywords:",
        "      HELP                Display this info.",
        "      SORT                Sort input array.",
        "      UNIQUE              Discard repeated integers (will also force sort)",
        "      VERBOSE             Verbosity, default is 0 (only error output)."
    };
    constexpr int ints2str_help_sz = sizeof(ints2str_help) / sizeof(string);
}

BaseGDL* redux::rdx_ints2str_gdl( EnvT* e ) {
    
    static int helpIx = e->KeywordIx("HELP");
    static int sortIx = e->KeywordIx("SORT");
    static int uniqIx = e->KeywordIx("UNIQUE");
    static int verbIx = e->KeywordIx("VERBOSE");
    
    if( e->KeywordSet( helpIx ) ) {
        e->Help( ints2str_help, ints2str_help_sz );
    }
    
    SizeT nParam RDX_UNUSED = e->NParam( 1 );
    BaseGDL* uints = e->GetPar(0);

    SizeT nEl = uints->N_Elements();
    if ( nEl < 1 ) {
        return new DStringGDL("");
    }
    
    bool sortSet = e->KeywordSet( sortIx );
    bool uniqSet = e->KeywordPresent( uniqIx );
    bool verbSet = e->KeywordPresent( verbIx );

    vector<uint64_t> tmp;
    switch( uints->Type() ) {
        case GDL_BYTE: {
            auto beg = reinterpret_cast<uint8_t*>(uints->DataAddr());
            std::copy( beg, beg + nEl , back_inserter(tmp) );
            break;
        }
        case GDL_INT: {
            auto beg = reinterpret_cast<int16_t*>(uints->DataAddr());
            std::copy( beg, beg + nEl , back_inserter(tmp) );
            break;
        }
        case GDL_UINT: {
            auto beg = reinterpret_cast<uint16_t*>(uints->DataAddr());
            std::copy( beg, beg + nEl , back_inserter(tmp) );
            break;
        }
        case GDL_LONG: {
            auto beg = reinterpret_cast<int32_t*>(uints->DataAddr());
            std::copy( beg, beg + nEl , back_inserter(tmp) );
            break;
        }
        case GDL_ULONG: {
            auto beg = reinterpret_cast<uint32_t*>(uints->DataAddr());
            std::copy( beg, beg + nEl , back_inserter(tmp) );
            break;
        }
        case GDL_LONG64: {
            auto beg = reinterpret_cast<int64_t*>(uints->DataAddr());
            std::copy( beg, beg + nEl , back_inserter(tmp) );
            break;
        }
        case GDL_ULONG64: {
            auto beg = reinterpret_cast<uint64_t*>(uints->DataAddr());
            std::copy( beg, beg + nEl , back_inserter(tmp) );
            break;
        }
        default: e->Throw("rdx_ints2str: input array must be of type BYTE, INT, LONG or LONG64.");
    }

    
    if( sortSet || uniqSet ) {
        std::sort( tmp.begin(), tmp.end() );
    }
    if( uniqSet ) {
        auto it = std::unique( tmp.begin(), tmp.end() );
        tmp.resize( std::distance( tmp.begin(), it ) );
    }

    string ret;
    try {
        ret = redux::util::uIntsToString( tmp );
    } catch( const exception & ex ) {
        e->Throw(string("rdx_ints2str: Error during call to redux::util::uIntsToString() : ") + ex.what());
    }
    
    return new DStringGDL(ret);

}


namespace {
    const string str2ints_help[]={
        "RDX_STR2INTS:       integer_array = rdx_str2ints(str, /KEYWORDS)",
        "   Accepted Keywords:",
        "      HELP                Display this info.",
        "      SORT                Sort output array.",
        "      UNIQUE              Discard repeated integers (will also force sort)",
        "      VERBOSE             Verbosity, default is 0 (only error output)."
    };
    constexpr int str2ints_help_sz = sizeof(str2ints_help) / sizeof(string);
}

BaseGDL* redux::rdx_str2ints_gdl( EnvT* e ) {
 
    static int helpIx = e->KeywordIx("HELP");
    static int sortIx = e->KeywordIx("SORT");
    static int uniqIx = e->KeywordIx("UNIQUE");
    static int verbIx = e->KeywordIx("VERBOSE");
    
    if( e->KeywordSet( helpIx ) ) {
        e->Help( ints2str_help, ints2str_help_sz );
    }
    
    SizeT nParam RDX_UNUSED = e->NParam( 1 );
    DString intList;
    e->AssureStringScalarPar( 0, intList );

    bool sortSet = e->KeywordSet( sortIx );
    bool uniqSet = e->KeywordPresent( uniqIx );
    bool verbSet = e->KeywordPresent( verbIx );

    vector<uint32_t> uints;
    try {
        uints = redux::util::stringToUInts<uint32_t>( intList );
    } catch( const exception& ex ) {
        e->Throw(string("rdx_str2ints: Error during call to redux::util::stringToUInts() : ") + ex.what());
    }

    if( sortSet || uniqSet ) {
        std::sort( uints.begin(), uints.end() );
    }
    if( uniqSet ) {
        auto it = std::unique( uints.begin(), uints.end() );
        uints.resize( std::distance( uints.begin(), it ) );
    }

    if ( uints.empty() ) {
        return new DUIntGDL(0);
    }
    
    DUIntGDL* ret = new DUIntGDL(dimension(uints.size()));
    for( SizeT i(0); i<uints.size(); ++i ) {
        (*ret)[i] = uints[i];
    }

    return ret;
    
}


void redux::rdx_gdl( EnvT* e ) {
    
    SizeT nParam RDX_UNUSED = e->NParam( 0 );

    static int vIx = e->KeywordIx("VERSION");
    bool get_version = e->KeywordPresent( vIx );

    if( get_version ) {
        DStringGDL* rdxVersion = new DStringGDL( 1 );
        e->SetKW( vIx, rdxVersion );
        *rdxVersion = getLongVersionString(false);
        return;
    }

}


extern "C" {

    BaseGDL* rdx_cacheget( EnvT* e ) {  return redux::rdx_cacheget_gdl(e); }
    BaseGDL* rdx_filetype( EnvT* e ) {  return redux::rdx_filetype_gdl(e); }
    BaseGDL* rdx_fillpix( EnvT* e ) {  return redux::rdx_fillpix_gdl(e); }
    BaseGDL* rdx_hasopencv( EnvT* e ) {  return redux::rdx_hasopencv_gdl(e); }
    BaseGDL* rdx_readdata( EnvT* e ) {  return redux::rdx_readdata_gdl(e); }
    BaseGDL* rdx_readhead( EnvT* e ) {  return redux::rdx_readhead_gdl(e); }
    BaseGDL* rdx_segment( EnvT* e ) {  return redux::rdx_segment_gdl(e); }
    BaseGDL* rdx_ints2str( EnvT* e ) {  return redux::rdx_ints2str_gdl(e); }
    BaseGDL* rdx_str2ints( EnvT* e ) {  return redux::rdx_str2ints_gdl(e); }

    void rdx( EnvT* e ) {  return redux::rdx_gdl(e); }
    void rdx_cache( EnvT* e ) {  return redux::rdx_cache_gdl(e); }
    void rdx_cacheclear( EnvT* e ) {  return redux::rdx_cacheclear_gdl(e); }
    void rdx_cachedel( EnvT* e ) {  return redux::rdx_cachedel_gdl(e); }
    void rdx_cacheinfo( EnvT* e ) {  return redux::rdx_cacheinfo_gdl(e); }
    
}
