#include "idlutil.hpp"

#include <redux/util/datautil.hpp>
#include <redux/util/stringutil.hpp>


#include <algorithm>

#include <boost/algorithm/string.hpp>

using namespace redux::util;
using namespace redux;

using namespace std;


typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_INT help;
    IDL_INT sort;
    IDL_INT unique;
    IDL_INT verbose;
} STR_TO_INT_KW;


// NOTE:  The keywords MUST be listed in alphabetical order !!
static IDL_KW_PAR str_to_int_kw_pars[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "HELP",      IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*)IDL_KW_OFFSETOF2(STR_TO_INT_KW,help) },
    { (char*) "SORT",      IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*)IDL_KW_OFFSETOF2(STR_TO_INT_KW,sort) },
    { (char*) "UNIQUE",    IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*)IDL_KW_OFFSETOF2(STR_TO_INT_KW,unique) },
    { (char*) "VERBOSE",   IDL_TYP_INT, 1, IDL_KW_ZERO, 0, (char*)IDL_KW_OFFSETOF2(STR_TO_INT_KW,verbose) },
    { NULL }
};


string str2ints_info( int lvl ) {
    
    string ret = "RDX_STR2INTS";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"       ");          // newline if lvl>1
        ret += "   Syntax:   arr = rdx_str2ints(str,/KEYWORDS)\n";
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      SORT                Sort output array.\n"
                    "      UNIQUE              Discard repeated integers (will also force sort)\n"
                    "      VERBOSE             Verbosity, default is 0 (only error output).\n";
        }
    } else ret += "\n";

    return ret;
    
}


IDL_VPTR str2ints(int argc, IDL_VPTR* argv, char* argk) {
 
    if (argc < 1) {
        return IDL_GettmpInt(0);
    }
    
    IDL_VPTR ret;
    vector<uint32_t> uints;
    IDL_VPTR str = argv[0];
    IDL_ENSURE_SIMPLE(str);
    IDL_ENSURE_STRING(str);
    
    STR_TO_INT_KW kw;
    kw.sort = 0;
    kw.verbose = 0;
    kw.unique = 0;
    (void) IDL_KWProcessByOffset (argc, argv, argk, str_to_int_kw_pars, (IDL_VPTR*) 0, 255, &kw);
    if( kw.help ) {
        cout << str2ints_info(2) << endl;
        return IDL_GettmpInt(0);
    }
    
 
    try {
        uints = redux::util::stringToUInts<uint32_t>( IDL_VarGetString(str) );
    } catch( const exception& e ) {
        cout << "rdx_str2ints: Error during call to redux::util::stringToUInts() : " << e.what() << endl;
    }

    if( kw.sort || kw.unique ) {
        std::sort( uints.begin(), uints.end() );
    }
    if( kw.unique ) {
        auto it = std::unique( uints.begin(), uints.end() );
        uints.resize( std::distance( uints.begin(), it ) );
    }

    if ( uints.empty() ) {
        return IDL_GettmpInt(0);
    }
    
    IDL_MEMINT dims[] = { (int)uints.size() };
    IDL_LONG* tmpData = (IDL_LONG*)IDL_MakeTempArray(IDL_TYP_LONG,1,dims,IDL_ARR_INI_NOP,&ret);
    std::copy( uints.begin(), uints.end(), tmpData );

    return ret;

}


string ints2str_info( int lvl ) {
    
    string ret = "RDX_INTS2STR";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"       ");          // newline if lvl>1
        ret += "   Syntax:   str = rdx_ints2str(arr,/KEYWORDS)\n";
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      SORT                Sort input array.\n"
                    "      UNIQUE              Discard repeated integers (will also force sort)\n"
                    "      VERBOSE             Verbosity, default is 0 (only error output).\n";
        }
    } else ret += "\n";

    return ret;
    
}


IDL_VPTR ints2str( int argc, IDL_VPTR* argv, char* argk ) {
    
    STR_TO_INT_KW kw;
    kw.sort = 0;
    kw.verbose = 0;
    kw.unique = 0;
    int nPlainArgs = IDL_KWProcessByOffset (argc, argv, argk, str_to_int_kw_pars, (IDL_VPTR*) 0, 255, &kw);
    
    string ret;
    if (nPlainArgs < 1) {
        return IDL_StrToSTRING( (char*)"" );
    }
    
    IDL_VPTR uints = argv[0];
    IDL_ENSURE_SIMPLE(uints);
    IDL_ENSURE_ARRAY(uints);
    
    if( kw.help ) {
        cout << ints2str_info(2) << endl;
        return IDL_StrToSTRING( (char*)"" );
    }
    
    vector<uint64_t> tmp;
    if( uints->type == IDL_TYP_BYTE ) {
        auto beg = uints->value.arr->data;
        std::copy( beg, beg + uints->value.arr->n_elts , back_inserter(tmp) );
    } else if (uints->type == IDL_TYP_INT ) {
        auto beg = reinterpret_cast<int16_t*>(uints->value.arr->data);
        std::copy( beg, beg + uints->value.arr->n_elts , back_inserter(tmp) );
    } else if (uints->type == IDL_TYP_LONG ) {
        auto beg = reinterpret_cast<int32_t*>(uints->value.arr->data);
        std::copy( beg, beg + uints->value.arr->n_elts , back_inserter(tmp) );
    } else  {
        cout << "rdx_ints2str: input array must be of type BYTE, INT or LONG." << endl;
        return IDL_GettmpInt (0);
    }
    
    if( kw.sort || kw.unique ) {
        std::sort( tmp.begin(), tmp.end() );
    }
    if( kw.unique ) {
        auto it = std::unique( tmp.begin(), tmp.end() );
        tmp.resize( std::distance( tmp.begin(), it ) );
    }

    try {
        ret = redux::util::uIntsToString( tmp );
    } catch( const exception & e ) {
        cout << "rdx_ints2str: Error during call to redux::util::uIntsToString() : " << e.what() << endl;
    }
    
    return IDL_StrToSTRING( (char*)ret.c_str() );

}


typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_VPTR arglist;
    IDL_INT help;
    IDL_STRING split_chars;
} MT_KW;


// NOTE:  The keywords MUST be listed in alphabetical order !!
static IDL_KW_PAR mt_kw_pars[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "ARG_LIST",    IDL_TYP_UNDEF,  1, IDL_KW_OUT|IDL_KW_ZERO, 0, (char*) IDL_KW_OFFSETOF2(MT_KW,arglist) },
    { (char*) "HELP",        IDL_TYP_INT,    1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(MT_KW,help) },
    { (char*) "SPLIT_CHARS", IDL_TYP_STRING, 1, IDL_KW_ZERO,            0, (char*)IDL_KW_OFFSETOF2(MT_KW,split_chars) },
    { NULL }
};


string maketemplate_info( int lvl ) {
    
    string ret = "RDX_MAKE_TEMPLATE";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"  ");          // newline if lvl>1
        ret += "   Syntax:   tpl = rdx_make_template(list,/KEYWORDS)\n";
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      ARG_LIST            (output) List containing the template arguments.\n"
                    "      SPLIT_CHARS         Which character(s) to split the list-items on. (default: '.')\n"
                    "   The first split-char will also be used when generating the template.\n";
        }
    } else ret += "\n";

    return ret;
    
}


IDL_VPTR maketemplate(int argc, IDL_VPTR* argv, char* argk) {
    
    MT_KW kw;
    int nPlainArgs = IDL_KWProcessByOffset (argc, argv, argk, mt_kw_pars, (IDL_VPTR*)0, 255, &kw);
    
    if (nPlainArgs < 1) {
        return IDL_StrToSTRING( (char*)"" );
    }
    
    IDL_VPTR strlist = argv[0];
    IDL_ENSURE_STRING( strlist )
    IDL_ENSURE_SIMPLE( strlist );
    
    if ( !(strlist->flags & IDL_V_ARR) ) {
        return strlist;
    }

    IDL_ARRAY* strarr = strlist->value.arr;
    IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>(strlist->value.arr->data);

    if ( strarr->n_elts == 1) {
        return IDL_StrToSTRING( strptr->s );
    }
    
    if ( kw.help ) {
        cout << maketemplate_info(2) << endl;
        return IDL_StrToSTRING( (char*)"" );
    }
    

    string split_chars = ".";
    if( kw.split_chars.slen ) split_chars = kw.split_chars.s;

    std::vector<std::string> stringList;
    for( int i=0; i<strarr->n_elts; ++i ) {
        stringList.push_back(strptr[i].s);
    }
    
    std::string tpl;
    vector< set<string> > arg_list;
    try {
        arg_list = redux::util::make_template( stringList, tpl, split_chars );
    } catch( const exception & e ) {
        cout << "rdx_make_template: Error during call to redux::util::make_template : " << e.what() << endl;
    }
        
    IDL_MEMINT nArgs = arg_list.size();
    if( kw.arglist ) {
        if( nArgs ) {
            IDL_VPTR arglist;
            IDL_MEMINT dims[] = { nArgs };    // x/y, im#, nMatches 
            IDL_MakeTempArray( IDL_TYP_STRING, 1, dims, IDL_ARR_INI_ZERO, &arglist );
            strptr = reinterpret_cast<IDL_STRING*>(arglist->value.arr->data);
            for( int j=0; j<nArgs; ++j ) {
                if( arg_list[j].size() > 1 ) {
                    string tmp = boost::join( arg_list[j], " " );
                    IDL_StrStore(strptr++,const_cast<char*>(tmp.c_str()));
                }
            }
            IDL_VarCopy( arglist, kw.arglist );
        } else {
            IDL_VarCopy( IDL_StrToSTRING( (char*)"" ), kw.arglist );
        }
    }
    
    return IDL_StrToSTRING( (char*)tpl.c_str() );

}



typedef struct {
    IDL_KW_RESULT_FIRST_FIELD; /* Must be first entry in structure */
    IDL_INT help;
    IDL_INT max_n;
} STR_REPLACE_KW;


// NOTE:  The keywords MUST be listed in alphabetical order !!
static IDL_KW_PAR str_replace_pars[] = {
    IDL_KW_FAST_SCAN,
    { (char*) "HELP",        IDL_TYP_INT,    1, IDL_KW_ZERO,            0, (char*) IDL_KW_OFFSETOF2(STR_REPLACE_KW, help) },
    { (char*) "N",           IDL_TYP_INT,    1,           0,            0, (char*) IDL_KW_OFFSETOF2(STR_REPLACE_KW, max_n) },
    { NULL }
};


string strreplace_info( int lvl ) {
    
    string ret = "RDX_STRREPLACE";
    if( lvl > 0 ) {
        ret += ((lvl > 1)?"\n":"     ");          // newline if lvl>1
        ret += "   Syntax:   result = rdx_strreplace(input,findstring,replacestring,/KEYWORDS)\n";
        if( lvl > 1 ) {
            ret +=  "   Accepted Keywords:\n"
                    "      HELP                Display this info.\n"
                    "      N                   How many occurances to replace (default=1)\n";
        }
    } else ret += "\n";

    return ret;
    
}


IDL_VPTR str_replace(int argc, IDL_VPTR* argv, char* argk) {
    
    STR_REPLACE_KW kw;
    kw.max_n = 1;
    int nPlainArgs = IDL_KWProcessByOffset (argc, argv, argk, str_replace_pars, (IDL_VPTR*)0, 255, &kw);
    
    if( nPlainArgs < 3 ) {
        return IDL_StrToSTRING( (char*)"" );
    }

    IDL_VPTR strlist = argv[0];
    IDL_ENSURE_STRING( strlist );
    IDL_ENSURE_SIMPLE( strlist );
    IDL_ENSURE_STRING( argv[1] );
    IDL_ENSURE_SIMPLE( argv[1] );
    IDL_ENSURE_STRING( argv[2] );
    IDL_ENSURE_SIMPLE( argv[2] );
    
    IDL_STRING findStr = argv[1]->value.str;
    IDL_STRING replaceStr = argv[2]->value.str;
    
    if ( !(strlist->flags & IDL_V_ARR) ) {
        IDL_STRING strptr = strlist->value.str;
        string tmp = replace_n( strptr.s, findStr.s, replaceStr.s, kw.max_n );
        return IDL_StrToSTRING( (char*)tmp.c_str() );
    }
    
    IDL_VPTR results;
    IDL_ARRAY* strarr = strlist->value.arr;
    IDL_STRING* strptr = reinterpret_cast<IDL_STRING*>(strarr->data);
    IDL_MEMINT dims[] = { strarr->n_elts };    // x/y, im#, nMatches 
    IDL_MakeTempArray( IDL_TYP_STRING, 1, dims, IDL_ARR_INI_ZERO, &results );
    IDL_STRING* resptr = reinterpret_cast<IDL_STRING*>(results->value.arr->data);
    for( int i=0; i<strarr->n_elts; ++i ) {
        string tmp = replace_n( strptr[i].s, findStr.s, replaceStr.s, kw.max_n );
        IDL_StrStore( &(resptr[i]), (char*)tmp.c_str());
    }
    
    return results;


}


namespace {
    static int dummy RDX_UNUSED =
    IdlContainer::registerRoutine( {{(IDL_SYSRTN_GENERIC)str2ints}, (char*)"RDX_STR2INTS", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, str2ints_info ) +
    IdlContainer::registerRoutine( {{(IDL_SYSRTN_GENERIC)ints2str}, (char*)"RDX_INTS2STR", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, ints2str_info ) +
    IdlContainer::registerRoutine( {{(IDL_SYSRTN_GENERIC)maketemplate}, (char*)"RDX_MAKE_TEMPLATE", 1, 1, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, maketemplate_info ) +
    IdlContainer::registerRoutine( {{(IDL_SYSRTN_GENERIC)str_replace}, (char*)"RDX_STRREPLACE", 3, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0 }, 1, strreplace_info );
}

