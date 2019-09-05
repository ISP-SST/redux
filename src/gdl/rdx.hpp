#ifndef REDUX_GDL_RDX_HPP
#define REDUX_GDL_RDX_HPP

#include <envt.hpp>

namespace redux {

    // GDL functions
    BaseGDL* rdx_cacheget_gdl( EnvT* );
    BaseGDL* rdx_filetype_gdl( EnvT* );
    BaseGDL* rdx_fillpix_gdl( EnvT* );
    BaseGDL* rdx_hasopencv_gdl( EnvT* );
    BaseGDL* rdx_readdata_gdl( EnvT* );
    BaseGDL* rdx_readhead_gdl( EnvT* );
    BaseGDL* rdx_segment_gdl( EnvT* );
    BaseGDL* rdx_ints2str_gdl( EnvT* );
    BaseGDL* rdx_str2ints_gdl( EnvT* );
    
    // GDL Procedures
    void rdx_cache_gdl( EnvT* );
    void rdx_cacheclear_gdl( EnvT* );
    void rdx_cachedel_gdl( EnvT* );
    void rdx_cacheinfo_gdl( EnvT* );
    void rdx_gdl( EnvT* );

}

#endif  // REDUX_GDL_RDX_HPP
