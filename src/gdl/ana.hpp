#ifndef REDUX_GDL_ANA_HPP
#define REDUX_GDL_ANA_HPP

#include <envt.hpp>

namespace redux {

    // GDL functions
    BaseGDL* rdx_fzhead_gdl( EnvT* );
    BaseGDL* rdx_f0_gdl( EnvT* );
    
    // GDL Procedures
    void rdx_fzread_gdl( EnvT* );
    void rdx_fzwrite_gdl( EnvT* );
    void rdx_fcwrite_gdl( EnvT* );

}


#endif  // REDUX_GDL_ANA_HPP
