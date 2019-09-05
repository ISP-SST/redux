; GDL MODULE RDX
; DESCRIPTION REDUX tools
; VERSION ${LIBREDUX_VERSION_MAJOR}.${LIBREDUX_VERSION_MINOR}.${LIBREDUX_VERSION_PATCH}
; SOURCE Redux Pipeline - ${LIBREDUX_COMMIT_INFO}
; BUILD_DATE ${LIBREDUX_COMMIT_TIME}
PRO unload_rdx_gdl, verbose=verbose

    lib_name = '${GDL_DLM_DIR}/librdx_gdl.so'
    IF ( keyword_set(verbose) ) THEN print, 'Unlinking DLM: ' + lib_name
    unlinkimage, lib_name, /FORCE
    
END


PRO load_rdx_gdl, verbose=verbose
    
    lib_name = '${GDL_DLM_DIR}/librdx_gdl.so'
    IF ( keyword_set(verbose) ) THEN print, 'Linking DLM: ' + lib_name
    linkimage, 'rdx', lib_name, MIN_ARGS=0, MAX_ARGS=1, KEYWORDS=['VERSION']
    linkimage, 'rdx_cache', lib_name, MIN_ARGS=2, MAX_ARGS=2
    linkimage, 'rdx_cacheclear', lib_name, MIN_ARGS=0, MAX_ARGS=0
    linkimage, 'rdx_cachedel', lib_name, MIN_ARGS=1, MAX_ARGS=1, KEYWORDS=['COUNT']
    linkimage, 'rdx_cacheinfo', lib_name, MIN_ARGS=0, MAX_ARGS=0
    linkimage, 'rdx_cacheget', lib_name, /FUNCT, MIN_ARGS=1, MAX_ARGS=1, KEYWORDS=['COUNT']
    linkimage, 'rdx_f0', lib_name, /FUNCT, MIN_ARGS=1, MAX_ARGS=1, KEYWORDS=['FITS','HEADER','HELP']
    linkimage, 'rdx_fzhead', lib_name, /FUNCT, MIN_ARGS=1, MAX_ARGS=1, KEYWORDS=['FITS','HELP']
    linkimage, 'rdx_fzread', lib_name, MIN_ARGS=2, MAX_ARGS=3, KEYWORDS=['FITS','HELP']
    linkimage, 'rdx_fcwrite', lib_name, MIN_ARGS=2, MAX_ARGS=3, KEYWORDS=['COMPRESS','HELP']
    linkimage, 'rdx_fzwrite', lib_name, MIN_ARGS=2, MAX_ARGS=3, KEYWORDS=['COMPRESS','HELP']
    linkimage, 'rdx_hasopencv', lib_name, /FUNCT, MIN_ARGS=0, MAX_ARGS=0
    linkimage, 'rdx_filetype', lib_name, /FUNCT, MIN_ARGS=1, MAX_ARGS=1
    linkimage, 'rdx_fillpix', lib_name, /FUNCT, MIN_ARGS=1, MAX_ARGS=1, KEYWORDS=['HELP','MASK','NTHREADS','THRESHOLD','VERBOSE']
    linkimage, 'rdx_readdata', lib_name, /FUNCT, MIN_ARGS=0, MAX_ARGS=1, KEYWORDS=['ALL','DATE_BEG','FRAMENUMBERS','HEADER','HELP','RAW','STATUS']
    linkimage, 'rdx_readhead', lib_name, /FUNCT, MIN_ARGS=0, MAX_ARGS=1, KEYWORDS=['ALL','DATE_BEG','FRAMENUMBERS','HELP','RAW','STATUS']
    linkimage, 'rdx_segment', lib_name, /FUNCT, MIN_ARGS=3, MAX_ARGS=4, KEYWORDS=['HELP','MOMFBD']
    linkimage, 'rdx_ints2str', lib_name, /FUNCT, MIN_ARGS=1, MAX_ARGS=1, KEYWORDS=['HELP','SORT','UNIQUE','VERBOSE']
    linkimage, 'rdx_str2ints', lib_name, /FUNCT, MIN_ARGS=1, MAX_ARGS=1, KEYWORDS=['HELP','SORT','UNIQUE','VERBOSE']
    
    ; legacy names, will call the same functions as the rdx_* versions.
    linkimage, 'f0', lib_name, 1, 'rdx_f0', MIN_ARGS=1, MAX_ARGS=1, KEYWORDS=['FITS','HEADER','HELP']
    linkimage, 'fzhead', lib_name, 1, 'rdx_fzhead', MIN_ARGS=1, MAX_ARGS=1, KEYWORDS=['FITS','HELP']
    linkimage, 'fzread', lib_name, 0, 'rdx_fzread', MIN_ARGS=2, MAX_ARGS=3, KEYWORDS=['FITS','HELP']
    linkimage, 'fcwrite', lib_name, 0, 'rdx_fcwrite', MIN_ARGS=2, MAX_ARGS=3, KEYWORDS=['COMPRESS','HELP']
    linkimage, 'fzwrite', lib_name, 0, 'rdx_fzwrite', MIN_ARGS=2, MAX_ARGS=3, KEYWORDS=['COMPRESS','HELP']
    
END



