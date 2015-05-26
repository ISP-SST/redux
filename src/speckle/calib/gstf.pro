pro plotbeta, betas, file, tit=tit
    if keyword_set(title) then title=tit else title = (strsplit(file, '_', /extract))[2]
    set_plot, 'PS' 
    device, /landscape, filename='betas_'+title+'.ps'
    plot, betas, yrange=[-1,1], xrange=[0,155], /xstyle, title = 'betas at '+title, psym=10
    plot, betas, yrange=[0,1], xrange=[0,30], /xstyle, title = 'betas at '+title, psym=10
    device, /close
    set_plot, 'X'
    plot, betas, yrange=[-1,1], xrange=[0,155], /xstyle,title = 'betas at '+title, psym=10
end

pro calcbeta, path, resbetas
    start:

    if n_elements(path) eq 1 then begin
    totfiles = dialog_pickfile(Filter=['*totalerror*.fits'], path=path, title='Select total wferror dataset(s)', /multiple)
    if totfiles[0] eq '' then begin 
        print, "No files choosed!"
        goto, start
    endif

    resfiles = dialog_pickfile(Filter=['*reserror*.fits'], path=path, title='Select residual wferror dataset(s)', /multiple)
    if resfiles[0] eq '' then begin 
        print, "No files choosed!"
        goto, start
    endif

    nrfiles = n_elements(totfiles)
    if not nrfiles eq n_elements(resfiles) then begin
        print, "Choose as many total as residual wferror files"
        goto, start
    endif

    endif else begin
        nrfiles = 1
        totfiles = path[0]
        resfiles = path[1]
    endelse
    

    for i=0, nrfiles - 1  do begin
        readwfs, totfiles[i], tot, htot
        readwfs, resfiles[i], res, hres 
;        sres = stddev(res, dim=2)
;        stot = stddev(tot, dim=2)
        vres = variance(res, dim=2)
        vtot = variance(tot, dim=2)
        zeros = where(vtot eq 0) ; if there is no correction set beta to 0
        vres[zeros] = 1.         ; so residual error 1
        vtot[zeros] = 1.         ; and total error 1
;        betas = 1 - sqrt(sres*sres/(stot*stot))
        betas = 1 - sqrt(vres/vtot)
        if i eq 0 then mbetas = betas[*,0]*0.
        mbetas += betas
    endfor

    mbetas /= nrfiles

    resbetas = mbetas

    plot, mbetas

    print, "total error files:"
    print, totfiles
    print, "residual error files:"
    print, resfiles
end

pro readwfs, file, data, header
    tmp = readfits(file, h)
    tmp = 0.

    nraxis = fix((strsplit(h[2], ' ', /extract))[2])
    bitpix = float((strsplit(h[1], ' ', /extract))[2])

    if nraxis eq 1 then begin 
        naxis1 = fix((strsplit(h[3], ' ', /extract))[2]) 
    endif else begin 
        if fix((strsplit(h[3], ' ', /extract))[2]) eq 1 then begin 
            nraxis = 1
            naxis1 = fix((strsplit(h[4], ' ', /extract))[2]) 
        endif else begin
            naxis1 = fix((strsplit(h[3], ' ', /extract))[2])
            naxis2 = fix((strsplit(h[4], ' ', /extract))[2])
        endelse
    endelse

    case bitpix of
        8: if nraxis eq 1 then data = bytarr(naxis1) else data = bytarr(naxis1, naxis2)
        16: if nraxis eq 1 then data = intarr(naxis1) else data = intarr(naxis1, naxis2)
        32: if nraxis eq 1 then data = lonarr(naxis1) else data = lonarr(naxis1, naxis2)
        64: if nraxis eq 1 then data = lon64arr(naxis1) else data = lon64arr(naxis1, naxis2)
        -64: if nraxis eq 1 then data = dblarr(naxis1) else data = dblarr(naxis1, naxis2)
        else: print, 'Unknown datatype'
        endcase

    openr, lun, file, /get_lun, /swap_if_big_endian
    skip_lun, lun, 2880
    readu, lun, data
    close, lun
    free_lun, lun

    header = h

end


pro readresult, file, ltf, stf, sr, freqs, alphas, eps=eps
    nlines = 0.
    nfreqs = 1
    nalphas = 1
    ncorrlev = 0
    line = ''
    zernikes = 0
    cl = 0

    openr, lun, file, /get_lun
    while not EOF(lun) do begin 
        readf, lun, line
        if nlines eq 0 then zernikes = fix((strsplit(line, ' ', /extract))[4])
        if line eq '# spatial frequency' then nfreqs = nlines
        if line eq '# seeing alpha' then begin
            nfreqs = nlines - nfreqs - 1
            nalphas = nlines
        endif
        if line eq '# correction levels' or line eq '# angles'then begin
            if line eq '# correction levels' then cl = 1
            nalphas = nlines - nalphas - 1
            ncorrlev = nlines
        endif
        if line eq '# functions' then begin
            if ncorrlev eq 0 then begin 
                nalphas = nlines - nalphas - 1
                ncorrlev = 1
            endif else begin
                ncorrlev = nlines - ncorrlev - 1 
            endelse
        endif
        nlines += 1
    endwhile

    close, lun
    free_lun, lun


    data = strarr(nlines)
    alphas = fltarr(nalphas)
    freqs = fltarr(nfreqs)
    angles = fltarr(ncorrlev)
    ltf = fltarr(nfreqs, 2, nalphas, ncorrlev)
    stf = fltarr(nfreqs, 2, nalphas, ncorrlev)
    sr = fltarr(nfreqs, 2, nalphas, ncorrlev)

    openr, lun, file, /get_lun
    for i=0, nlines-1 do begin
        readf, lun, line
        data[i] = line
    endfor
    close, lun
    free_lun, lun

    offset = 2
    freqs[*] = float(data[offset : offset + nfreqs - 1])
    offset += nfreqs + 1
    alphas[*] = float(data[offset : offset + nalphas - 1 ])
    offset += nalphas + 1

    if ncorrlev gt 1 then begin
        for i=0, ncorrlev-1 do begin
            angles[i] = float(data[offset + i])
        endfor
        offset += ncorrlev + 2
        for k = 0, ncorrlev - 1 do begin
            for i=0, nalphas - 1 do begin
                for j=0, nfreqs - 1 do begin
                    tmp = strsplit(data[offset + i*nfreqs  + j + k*nfreqs*nalphas], ' ', /extract)
                    ltf[j,0, i, k] = float(tmp[0])
                    ltf[j,1, i, k] = float(tmp[1])
                    stf[j,0, i, k] = float(tmp[2])
                    stf[j,1, i, k] = float(tmp[3])
                    sr[j,0, i, k] = float(tmp[4])
                    sr[j,1, i, k] = float(tmp[5])
                endfor
                offset += 1
            endfor
            if cl then offset += 1
        endfor
   endif else begin
        offset += 1
        for i=0, nalphas - 1 do begin
            for j=0, nfreqs - 1 do begin
                tmp = strsplit(data[offset + i*nfreqs  + j ], ' ', /extract)
                ltf[j,0, i] = float(tmp[0])
                ltf[j,1, i] = float(tmp[1])
                stf[j,0, i] = float(tmp[2])
                stf[j,1, i] = float(tmp[3])
                sr[j,0, i] = float(tmp[4])
                sr[j,1, i] = float(tmp[5])
            endfor
            offset += 1
        endfor
    endelse

    window, 0, title="Speckle transfer function"
    plot, freqs, stf[*,0,0], title="Speckle transfer function"
    for i=0, nalphas-1 do oplot, freqs, stf[*,0,i]
    window, 1, title="Long exposure transfer function"
    plot, freqs, ltf[*,0,0], title="Long exposure transfer function"
    for i=0, nalphas-1 do oplot, freqs, ltf[*,0,i]
    window, 2, title="Spectral ratio"
    plot, freqs, sr[*,0,0], title="Spectral ratio"
    for i=0, nalphas-1 do oplot, freqs, sr[*,0,i]

    if keyword_set(eps) then begin 
       set_plot, 'ps'  
       device, filename="stf.ps", /encaps
       plot, freqs, stf[*,0,0], title="Speckle transfer function"
       for i=0, nalphas-1 do oplot, freqs, stf[*,0,i]
       device, /close

       device, filename="ltf.ps", /encaps
       plot, freqs, ltf[*,0,0], title="Long exposure transfer function"
       for i=0, nalphas-1 do oplot, freqs, ltf[*,0,i]
       device, /close

       device, filename="sr.ps", /encaps
       plot, freqs, sr[*,0,0], title="Spectral ratio"
       for i=0, nalphas-1 do oplot, freqs, sr[*,0,i]
       device, /close
    set_plot, 'X'
    endif
end

