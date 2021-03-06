; uses histograms and delta_z <= dz_coeff*0.005*(1+z)

PRO latefrac_byField_stats_FD, outputFormat
	input_dm = 1.0
	string_dm = strtrim(strcompress(string(input_dm, format='(f20.1)')),1)

;	IF IPtype EQ 1 THEN BEGIN
;		IPsuffix = 'SF'
;		IPtypeLabel = 'Late-type'
;	ENDIF ELSE BEGIN
;		IPsuffix = 'Q'
;		IPtypeLabel = 'Early-type'
;	ENDELSE

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/latefrac_byField_norm_diff_' + string_dm + 'dex_FD', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	!P.CHARSIZE=1.1

	dz_coeff = 2.
	zmin = 0.2
	zmax = 1.0

	n_annuli 	= 60
	dRproj		= 0.25
	printEvery	= 100

	IPdataPath 	= '~/field-dependent_SFR_vs_mass/'
	zerodInputFile 	= 'zerodSFQ_IP_dz' + strtrim(strcompress(string(dz_coeff, format='(f20.1)')),1) + '_' + string_dm + 'dex.fits'
	zerodInput 	= mrdfits(IPdataPath + zerodInputFile, 1)
	fields 		= zerodInput[uniq(zerodInput.field, sort(zerodInput.field))].field

	Rproj_array = float(dRproj)*findgen(n_annuli+1)

	xmin = 0
        xmax = n_annuli*dRproj
        xr = xmax-xmin

        ymin = -0.1
        ymax = 0.2
        yr = ymax-ymin

	!P.MULTI = 0
	!P.CHARSIZE = 1.25

	colors = ['Green', 'Cyan', 'Blue', 'Magenta', 'Red']

	PLOT, findgen(10), findgen(10), xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Normalized % diff. in late-type frac between LT/ET IPs', /NODATA, $
		title = 'Norm. conformity signal strengh by field ' + textoidl('(M13+1.0dex; \Deltaz=2.0)')

	LEGEND, [strcompress(fields)], LINESTYLE=[0,0,0,0,0], color=cgColor(colors), BOX=0, /TOP, /RIGHT
;	LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], BOX=0, /BOTTOM, /LEFT
;	XYOUTS, xmin+0.05*xr, ymin+0.05*yr, '0.2 < z < 1.0', ALIGNMENT=0.0
;	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, IPtypeLabel + ' IPs', ALIGNMENT=1.0

	masterArray = ['Mpc', string(Rproj_array[1:*])]
;	PRINT, 'Length of masterArray: ', n_elements(masterArray)

	FOR i=0,n_elements(fields)-1 DO BEGIN
		PRINT, 'Field: ' + strtrim(fields[i],1)

		get_latefrac, zerodInput, string_dm, dz_coeff, zmin, zmax, n_annuli, dRproj, printEvery, n_tot, n_late, fields[i], 1
		fracIPSF = n_late/n_tot

		get_latefrac, zerodInput, string_dm, dz_coeff, zmin, zmax, n_annuli, dRproj, printEvery, n_tot, n_late, fields[i], 0
		fracIPQ = n_late/n_tot

;		frac	= n_late/n_tot
;		dfrac	= poissonError(n_late, n_tot)

		fieldArray = [strtrim(fields[i],1), (fracIPSF - fracIPQ)/((fracIPSF + fracIPQ)/2.)]

		OPLOT, Rproj_array[1:*], fieldArray[1:*], LINESTYLE=0, color=cgColor(colors[i])
;		ERRPLOT, Rproj_array+0.5*dRproj-0.025*i, frac-dfrac, frac+dfrac, color=cgColor(colors[i])

;		PRINT, 'Length of fieldArray: ', n_elements(fieldArray)

		masterArray = [[masterArray], [fieldArray]]
	ENDFOR

	PRINT, TRANSPOSE(masterArray)

	SAVE, masterArray, FILE='~/field-dependent_SFR_vs_mass/latefrac_byField_data_FD.sav'
	
	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

PRO get_latefrac, zerodInput, string_dm, dz_coeff, zmin, zmax, n_annuli, dRproj, printEvery, n_tot, n_late, field, IPtype
	data = zerodInput[where((zerodInput.zprimus GE zmin) AND (zerodInput.zprimus LE zmax) AND (zerodInput.field EQ field))]

;	PRINT, 'total galaxies: ', n_elements(zerodInput), $
;	  ' zmin:', zmin, $
;	  ' zmax:', zmax, $
;	  ' galaxies within z range:', n_elements(data)

	Rproj_array = float(dRproj)*findgen(n_annuli+1)

	data_isIP = data[where(data.IP eq 1)]

	n_tot 	= 0.0*findgen(n_annuli)
	n_late 	= 0.0*findgen(n_annuli)

	outputIndex	= 0L

	IP 	= data_isIP[where(data_isIP.SFQ eq IPtype)]

	subtotal_neighbors	= 0L*indgen(n_annuli)
	FOR i=0,(n_elements(IP)-1) DO BEGIN 
;	FOR i=0,99 DO BEGIN 
		currentIP 	= IP[i]

		; keep galaxies within dz=dz+coeff*0.005*(1+z) of IP candidate (and in same field) 
		dz_slice_all 	= data[ where( (data.field eq currentIP.field) AND (data.objname ne currentIP.objname) AND $
			(ABS(data.zprimus - currentIP.zprimus) LE dz_coeff*0.005*(1.+currentIP.zprimus)) ) ]
		dz_slice_all_dists = SQRT( (currentIP.xprop - dz_slice_all.xprop)^2 + (currentIP.yprop - dz_slice_all.yprop)^2 )
	
		FOR n=0,n_annuli-1 DO BEGIN
			currentMinR = dRproj*float(n)
			currentMaxR = dRproj*float(n+1)
			annulus_dists = dz_slice_all_dists[where( (dz_slice_all_dists GE currentMinR) AND (dz_slice_all_dists LE currentMaxR), /NULL )]
;			subtotal_annuli_dists = dz_slice_all_dists[where( (dz_slice_all_dists LE currentMaxR), /NULL )]
;			PRINT, currentMinR, ' ', currentMaxR, ' ', n_elements(annulus_dists), ' ', n_elements(subtotal_annuli_dists)
			subtotal_neighbors[n] += n_elements(dz_slice_all_dists[where( (dz_slice_all_dists LE currentMaxR), /NULL) ])
;			subtotal_neighbors[n] += n_elements(annulus_dists)
		ENDFOR
;		PRINT, strcompress(outputIndex), ') ', subtotal_neighbors

		; add stats for current IP to n_tot
		dR_array_tot	= histogram(dz_slice_all_dists, bin=dRproj, min=0., nbins=n_annuli)
		n_tot		= n_tot + dR_array_tot[0:n_annuli-1]
		
		dz_slice_late	= dz_slice_all[where(dz_slice_all.SFQ eq 1)]
		dz_slice_late_dists = ((currentIP.xprop - dz_slice_late.xprop)^2 + $
				   (currentIP.yprop - dz_slice_late.yprop)^2)^(0.5)
		dR_array_late	= histogram(dz_slice_late_dists, bin=dRproj, min=0., nbins=n_annuli)
		n_late		= n_late + dR_array_late[0:n_annuli-1]

		IF outputIndex mod printEvery eq 0 THEN PRINT, outputIndex
		outputIndex += 1
	ENDFOR
END
