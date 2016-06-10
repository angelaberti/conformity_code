; uses histograms and delta_z <= dz_coeff*0.005*(1+z)

PRO latefrac_byField_stats_faster_FD, outputFormat
	input_dm = 1.0
	string_dm = strtrim(strcompress(string(input_dm, format='(f20.1)')),1)

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

	RESTORE, '~/field-dependent_SFR_vs_mass/latefrac_byField_data_FD.sav'
	masterArray = TRANSPOSE(masterArray)

	FOR i=0,n_elements(fields)-1 DO BEGIN
		PRINT, 'Field: ' + strtrim(fields[i],1)
		OPLOT, Rproj_array[1:*], masterArray[i+1,1:*], LINESTYLE=0, color=cgColor(colors[i])
	ENDFOR

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
