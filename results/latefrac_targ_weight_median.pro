; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

PRO latefrac_targ_weight_median, input_dm, outputFormat;, zmin, zmax;, dz_coeff, printEvery

        string_dm = strtrim(strcompress(string(input_dm, format='(f20.1)')),1)

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/latefracplot_targ_weight_median' + string_dm + 'dex_zall', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

        IPdataPath = '~/conformity/results/variable_mass+' + string_dm + 'dex/IP_data/'
        zerodInputFile = 'zerodSFQ_IP_dz' + strtrim(strcompress(string(dz_coeff, format='(f20.1)')),1) + '.fits'

;	PRINT, 'Input data: ', zerodInputFile

        zerodInput = mrdfits(IPdataPath + zerodInputFile, 1)

	Rmax 		= 15.
	dRproj		= 1.0
	n_annuli 	= Ceil(Rmax/dRproj)
	printEvery	= 100

	data = zerodInput[where((zerodInput.zprimus GE zmin) and (zerodInput.zprimus LE zmax))]
	data = data[where(data.targ_weight GE 1.)]

	!p.multi=0
	!p.charsize=1.25

	Rproj_array = float(dRproj)*(findgen(n_annuli+1)+0.5)

	IF (string(outputFormat) EQ 'ps') THEN (medColor  = cgColor('Black')) ELSE (medColor = cgColor('White'))
	meanColor = cgColor('Blue')
	colors = [meanColor, meanColor, medColor, medColor]

	xmin = 0
	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = 0.5
	ymax = 1
	yr = ymax-ymin

	get_targ_weight, data, dz_coeff, 1, IP_latefrac_mean_array, IP_latefrac_median_array, n_annuli, dRproj, printEvery

	PLOT, Rproj_array, IP_latefrac_mean_array, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', $
	  title = 'M13 Mass Limit + ' + string_dm + ' (0.0) dex for SF (Q) IP ' + textoidl('(\Deltaz=2.0)'), /NODATA
	OPLOT, Rproj_array, IP_latefrac_mean_array, LINESTYLE=0, color=meanColor
	LEGEND, ['LT IP mean', 'ET IP mean', 'LT IP median', 'ET IP median'], LINESTYLE=[0,2,0,2], color=colors, BOX=0, /BOTTOM, /LEFT
	OPLOT, Rproj_array, IP_latefrac_median_array, LINESTYLE=0, color=medColor

	get_targ_weight, data, dz_coeff, 0, IP_latefrac_mean_array, IP_latefrac_median_array, n_annuli, dRproj, printEvery

	OPLOT, Rproj_array, IP_latefrac_mean_array, LINESTYLE=2, color=meanColor
	OPLOT, Rproj_array, IP_latefrac_median_array, LINESTYLE=2, color=medColor

	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, '0.2 < z < 1.0', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

PRO get_targ_weight, data, dz_coeff, IPSFstatus, IP_latefrac_mean_array, IP_latefrac_median_array, n_annuli, dRproj, printEvery
	data_isIP = data[where(data.IP EQ 1)]

	IP_latefrac_array = []

        IP = data_isIP[where(data_isIP.SFQ EQ IPSFstatus)]
	IF (IPSFstatus EQ 1) THEN (IPtype = 'Star-forming') ELSE (IPtype = 'Quiescent')
	PRINT, IPtype + ' IP galaxies: ' + strtrim(n_elements(IP),1)

	outputIndex = 0L

	FOR i=0,(n_elements(IP)-1) DO BEGIN
;	FOR i=0,99 DO BEGIN
                currentIP = IP[i]

                ; keep galaxies within dz=dz+coeff*0.005*(1+z) of IP candidate (and in same field)
                all_dz_neigh = data[ where( (data.field EQ currentIP.field) AND (data.objname NE currentIP.objname) AND $
                        (ABS(data.zprimus - currentIP.zprimus) le dz_coeff*0.005*(1.+currentIP.zprimus)) ) ]
		all_dz_neigh_dists = [SQRT((all_dz_neigh.xprop-currentIP.xprop)^2 + (all_dz_neigh.yprop-currentIP.yprop)^2)]

		currentIP_latefrac_array = []

		; GET LATE FRACTION FOR EACH ANNULUS
                FOR n=0,n_annuli-1 DO BEGIN
			currentRmin = dRproj*float(n)
			currentRmax = dRproj*float(n+1)
			
;			PRINT, currentRmin, ' ', currentRmax

			annulus_neigh = all_dz_neigh[where( (all_dz_neigh_dists GE currentRmin) AND (all_dz_neigh_dists LE currentRmax), /NULL )]
;			IF annulus_neigh NE !NULL THEN PRINT, 'Annulus neighbors SF status: ', annulus_neigh.SFQ

			; DEFAULT
			annulus_latefrac = -99

			; NUMBER OF GALAXIES IN ANNULUS > 0
			IF (annulus_neigh NE !NULL) THEN BEGIN
				annulus_neigh_SF = annulus_neigh[where(annulus_neigh.SFQ EQ 1, /NULL)]
;				PRINT, 'SF neighbors in annulus: ', n_elements(annulus_neigh_SF)
;				PRINT, 'Q neighbors in annulus:  ', n_elements(annulus_neigh[where(annulus_neigh.SFQ EQ 0, /NULL)])
;				PRINT, 'Annulus late-type fraction: ', float(n_elements(annulus_neigh_SF))/n_elements(annulus_neigh)

;				PRINT, 'annulus_neigh_SF.targ_weight: ', annulus_neigh_SF.targ_weight
;				PRINT, 'annulus_neigh.targ_weight: ', annulus_neigh.targ_weight
				IF (annulus_neigh_SF EQ !NULL) THEN (annulus_latefrac = 0.) ELSE $
				  annulus_latefrac = TOTAL(annulus_neigh_SF.targ_weight)/TOTAL(annulus_neigh.targ_weight)				
;				PRINT, 'Annulus WEIGHTED late-type fraction: ', annulus_latefrac
;				PRINT, ''
			ENDIF

			; ADD ANNULUS LATE FRACTION TO ARRAY FOR SPECIFIC IP
			currentIP_latefrac_array = [currentIP_latefrac_array, annulus_latefrac]
		ENDFOR

		IP_latefrac_array = [[IP_latefrac_array], [currentIP_latefrac_array]]

                IF outputIndex mod printEvery EQ 0 THEN PRINT, outputIndex
		outputIndex += 1
        ENDFOR

	IP_latefrac_mean_array   = []
	IP_latefrac_median_array = []	

	; ARRAY[ COL, ROW ]
	FOR n=0,n_annuli-1 DO BEGIN
		IP_latefrac_array_col 	  = TRANSPOSE([IP_latefrac_array[n,*]])
;		PRINT, 'IP_latefrac_array_col: ', IP_latefrac_array_col

		IP_latefrac_array_col_pos = IP_latefrac_array_col[where(IP_latefrac_array_col GE 0., /NULL)]
;		PRINT, 'IP_latefrac_array_col_pos: ', IP_latefrac_array_col_pos

		IP_latefrac_mean_array    = [IP_latefrac_mean_array, MEAN(IP_latefrac_array_col_pos)]
		IP_latefrac_median_array  = [IP_latefrac_median_array, MEDIAN(IP_latefrac_array_col_pos)]
	ENDFOR

	PRINT, 'annulus targ_weight means:   ', IP_latefrac_mean_array
	PRINT, 'annulus targ_weight medians: ', IP_latefrac_median_array
END
