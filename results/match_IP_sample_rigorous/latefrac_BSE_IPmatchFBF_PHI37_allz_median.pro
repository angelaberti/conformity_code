; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; latefrac:	calculates late-type fractions
; BSE:		uses bootstrap error method (doesn't apply when median is used)
; IPmatchFBF:	SF and Q IP populations have been selected to have same mass and redshift distribution WITHIN EACH FIELD
; PHI3.7:	IP matching procedure uses high mass cutoff threshold based on Moustakas13 stellar mass function results
; allz: 	late-type fraction computed for 0.2<z<1.0
; median:	median of all late-type fractions for each annulus is taken

PRO latefrac_BSE_IPmatchFBF_PHI37_allz_median, outputFormat;, zmin, zmax;, dz_coeff, printEvery
	PHI = 3.7
	stringPHI = strtrim(string(PHI, format='(f20.1)'),2)

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	; read in data files
;	dataIPallz  = MRDFITS('~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)
	dataIPallz  = MRDFITS('~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
	dataAllallz = MRDFITS('~/results/zerodSFQ_all_cart.fits', 1)

	Rmax 		= 5.
	dRproj		= .25
	n_annuli 	= CEIL(Rmax/dRproj)
	printEvery	= 100

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 0.50 THEN string_dR = '_dR500kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

	; eliminate data with targ_weight < 1
        dataAllallz = dataAllallz[where(dataAllallz.targ_weight GE 1.)]
        dataIPallz  = dataIPallz[where(dataIPallz.targ_weight GE 1.)]

	; divide IP population into 3 redshift bins with equal numbers of IPs
        orderedRedshifts = dataIPallz[SORT(dataIPallz.zprimus)].zprimus

        lower_index = CEIL(n_elements(dataIPallz)/3.)
        upper_index = 2*FLOOR(n_elements(dataIPallz)/3.)

        lower_z = orderedRedshifts[lower_index]
        upper_z = orderedRedshifts[upper_index]

        zArray = [0.2, lower_z, upper_z, 1.0]
        PRINT, stringPHI, ' ', zArray

;	dataAll = dataAll[where(dataAll.targ_weight GE 1.)]
;	dataIP  = dataIP[where(dataIP.targ_weight GE 1.)]

	!p.multi=0
	!p.charsize=1.25

	Rproj_array = float(dRproj)*(findgen(n_annuli+1))	
;	PRINT, Rproj_array

	xmin = 0
	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = 0.5
	ymax = 1
	yr = ymax-ymin

;  FOR i=0,2 DO BEGIN
;	z_low   = zArray[i]
	z_low = 0.2
;	z_high  = zArray[i+1]
	z_high = 1.0

;	zsuffix = 'T' + strtrim(i+1,1)
	zsuffix = 'allz'
	zlabel  = strtrim(string(z_low,format='(f20.2)'),1) + ' < z < ' + strtrim(string(z_high,format='(f20.2)'),1)

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsampleFBF/latefracplot' + string_dR + '_IPmatchFBF_PHI37_allz_median', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

;	dataIP  = dataIPallz[where( (dataIPallz.zprimus GE z_low) AND (dataIPallz.zprimus LE z_high) )]
;	dataAll = dataAllallz[where( (dataAllallz.zprimus GE z_low) AND (dataAllallz.zprimus LE z_high) )]
	dataIP  = dataIPallz
	dataAll = dataAllallz

;	RESTORE, 'IPSF_latefrac_array' + string_dR + '.sav'
;	RESTORE, 'IPQ_latefrac_array' + string_dR + '.sav'

	get_targ_weight, dataIP, dataAll, dz_coeff, 1, IP_latefrac_median_array, n_annuli, dRproj, printEvery, string_dR
	IPSF_latefrac_median_array = IP_latefrac_median_array

;	PLOT, Rproj_array+0.5*dRproj, float(neigh_late_grand_total)/neigh_all_grand_total, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Median Late-type Fraction', $
	PLOT, Rproj_array+0.5*dRproj, IP_latefrac_median_array, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', $
	  title = 'Field-by-Field Matched SF & Q IP Samples ' + textoidl('(\Deltaz=2.0)')
	LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], BOX=0, /BOTTOM, /LEFT

	get_targ_weight, dataIP, dataAll, dz_coeff, 0, IP_latefrac_median_array, n_annuli, dRproj, printEvery, string_dR
	IPQ_latefrac_median_array = IP_latefrac_median_array

;	OPLOT, Rproj_array+0.5*dRproj, frac_IPQ, LINESTYLE=2
	OPLOT, Rproj_array+0.5*dRproj, IP_latefrac_median_array, LINESTYLE=2

	XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'SF IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 1)])), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Q IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 0)])), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, zlabel, ALIGNMENT=1.0

	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'

	templateRow     = create_struct('Rmin', 0.0, 'Rmax', 0.0, 'IPSF_latefrac_median', 0.0, 'IPQ_latefrac_median', 0.0)
	outputStruct    = replicate(templateRow, n_annuli)

	FOR j=0,n_annuli-1 DO BEGIN
		newRow  = create_struct('Rmin', Rproj_array[j], 'Rmax', Rproj_array[j+1], $
		  'IPSF_latefrac_median', IPSF_latefrac_median_array[j], 'IPQ_latefrac_median', IPQ_latefrac_median_array[j])
		outputStruct[j] = newRow
	ENDFOR

	PRINT, Rproj_array
	PRINT, (IPSF_latefrac_median_array - IPQ_latefrac_median_array)/((IPSF_latefrac_median_array + IPQ_latefrac_median_array)/2)

	MWRFITS, outputStruct, '~/results/match_IP_sample_rigorous/latefrac_' + zsuffix + '_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_median.fits', /create
;  ENDFOR
END

PRO get_targ_weight, dataIP, dataAll, dz_coeff, IPSFstatus, IP_latefrac_median_array, n_annuli, dRproj, printEvery, string_dR
	data_isIP = dataIP
	data = dataAll
	
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
                        (ABS(data.zprimus - currentIP.zprimus) LE dz_coeff*0.005*(1.+currentIP.zprimus)), /NULL ) ]
;		PRINT, 'all_dz_neigh_dists: ', n_elements(all_dz_neigh_dists)
		
		IF (all_dz_neigh NE !NULL) THEN (all_dz_neigh_dists = [SQRT((all_dz_neigh.xprop-currentIP.xprop)^2 + (all_dz_neigh.yprop-currentIP.yprop)^2)]) $
		  ELSE all_dz_neigh_dists = []
		
		currentIP_latefrac_array = []		

		; GET LATE FRACTION FOR EACH ANNULUS
                FOR n=0,n_annuli-1 DO BEGIN
;			PRINT, 'dRproj: ', dRproj
		
			currentRmin = dRproj*float(n)
			currentRmax = dRproj*float(n+1)

;			PRINT, currentRmin, currentRmax
			
			IF (all_dz_neigh NE !NULL) THEN (annulus_neigh = all_dz_neigh[where((all_dz_neigh_dists GE currentRmin) AND (all_dz_neigh_dists LE currentRmax), /NULL)]) $
			  ELSE (annulus_neigh = [])
			IF (annulus_neigh NE !NULL) THEN (annulus_neigh_SF = annulus_neigh[where(annulus_neigh.SFQ EQ 1, /NULL)]) $
			  ELSE (annulus_neigh_SF = [])
;			PRINT, 'annulus_neigh: ', n_elements(annulus_neigh)
;			PRINT, 'annulus_neigh_SF: ', n_elements(annulus_neigh_SF)

			; DEFAULT
			annulus_latefrac = -99

			; NUMBER OF GALAXIES IN ANNULUS > 0
			IF (annulus_neigh NE !NULL) THEN BEGIN
				IF (annulus_neigh_SF EQ !NULL) THEN (annulus_latefrac = 0.) ELSE $
				  annulus_latefrac = TOTAL(annulus_neigh_SF.targ_weight)/TOTAL(annulus_neigh.targ_weight)
			ENDIF

			; ADD ANNULUS LATE FRACTION TO ARRAY FOR SPECIFIC IP
			currentIP_latefrac_array = [currentIP_latefrac_array, annulus_latefrac]
		ENDFOR

		IP_latefrac_array = [[IP_latefrac_array], [currentIP_latefrac_array]]

		IF (IPSFstatus EQ 1) THEN BEGIN
			SAVE, IP_latefrac_array, FILE='IPSF_latefrac_array' + string_dR + '.sav'
		ENDIF ELSE BEGIN
			SAVE, IP_latefrac_array, FILE='IPQ_latefrac_array' + string_dR + '.sav'
		ENDELSE

		IF outputIndex mod printEvery EQ 0 THEN PRINT, outputIndex
		outputIndex += 1
        ENDFOR

	IP_latefrac_median_array = []
	
	; ARRAY[ COL, ROW ]
	FOR n=0,n_annuli-1 DO BEGIN
		IP_latefrac_array_col     = TRANSPOSE([IP_latefrac_array[n,*]])
;		PRINT, 'IP_latefrac_array_col: ', IP_latefrac_array_col

		IP_latefrac_array_col_pos = IP_latefrac_array_col[where(IP_latefrac_array_col GE 0., /NULL)]
;		PRINT, 'IP_latefrac_array_col_pos: ', IP_latefrac_array_col_pos

;		IP_latefrac_mean_array    = [IP_latefrac_mean_array, MEAN(IP_latefrac_array_col_pos)]
		IP_latefrac_median_array  = [IP_latefrac_median_array, MEDIAN(IP_latefrac_array_col_pos)]
	ENDFOR

;	PRINT, 'annulus targ_weight medians: ', IP_latefrac_median_array
END
