PRO latefrac_BSE_IPmatchFBF_PHI37_allz_SFRquartiles, outputFormat;, zmin, zmax;, dz_coeff, printEvery
	PHI = 3.7
	stringPHI = strtrim(string(PHI, format='(f20.1)'),2)

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	; read in data files
	allDataAllz = MRDFITS('~/results/zerodSFQ_all_cart.fits', 1)
	IPdataAllz  = MRDFITS('~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
	IPdataAllz  = IPdataAllz[WHERE(IPdataAllz.IP EQ 1)]

	Rmax 		= 10.
	dRproj		= 1.
	n_annuli 	= Ceil(Rmax/dRproj)
	printEvery	= 100

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 0.50 THEN string_dR = '_dR500kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

; eliminate data with targ_weight < 1
        allDataAllz     = allDataAllz[WHERE(allDataAllz.targ_weight GE 1.)]
        dd              = IPdataAllz[WHERE(IPdataAllz.targ_weight GE 1.)]
       	dd_unique	= dd[UNIQ(dd.objname, SORT(dd.objname))]
        PRINT, 'All IP: ', n_elements(dd)
        PRINT, 'Unique IP: ', n_elements(dd_unique)

; USE ALL IPs TO DETERMINE QUARTILES
        SFRperM_med     = MEDIAN(dd.SFR-dd.mstar)
        dd_lowHalf	= dd[WHERE((dd.SFR-dd.mstar) LE SFRperM_med)]
        dd_highHalf     = dd[WHERE((dd.SFR-dd.mstar) GE SFRperM_med)]

        SFRperM_low     = MEDIAN(dd_lowHalf.SFR-dd_lowHalf.mstar)
        SFRperM_high    = MEDIAN(dd_highHalf.SFR-dd_highHalf.mstar)

        PRINT, SFRperM_low, SFRperM_med, SFRperM_high

; LOWEST
        IP_QRT1 = dd[WHERE( (dd.SFR-dd.mstar) LT SFRperM_low )]

        IP_QRT2 = dd[WHERE( ((dd.SFR-dd.mstar) GE SFRperM_low) AND ((dd.SFR-dd.mstar) LT SFRperM_med) )]
        IP_QRT3 = dd[WHERE( ((dd.SFR-dd.mstar) GE SFRperM_med) AND ((dd.SFR-dd.mstar) LT SFRperM_high) )]
; HIGEST
        IP_QRT4 = dd[WHERE( (dd.SFR-dd.mstar) GE SFRperM_high )]

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

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/latex/figures/latefracplot_BSE_SFRquartiles', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

  colors = ['red','orange','cyan','blue']
  dataAll = allDataAllz
  FOR i=0,3 DO BEGIN
	IF (i EQ 0) THEN dataIP=IP_QRT1
	IF (i EQ 1) THEN dataIP=IP_QRT2
	IF (i EQ 2) THEN dataIP=IP_QRT3
	IF (i EQ 3) THEN dataIP=IP_QRT4

	PRINT, 'Quartile ', STRTRIM(i+1,2), ': ', median(dataIP.zprimus), median(dataIP.mstar)

	get_targ_weight, dataIP, dataAll, dz_coeff, neigh_all_grand_total, neigh_late_grand_total, error_array, n_annuli, dRproj, printEvery

	n_tot_IP  = neigh_all_grand_total
	n_late_IP = neigh_late_grand_total
	errors_IP = error_array
	frac_IP   = n_late_IP/n_tot_IP

    IF (i EQ 0) THEN $
	PLOT, Rproj_array+0.5*dRproj, float(neigh_late_grand_total)/neigh_all_grand_total, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', /NODATA
;	PLOT, dd.mstar, dd.SFR, xtitle='stellar mass', ytitle='SFR', /NODATA
;	LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], BOX=0, /BOTTOM, /LEFT

	OPLOT, Rproj_array+0.5*dRproj, frac_IP, LINESTYLE=0, color=cgColor(colors[i])
	ERRPLOT, Rproj_array+0.5*dRproj, frac_IP-errors_IP, frac_IP+errors_IP, color=cgColor(colors[i])
;	OPLOT, dataIP.mstar, dataIP.SFR, PSYM=3, color=cgColor(colors[i])

        templateRow     = create_struct('Rmin', 0.0, 'Rmax', 0.0, 'n_tot_IP', 0.0, 'n_late_IP', 0.0, 'errors_IP', 0.0)
        outputStruct    = replicate(templateRow, n_annuli)

        FOR j=0,n_annuli-1 DO BEGIN
                newRow  = create_struct('Rmin', Rproj_array[j], 'Rmax', Rproj_array[j+1], $
                  'n_tot_IP', n_tot_IP[j], 'n_late_IP', n_late_IP[j], 'errors_IP', errors_IP[j])
                outputStruct[j] = newRow
        ENDFOR

	MWRFITS, outputStruct, '~/results/match_IP_sample_rigorous/jackknife_error/SFR_quartiles/latefrac_BSE_quart'+STRTRIM(i+1,2)+'.fits', /CREATE
  ENDFOR
	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

PRO get_targ_weight, dataIP, dataAll, dz_coeff, neigh_all_grand_total, neigh_late_grand_total, error_array, n_annuli, dRproj, printEvery
	IP = dataIP
	data = dataAll
	
	allIP_total_array = [] ; list of all (weighted) neighbors around all IPs by annulus
	allIP_late_array  = [] ; list of late-type (weighted) neighbors around all IPs by annulus

	PRINT,'IP galaxies: ' + strtrim(n_elements(IP),1)

	outputIndex = 0L

;	FOR i=0,(n_elements(IP)-1) DO BEGIN
	FOR i=0,9 DO BEGIN
                currentIP = IP[i]

                ; keep galaxies within dz=dz+coeff*0.005*(1+z) of IP candidate (and in same field)
                all_dz_neigh = data[ where( (data.field EQ currentIP.field) AND (data.objname NE currentIP.objname) AND $
                        (ABS(data.zprimus - currentIP.zprimus) le dz_coeff*0.005*(1.+currentIP.zprimus)) ) ]
		all_dz_neigh_dists = [SQRT((all_dz_neigh.xprop-currentIP.xprop)^2 + (all_dz_neigh.yprop-currentIP.yprop)^2)]

		currentIP_total_array = [] ; weight total of all neighbors by annulus for current IP
		currentIP_late_array  = [] ; weight total of late-type neighbors by annulus for current IP

		; GET LATE FRACTION FOR EACH ANNULUS
                FOR n=0,n_annuli-1 DO BEGIN
			currentRmin = dRproj*float(n)
			currentRmax = dRproj*float(n+1)
			
			annulus_neigh = all_dz_neigh[where((all_dz_neigh_dists GE currentRmin) AND (all_dz_neigh_dists LE currentRmax), /NULL)]
			IF (annulus_neigh NE !NULL) THEN $
				annulus_neigh_SF = annulus_neigh[where(annulus_neigh.SFQ EQ 1, /NULL)] ELSE $
				annulus_neigh_SF = []

			IF (annulus_neigh NE !NULL) THEN $
				currentIP_annulus_total = TOTAL(annulus_neigh.targ_weight) ELSE $
				currentIP_annulus_total = 0
			IF (annulus_neigh_SF NE !NULL) THEN $
				currentIP_annulus_late = TOTAL(annulus_neigh_SF.targ_weight) ELSE $
				currentIP_annulus_late = 0
			
			currentIP_total_array = [currentIP_total_array, currentIP_annulus_total]
			currentIP_late_array  = [currentIP_late_array, currentIP_annulus_late]
		ENDFOR

		allIP_total_array = [[allIP_total_array], [currentIP_total_array]]
		allIP_late_array  = [[allIP_late_array], [currentIP_late_array]]

;		IF outputIndex mod printEvery EQ 0 THEN PRINT, outputIndex
		outputIndex += 1
        ENDFOR

	; initialize weighted total and error on late-type fraction for each annulus
	neigh_all_grand_total  = []
	neigh_late_grand_total = []
	
	error_array = []

	; ARRAY[ COL, ROW ]
	FOR n=0,n_annuli-1 DO BEGIN
		neigh_all_grand_total_annulus  = TOTAL(TRANSPOSE([allIP_total_array[n,*]]))
		neigh_late_grand_total_annulus = TOTAL(TRANSPOSE([allIP_late_array[n,*]]))
		neigh_all_grand_total 	= [neigh_all_grand_total, neigh_all_grand_total_annulus]
		neigh_late_grand_total	= [neigh_late_grand_total, neigh_late_grand_total_annulus]

		; WANT TO COMPUTE THE LATE-TYPE FRACTION AT EACH ANNULUS 200 TIMES USING A RANDOM 90% (W/REPLACEMENT) OF THE IPS EACH TIME
		annulus_frac_200 = []
		FOR s=1,200 DO BEGIN
			seed = s
			tt_100	= TRANSPOSE(allIP_total_array[n,*])
			random_indices = ROUND( RANDOMU( seed, ROUND( 0.9*n_elements(tt_100) ) )*n_elements(tt_100) )
			ll_100	= TRANSPOSE(allIP_late_array[n,*])
			tt_90	= tt_100[random_indices]
			ll_90	= ll_100[random_indices] ; USE SAME INDICES AS FOR tt_100

			annulus_frac	 = TOTAL(ll_90)/TOTAL(tt_90)
			annulus_frac_200 = [annulus_frac_200, annulus_frac]
			annulus_error	 = STDDEV(annulus_frac_200)
;			PRINT, annulus_frac, annulus_error
		ENDFOR

		error_array = [error_array, annulus_error]
	ENDFOR

;	PRINT, 'all neighbor weighted totals: ', neigh_all_grand_total
;	PRINT, 'late neighbor weighted totals: ', neigh_late_grand_total
;	PRINT, 'errors: ', error_array
END
