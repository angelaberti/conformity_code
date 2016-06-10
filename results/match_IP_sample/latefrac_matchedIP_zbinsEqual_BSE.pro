; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

PRO latefrac_matchedIP_zbinsEqual_BSE, outputFormat;, zmin, zmax;, dz_coeff, printEvery
	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

        dataIPallz  = MRDFITS('~/conformity/results/match_IP_sample/matchedIPsample.fits', 1)
	dataAllallz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

	Rmax 		= 15.
	dRproj		= 0.25
	n_annuli 	= Ceil(Rmax/dRproj)
	printEvery	= 100

	dataAllallz = dataAllallz[where(dataAllallz.targ_weight GE 1.)]
	dataIPallz  = dataIPallz[where(dataIPallz.targ_weight GE 1.)]

	med_z = median(dataIPallz.zprimus)
	dataLower = dataIPallz[where(dataIPallz.zprimus LE med_z)]
	dataUpper = dataIPallz[where(dataIPallz.zprimus GE med_z)]

	med_z_lower = median(dataLower.zprimus)
	med_z_upper = median(dataUpper.zprimus)

	zArray = [0.2, med_z_lower, med_z, med_z_upper, 1.0]
	PRINT, zArray

	!p.multi=0
	!p.charsize=1.1

	Rproj_array = float(dRproj)*(findgen(n_annuli+1)+0.5)

	xmin = 0
	xmax = 5
;	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = 0.5
	ymax = 1
	yr = ymax-ymin

	FOR i=0,3 DO BEGIN
		z_low	= zArray[i]
		z_high	= zArray[i+1]	

;		zsuffix = strtrim(string(0.2+0.2*i,format='(f20.1)'),1) + '_' + strtrim(string(0.2+0.2*(i+1),format='(f20.1)'),1)
		zsuffix = 'Q' + strtrim(i+1,1)
;        	zlabel  = strtrim(string(0.2+0.2*i,format='(f20.1)'),1) + ' < z < ' + strtrim(string(0.2+0.2*(i+1),format='(f20.1)'),1)
        	zlabel  = strtrim(string(z_low,format='(f20.2)'),1) + ' < z < ' + strtrim(string(z_high,format='(f20.2)'),1)
		
;		PRINT, zlabel

        	IF (string(outputFormat) eq 'ps') THEN BEGIN
        	        PS_OPEN, '~/figures/latefracplot_matchedIP_z' + zsuffix + '_BSE', THICK=3, /ENCAP
        	        DEVICE, /INCH, XS=8, YS=6
        	ENDIF ELSE BEGIN
        	        SET_PLOT, 'X'
        	        THICK=1
        	ENDELSE

        	dataIP  = dataIPallz[where( (dataIPallz.zprimus GE z_low) AND (dataIPallz.zprimus LE z_high) )]
        	dataAll = dataAllallz[where( (dataAllallz.zprimus GE z_low) AND (dataAllallz.zprimus LE z_high) )]

		get_targ_weight, dataIP, dataAll, dz_coeff, 1, neigh_all_grand_total, neigh_late_grand_total, error_array, n_annuli, dRproj, printEvery

		n_tot_IPSF  = neigh_all_grand_total
		n_late_IPSF = neigh_late_grand_total
		errors_IPSF = error_array
		frac_IPSF   = n_late_IPSF/n_tot_IPSF

		PLOT, Rproj_array, float(neigh_late_grand_total)/neigh_all_grand_total, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', $
		  title = 'Matched SF and Q IP Samples ' + textoidl('(\Deltaz=2.0)'), /NODATA
		LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], BOX=0, /BOTTOM, /LEFT

		OPLOT, Rproj_array, frac_IPSF, LINESTYLE=0
		ERRPLOT, Rproj_array, frac_IPSF-errors_IPSF, frac_IPSF+errors_IPSF

		get_targ_weight, dataIP, dataAll, dz_coeff, 0, neigh_all_grand_total, neigh_late_grand_total, error_array, n_annuli, dRproj, printEvery

		n_tot_IPQ  = neigh_all_grand_total
		n_late_IPQ = neigh_late_grand_total
		errors_IPQ = error_array
		frac_IPQ   = n_late_IPQ/n_tot_IPQ

		OPLOT, Rproj_array, frac_IPQ, LINESTYLE=2
		ERRPLOT, Rproj_array, frac_IPQ-errors_IPQ, frac_IPQ+errors_IPQ

		XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'SF IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 1)])), ALIGNMENT=1.0
		XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Q IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 0)])), ALIGNMENT=1.0
		XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
		XYOUTS, xmin+0.95*xr, ymin+0.05*yr, zlabel, ALIGNMENT=1.0

		XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0

		IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
		SET_PLOT, 'X'

		templateRow     = create_struct('Rmin', 0.0, 'Rmax', 0.0, 'n_tot_IPSF', 0.0, 'n_late_IPSF', 0.0, 'errors_IPSF', 0.0, 'n_tot_IPQ', 0.0, 'n_late_IPQ', 0.0, 'errors_IPQ', 0.0)
		outputStruct    = replicate(templateRow, n_annuli)

		FOR j=0,n_annuli-1 DO BEGIN
			newRow  = create_struct('Rmin', Rproj_array[j], 'Rmax', Rproj_array[j+1], $
			  'n_tot_IPSF', n_tot_IPSF[j], 'n_late_IPSF', n_late_IPSF[j], 'errors_IPSF', errors_IPSF[j], 'n_tot_IPQ', n_tot_IPQ[j], 'n_late_IPQ', n_late_IPQ[j], 'errors_IPQ', errors_IPQ[j])
			outputStruct[j] = newRow
		ENDFOR

		MWRFITS, outputStruct, '~/conformity/results/match_IP_sample/latefrac_' + zsuffix + '_targ_weight_matchedIPsample_BSE.fits', /create
	ENDFOR
END

PRO get_targ_weight, dataIP, dataAll, dz_coeff, IPSFstatus, neigh_all_grand_total, neigh_late_grand_total, error_array, n_annuli, dRproj, printEvery
	data_isIP = dataIP
	data = dataAll
	
	allIP_total_array = [] ; list of all (weighted) neighbors around all IPs by annulus
	allIP_late_array  = [] ; list of late-type (weighted) neighbors around all IPs by annulus

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

		; WANT TO COMPUTE THE LATE-TYPE FRACTION AT EACH ANNULUS 200 TIMES USING A RANDOM 90% OF THE IPS EACH TIME
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
