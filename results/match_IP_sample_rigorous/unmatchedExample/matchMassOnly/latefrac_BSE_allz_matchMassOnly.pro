; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; latefrac:	calculates late-type fractions
; BSE:		uses bootstrap error method

PRO latefrac_BSE_allz_matchMassOnly, outputFormat;, zmin, zmax;, dz_coeff, printEvery
  IPdatafiles = ['~/conformity/results/match_IP_sample_rigorous/unmatchedExample/matchMassOnly/matchedIPsampleMassOnlyFBF_PHI3.7.fits']
;	'~/conformity/results/default_parameters/IP_data/zerodSFQ_IP_dz2.0_dm0.5.fits', $
;	'~/conformity/results/single_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_singleMass.fits', $
;	'~/conformity/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits', $
;	'~/conformity/results/conservative+0.3dex/IP_data/zerodSFQ_IP_dz2.0_dm0.3.fits']
;	'~/conformity/results/variable_mass+1.0dex/IP_data/zerodSFQ_IP_dz2.0.fits']
;	'~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits']

;  titles = ['Default Parameters', 'Variable Mass Cutoff (+1.0 dex for SF)']
;  tags = ['default', 'var1.0dex'] ; DONE
   tags = ['matchMassOnly'] ; DONE
;  tags = ['conservative', 'conservative0.3dex']

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	Rmax 		= 15.
	dRproj		= 1.
	n_annuli 	= CEIL(Rmax/dRproj)
	printEvery	= 100

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

	ERASE
	!P.MULTI=0
	!P.CHARSIZE=1.5

	Rproj_array = float(dRproj)*(findgen(n_annuli+1))	
;	PRINT, Rproj_array

	xmin = 0
	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = 0.5
	ymax = 1
	yr = ymax-ymin

	z_low = 0.2
	z_high = 1.0

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsampleFBF/latefracplot_unmatchedIPsampleCompare', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=10
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	dataAll_allz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)
	; eliminate data with targ_weight < 1
        dataAll_allz = dataAll_allz[where(dataAll_allz.targ_weight GE 1.)]

  FOR i=0,n_elements(IPdatafiles)-1 DO BEGIN
	PRINT, IPdatafiles[i]
	IPstats, IPdatafiles[i]

	dataIP_allz  = MRDFITS(IPdatafiles[i], 1)

	; eliminate data with targ_weight < 1
        dataIP_allz  = dataIP_allz[where((dataIP_allz.targ_weight GE 1.) AND (dataIP_allz.IP EQ 1))]

	dataAll = dataAll_allz
	dataIP = dataIP_allz

	zlabel = 'z=[' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']'

	get_targ_weight, dataIP, dataAll, dz_coeff, 1, neigh_all_grand_total, neigh_late_grand_total, error_array, n_annuli, dRproj, printEvery

	n_tot_IPSF  = neigh_all_grand_total
	n_late_IPSF = neigh_late_grand_total
	errors_IPSF = error_array
	frac_IPSF   = n_late_IPSF/n_tot_IPSF

	PLOT, Rproj_array+0.5*dRproj, float(neigh_late_grand_total)/neigh_all_grand_total, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', $
	  /NODATA
;	  title = 'Field-by-Field Matched SF & Q IP Samples ' + textoidl('(\Deltaz=2.0)'), /NODATA
;	LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], BOX=0, /TOP, /RIGHT
;	XYOUTS, xmin+0.5*xr, ymin+0.925*yr, titles[i], ALIGNMENT=0.5

	OPLOT, Rproj_array+0.5*dRproj, frac_IPSF, LINESTYLE=0
	ERRPLOT, Rproj_array+0.5*dRproj, frac_IPSF-errors_IPSF, frac_IPSF+errors_IPSF

	get_targ_weight, dataIP, dataAll, dz_coeff, 0, neigh_all_grand_total, neigh_late_grand_total, error_array, n_annuli, dRproj, printEvery

	n_tot_IPQ  = neigh_all_grand_total
	n_late_IPQ = neigh_late_grand_total
	errors_IPQ = error_array
	frac_IPQ   = n_late_IPQ/n_tot_IPQ

	OPLOT, Rproj_array+0.5*dRproj, frac_IPQ, LINESTYLE=2
	ERRPLOT, Rproj_array+0.5*dRproj, frac_IPQ-errors_IPQ, frac_IPQ+errors_IPQ

;	XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'SF IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 1)])), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Q IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 0)])), ALIGNMENT=1.0

	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, textoidl('Median M_{*} = ') + decimal(MEDIAN(dataIP[WHERE(dataIP.SFQ EQ 1)].mstar),2) + ' (SF) ' + $
		decimal(MEDIAN(dataIP[WHERE(dataIP.SFQ EQ 0)].mstar),2) + ' (Q)', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, textoidl('Median z = ') + decimal(MEDIAN(dataIP[WHERE(dataIP.SFQ EQ 1)].zprimus),2) + ' (SF) ' + $
		decimal(MEDIAN(dataIP[WHERE(dataIP.SFQ EQ 0)].zprimus),2) + ' (Q)', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
	XYOUTS, xmin+0.05*xr, ymin+0.05*yr, zlabel, ALIGNMENT=0.0

;	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'

	templateRow     = create_struct('Rmin', 0.0, 'Rmax', 0.0, 'n_tot_IPSF', 0.0, 'n_late_IPSF', 0.0, 'errors_IPSF', 0.0, 'n_tot_IPQ', 0.0, 'n_late_IPQ', 0.0, 'errors_IPQ', 0.0)
	outputStruct    = replicate(templateRow, n_annuli)

	FOR j=0,n_annuli-1 DO BEGIN
		newRow  = create_struct('Rmin', Rproj_array[j], 'Rmax', Rproj_array[j+1], $
		  'n_tot_IPSF', n_tot_IPSF[j], 'n_late_IPSF', n_late_IPSF[j], 'errors_IPSF', errors_IPSF[j], 'n_tot_IPQ', n_tot_IPQ[j], 'n_late_IPQ', n_late_IPQ[j], 'errors_IPQ', errors_IPQ[j])
		outputStruct[j] = newRow
	ENDFOR

	MWRFITS, outputStruct, '~/conformity/results/match_IP_sample_rigorous/unmatchedExample/latefrac_' + tags[i] + '_allz' + string_dR + '_BSE.fits', /CREATE
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
