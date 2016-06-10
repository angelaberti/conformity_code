; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; latefrac:	calculates late-type fractions
; BSE:		uses bootstrap error method
; IPmatchFBF:	SF and Q IP populations have been selected to have same mass and redshift distribution WITHIN EACH FIELD
; PHI3.7:	IP matching procedure uses high mass cutoff threshold based on Moustakas13 stellar mass function results

PRO latefrac_BSE_allz, set
	PHI = 3.7
	stringPHI = strtrim(string(PHI, format='(f20.1)'),2)

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	Rmax 		= 15.
	printEvery	= 10

	dRproj=1.
	
        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 0.50 THEN string_dR = '_dR500kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

; Matched M* and Redshift
; DONE
  IF (set EQ 0) THEN BEGIN
	dataIPallz  = MRDFITS('~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
	zsuffix = 'allz'
	outputPath = '~/results/match_IP_sample_rigorous/unmatchedExample/latefrac_' + zsuffix + '_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_BSE.fits'
  ENDIF

; No Matching
; 
  IF (set EQ 1) THEN BEGIN
	dataIPallz  = MRDFITS('~/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits', 1)
	outputPath = '~/results/match_IP_sample_rigorous/unmatchedExample/latefrac_conservative_allz_dR1Mpc_BSE.fits'
;	dataIPallz  = MRDFITS('~/results/match_IP_sample_rigorous/unmatchedExample/noMatching/unmatchedSampleFBF_PHI3.7.fits', 1)
;	outputPath = '~/results/match_IP_sample_rigorous/unmatchedExample/latefrac_noMatching_allz_dR1Mpc_BSE.fits'
  ENDIF

; Matched Redshift Only
; DONE
  IF (set EQ 2) THEN BEGIN
	dataIPallz  = MRDFITS('~/results/match_IP_sample_rigorous/unmatchedExample/matchRedshiftOnly/matchedIPsampleRedshiftOnlyFBF_PHI3.7.fits', 1)
	outputPath = '~/results/match_IP_sample_rigorous/unmatchedExample/latefrac_matchRedshiftOnly_allz_dR1Mpc_BSE.fits'
;	dataIPallz  = MRDFITS('~/results/single_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_singleMass.fits', 1)
;	outputPath = '~/results/match_IP_sample_rigorous/unmatchedExample/latefrac_singleMass_allz_dR1Mpc_BSE.fits'
  ENDIF

; Matched M* Only
;
  IF (set EQ 3) THEN BEGIN
	dataIPallz  = MRDFITS('~/results/match_IP_sample_rigorous/unmatchedExample/matchMassOnly/matchedIPsampleMassOnlyFBF_PHI3.7.fits', 1)
	outputPath = '~/results/match_IP_sample_rigorous/unmatchedExample/latefrac_matchMassOnly_allz_dR1Mpc_BSE.fits'
  ENDIF

	dataAllallz = MRDFITS('~/results/zerodSFQ_all_cart.fits', 1)

; eliminate data with targ_weight < 1
        dataAllallz = dataAllallz[where(dataAllallz.targ_weight GE 1.)]
        dataIPallz  = dataIPallz[where(dataIPallz.targ_weight GE 1.)]

	!P.MULTI=0
	!P.CHARSIZE=1.25

	Rproj_array = FINDGEN(Rmax+1)
;	Rproj_array = [0.0, 0.5, 1.0+FINDGEN(Rmax)]
	Rplot_array = FINDGEN(Rmax) + 0.5*dRproj
;	Rplot_array = Rproj_array[0:n_elements(Rproj_array)-2] + 0.5
;	Rplot_array = [0.25, 0.75, 1.5+FINDGEN(Rmax-1)]

	PRINT, 'Rproj_array: ', Rproj_array
	PRINT, 'Rplot_array: ', Rplot_array

	xmin = 0
	xmax = 15
	xr = xmax-xmin

	ymin = 0.6
	ymax = 0.88
	yr = ymax-ymin

	z_low = 0.2
	z_high = 1.0

	zsuffix = 'allz'
	zlabel  = strtrim(string(z_low,format='(f20.2)'),1) + ' < z < ' + strtrim(string(z_high,format='(f20.2)'),1)

        SET_PLOT, 'X'

	dataIP  = dataIPallz[WHERE(dataIPallz.IP EQ 1)]
	dataAll = dataAllallz

  get_targ_weight, dataIP, dataAll, dz_coeff, 1, neigh_all_grand_total, neigh_late_grand_total, error_array, n_annuli, Rproj_array, printEvery

	n_tot_IPSF  = neigh_all_grand_total
	n_late_IPSF = neigh_late_grand_total
	errors_IPSF = error_array
	frac_IPSF   = n_late_IPSF/n_tot_IPSF

	PLOT, Rplot_array, FLOAT(neigh_late_grand_total)/neigh_all_grand_total, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', /NODATA
;	  title = 'Field-by-Field Matched SF & Q IP Samples ' + textoidl('(\Deltaz=2.0)'), /NODATA
	LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], BOX=0, /BOTTOM, /LEFT

	OPLOT, Rplot_array, frac_IPSF, LINESTYLE=0
	ERRPLOT, Rplot_array, frac_IPSF-errors_IPSF, frac_IPSF+errors_IPSF

  get_targ_weight, dataIP, dataAll, dz_coeff, 0, neigh_all_grand_total, neigh_late_grand_total, error_array, n_annuli, Rproj_array, printEvery

	n_tot_IPQ  = neigh_all_grand_total
	n_late_IPQ = neigh_late_grand_total
	errors_IPQ = error_array
	frac_IPQ   = n_late_IPQ/n_tot_IPQ

	OPLOT, Rplot_array, frac_IPQ, LINESTYLE=2
	ERRPLOT, Rplot_array, frac_IPQ-errors_IPQ, frac_IPQ+errors_IPQ

        templateRow     = create_struct('Rmin', 0.0, 'Rmax', 0.0, 'n_tot_IPSF', 0.0, 'n_late_IPSF', 0.0, 'errors_IPSF', 0.0, 'n_tot_IPQ', 0.0, 'n_late_IPQ', 0.0, 'errors_IPQ', 0.0)
        outputStruct    = replicate(templateRow, n_elements(Rplot_array))

        FOR j=0,n_elements(Rproj_array)-2 DO BEGIN
                newRow  = create_struct('Rmin', Rproj_array[j], 'Rmax', Rproj_array[j+1], $
                  'n_tot_IPSF', n_tot_IPSF[j], 'n_late_IPSF', n_late_IPSF[j], 'errors_IPSF', errors_IPSF[j], 'n_tot_IPQ', n_tot_IPQ[j], 'n_late_IPQ', n_late_IPQ[j], 'errors_IPQ', errors_IPQ[j])
                outputStruct[j] = newRow
;		PRINT, newRow
        ENDFOR
	
	PRINT, outputStruct

;	MWRFITS, outputStruct, '~/results/match_IP_sample_rigorous/latefrac_' + zsuffix + '_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_BSE.fits', /create
	MWRFITS, outputStruct, outputPath, /CREATE
END

PRO get_targ_weight, dataIP, dataAll, dz_coeff, IPSFstatus, neigh_all_grand_total, neigh_late_grand_total, error_array, n_annuli, Rproj_array, printEvery
	data_isIP = dataIP
	data = dataAll
	
	allIP_total_array = [] ; list of all (weighted) neighbors around all IPs by annulus
	allIP_late_array  = [] ; list of late-type (weighted) neighbors around all IPs by annulus

        IP = data_isIP[where(data_isIP.SFQ EQ IPSFstatus)]
	IF (IPSFstatus EQ 1) THEN (IPtype = 'Star-forming') ELSE (IPtype = 'Quiescent')
	PRINT, IPtype + ' IP galaxies: ' + strtrim(n_elements(IP),1)

	outputIndex = 0L

	FOR i=0,(n_elements(IP)-1) DO BEGIN
;	FOR i=0,199 DO BEGIN
                currentIP = IP[i]

                ; keep galaxies within dz=dz+coeff*0.005*(1+z) of IP candidate (and in same field)
                all_dz_neigh = data[ where( (data.field EQ currentIP.field) AND (data.objname NE currentIP.objname) AND $
                        (ABS(data.zprimus - currentIP.zprimus) le dz_coeff*0.005*(1.+currentIP.zprimus)) ) ]
		all_dz_neigh_dists = [SQRT((all_dz_neigh.xprop-currentIP.xprop)^2 + (all_dz_neigh.yprop-currentIP.yprop)^2)]

		currentIP_total_array = [] ; weight total of all neighbors by annulus for current IP
		currentIP_late_array  = [] ; weight total of late-type neighbors by annulus for current IP

		; GET LATE FRACTION FOR EACH ANNULUS
                FOR n=0,n_elements(Rproj_array)-2 DO BEGIN
			currentRmin = Rproj_array[n]
			currentRmax = Rproj_array[n+1]
			
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

		IF outputIndex mod printEvery EQ 0 THEN PRINT, outputIndex

		outputIndex += 1
        ENDFOR

	; initialize weighted total and error on late-type fraction for each annulus
	neigh_all_grand_total  = []
	neigh_late_grand_total = []
	
	error_array = []

	; ARRAY[ COL, ROW ]
	FOR n=0,n_elements(Rproj_array)-2 DO BEGIN
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

	PRINT, 'all neighbor weighted totals: ', neigh_all_grand_total
	PRINT, 'late neighbor weighted totals: ', neigh_late_grand_total
	PRINT, 'errors: ', error_array
END
