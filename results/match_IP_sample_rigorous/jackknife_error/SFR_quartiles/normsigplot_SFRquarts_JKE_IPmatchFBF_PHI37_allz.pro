; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; normsigplot:	plots normalized conformity signal
; JKE:		uses jackknife error method
; IPmatchFBF:	SF and Q IP populations have been selected to have same mass and redshift distribution WITHIN EACH FIELD
; PHI3.7:	IP matching procedure uses high mass cutoff threshold based on Moustakas13 stellar mass function results
; binsEqualThirds: late-type fractions are computed for three redshift ranges with equal # IP

PRO normsigplot_SFRquarts_JKE_IPmatchFBF_PHI37_allz, outputFormat;, zmin, zmax;, dz_coeff, printEvery
	!P.FONT=0

	COMMON SHARE, xmin, xmax, xr, ymin, ymax, yr, Rmax, plotMax, dRproj, n_annuli, Rproj_array, z_low, z_high, zsuffix, zlabel
	COScolor = 'red'
	IF (outputFormat EQ 'ps') THEN (REGcolor = 'black') ELSE (REGcolor = 'white')
	PHI = 3.7
	stringPHI = strtrim(string(PHI, format='(f20.1)'),2)

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	Rmax 		= 15.
	plotMax 	= Rmax-1
	dRproj		= 1.
	n_annuli 	= CEIL(Rmax/dRproj)
	Rproj_array 	= FLOAT(dRproj)*(findgen(n_annuli+1))	

;	Rproj_array = [0.0, 0.5, 1.0+FINDGEN(Rmax)]
;	Rplot_array = [0.25, 0.75, 1.5+FINDGEN(Rmax-1)]
	Rplot_array = Rproj_array[0:n_elements(Rproj_array)-2]+0.5

	xmin = 0
	xmax = 15
	xr = xmax-xmin

	ymin = -0.04
	ymax = 0.11
	yr = ymax-ymin
	
	z_low = zmin
	z_high = zmax
	zsuffix = 'allz'
	zlabel  = 'z = [' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']'

	; read in data files
;	dataIPallz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/SFR_quartiles/matchedIPsampleFBF.fits', 1)
	dataIPallz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
	dataAllallz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 0.50 THEN string_dR = '_dR500kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

	; eliminate data with targ_weight < 1
        dataAll = dataAllallz[where(dataAllallz.targ_weight GE 1.)]
        dataIP  = dataIPallz[where(dataIPallz.targ_weight GE 1.)]

	; divide IP population into 3 redshift bins with equal numbers of IPs
        orderedRedshifts = dataIPallz[SORT(dataIPallz.zprimus)].zprimus

        lower_index = CEIL(n_elements(dataIPallz)/3.)
        upper_index = 2*FLOOR(n_elements(dataIPallz)/3.)

        lower_z = orderedRedshifts[lower_index]
        upper_z = orderedRedshifts[upper_index]

        zArray = [0.2, lower_z, upper_z, 1.0]
        PRINT, stringPHI, ' ', zArray

	!P.MULTI=0
	!P.CHARSIZE=1.5

	colors = ['blue', 'green', 'red']

	IF (string(outputFormat) eq 'ps') THEN BEGIN
;		PS_OPEN, '~/latex/figures/normsigplot_JKE_IPmatchFBF_PHI37_allz', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=6, YS=5, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	PLOT, [xmin,xmax], [0,0], xrange=[xmin,xmax], yrange=[ymin,ymax], LINESTYLE=1, xtitle='Projected Radius (Mpc)', ytitle='Normalized Conformity Signal', THICK=2
;	  title = 'Field-by-Field Matched SF & Q IP Samples ' + textoidl('(\Deltaz=2.0)')

;	data = MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/SFR_quartiles/normsig_' + zsuffix + '_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_JKE.fits', 1)
;	data = MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/SFR_quartiles/normsig_matchedIPsampleFBF_allIPquarts.fits', 1)
	data = MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/SFR_quartiles/normsig_matchedIPsampleFBF_highLowQuarts.fits', 1)
	makePlots, data, cgcolor(REGcolor), 0.10

;  IF (COSMOScomp EQ 1) THEN LEGEND, ['Late-type IP', 'Early-type IP', 'COSMOS field excluded'], LINESTYLE=[0,2,0], colors=cgcolor([REGcolor,REGcolor,COScolor]), BOX=0, /TOP, /RIGHT
;  IF (COSMOScomp NE 1) THEN BEGIN
;	LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], BOX=0, /TOP, /RIGHT
;	XYOUTS, xmin+0.95*xr, ymin+0.25*yr, 'SF IP: ' + bigint(n_elements(UNIQ(dataIP[where(dataIP.SFQ EQ 1)].objname))), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'Q IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 0)])), ALIGNMENT=1.0
;  ENDIF
;	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, zlabel, ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

PRO makePlots, data, color, sigmaHeight
	COMMON SHARE

	normsig	= data.normsig
	normsig_errors  = data.normsig_errors

	n_tot_IPhigh  = data.n_tot_IPhigh
	n_late_IPhigh = data.n_late_IPhigh
;	errors_IPhigh = data.errors_IPhigh
	n_tot_IPlow   = data.n_tot_IPlow
	n_late_IPlow  = data.n_late_IPlow
;	errors_IPlow  = data.errors_IPlow

;	frac_IPhigh   = n_late_IPhigh/n_tot_IPhigh
;	frac_IPlow    = n_late_IPlow/n_tot_IPlow

;	OPLOT, Rplot_array[0:plotMax], normsig[0:plotMax], color=color
;	ERRPLOT, Rplot_array[0:plotMax], normsig-normsig_errors, normsig+normsig_errors, color=color

	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, normsig[0:plotMax], color=color
	ERRPLOT, Rproj_array[0:plotMax]+0.5*dRproj, normsig-normsig_errors, normsig+normsig_errors, color=color

;	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, frac_IPlow[0:plotMax], LINESTYLE=2, color=color
;	ERRPLOT, Rproj_array+0.5*dRproj, frac_IPlow-errors_IPlow, frac_IPlow+errors_IPlow, color=color

;	sigmas = textoidl('\sigma_{0-1}=') + getSigmaRange_JKE(data, 0, 0) $
;		+ textoidl('; \sigma_{1-2}=') + getSigmaRange_JKE(data, 1, 1) $
;		+ textoidl('; \sigma_{1-3}=') + getSigmaRange_JKE(data, 1, 2) $
;		+ textoidl('; \sigma_{3-5}=') + getSigmaRange_JKE(data, 3, 4)

;	XYOUTS, xmin+0.95*xr, ymin+sigmaHeight*yr, sigmas, ALIGNMENT=1.0, color=color, charsize=1.5
END
