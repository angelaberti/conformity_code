; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; normsigplot:	plots normalized conformity signal
; JKE:		uses jackknife error method
; IPmatchFBF:	SF and Q IP populations have been selected to have same mass and redshift distribution WITHIN EACH FIELD
; PHI3.7:	IP matching procedure uses high mass cutoff threshold based on Moustakas13 stellar mass function results
; binsEqualThirds: late-type fractions are computed for three redshift ranges with equal # IP

PRO normsigplot_JKE_IPmatchFBFnoCosmos_PHI37_allz, outputFormat;, zmin, zmax;, dz_coeff, printEvery
	!P.FONT=0

	COMMON SHARE, xmin, xmax, xr, ymin, ymax, yr, Rmax, plotMax, dRproj, n_annuli, Rproj_array, z_low, z_high, zsuffix, zlabel
	COScolor = 'red'
	IF (outputFormat EQ 'ps') THEN (REGcolor = 'black') ELSE (REGcolor = 'white')
	PHI = 3.7
	stringPHI = strtrim(string(PHI, format='(f20.1)'),2)

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	Rmax 		= 5.
	plotMax 	= Rmax-1
	dRproj		= 1.
	n_annuli 	= CEIL(Rmax/dRproj)
	Rproj_array 	= FLOAT(dRproj)*(findgen(n_annuli+1))	

	xmin = 0
	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = -0.04
	ymax = 0.11
	yr = ymax-ymin
	
	z_low = zmin
	z_high = zmax
	zsuffix = 'allz'
	zlabel  = 'z = [' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']'

	; read in data files
;	dataIPallz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)
	dataIPallz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
	dataAllallz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 0.50 THEN string_dR = '_dR500kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

	; eliminate data with targ_weight < 1
        dataAll = dataAllallz[where(dataAllallz.targ_weight GE 1.)]
        dataIP  = dataIPallz[where( (dataIPallz.targ_weight GE 1.) AND (dataIPallz.field NE 'cosmos    ') )]

	; divide IP population into 3 redshift bins with equal numbers of IPs
        orderedRedshifts = dataIPallz[SORT(dataIPallz.zprimus)].zprimus

        lower_index = CEIL(n_elements(dataIPallz)/3.)
        upper_index = 2*FLOOR(n_elements(dataIPallz)/3.)

        lower_z = orderedRedshifts[lower_index]
        upper_z = orderedRedshifts[upper_index]

        zArray = [0.2, lower_z, upper_z, 1.0]
        PRINT, stringPHI, ' ', zArray

	!p.multi=0
	!p.charsize=1.75

	colors = ['blue', 'green', 'red']

	IF (string(outputFormat) eq 'ps') THEN BEGIN
		PS_OPEN, '~/latex/figures/normsigplot_JKE_IPmatchFBFnoCosmos_PHI37_allz', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	PLOT, [xmin,xmax], [0,0], xrange=[xmin,xmax], yrange=[ymin,ymax], LINESTYLE=1, THICK=2, xtitle='Projected Radius (Mpc)', ytitle='Normalized Conformity Signal' 
;	  title = 'Field-by-Field Matched SF & Q IP Samples ' + textoidl('(\Deltaz=2.0)')

	data = MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/noCosmos/normsig_' + zsuffix + '_targ_weight_IPmatchFBFnoCosmos_PHI' + stringPHI + string_dR + '_JKE.fits', 1)
	makePlots, data, cgcolor(REGcolor), 0.10

;  IF (COSMOScomp EQ 1) THEN LEGEND, ['Late-type IP', 'Early-type IP', 'COSMOS field excluded'], LINESTYLE=[0,2,0], colors=cgcolor([REGcolor,REGcolor,COScolor]), BOX=0, /TOP, /RIGHT
;  IF (COSMOScomp NE 1) THEN BEGIN
;	LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], BOX=0, /TOP, /RIGHT
;	XYOUTS, xmin+0.95*xr, ymin+0.25*yr, 'SF IP: ' + bigint(n_elements(UNIQ(dataIP[where(dataIP.SFQ EQ 1)].objname))), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'Q IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 0)])), ALIGNMENT=1.0
;  ENDIF
;	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
	XYOUTS, xmin+0.05*xr, ymin+0.05*yr, 'COSMOS field excluded', ALIGNMENT=0.0
	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, zlabel, ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.925*yr, 'COSMOS field excluded', ALIGNMENT=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

PRO makePlots, data, color, sigmaHeight
	COMMON SHARE

	normsig	= data.normsig
	normsig_errors  = data.normsig_errors

	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPSF = data.n_late_IPSF
;	errors_IPSF = data.errors_IPSF
	n_tot_IPQ   = data.n_tot_IPQ
	n_late_IPQ  = data.n_late_IPQ
;	errors_IPQ  = data.errors_IPQ

;	frac_IPSF   = n_late_IPSF/n_tot_IPSF
;	frac_IPQ    = n_late_IPQ/n_tot_IPQ

	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, normsig[0:plotMax], color=color
	ERRPLOT, Rproj_array+0.5*dRproj, normsig-normsig_errors, normsig+normsig_errors, color=color

;	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, frac_IPQ[0:plotMax], LINESTYLE=2, color=color
;	ERRPLOT, Rproj_array+0.5*dRproj, frac_IPQ-errors_IPQ, frac_IPQ+errors_IPQ, color=color

	sigmas = textoidl('\sigma_{0-1}=') + getSigmaRange_JKE(data, 0, 0) $
		+ textoidl('; \sigma_{1-2}=') + getSigmaRange_JKE(data, 1, 1) $
		+ textoidl('; \sigma_{1-3}=') + getSigmaRange_JKE(data, 1, 2) $
		+ textoidl('; \sigma_{3-5}=') + getSigmaRange_JKE(data, 3, 4)

;	XYOUTS, xmin+0.95*xr, ymin+sigmaHeight*yr, sigmas, ALIGNMENT=1.0, color=color, charsize=1.5
END

