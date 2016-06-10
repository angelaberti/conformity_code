; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; normsigplot:	plots normalized conformity signal
; JKE:		uses jackknife error method
; IPmatchFBF:	SF and Q IP populations have been selected to have same mass and redshift distribution WITHIN EACH FIELD
; PHI37:	IP matching procedure uses high mass cutoff threshold based on Moustakas13 stellar mass function results
; binsEqualThirds: late-type fractions are computed for three redshift ranges with equal # IP

PRO normsigplot_JKE_IPmatchFBF_PHI37_massThirdsEqualIP, outputFormat;, zmin, zmax;, dz_coeff, printEvery
	!P.FONT=0

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
	plotMax 	= Rmax-1
	dRproj		= 1.
	n_annuli 	= Ceil(Rmax/dRproj)
	printEvery	= 100

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

	; eliminate data with targ_weight < 1
        dataAllallz = dataAllallz[where(dataAllallz.targ_weight GE 1.)]
        dataIPallz  = dataIPallz[where(dataIPallz.targ_weight GE 1.)]

	; divide IP population into 3 mass bins with equal numbers of IPs
	lower_index = CEIL(n_elements(dataIPallz)/3.)
	upper_index = 2*FLOOR(n_elements(dataIPallz)/3.)
	orderedMasses = dataIPallz[SORT(dataIPallz.mstar)].mstar
	lowerMass = orderedMasses(lower_index)
	upperMass = orderedMasses(upper_index)
	massArray = [min(dataIPallz.mstar), lowerMass, upperMass, max(dataIPallz.mstar)]
        PRINT, massArray

	!P.MULTI=0
	!P.CHARSIZE=1.75

	Rproj_array = float(dRproj)*(findgen(n_annuli+1))	
;	PRINT, Rproj_array

	xmin = 0
	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = -0.04
	ymax = 0.11
	yr = ymax-ymin
	
;	colors = ['blue', 'green', 'red']
	colors = ['purple', 'magenta', 'orange']
	ls = [0,2,3]
	legend = []

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/latex/figures/normsigplot_JKE_IPmatchFBF_PHI37_massThirdsEqualIP', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	deltaR = [-0.05, 0, 0.05]

  FOR i=0,2 DO BEGIN
	mass_low   = massArray[i]
	mass_high  = massArray[i+1]

	massSuffix = 'M' + strtrim(i+1,1) + 'equalIP'
	IF (i EQ 0) THEN $
;	  masslabel  = textoidl('log (M_{stellar}/M_{\odot}) = [') + decimal(mass_low,2) + ', ' + decimal(mass_high,2) + ']' $
	  masslabel  = textoidl('log ( M_{stellar} / M') + sunsymbol() + ' ) = [' + decimal(mass_low,2) + ', ' + decimal(mass_high,2) + ']' $
	ELSE $
	  masslabel  = '[' + decimal(mass_low,2) + ', ' + decimal(mass_high,2) + ']'

	dataIP  = dataIPallz[where( (dataIPallz.mstar GE mass_low) AND (dataIPallz.mstar LE mass_high) )]
;	dataAll = dataAllallz[where( (dataAllallz.mstar GE mass_low) AND (dataAllallz.mstar LE mass_high) )]
	dataAll = dataAllallz

	data = MRDFITS('~/results/match_IP_sample_rigorous/jackknife_error/massThirds/normsig_' + massSuffix + '_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_JKE.fits', 1)

	normsig		= data.normsig
	normsig_errors	= data.normsig_errors
	
	n_tot_IPSF  	= data.n_tot_IPSF
	n_late_IPSF	= data.n_late_IPSF
	n_tot_IPQ  	= data.n_tot_IPQ
	n_late_IPQ	= data.n_late_IPQ

  IF (i EQ 0) THEN BEGIN
	PLOT, [xmin,xmax], [0,0], xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Normalized Conformity Signal', LINESTYLE=1, THICK=2
;	  title = 'Field-by-Field Matched SF & Q IP Samples ' + textoidl('(\Deltaz=2.0)')
  ENDIF

	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, normsig[0:plotMax], color=cgcolor(colors[i]), LINESTYLE=ls[i]
	ERRPLOT, Rproj_array+0.5*dRproj+deltaR[i], normsig-normsig_errors, normsig+normsig_errors, color=cgcolor(colors[i])

	legend = [legend, masslabel]; $
;		+ textoidl('; \sigma_{0-1}=') + getSigmaRange_JKE(data, 0, 0) $
;		+ textoidl('; \sigma_{1-2}=') + getSigmaRange_JKE(data, 1, 1) $
;		+ textoidl('; \sigma_{1-3}=') + getSigmaRange_JKE(data, 1, 2) $
;		+ textoidl('; \sigma_{3-5}=') + getSigmaRange_JKE(data, 3, 4)]

  ENDFOR

	LEGEND, legend, LINESTYLE=ls, COLOR=cgColor(colors), BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=1.5

;	XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'SF IP: ' + bigint(n_elements(UNIQ(dataIP[where(dataIP.SFQ EQ 1)].objname))), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Q IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 0)])), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, masslabel, ALIGNMENT=1.0

;	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
