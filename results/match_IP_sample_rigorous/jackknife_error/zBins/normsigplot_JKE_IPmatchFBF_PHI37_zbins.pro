; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; normsigplot:	plots normalized conformity signal
; JKE:		uses jackknife error method
; IPmatchFBF:	SF and Q IP populations have been selected to have same mass and redshift distribution WITHIN EACH FIELD
; PHI37:	IP matching procedure uses high mass cutoff threshold based on Moustakas13 stellar mass function results

PRO normsigplot_JKE_IPmatchFBF_PHI37_zbins, n_bins, outputFormat
	!P.FONT=0

	PHI = 3.7
	stringPHI = strtrim(string(PHI, format='(f20.1)'),2)

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	; read in data files
;	IPdataAllz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)
	IPdataAllz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
	dataAllallz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

	Rmax 		= 5.
	dRproj		= 1.
	n_annuli 	= CEIL(Rmax/dRproj)
	plotMax 	= Rmax-1

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

; eliminate data with targ_weight < 1
        dataAllallz = dataAllallz[where(dataAllallz.targ_weight GE 1.)]
        IPdataAllz  = IPdataAllz[where(IPdataAllz.targ_weight GE 1.)]

	!P.MULTI=0
	!P.CHARSIZE=1.75

	Rproj_array = FINDGEN(Rmax+1)
	Rplot_array = FINDGEN(Rmax) + 0.5*dRproj
;	Rproj_array = [0.0, 0.5, 1.0+FINDGEN(Rmax)]
;	Rplot_array = [0.25, 0.75, 1.5+FINDGEN(Rmax-1)]
;	n_annuli = n_elements(Rplot_array)	

	xmin = 0
	xmax = 5
	xr = xmax-xmin

	ymin = -0.04
	ymax = 0.11
	yr = ymax-ymin

	colors	= ['purple', 'magenta', 'orange']
	ls	= [0,2,3]
	legend	= []

	deltaR	= [-0.05, 0, 0.05]

	IF (n_bins EQ 2) THEN BEGIN
		zArray = [0.2, MEDIAN(IPdataAllz.zprimus), 1.0]
		zSuffixArray = ['H1', 'H2']
		figureSuffix = 'zbinsEqualHalves'

		colors	= colors[0:n_elements(colors)-2]
		ls	= ls[0:n_elements(ls)-2]

		deltaR	= deltaR[1:n_elements(deltaR)-1]
	ENDIF ELSE BEGIN ; n_bins=3
		orderedRedshifts = IPdataAllz[SORT(IPdataAllz.zprimus)].zprimus

		lower_index = CEIL(n_elements(IPdataAllz)/3.)
		upper_index = 2*FLOOR(n_elements(IPdataAllz)/3.)

		zArray = [0.2, orderedRedshifts[lower_index], orderedRedshifts[upper_index], 1.0]
		zSuffixArray = ['T1', 'T2', 'T3']
		figureSuffix = 'zbinsEqualThirds'
	ENDELSE

	PRINT, 'N bins: ', n_bins
	PRINT, 'z array: ', zArray

	IF (string(outputFormat) EQ 'ps') THEN BEGIN
		PS_OPEN, '~/latex/figures/normsigplot_JKE_IPmatchFBF_PHI37_' + figureSuffix, THICK=5, /ENCAP
		DEVICE, /INCH, XS=8, YS=6, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
		SET_PLOT, 'X'
		THICK=1
	ENDELSE

  FOR n=0,n_elements(zArray)-2 DO BEGIN
	z_low   = zArray[n]
	z_high  = zArray[n+1]

	IF (n EQ 0) THEN ( zlabel = 'z=[' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']' ) $
	  ELSE ( zlabel = '[' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']' )

	dataIP  = IPdataAllz[where( (IPdataAllz.zprimus GE z_low) AND (IPdataAllz.zprimus LE z_high) )]
	dataAll = dataAllallz[where( (dataAllallz.zprimus GE z_low) AND (dataAllallz.zprimus LE z_high) )]

	dataIPSF = dataIP[where(dataIP.SFQ EQ 1)]
	dataIPQ  = dataIP[where(dataIP.SFQ EQ 0)]

	PRINT, zlabel
	PRINT, 'SF IP: ', n_elements(UNIQ(dataIPSF.objname))
	PRINT, 'Q IP:  ', n_elements(dataIPQ)

	data = MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/zBins/normsig_' + zSuffixArray[n] + '_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_JKE.fits', 1)

	normsig		= data.normsig
	normsig_errors	= data.normsig_errors
	
	n_tot_IPSF  	= data.n_tot_IPSF
	n_late_IPSF	= data.n_late_IPSF
	n_tot_IPQ  	= data.n_tot_IPQ
	n_late_IPQ	= data.n_late_IPQ

  IF (n EQ 0) THEN BEGIN
	PLOT, [xmin,xmax], [0,0], xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Normalized Conformity Signal', LINESTYLE=1, THICK=2
;	  title = 'Field-by-Field Matched SF & Q IP Samples ' + textoidl('(\Deltaz=2.0)'), LINESTYLE=2
  ENDIF

	OPLOT, Rplot_array[0:plotMax], normsig[0:plotMax], color=cgcolor(colors[n]), LINESTYLE=ls[n]
	ERRPLOT, Rplot_array+deltaR[n], normsig-normsig_errors, normsig+normsig_errors, color=cgcolor(colors[n])

	legend = [legend, zlabel]; $
;	legend = [legend, zlabel, $
;		+ textoidl('; \sigma_{0-1}=') + getSigmaRange_JKE(data, 0, 0) $
;		+ textoidl('; \sigma_{1-2}=') + getSigmaRange_JKE(data, 1, 1) $
;		+ textoidl('; \sigma_{1-3}=') + getSigmaRange_JKE(data, 1, 2) $
;		+ textoidl('; \sigma_{3-5}=') + getSigmaRange_JKE(data, 3, 4)]

  ENDFOR

	LEGEND, legend, LINESTYLE=ls, COLOR=cgColor(colors), BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=1.5
;	LEGEND, legend;, LINESTYLE=ls, COLOR=cgColor(colors), BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=1.5

;	XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'SF IP: ' + bigint(n_elements(UNIQ(dataIP[where(dataIP.SFQ EQ 1)].objname))), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Q IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 0)])), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, zlabel, ALIGNMENT=1.0

;	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
