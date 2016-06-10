; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; normsigplot:	plots normalized conformity signal
; JKE:		uses jackknife error method
; IPmatchFBF:	SF and Q IP populations have been selected to have same mass and redshift distribution WITHIN EACH FIELD
; PHI3.7:	IP matching procedure uses high mass cutoff threshold based on Moustakas13 stellar mass function results
; binsEqualThirds: late-type fractions are computed for three redshift ranges with equal # IP

PRO normsigplot_allz_errorCompare, outputFormat
	!P.FONT=0

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

	Rplot_array = Rproj_array[0:n_elements(Rproj_array)-2]+0.5

	xmin = 0
	xmax = 5
	xr = xmax-xmin

	ymin = -0.02
	ymax = 0.09
	yr = ymax-ymin
	
	z_low = zmin
	z_high = zmax
	zsuffix = 'allz'
	zlabel  = 'z = [' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']'

	; read in data files
	dataIPallz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
	dataAllallz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 0.50 THEN string_dR = '_dR500kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

; eliminate data with targ_weight < 1
        dataAll = dataAllallz[where(dataAllallz.targ_weight GE 1.)]
        dataIP  = dataIPallz[where(dataIPallz.targ_weight GE 1.)]

	!P.MULTI=0
	!P.CHARSIZE=1.5

	colors = ['red', 'green', 'red']

	IF (string(outputFormat) EQ 'ps') THEN BEGIN
		PS_OPEN, '~/latex/figures/normsigplot_allz_errorCompare', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=6, YS=5, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

;	PLOT, [xmin,xmax], [0,0], xrange=[xmin,xmax], yrange=[ymin,ymax], LINESTYLE=1, xtitle='Projected Radius (Mpc)', ytitle='Signal (%)', THICK=2
	PLOT, [xmin,xmax], [0,0], xrange=[xmin,xmax], yrange=[ymin,ymax], LINESTYLE=1, xtitle=textoidl('R_{proj} (Mpc)'), ytitle=textoidl('\xi_{norm}'), THICK=2

	data	= MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/normsig_' + zsuffix + '_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_JKE.fits', 1)
	normsig		= data.normsig
	normsig_errors  = data.normsig_errors ; jackknife errors

	dataBSE = MRDFITS('~/conformity/results/match_IP_sample_rigorous/latefrac_allz_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits',1)
	n_tot_IPSF 	= dataBSE.n_tot_IPSF
	n_late_IPSF	= dataBSE.n_late_IPSF
	errors_IPSF	= dataBSE.errors_IPSF
	n_tot_IPQ	= dataBSE.n_tot_IPQ
	n_late_IPQ	= dataBSE.n_late_IPQ
	errors_IPQ	= dataBSE.errors_IPQ

	sigmaBSE  	= ABS(normsig)/SQRT((errors_IPSF)^2+(errors_IPQ)^2)
	errorBSE	= ABS(normsig/sigmaBSE)

	w = 0.1
  FOR i=0,4 DO BEGIN
	POLYFILL, $
	[[0.5+Rproj_array[i]-w/2,normsig[i]-errorBSE[i]],$
	[0.5+Rproj_array[i]+w/2,normsig[i]-errorBSE[i]],$
	[0.5+Rproj_array[i]+w/2,normsig[i]+errorBSE[i]],$
	[0.5+Rproj_array[i]-w/2,normsig[i]+errorBSE[i]]], $
	 COLOR=cgColor('RED3')
  ENDFOR
;	ERRPLOT, Rproj_array[0:plotMax]+0.5*dRproj, normsig-errorBSE, normsig+errorBSE, COLOR=cgColor('RED3'), THICK=16
	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, normsig[0:plotMax]
	ERRPLOT, Rproj_array[0:plotMax]+0.5*dRproj, normsig-normsig_errors, normsig+normsig_errors
	USERSYM, 0.5*[0,1,2,1,1,0,1,2], $
		0.25+1.25*[-1,-1,-1,-1,1,1,1,1], THICK=5
;	LEGEND, ['Bootstrap Errors', 'Jackknife Errors'], PSYM=[8,8], COLOR=[cgColor('RED3'), cgColor('black')], BOX=0, NUMBER=2, /TOP, /RIGHT, PSPACING=1

	XW=0.05
	YW=0.0025
	XO=[4.375,4.625]
	YO=[0.07,0.06]

	XYOUTS, 4.25, YO[0]-YW, 'Bootstrap Errors', ALIGNMENT=1.0
	XYOUTS, 4.25, YO[1]-YW, 'Jackknife Errors', ALIGNMENT=1.0

	POLYFILL, TRANSPOSE([[XO[0]+XW*[0,2,2,0]], $
		[YO[0]+YW*[-1,-1,1,1]]]), COLOR=cgColor('RED3')
	POLYFILL, TRANSPOSE([[XO[1]+XW*[0,2,2,0]], $
		[YO[0]+YW*[-1,-1,1,1]]]), COLOR=cgColor('RED3')
	OPLOT, XO[0]+XW*[0,1,2,1,1,0,1,2], $
		YO[1]+YW*[-1,-1,-1,-1,1,1,1,1], PSYM=0
	OPLOT, XO[1]+XW*[0,1,2,1,1,0,1,2], $
		YO[1]+YW*[-1,-1,-1,-1,1,1,1,1], PSYM=0
	

	XYOUTS, xmin+0.05*xr, ymin+0.05*yr, 'z = [0.20, 1.00]', ALIGNMENT=0.0
	XYOUTS, xmin+0.05*xr, ymin+0.90*yr, 'Matched Sample', ALIGNMENT=0.0

	IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
