PRO normsigplot_zbins_1halo, outputFormat
	PHI = 3.7
	stringPHI = strtrim(string(PHI, format='(f20.1)'),2)

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	Rmax 		= 1.
	printEvery	= 100

	dRproj=1.
	
        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 0.50 THEN string_dR = '_dR500kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

; Matched M* and Redshift
	dataIPallz  = MRDFITS('~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
	dataAllallz = MRDFITS('~/results/zerodSFQ_all_cart.fits', 1)

; eliminate data with targ_weight < 1
        dataAllallz = dataAllallz[where(dataAllallz.targ_weight GE 1.)]
        dataIPallz  = dataIPallz[where(dataIPallz.targ_weight GE 1.)]

	!P.MULTI=0
	!P.CHARSIZE=1.25
	charsz=1

	Rproj_array = FINDGEN(Rmax+1)
	Rplot_array = FINDGEN(Rmax) + 0.5*dRproj
;	Rproj_array = [0.0, 0.5, 1.0+FINDGEN(Rmax)]
;	Rplot_array = [0.25, 0.75, 1.5+FINDGEN(Rmax-1)]

	PRINT, 'Rproj_array: ', Rproj_array
	PRINT, 'Rplot_array: ', Rplot_array

	xmin = 0
	xmax = 5
	xr = xmax-xmin

	ymin = -0.18
	ymax = 0.12
	yr = ymax-ymin

	z_low = 0.2
	z_high = 1.0

	zsuffix = 'allz'
	zlabel = 'z=[' + decimal(z_low,1) + ', ' + decimal(z_high,1) + ']'

	!P.FONT=0

	IF (string(outputFormat) EQ 'ps') THEN BEGIN
	        PS_OPEN, '~/latex/figures/normsigplot_XMM-SXDS_zbins_1halo', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=6, YS=5;, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	fields = dataAllallz[uniq(dataAllallz.field, sort(dataAllallz.field))].field

	COLORS	= ['ORG5', 'PUR5', 'YGB5', 'GRN5', 'BLK5']

	PSYMS	= [-4,-5,-6,-1,-2]

  FOR h=0,1 DO BEGIN
	datafile = '~/results/match_IP_sample_rigorous/fields/XMM-SXDS/latefrac_H'+STRTRIM(h,2)+'_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_BSE.fits'
	data = MRDFITS(datafile, 1)

	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPSF = data.n_late_IPSF
	errors_IPSF = data.errors_IPSF
	frac_IPSF   = n_late_IPSF/n_tot_IPSF

	n_tot_IPQ  = data.n_tot_IPQ
	n_late_IPQ = data.n_late_IPQ
	errors_IPQ = data.errors_IPQ
	frac_IPQ   = n_late_IPQ/n_tot_IPQ

	normsig = (frac_IPSF-frac_IPQ)/((frac_IPSF+frac_IPQ)/2)
	sigma	= ABS(normsig)/SQRT((errors_IPSF)^2+(errors_IPQ)^2)
	error	= ABS(normsig/sigma)

	PRINT, 'signal ', normsig[0:4]
	PRINT, 'error  ', error[0:4]
	PRINT, 'sigma  ', sigma[0:4]

	IF (h EQ 0) THEN BEGIN
		PLOT, [0,4], [0,0], xrange=[0,4], yrange=[ymin,ymax], xtickformat='(A1)', ytitle='Signal (%)', LINESTYLE=1, THICK=2, XMINOR=1, XTICKS=1
	ENDIF

	OPLOT, [2*h]+1, normsig[0:n_elements(Rplot_array)-1], PSYM=6, COLOR=cgColor(COLORS[0]), SYMSIZE=1, THICK=10
	ERRPLOT, [2*h]+1, normsig+error, normsig-error, COLOR=cgColor(COLORS[0])
  IF (h EQ 0) THEN $
	XYOUTS, [2*h]+1, normsig[0:n_elements(Rplot_array)-1]-error-0.01, textoidl('0.2<z<z_{med}'), ALIGNMENT=0.5, CHARSIZE=charsz $
  ELSE $
	XYOUTS, [2*h]+1, normsig[0:n_elements(Rplot_array)-1]-error-0.01, textoidl('z_{med}<z<1.0'), ALIGNMENT=0.5, CHARSIZE=charsz 
  ENDFOR

  FOR q=0,3 DO BEGIN
	datafile = '~/results/match_IP_sample_rigorous/fields/XMM-SXDS/latefrac_Q'+STRTRIM(q,2)+'_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_BSE.fits'
	data = MRDFITS(datafile, 1)

	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPSF = data.n_late_IPSF
	errors_IPSF = data.errors_IPSF
	frac_IPSF   = n_late_IPSF/n_tot_IPSF

	n_tot_IPQ  = data.n_tot_IPQ
	n_late_IPQ = data.n_late_IPQ
	errors_IPQ = data.errors_IPQ
	frac_IPQ   = n_late_IPQ/n_tot_IPQ

	normsig = (frac_IPSF-frac_IPQ)/((frac_IPSF+frac_IPQ)/2)
	sigma	= ABS(normsig)/SQRT((errors_IPSF)^2+(errors_IPQ)^2)
	error	= ABS(normsig/sigma)

	OPLOT, [q]+0.5, normsig[0:n_elements(Rplot_array)-1], PSYM=6, COLOR=cgColor(COLORS[1]), SYMSIZE=1, THICK=10
	ERRPLOT, [q]+0.5, normsig+error, normsig-error, COLOR=cgColor(COLORS[1])

	zlow  = decimal(0.2+0.2*q,1)
	zhigh = decimal(0.2+0.2*(q+1),1)

	XYOUTS, [q]+0.5, normsig[0:n_elements(Rplot_array)-1]+error+0.005, zlow+'<z<'+zhigh, ALIGNMENT=0.5, CHARSIZE=charsz
  ENDFOR

	XYOUTS, xmin+0.05*xr, ymin+0.10*yr, 'XMM-SXDS', ALIGNMENT=0.0
	XYOUTS, xmin+0.05*xr, ymin+0.05*yr, textoidl('0 < R_{proj} < 1 Mpc'), ALIGNMENT=0.0

	IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
