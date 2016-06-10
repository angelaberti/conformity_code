; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; latefrac:	calculates late-type fractions
; BSE:		uses bootstrap error method
; IPmatchFBF:	SF and Q IP populations have been selected to have same mass and redshift distribution WITHIN EACH FIELD
; PHI3.7:	IP matching procedure uses high mass cutoff threshold based on Moustakas13 stellar mass function results
; binsEqualThirds: late-type fractions are computed for three redshift ranges with equal # IP

PRO normsigplot_allz_byField_1halo, outputFormat
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
	dataIPallz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
	dataAllallz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

; eliminate data with targ_weight < 1
        dataAllallz = dataAllallz[where(dataAllallz.targ_weight GE 1.)]
        dataIPallz  = dataIPallz[where(dataIPallz.targ_weight GE 1.)]

	!P.MULTI=0
	!P.CHARSIZE=1.5
	charsz=1.25

	Rproj_array = FINDGEN(Rmax+1)
	Rplot_array = FINDGEN(Rmax) + 0.5*dRproj
;	Rproj_array = [0.0, 0.5, 1.0+FINDGEN(Rmax)]
;	Rplot_array = [0.25, 0.75, 1.5+FINDGEN(Rmax-1)]

	PRINT, 'Rproj_array: ', Rproj_array
	PRINT, 'Rplot_array: ', Rplot_array

	xmin = 0
	xmax = 5
	xr = xmax-xmin

	ymin = -0.02
	ymax = 0.15
	yr = ymax-ymin

	z_low = 0.2
	z_high = 1.0

	zsuffix = 'allz'
	zlabel = 'z=[' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']'

	!P.FONT=0

	IF (string(outputFormat) EQ 'ps') THEN BEGIN
	        PS_OPEN, '~/latex/figures/normsigplot_byField_1halo', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=6, YS=5, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	fields = dataAllallz[uniq(dataAllallz.field, sort(dataAllallz.field))].field

	COLORS	= ['ORG5', 'PUR5', 'YGB5', 'GRN5', 'BLK5']

	PSYMS	= [-4,-5,-6,-1,-2]

  FOR f=0,n_elements(fields)-1 DO BEGIN
	IF (fields[f] EQ 'cdfs      ') THEN (stringField = 'CDFS')
	IF (fields[f] EQ 'cfhtls_xmm') THEN (stringField = 'XMM-CFHTLS')
	IF (fields[f] EQ 'cosmos    ') THEN (stringField = 'COSMOS')
	IF (fields[f] EQ 'es1       ') THEN (stringField = 'ES1')
	IF (fields[f] EQ 'xmm       ') THEN (stringField = 'XMM-SXDS')

	datafile = '~/conformity/results/match_IP_sample_rigorous/fields/latefrac_allz_' + STRTRIM(fields[f],2) + '_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_BSE.fits'
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

	PRINT, fields[f]
	PRINT, 'signal ', normsig[0:4]
	PRINT, 'error  ', error[0:4]
	PRINT, 'sigma  ', sigma[0:4]

	normsigAll = 0.053
	errorAll = 0.015
	IF (f EQ 0) THEN BEGIN
;		PLOT, [0,6.25], [0,0], xrange=[xmin,6.25], yrange=[ymin,ymax], xtickformat='(A1)', ytitle='Signal (%)', LINESTYLE=1, THICK=2, XMINOR=1, XTICKS=1
		PLOT, [0,6.25], [0,0], xrange=[xmin,6.25], yrange=[ymin,ymax], xtickformat='(A1)', ytitle=textoidl('\xi_{norm}'), LINESTYLE=1, THICK=2, XMINOR=1, XTICKS=1
		OPLOT,  [5.5], [normsigAll], PSYM=6, SYMSIZE=1, THICK=10
		ERRPLOT, [5.5], [normsigAll+errorAll], [normsigAll-errorAll]
		XYOUTS, 5.5, normsigAll+errorAll+0.005, 'All Fields', ALIGNMENT=0.5, CHARSIZE=0.9*charsz
	ENDIF

	OPLOT, [f]+0.5, normsig[0:n_elements(Rplot_array)-1], PSYM=6, COLOR=cgColor(COLORS[f]), SYMSIZE=1, THICK=10
	ERRPLOT, [f]+0.5, normsig+error, normsig-error, COLOR=cgColor(COLORS[f])
	IF (f EQ 0) OR (f EQ 2) OR (f EQ 3) THEN $
		XYOUTS, f+0.5, normsig[0]-(error[0]+0.01), stringField, ALIGNMENT=0.5, CHARSIZE=0.9*charsz $
	ELSE XYOUTS, f+0.5, normsig[0]+error[0]+0.005, stringField, ALIGNMENT=0.5, CHARSIZE=0.9*charsz
  ENDFOR
	XYOUTS, xmin+0.05*xr, ymin+0.90*yr, textoidl('0 < R_{proj} < 1 Mpc'), ALIGNMENT=0.0
	XYOUTS, xmin+0.05*xr, ymin+0.05*yr, 'z = [0.20, 1.00]', ALIGNMENT=0.0

	IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
