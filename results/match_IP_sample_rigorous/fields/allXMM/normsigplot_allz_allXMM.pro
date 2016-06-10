; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; latefrac:	calculates late-type fractions
; BSE:		uses bootstrap error method
; IPmatchFBF:	SF and Q IP populations have been selected to have same mass and redshift distribution WITHIN EACH FIELD
; PHI3.7:	IP matching procedure uses high mass cutoff threshold based on Moustakas13 stellar mass function results
; binsEqualThirds: late-type fractions are computed for three redshift ranges with equal # IP

PRO normsigplot_allz_allXMM, outputFormat
	PHI = 3.7
	stringPHI = strtrim(string(PHI, format='(f20.1)'),2)

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	Rmax 		= 5.
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
	ymax = 0.09
	yr = ymax-ymin

	z_low = 0.2
	z_high = 1.0

	zsuffix = 'allz'
	zlabel = 'z=[' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']'

	!P.FONT=0

	IF (string(outputFormat) EQ 'ps') THEN BEGIN
	        PS_OPEN, '~/latex/figures/normsigplot_byField', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=6, YS=5, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	fields = dataAllallz[uniq(dataAllallz.field, sort(dataAllallz.field))].field
;	PRINT, fields

;	COLORS	= ['dark_red','magenta','purple','blue','dark_green']
	COLORS	= ['ORG5', 'PUR5', 'YGB5', 'GRN5', 'BLK5']

	PSYMS	= [-4,-5,-6,-1,-2]
;	PSYMS	= [0,3,4,5,8]
;	PSYMS	= [cgsymcat(14), $
;		cgsymcat(15), $
;		cgsymcat(16), $
;		cgsymcat(17), $
;		cgsymcat(46)]

;	PSYMS = [14,15,16,17,46]

;	PRINT, PSYMS
	legSyms = []

  FOR f=0,0 DO BEGIN
;	IF (fields[f] EQ 'cdfs      ') THEN (stringField = 'CDFS')
;	IF (fields[f] EQ 'cfhtls_xmm') THEN (stringField = 'CFHTLS-XMM')
;	IF (fields[f] EQ 'cosmos    ') THEN (stringField = 'COSMOS')
;	IF (fields[f] EQ 'es1       ') THEN (stringField = 'ES1')
;	IF (fields[f] EQ 'xmm       ') THEN (stringField = 'XMM')

	stringField = 'All XMM'

;	datafile = '~/conformity/results/match_IP_sample_rigorous/fields/latefrac_allz_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_BSE_halo12_' + stringField + '.fits'
	datafile = '~/conformity/results/match_IP_sample_rigorous/fields/allXMM/latefrac_allz_allXMM_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_BSE.fits'
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

	IF (f EQ 0) THEN BEGIN
		PLOT, [xmin,xmax], [0,0], xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Signal (%)', LINESTYLE=1, THICK=2
	ENDIF

	PP = cgsymcat(PSYMS[f])

	OPLOT, Rplot_array+0.05*(2-f), normsig[0:n_elements(Rplot_array)-1], PSYM=0, COLOR=cgColor(COLORS[f])
;	OPLOT, Rplot_array, normsig[0:n_elements(Rplot_array)-1], PSYM=CGSYMCAT(PSYMS[f]), COLOR=cgColor(COLORS[f]), SYMSIZE=2
	OPLOT, Rplot_array+0.05*(2-f), normsig[0:n_elements(Rplot_array)-1], PSYM=PSYMS[f], COLOR=cgColor(COLORS[f]), SYMSIZE=1
	ERRPLOT, Rplot_array+0.05*(2-f), normsig+error, normsig-error, COLOR=cgColor(COLORS[f])
  ENDFOR

;	LEGEND, ['CDFS', 'CFHTLS-XMM', 'COSMOS', 'ES1', 'XMM'], PSYM=PSYMS, COLOR=cgColor(COLORS), BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3, CHARSIZE=charsz
;	AL_LEGEND, ['CDFS', 'CFHTLS-XMM', 'COSMOS', 'ES1', 'XMM'], PSYM=PSYMS, BOX=0, /TOP, /RIGHT, SYMSIZE=2

	IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'

END
