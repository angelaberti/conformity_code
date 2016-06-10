; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; latefracplot:	plots late-type fractions
; BSE:		uses bootstrap error method
; IPmatchFBF:	SF and Q IP populations have been selected to have same mass and redshift distribution WITHIN EACH FIELD
; PHI3.7:	IP matching procedure uses high mass cutoff threshold based on Moustakas13 stellar mass function results
; binsEqualThirds: late-type fractions are computed for three redshift ranges with equal # IP

PRO latefracplot_BSE_IPmatchFBFnoCosmos_PHI37_allz, outputFormat;, zmin, zmax;, dz_coeff, printEvery
	PHI = 3.7
	stringPHI = strtrim(string(PHI, format='(f20.1)'),2)

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	; read in data files
;	dataIPallz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)
	dataIPallz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
	dataAllallz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

	Rmax 		= 5.
	plotMax 	= Rmax-1
	dRproj		= 1.
	n_annuli 	= Ceil(Rmax/dRproj)
	printEvery	= 100

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 0.50 THEN string_dR = '_dR500kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

	; eliminate data with targ_weight < 1
        dataAllallz = dataAllallz[where(dataAllallz.targ_weight GE 1.)]
        dataIPallz  = dataIPallz[where(dataIPallz.targ_weight GE 1.)]

	; divide IP population into 3 redshift bins with equal numbers of IPs
        orderedRedshifts = dataIPallz[SORT(dataIPallz.zprimus)].zprimus

        lower_index = CEIL(n_elements(dataIPallz)/3.)
        upper_index = 2*FLOOR(n_elements(dataIPallz)/3.)

        lower_z = orderedRedshifts[lower_index]
        upper_z = orderedRedshifts[upper_index]

        zArray = [0.2, lower_z, upper_z, 1.0]
        PRINT, stringPHI, ' ', zArray

	!p.multi=0
	!p.charsize=1.25

	Rproj_array = float(dRproj)*(findgen(n_annuli+1))	

	xmin = 0
	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = 0.65
	ymax = 0.85
	yr = ymax-ymin
	
	colors = ['blue', 'green', 'red']

	legend = []

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsampleFBF/latefracplot_BSE_IPmatchFBFnoCosmos_PHI3.7_allz', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

;  FOR i=0,2 DO BEGIN
;	z_low   = zArray[i]
	z_low = 0.2
;	z_high  = zArray[i+1]
	z_high = 1.0

;	zsuffix = 'T' + strtrim(i+1,1)
	zsuffix = 'allz'
	zlabel  = 'z=[' + decimal(z_low,2) + '; ' + decimal(z_high,2) + ']'

;	dataIP  = dataIPallz[where( (dataIPallz.zprimus GE z_low) AND (dataIPallz.zprimus LE z_high) )]
;	dataAll = dataAllallz[where( (dataAllallz.zprimus GE z_low) AND (dataAllallz.zprimus LE z_high) )]
	dataIP  = dataIPallz
	dataAll = dataAllallz
	dataIP  = dataIP[where(dataIP.field NE 'cosmos    ')]
	dataAll = dataAll[where(dataAll.field NE 'cosmos    ')]

	data = MRDFITS('~/conformity/results/match_IP_sample_rigorous/noCosmos/latefrac_' + zsuffix + '_targ_weight_IPmatchFBFnoCosmos_PHI' + stringPHI + string_dR + '_BSE.fits', 1)

	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPSF = data.n_late_IPSF
	errors_IPSF = data.errors_IPSF
	n_tot_IPQ  = data.n_tot_IPQ
	n_late_IPQ = data.n_late_IPQ
	errors_IPQ = data.errors_IPQ

	frac_IPSF   = n_late_IPSF/n_tot_IPSF
	frac_IPQ   = n_late_IPQ/n_tot_IPQ
;  IF (i EQ 0) THEN BEGIN
	PLOT, [1], [1], xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', $
	  title = 'Field-by-Field Matched SF & Q IP Samples ' + textoidl('(\Deltaz=2.0)'), /NODATA
	LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], BOX=0, /TOP, /RIGHT
;  ENDIF
	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, frac_IPSF[0:plotMax], LINESTYLE=0;, color=cgcolor(colors[i])
	ERRPLOT, Rproj_array+0.5*dRproj, frac_IPSF-errors_IPSF, frac_IPSF+errors_IPSF;, color=cgcolor(colors[i])

	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, frac_IPQ[0:plotMax], LINESTYLE=2;, color=cgcolor(colors[i])
	ERRPLOT, Rproj_array+0.5*dRproj, frac_IPQ-errors_IPQ, frac_IPQ+errors_IPQ;, color=cgcolor(colors[i])

;	legend = [legend, zlabel $
	sigmas = textoidl('\sigma_{0-1}=') + getSigmaRange(data, 0, 0) $
		+ textoidl('; \sigma_{1-2}=') + getSigmaRange(data, 1, 1) $
		+ textoidl('; \sigma_{1-3}=') + getSigmaRange(data, 1, 2) $
		+ textoidl('; \sigma_{3-5}=') + getSigmaRange(data, 3, 4)
;  ENDFOR

;	LEGEND, legend, LINESTYLE=[0,0,0], COLOR=cgColor(colors), BOX=0, /BOTTOM, /RIGHT

	XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'SF IP: ' + bigint(n_elements(UNIQ(dataIP[where(dataIP.SFQ EQ 1)].objname))), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Q IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 0)])), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
	XYOUTS, xmin+0.05*xr, ymin+0.875*yr, zlabel, ALIGNMENT=0.0

	XYOUTS, xmin+0.05*xr, ymin+0.925*yr, 'COSMOS field excluded', ALIGNMENT=0.0
	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, sigmas, ALIGNMENT=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
