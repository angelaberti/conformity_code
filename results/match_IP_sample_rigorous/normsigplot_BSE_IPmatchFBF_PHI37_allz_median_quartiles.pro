; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; normsigplot:	plots normalized conformity signal
; BSE:		uses bootstrap error method (N/A for median, etc.)
; IPmatchFBF:	SF and Q IP populations have been selected to have same mass and redshift distribution WITHIN EACH FIELD
; PHI37:	IP matching procedure uses high mass cutoff threshold based on Moustakas13 stellar mass function results
; binsEqualThirds: late-type fractions are computed for three redshift ranges with equal # IP

PRO normsigplot_BSE_IPmatchFBF_PHI37_allz_median_quartiles, outputFormat
	COMMON SHARE, xmin, xmax, xr, ymin, ymax, yr, Rmax, plotMax, dRproj, n_annuli, Rproj_array, z_low, z_high, zsuffix, zlabel, P50color, P25color, P75color, oneSigColor, meanColor

	P25color = cgColor('Blue')
	P75color = P25color
	oneSigColor = cgColor('Magenta')
	meanColor = cgColor('Red')

	IF (outputFormat EQ 'ps') THEN (P50color = cgColor('Black')) ELSE (P50color = cgColor('White'))
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

	IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
	IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

	xmin = 0
	xmax = n_annuli*dRproj
	xr   = xmax-xmin

	ymin = -0.1
	ymax = 0.15
	yr   = ymax-ymin
	
	z_low   = zmin
	z_high  = zmax
	zsuffix = 'allz'
	zlabel  = 'z=[' + decimal(z_low,2) + '; ' + decimal(z_high,2) + ']'

	; read in data files
	dataIPallz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
;	dataAllallz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

	; eliminate data with targ_weight < 1
;	dataAll = dataAllallz[where(dataAllallz.targ_weight GE 1.)]
        dataIP  = dataIPallz[where(dataIPallz.targ_weight GE 1.)]

	!p.multi=0
	!p.charsize=1.1

	IF (string(outputFormat) eq 'ps') THEN BEGIN
		PS_OPEN, '~/figures/matchedIPsampleFBF/normsigplot_BSE_IPmatchFBF_PHI37_allz_median_quartiles', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
;		THICK=1
	ENDELSE

	PLOT, [xmin,xmax], [0,0], xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Normalized % Difference', $
	  title = 'Field-by-Field Matched SF & Q IP Samples ' + textoidl('(\Deltaz=2.0)'), LINESTYLE=2

	RESTORE, 'IPQ_latefrac_array_dR1Mpc.sav'
	IPQarray  = IP_latefrac_array
	plotQuartiles, IPQarray, P25IPQ, P50IPQ, P75IPQ, meanIPQ, oneSigmaLowerIPQ, oneSigmaUpperIPQ

	RESTORE, 'IPSF_latefrac_array_dR1Mpc.sav'
	IPSFarray = IP_latefrac_array
	plotQuartiles, IPSFarray, P25IPSF, P50IPSF, P75IPSF, meanIPSF, oneSigmaLowerIPSF, oneSigmaUpperIPSF

	ls = 0

	P50array = NORMSIG(P50IPSF, P50IPQ)
	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, P50array[0:plotMax], LINESTYLE=ls, color=P50color
	meanArray = NORMSIG(meanIPSF, meanIPQ)
	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, meanArray[0:plotMax], LINESTYLE=ls, color=meanColor

	P25array = NORMSIG(P25IPSF, P25IPQ)
	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, P25array[0:plotMax], LINESTYLE=ls, color=P25Color
	P75array = NORMSIG(P75IPSF, P75IPQ)
	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, P75array[0:plotMax], LINESTYLE=ls, color=P75color

	oneSigmaLowerArray = NORMSIG(oneSigmaLowerIPSF, oneSigmaLowerIPQ)
	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, oneSigmaLowerArray[0:plotMax], LINESTYLE=ls, color=oneSigColor
	oneSigmaUpperArray = NORMSIG(oneSigmaUpperIPSF, oneSigmaUpperIPQ)
	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, oneSigmaUpperArray[0:plotMax], LINESTYLE=ls, color=oneSigColor


;	LEGEND, ['Late-type IP median', 'Early-type IP median', 'mean', '25 & 75 percentile', textoidl('1\sigma')], LINESTYLE=[0,2,0,0,0], $
;		color=[P50color, P50color, meanColor, P25color, oneSigColor], BOX=0, /BOTTOM, /RIGHT

;	XYOUTS, xmin+0.05*xr, ymin+0.15*yr, 'SF IP: ' + bigint(n_elements(UNIQ(dataIP[where(dataIP.SFQ EQ 1)].objname))), ALIGNMENT=0.0
;	XYOUTS, xmin+0.05*xr, ymin+0.10*yr, 'Q IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 0)])), ALIGNMENT=0.0
	XYOUTS, xmin+0.05*xr, ymin+0.05*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=0.0
;	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Median late-type fraction', ALIGNMENT=1.0

	XYOUTS, xmin+0.95*xr, ymin+0.925*yr, zlabel, ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

PRO plotQuartiles, data, P25array, P50array, P75array, meanArray, oneSigmaLowerArray, oneSigmaUpperArray
	COMMON SHARE

;	[ROW,COLUMN]

	P25array = []
	P50array = []
	P75array = []
	oneSigmaLowerArray = []
	oneSigmaUpperArray = []
	meanArray = []
		
	FOR i=0,n_elements(data[*,0])-1 DO BEGIN
		annulusCol = data[i,*]
		annulusCol = annulusCol[where(annulusCol GE 0)]
		stddev = STDDEV(annulusCol)
		PRINT, stddev

		P50 = MEDIAN(annulusCol)
		mean = MEAN(annulusCol)

		lowerHalf = annulusCol[where(annulusCol LE P50)]
		P25 = MEDIAN(lowerHalf)
		upperHalf = annulusCol[where(annulusCol GE P50)]
		P75 = MEDIAN(upperHalf)

		lowerHalfOrdered = lowerHalf[SORT(lowerHalf)]
		upperHalfOrdered = upperHalf[SORT(upperHalf)]

;		oneSigmaLower = lowerHalfOrdered[CEIL((0.6827/2)*n_elements(lowerHalf))]
;		oneSigmaUpper = upperHalfOrdered[CEIL((2*0.3173)*n_elements(upperHalf))]

		P25array = [P25array, P25]
		P50array = [P50array, P50]
		P75array = [P75array, P75]
		oneSigmaLowerArray = [oneSigmaLowerArray, mean-stddev]
		oneSigmaUpperArray = [oneSigmaUpperArray, mean+stddev]
		meanArray = [meanArray, mean]
	ENDFOR
END

FUNCTION NORMSIG, SFdata, Qdata
	RETURN, (SFdata-Qdata)/((SFdata+Qdata)/2.)
END
