; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; latefracplot:	plots late-type fractions
; BSE:		uses bootstrap error method
; IPmatchFBF:	SF and Q IP populations have been selected to have same mass and redshift distribution WITHIN EACH FIELD
; PHI3.7:	IP matching procedure uses high mass cutoff threshold based on Moustakas13 stellar mass function results
; binsEqualThirds: late-type fractions are computed for three redshift ranges with equal # IP

PRO latefracplot_BSE_IPmatchFBF_PHI37_allz_median_quartiles, outputFormat
	!P.FONT=0

	P25color	= cgColor('gray')
	P75color	= P25color
	meanSFcolor	= cgColor('blue')
	meanQcolor	= cgColor('red')
	medianSFcolor	= cgColor('purple')
	medianQcolor	= cgColor('magenta')

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

	ymin = 0.63
	ymax = 1.0
	yr   = ymax-ymin
	
	z_low   = zmin
	z_high  = zmax
	zsuffix = 'allz'
	zlabel  = 'z = [' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']'

	; read in data files
	dataIPallz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
;	dataAllallz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

	; eliminate data with targ_weight < 1
;	dataAll = dataAllallz[where(dataAllallz.targ_weight GE 1.)]
        dataIP  = dataIPallz[where(dataIPallz.targ_weight GE 1.)]

	!p.multi=0
	!p.charsize=1.5

	IF (string(outputFormat) EQ 'ps') THEN BEGIN
		PS_OPEN, '~/latex/figures/latefracplot_BSE_IPmatchFBF_PHI37_allz_median_quartiles', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=6, YS=5, DECOMPOSED=1, COLOR=1, BITS_PER_PIXEL=8, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

;	PLOT, [1], [1], xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', /NODATA
	PLOT, [1], [1], xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle=textoidl('R_{proj} (Mpc)'), ytitle=textoidl('f_{late}'), /NODATA
;	  title = 'Field-by-Field Matched SF & Q IP Samples ' + textoidl('(\Deltaz=2.0)'), /NODATA

	medianSFls	= 3
	medianQls	= 1
	meanSFls	= 0
	meanQls		= 2

	RESTORE, 'IPQ_latefrac_array_dR1Mpc.sav'
	IPQarray  = IP_latefrac_array

	RESTORE, 'IPSF_latefrac_array_dR1Mpc.sav'
	IPSFarray = IP_latefrac_array

;	[ COL, ROW ]
	allIParray = []
	FOR c=0,n_elements(IPSFarray[*,0])-1 DO BEGIN	
		col = [TRANSPOSE(IPSFarray[c,*]), TRANSPOSE(IPQarray[c,*])]
		allIParray = [[allIParray], [col]]
	ENDFOR

	plotQuartiles, TRANSPOSE(allIParray), P25color, 0, 0, 0, 1, xmin, xmax, xr, ymin, ymax, yr, Rmax, plotMax, dRproj, n_annuli, Rproj_array, z_low, z_high, zsuffix, zlabel

	plotQuartiles, IPQarray, meanQcolor, meanQls, 1, 0, 0, xmin, xmax, xr, ymin, ymax, yr, Rmax, plotMax, dRproj, n_annuli, Rproj_array, z_low, z_high, zsuffix, zlabel
	plotQuartiles, IPQarray, medianQcolor, medianQls, 0, 1, 0, xmin, xmax, xr, ymin, ymax, yr, Rmax, plotMax, dRproj, n_annuli, Rproj_array, z_low, z_high, zsuffix, zlabel

	plotQuartiles, IPSFarray, meanSFcolor, meanSFls, 1, 0, 0, xmin, xmax, xr, ymin, ymax, yr, Rmax, plotMax, dRproj, n_annuli, Rproj_array, z_low, z_high, zsuffix, zlabel
	plotQuartiles, IPSFarray, medianSFcolor, medianSFls, 0, 1, 0, xmin, xmax, xr, ymin, ymax, yr, Rmax, plotMax, dRproj, n_annuli, Rproj_array, z_low, z_high, zsuffix, zlabel
	
	LEGEND, ['SF IP median', 'Q IP median', 'SF IP mean', 'Q IP mean'], $
		LINESTYLE=[medianSFls, medianQls, meanSFls, meanQls], $
		COLOR=[medianSFcolor, medianQcolor, meanSFcolor, meanQcolor], BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3

;	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'SF IP: ' + bigint(n_elements(UNIQ(dataIP[where(dataIP.SFQ EQ 1)].objname))), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Q IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 0)])), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Median late-type fraction', ALIGNMENT=1.0

	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, zlabel, ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0

	IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

PRO plotQuartiles, data, col, ls, av, med, quartiles, xmin, xmax, xr, ymin, ymax, yr, Rmax, plotMax, dRproj, n_annuli, Rproj_array, z_low, z_high, zsuffix, zlabel
;	P25color	= cgColor('gray')
;	P75color	= P25color
;	meanSFcolor	= cgColor('blue')
;	meanQcolor	= cgColor('red')
;	medianSFcolor	= cgColor('purple')
;	medianQcolor	= cgColor('magenta')

;	[ROW,COLUMN]

	P25array = []
	P50array = []
	P75array = []
;	oneSigmaLowerArray = []
;	oneSigmaUpperArray = []
	meanArray = []
		
	FOR i=0,n_elements(data[*,0])-1 DO BEGIN
		annulusCol = data[i,*]
		annulusCol = annulusCol[where(annulusCol GE 0)]
		stddev = STDDEV(annulusCol)
;		PRINT, stddev

		P50 = MEDIAN(annulusCol)
		mean = MEAN(annulusCol)

		lowerHalf = annulusCol[where(annulusCol LE P50)]
		P25 = MEDIAN(lowerHalf)
		upperHalf = annulusCol[where(annulusCol GE P50)]
		P75 = MEDIAN(upperHalf)

		lowerHalfOrdered = lowerHalf[SORT(lowerHalf)]
		upperHalfOrdered = upperHalf[SORT(upperHalf)]

		P25array = [P25array, P25]
		P50array = [P50array, P50]
		P75array = [P75array, P75]
		meanArray = [meanArray, mean]
	ENDFOR

	IF (med EQ 1) THEN $
	OPLOT, Rproj_array[1:plotMax]+0.5*dRproj, P50array[1:plotMax], LINESTYLE=ls, color=col
	IF (av EQ 1) THEN $
	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, meanArray[0:plotMax], LINESTYLE=ls, color=col

	IF (quartiles EQ 1) THEN BEGIN $
	  POLYFILL, [REVERSE(Rproj_array[1:plotMax]+0.5*dRproj), Rproj_array[1:plotMax]+0.5*dRproj], $
	  [REVERSE(P25array[1:plotMax]), P75array[1:plotMax]], COLOR=cgColor('lightgray')
	ENDIF
END
