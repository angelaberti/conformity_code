; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; latefracplot:	plots late-type fractions
; BSE:		uses bootstrap error method
; IPmatchFBF:	SF and Q IP populations have been selected to have same mass and redshift distribution WITHIN EACH FIELD
; PHI3.7:	IP matching procedure uses high mass cutoff threshold based on Moustakas13 stellar mass function results
; binsEqualThirds: late-type fractions are computed for three redshift ranges with equal # IP

PRO latefracplot_BSE_IPmatchFBF_COSMOS_PHI37_allz, COSMOScomp, outputFormat;, zmin, zmax;, dz_coeff, printEvery
	!P.FONT=0

	COMMON SHARE, xmin, xmax, xr, ymin, ymax, yr, Rmax, plotMax, dRproj, n_annuli, Rproj_array, z_low, z_high, zsuffix, zlabel
	SFcolor = cgcolor('blue')
	Qcolor  = cgcolor('red')
	COScolor = 'red'
	IF (outputFormat EQ 'ps') THEN (REGcolor = 'black') ELSE (REGcolor = 'white')
	PHI = 3.7
	stringPHI = strtrim(string(PHI, format='(f20.1)'),2)

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	Rmax 		= 15.
	plotMax 	= Rmax-1
	dRproj		= 1.
	n_annuli 	= Ceil(Rmax/dRproj)
	Rproj_array 	= FLOAT(dRproj)*(findgen(n_annuli+1))	

	xmin = 0
	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = 0.71
	ymax = 0.88
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
        dataIP  = dataIPallz[where(dataIPallz.targ_weight GE 1.)]

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
	  IF (COSMOScomp EQ 1) THEN $
	        PS_OPEN, '~/figures/matchedIPsampleFBF/latefracplot_COSMOScomp_BSE_IPmatchFBF_PHI37_allz', THICK=5, /ENCAP $
;	  ELSE PS_OPEN, '~/figures/matchedIPsampleFBF/latefracplot_BSE_IPmatchFBF_PHI37_allz', THICK=5, /ENCAP
	  ELSE PS_OPEN, '~/latex/figures/latefracplot_BSE_IPmatchFBF_COSMOS_PHI37_allz', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	PLOT, [1], [1], xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', /NODATA ;, $ 
;	  title = 'Field-by-Field Matched SF & Q IP Samples ' + textoidl('(\Deltaz=2.0)')

	;latefrac_allz_cosmos_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits
	data = MRDFITS('~/conformity/results/match_IP_sample_rigorous/fields/latefrac_' + zsuffix + '_cosmos_targ_weight_IPmatchFBF_PHI' + stringPHI + string_dR + '_BSE.fits', 1)
	makePlots, data, cgcolor(REGcolor), 0.10
  IF (COSMOScomp EQ 1) THEN BEGIN
	data = MRDFITS('~/conformity/results/match_IP_sample_rigorous/noCosmos/latefrac_' + zsuffix + '_targ_weight_IPmatchFBFnoCosmos_PHI' + stringPHI + string_dR + '_BSE.fits', 1)
	makePlots, data, cgColor(COScolor), 0.05
  ENDIF

  IF (COSMOScomp EQ 1) THEN LEGEND, ['Late-type IP', 'Early-type IP', 'COSMOS field excluded'], LINESTYLE=[0,2,0], colors=cgcolor([REGcolor,REGcolor,COScolor]), BOX=0, /TOP, /RIGHT
  IF (COSMOScomp NE 1) THEN BEGIN
	LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], COLOR=[cgcolor('blue'), cgcolor('red')], BOX=0, /TOP, /RIGHT, NUMBER=2, PSPACING=3
;	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'SF IP: ' + bigint(n_elements(UNIQ(dataIP[where(dataIP.SFQ EQ 1)].objname))), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Q IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 0)])), ALIGNMENT=1.0
  ENDIF
;	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, zlabel, ALIGNMENT=1.0
	XYOUTS, xmin+0.05*xr, ymin+0.05*yr, 'COSMOS field only', ALIGNMENT=0.0
;	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

PRO makePlots, data, color, sigmaHeight
	COMMON SHARE

	SFcolor = cgcolor('blue')
	Qcolor  = cgcolor('red')

;	dataIP  = dataIPallz
;	dataAll = dataAllallz

;	frac_IPSF = data.IPSF_latefrac_median_array
;	frac_IPQ  = data.IPQ_latefrac_median_array

	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPSF = data.n_late_IPSF
	errors_IPSF = data.errors_IPSF
	n_tot_IPQ  = data.n_tot_IPQ
	n_late_IPQ = data.n_late_IPQ
	errors_IPQ = data.errors_IPQ

	frac_IPSF   = n_late_IPSF/n_tot_IPSF
	frac_IPQ   = n_late_IPQ/n_tot_IPQ

	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, frac_IPSF[0:plotMax], LINESTYLE=ls, color=SFcolor
	ERRPLOT, Rproj_array+0.5*dRproj, frac_IPSF-errors_IPSF, frac_IPSF+errors_IPSF, color=SFcolor

	OPLOT, Rproj_array[0:plotMax]+0.5*dRproj, frac_IPQ[0:plotMax], LINESTYLE=2, color=Qcolor
	ERRPLOT, Rproj_array+0.5*dRproj, frac_IPQ-errors_IPQ, frac_IPQ+errors_IPQ, color=Qcolor

;	sigmas = textoidl('\sigma_{0-1}=') + getSigmaRange(data, 0, 0) $
;		+ textoidl('; \sigma_{1-2}=') + getSigmaRange(data, 1, 1) $
;		+ textoidl('; \sigma_{1-3}=') + getSigmaRange(data, 1, 2) $
;		+ textoidl('; \sigma_{3-5}=') + getSigmaRange(data, 3, 4)

;	XYOUTS, xmin+0.95*xr, ymin+sigmaHeight*yr, sigmas, ALIGNMENT=1.0, color=color, charsize=1.5
END
