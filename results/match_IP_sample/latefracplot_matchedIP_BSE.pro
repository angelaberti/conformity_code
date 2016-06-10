; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

PRO latefracplot_matchedIP_BSE, outputFormat;, zmin, zmax;, dz_coeff, printEvery
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/latefracplot_matchedIP_zall_BSE', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

;	dataIP  = MRDFITS('~/conformity/results/match_IP_sample/matchedIPsample.fits', 1)
;	dataAll = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

	Rmax 		= 15.
	dRproj		= 0.25
	n_annuli 	= Ceil(Rmax/dRproj)
	printEvery	= 100

;	dataAll = dataAll[where(dataAll.targ_weight GE 1.)]
;	dataIP  = dataIP[where(dataIP.targ_weight GE 1.)]

	!p.multi=0
	!p.charsize=1.1

	Rproj_array = float(dRproj)*(findgen(n_annuli+1)+0.5)

	xmin = 0
	xmax = 5	
;	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = 0.5
	ymax = 1
	yr = ymax-ymin

	dataIP	= MRDFITS('~/conformity/results/match_IP_sample/matchedIPsample.fits',1)
	data	= MRDFITS('~/conformity/results/match_IP_sample/latefrac_' + strtrim(strcompress(string(zmin, format='(f20.1)')), 1) + '_' $
		+ strtrim(strcompress(string(zmax, format='(f20.1)')) ,1) + '_targ_weight_matchedIPsample_BSE.fits',1)

	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPSF = data.n_late_IPSF
	errors_IPSF = data.errors_IPSF
	frac_IPSF   = n_late_IPSF/n_tot_IPSF

	n_tot_IPQ  = data.n_tot_IPQ
	n_late_IPQ = data.n_late_IPQ
	errors_IPQ = data.errors_IPQ
	frac_IPQ   = n_late_IPQ/n_tot_IPQ

	PLOT, findgen(10), findgen(10), xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', $
	  title = 'Matched SF and Q IP Samples ' + textoidl('(\Deltaz=2.0)'), /NODATA
	LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], BOX=0, /BOTTOM, /LEFT
	OPLOT, Rproj_array, frac_IPSF, LINESTYLE=0
	ERRPLOT, Rproj_array, frac_IPSF-errors_IPSF, frac_IPSF+errors_IPSF
	OPLOT, Rproj_array, frac_IPQ, LINESTYLE=2
	ERRPLOT, Rproj_array, frac_IPQ-errors_IPQ, frac_IPQ+errors_IPQ

	XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'SF IP: ' + bigint(n_elements(dataIP[where((dataIP.SFQ EQ 1) AND (dataIP.targ_weight GE 1.))])) + $
		', Med M* ' + decimal(median(dataIP[where(dataIP.SFQ EQ 1)].mstar),2) + $
		', Med z ' + decimal(median(dataIP[where(dataIP.SFQ EQ 1)].zprimus),2), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Q IP: ' + bigint(n_elements(dataIP[where((dataIP.SFQ EQ 0) AND (dataIP.targ_weight GE 1.))])) + $
		', Med M* ' + decimal(median(dataIP[where(dataIP.SFQ EQ 0)].mstar),2) + $
		', Med z ' + decimal(median(dataIP[where(dataIP.SFQ EQ 0)].zprimus),2), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'SF IP: ' + bigint(n_elements(dataIP[where((dataIP.SFQ EQ 1) AND (dataIP.targ_weight GE 1.))])), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Q IP: ' + bigint(n_elements(dataIP[where((dataIP.SFQ EQ 0) AND (dataIP.targ_weight GE 1.))])), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, '0.2 < z < 1.0', ALIGNMENT=1.0

	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
