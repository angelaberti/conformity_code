
PRO norm_diff_plot_matchedIP_byField_highzThird, outputFormat
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/norm_diff_plot_matchedIP_byField_highzThird', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	!P.CHARSIZE=1.2
	
	xmin = 0
	xmax = 5
	xr = xmax - xmin

	ymin = -0.3
	ymax = 0.35
	yr = ymax - ymin

	data = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)
	fields = data[uniq(data.field, sort(data.field))].field	

	PLOT, findgen(10), findgen(10), xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle=textoidl('Normalized % Difference'), $
		title='Normalized % diff. between LT- and ET-IP late-type fraction', /NODATA
	
	colors = [cgColor('Red'), cgColor('Magenta'), cgColor('Blue'), cgColor('Cyan'), cgColor('Green')]
;	colors = [cgColor('Yellow'), cgColor('Magenta'), cgColor('Blue'), cgColor('Cyan'), cgColor('Green')]

	legend = []
	zmin = 0.699939

	FOR i=0,4 DO BEGIN
		field = fields[i]
		stringField = strtrim(strcompress(field),2)

		dataIP = MRDFITS('~/conformity/results/match_IP_sample/matchedIPsample.fits',1)
		dataIPfield = dataIP[where((dataIP.field EQ field) AND (dataIP.targ_weight GE 1.) AND (dataIP.zprimus GE zmin))]

		legend = [legend, stringField + '; Med M*, z: ' + $
			decimal(median(dataIPfield[where(dataIPfield.SFQ EQ 1)].mstar),2) + ', ' + $
			decimal(median(dataIPfield[where(dataIPfield.SFQ EQ 1)].zprimus),2) + ' (SF) ' + $
			decimal(median(dataIPfield[where(dataIPfield.SFQ EQ 0)].mstar),2) + ', ' + $
			decimal(median(dataIPfield[where(dataIPfield.SFQ EQ 0)].zprimus),2) + ' (Q)']

		zrange = '0.70_1.00'
		data = MRDFITS('~/conformity/results/match_IP_sample/byField/latefrac_' + zrange + '_targ_weight_matchedIPsample_' + stringField + '_dR1Mpc_BSE.fits', 1)

		Rmin		= data.rmin
		Rmax		= data.rmax
		total_IPSF	= data.n_tot_IPSF
		total_IPQ	= data.n_tot_IPQ
		late_IPSF	= data.n_late_IPSF
		late_IPQ	= data.n_late_IPQ

		frac_IPSF	= float(late_IPSF)/total_IPSF
		frac_IPQ	= float(late_IPQ)/total_IPQ

		norm_diff	= (frac_IPSF - frac_IPQ)/((frac_IPSF + frac_IPQ)/2.)

		ls = 0
		OPLOT, Rmin, norm_diff, LINESTYLE=ls, COLOR=colors[i]
	ENDFOR
	
	XYOUTS, xmin+0.05*xr, ymin+0.90*yr, 'Matched SF and Q IP Samples', ALIGNMENT=0.0
	XYOUTS, xmin+0.05*xr, ymin+0.85*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=0.0
	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, decimal(zmin,2) + ' < z < 1.0', ALIGNMENT=1.0

	LEGEND, legend, LINESTYLE=[0,0,0,0,0], COLOR=colors, BOX=0, /BOTTOM, /RIGHT

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
