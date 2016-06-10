PRO norm_diff_plot_IPmatchFBFonlyCOSMOS_PHI37, outputFormat
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsampleFBF/normsigplot_IPmatchFBFonlyCOSMOS_PHI37_allz', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	ERASE
	!P.MULTI=1
	!P.CHARSIZE=1.2
	
	xmin = 0
	xmax = 15
	xr = xmax - xmin

	ymin = -0.1
	ymax = 0.1
	yr = ymax - ymin

	PLOT, [xmin,xmax], [0,0], xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle=textoidl('Normalized % Difference'), $
		title='Normalized % diff. between LT- and ET-IP late-type fraction', LINESTYLE=2
	
	colors = [cgColor('Magenta'), cgColor('Purple'), cgColor('Blue')]

	zrange = 'allz_cosmos'
	data = MRDFITS('~/results/match_IP_sample_rigorous/fields/latefrac_' + zrange + '_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits', 1)

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
	OPLOT, Rmin+0.5*(Rmax[1]-Rmax[0]), norm_diff, LINESTYLE=ls;, COLOR=colors[i]
	
	XYOUTS, xmin+0.05*xr, ymin+0.10*yr, 'Matched SF and Q IP Samples', ALIGNMENT=0.0
	XYOUTS, xmin+0.05*xr, ymin+0.05*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=0.0
	XYOUTS, xmin+0.05*xr, ymin+0.925*yr, 'COSMOS field only', ALIGNMENT=0.0
	XYOUTS, xmin+0.05*xr, ymin+0.875*yr, 'z = [0.20, 1.00]', ALIGNMENT=0.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'

END
