PRO norm_diff_plot_IPmatchFBF_PHI37_zbinsEqualThirds, dR, COSMOScompare, outputFormat
	string_dR = 'dR1Mpc'
	datapath = '~/conformity/results/match_IP_sample_rigorous/'

	IF (dR EQ 250.) THEN BEGIN
		string_dR = 'dR250kpc'
		datapath = '~/conformity/results/match_IP_sample_rigorous/dR_250kpc/'
	ENDIF

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	  IF (COSMOScompare EQ 1) THEN $
		PS_OPEN, '~/figures/matchedIPsampleFBF/normDiffPlot_COSMOScomp_IPmatchFBF_PHI37_zbinsEqualThirds_' + string_dR, THICK=3, /ENCAP $
	  ELSE PS_OPEN, '~/figures/matchedIPsampleFBF/normDiffPlot_IPmatchFBF_PHI37_zbinsEqualThirds_' + string_dR, THICK=3, /ENCAP
		DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE
		
	!P.CHARSIZE=1.2
	
	xmin = 0
	xmax = 5
	plotMax = xmax-1
	xr = xmax - xmin

	ymin = -0.05
	ymax = 0.1
	yr = ymax - ymin

	zArray = [0.2, 0.482589, 0.684636, 1.0]

	PLOT, findgen(10), findgen(10), xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle=textoidl('Normalized % Difference'), $
		title='Normalized % diff. between LT- and ET-IP late-type fraction', /NODATA
	
	colors = [cgColor('Blue'), cgColor('Green'), cgColor('Red')]
;	colors = [cgColor('Black'), cgColor('Magenta'), cgColor('Blue'), cgColor('Cyan'), cgColor('Green')]
;	colors = [cgColor('Yellow'), cgColor('Magenta'), cgColor('Blue'), cgColor('Cyan'), cgColor('Green')]

	legend = []
	FOR i=0,2 DO legend = [legend, 'z = [' + decimal(zArray[i],2) + ', ' + decimal(zArray[i+1],2) + ']']

	FOR i=0,2 DO BEGIN
		zrange = 'T' + strtrim(string(i+1, format='(i20)'),2)

		data = MRDFITS(datapath + 'latefrac_' + zrange + '_targ_weight_IPmatchFBF_PHI3.7_' + string_dR + '_BSE.fits', 1)

		Rmin		= data.rmin
		Rmax		= data.rmax
		total_IPSF	= data.n_tot_IPSF
		total_IPQ	= data.n_tot_IPQ
		late_IPSF	= data.n_late_IPSF
		late_IPQ	= data.n_late_IPQ

		frac_IPSF	= float(late_IPSF)/total_IPSF
		frac_IPQ	= float(late_IPQ)/total_IPQ
		errors_IPSF	= data.errors_IPSF
		errors_IPQ	= data.errors_IPQ

		norm_diff	= (frac_IPSF - frac_IPQ)/((frac_IPSF + frac_IPQ)/2.)

		norm_diff_err	= norm_diff*SQRT((errors_IPSF/frac_IPSF)^2 + (errors_IPQ/frac_IPQ)^2)

;		IF i EQ 0 THEN ls = 2 ELSE ls = 0
		ls = 0
		OPLOT, Rmin+0.5*(Rmax[1]-Rmax[0]), norm_diff[0:plotMax], LINESTYLE=ls, COLOR=colors[i]
		ERRPLOT, Rmin+0.5*(Rmax[1]-Rmax[0]), norm_diff+norm_diff_err, norm_diff-norm_diff_err, COLOR=colors[i]
;		PRINT, Rmin, Rmax
	ENDFOR
  IF (COSMOScompare EQ 1) THEN BEGIN
	FOR i=0,2 DO BEGIN
		zrange = 'T' + strtrim(string(i+1, format='(i20)'),2)

		data = MRDFITS(datapath + 'noCosmos/latefrac_' + zrange + '_targ_weight_IPmatchFBFnoCosmos_PHI3.7_' + string_dR + '_BSE.fits', 1)

		Rmin		= data.rmin
		Rmax		= data.rmax
		total_IPSF	= data.n_tot_IPSF
		total_IPQ	= data.n_tot_IPQ
		late_IPSF	= data.n_late_IPSF
		late_IPQ	= data.n_late_IPQ

		frac_IPSF	= float(late_IPSF)/total_IPSF
		frac_IPQ	= float(late_IPQ)/total_IPQ
		errors_IPSF	= data.errors_IPSF
		errors_IPQ	= data.errors_IPQ

		norm_diff	= (frac_IPSF - frac_IPQ)/((frac_IPSF + frac_IPQ)/2.)

		norm_diff_err	= norm_diff*SQRT((errors_IPSF/frac_IPSF)^2 + (errors_IPQ/frac_IPQ)^2)

		ls = 3
		OPLOT, Rmin+0.5*(Rmax[1]-Rmax[0]), norm_diff[0:plotMax], LINESTYLE=ls, COLOR=colors[i]
		ERRPLOT, Rmin+0.5*(Rmax[1]-Rmax[0]), norm_diff+norm_diff_err, norm_diff-norm_diff_err, COLOR=colors[i]
	ENDFOR
	LEGEND, ['All fields', 'COSMOS excluded'], LINESTYLE=[0,3], BOX=0, /TOP, /RIGHT
  ENDIF	
	XYOUTS, xmin+0.05*xr, ymin+0.10*yr, 'Matched SF and Q IP samples', ALIGNMENT=0.0
	XYOUTS, xmin+0.05*xr, ymin+0.05*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=0.0
	LEGEND, legend, LINESTYLE=[0,0,0], COLOR=colors, BOX=0, /BOTTOM, /RIGHT
;	LEGEND, ['Unweighted totals (default)', 'Weighted by TARG_WEIGHT', 'Annulus median of all IPs'], COLOR=colors, LINESTYLE=[0,0,0], BOX=0, /TOP, /RIGHT

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'

END
