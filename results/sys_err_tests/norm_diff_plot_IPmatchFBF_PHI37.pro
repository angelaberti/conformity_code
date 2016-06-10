PRO norm_diff_plot_IPmatchFBF_PHI37, outputFormat
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsampleFBF/norm_diff_plot_IPmatchFBF_PHI37', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	!P.CHARSIZE=1.2
	
	xmin = 0
	xmax = 5
	xr = xmax - xmin

	ymin = -0.2
	ymax = 0.2
	yr = ymax - ymin

	zArray = [0.2, 0.482589, 0.684636, 1.0]

	PLOT, findgen(10), findgen(10), xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle=textoidl('Normalized % Difference'), $
		title='Normalized % diff. between LT- and ET-IP late-type fraction', /NODATA
	
	colors = [cgColor('Magenta'), cgColor('Purple'), cgColor('Blue')]
;	colors = [cgColor('Black'), cgColor('Magenta'), cgColor('Blue'), cgColor('Cyan'), cgColor('Green')]
;	colors = [cgColor('Yellow'), cgColor('Magenta'), cgColor('Blue'), cgColor('Cyan'), cgColor('Green')]

;	legend = ['z = [0.2, 1.0]']
;	FOR i=0,3 DO legend = [legend, 'z = [' + strtrim(string(0.2 + 0.2*i, format='(f20.1)'),1) + ', ' + strtrim(string(0.2 + 0.2*(i+1), format='(f20.1)'),1) + ']']

	FOR i=0,2 DO BEGIN
		zrange = 'T' + strtrim(string(i+1, format='(i20)'),2)
;		IF i EQ 0 THEN zrange = '0.2_1.0' $
;		ELSE zrange = strtrim(string(0.2 + 0.2*(i-1), format='(f20.1)'),1) + '_' + strtrim(string(0.2 + 0.2*(i), format='(f20.1)'),1)
		data = MRDFITS('~/results/match_IP_sample_rigorous/latefrac_' + zrange + '_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits', 1)

		Rmin		= data.rmin
		Rmax		= data.rmax
		total_IPSF	= data.n_tot_IPSF
		total_IPQ	= data.n_tot_IPQ
		late_IPSF	= data.n_late_IPSF
		late_IPQ	= data.n_late_IPQ

		frac_IPSF	= float(late_IPSF)/total_IPSF
		frac_IPQ	= float(late_IPQ)/total_IPQ

		norm_diff	= (frac_IPSF - frac_IPQ)/((frac_IPSF + frac_IPQ)/2.)

;		IF i EQ 0 THEN ls = 2 ELSE ls = 0
		ls = 0
		OPLOT, Rmin+0.5*(Rmax[1]-Rmax[0]), norm_diff, LINESTYLE=ls, COLOR=colors[i]
;		PRINT, Rmin, Rmax
	ENDFOR
	
	XYOUTS, xmin+0.05*xr, ymin+0.10*yr, 'Matched SF and Q IP Samples', ALIGNMENT=0.0
	XYOUTS, xmin+0.05*xr, ymin+0.05*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=0.0
	LEGEND, ['','',''], LINESTYLE=[0,0,0], COLOR=colors, BOX=0, /BOTTOM, /RIGHT
;	LEGEND, ['Unweighted totals (default)', 'Weighted by TARG_WEIGHT', 'Annulus median of all IPs'], COLOR=colors, LINESTYLE=[0,0,0], BOX=0, /TOP, /RIGHT

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'

END
