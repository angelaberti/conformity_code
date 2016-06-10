; 'type' can be (1) 'equalIP' (2) 'equalRange' (3) 'halfDex'

PRO norm_diff_plot_IPmatchFBFnoCosmos_PHI37_massThirds, type, outputFormat
	string_dR = 'dR1Mpc'
	datapath = '~/conformity/results/match_IP_sample_rigorous/noCosmos/massThirds/'

;	IF (dR EQ 250.) THEN BEGIN
;		string_dR = 'dR250kpc'
;		datapath = '~/conformity/results/match_IP_sample_rigorous/dR_250kpc/'
;	ENDIF

	!P.CHARSIZE=1.1
	
	xmin = 0
	xmax = 5
	plotMax = xmax-1
	xr = xmax - xmin

	ymin = -0.15
	ymax = 0.15
	yr = ymax - ymin

	IF (type EQ 'halfDex') THEN BEGIN
		massArray = [10., 10.5, 11.0, 11.5]
		typeLabel = '0.5 dex bins'
	ENDIF
	IF (type EQ 'equalRange') THEN BEGIN
		massArray = [9.1278, 9.86134, 10.5949, 11.3284]
		typeLabel = 'Equal mass range per bin'
	ENDIF
	IF (type EQ 'equalIP') THEN BEGIN
		massArray = [9.1278, 10.6679, 10.953, 11.3284]
		typeLabel = 'Equal number of IP per bin'
	ENDIF

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsampleFBF/norm_diff_plot_IPmatchFBFnoCosmos_PHI37_massThirds' + type + '_' + string_dR, THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	PLOT, [0,xmax], [0,0], xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle=textoidl('Normalized % Difference'), $
		title='Normalized % diff. between LT- and ET-IP late-type fraction', LINESTYLE=2
	
	colors = [cgColor('Purple'), cgColor('Magenta'), cgColor('Blue')]

	legend = []
	FOR i=0,2 DO BEGIN
		massRange = 'M' + strtrim(string(i+1, format='(i20)'),2)
		masslabel = textoidl('M_{*} = [') + decimal(massArray[i],2) + ', ' + decimal(massArray[i+1],2) + ']'
		data = MRDFITS(datapath + 'latefrac_' + massRange + type + '_targ_weight_IPmatchFBFnoCosmos_PHI3.7_' + string_dR + '_BSE.fits', 1)

		legend = [legend, masslabel $
			+ textoidl('; \sigma_{0-1}=') + getSigmaRange(data, 0, 0) $
			+ textoidl('; \sigma_{1-2}=') + getSigmaRange(data, 1, 1) $
			+ textoidl('; \sigma_{1-3}=') + getSigmaRange(data, 1, 2) $
			+ textoidl('; \sigma_{3-5}=') + getSigmaRange(data, 3, 4)]

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
		OPLOT, Rmin[0:plotMax]+0.5*(Rmax[1]-Rmax[0]), norm_diff[0:plotMax], LINESTYLE=ls, COLOR=colors[i]
		ERRPLOT, Rmin+0.5*(Rmax[1]-Rmax[0]), norm_diff+norm_diff_err, norm_diff-norm_diff_err, COLOR=colors[i]
;		PRINT, Rmin, Rmax
	ENDFOR
	
;	XYOUTS, xmin+0.05*xr, ymin+0.10*yr, 'FBF matched SF and Q IP samples', ALIGNMENT=0.0
;	XYOUTS, xmin+0.05*xr, ymin+0.05*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=0.0
	XYOUTS, xmin+0.05*xr, ymin+0.925*yr, 'COSMOS field excluded', ALIGNMENT=0.0
	XYOUTS, xmin+0.95*xr, ymin+0.85*yr, '0.2 < z < 1.0', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, typeLabel, ALIGNMENT=1.0
	LEGEND, legend, LINESTYLE=[0,0,0], COLOR=colors, BOX=0, THICK=3, /BOTTOM, /RIGHT
;	LEGEND, ['Unweighted totals (default)', 'Weighted by TARG_WEIGHT', 'Annulus median of all IPs'], COLOR=colors, LINESTYLE=[0,0,0], BOX=0, /TOP, /RIGHT

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'

END
