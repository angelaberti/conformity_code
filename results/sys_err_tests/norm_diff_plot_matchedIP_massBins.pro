PRO norm_diff_plot_matchedIP_massBins, outputFormat
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/norm_diff_plot_matchedIP_massBins', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	!P.CHARSIZE=1.2
	
	xmin = 0
	xmax = 5
	xr = xmax - xmin

	ymin = -0.5
	ymax = 0.5
	yr = ymax - ymin

	PLOT, findgen(10), findgen(10), xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle=textoidl('Normalized % Difference'), $
		title='Normalized % diff. between LT- and ET-IP late-type fraction', /NODATA
	
	colors = [cgColor('Magenta'), cgColor('Blue'), cgColor('Green')]
;	colors = [cgColor('Black'), cgColor('Magenta'), cgColor('Blue'), cgColor('Green')]
;	colors = [cgColor('Yellow'), cgColor('Magenta'), cgColor('Blue'), cgColor('Cyan'), cgColor('Green')]

	legend = []

	massArray = [9.14254, 10.6708, 10.9550, 11.3974]

        dataIPallm  = MRDFITS('~/conformity/results/match_IP_sample/matchedIPsample.fits', 1)
        dataAllallm = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

        dataAllallm = dataAllallm[where(dataAllallm.targ_weight GE 1.)]
        dataIPallm  = dataIPallm[where(dataIPallm.targ_weight GE 1.)]

	FOR i=0,n_elements(massArray)-2 DO BEGIN
                mass_low  = massArray[i]
                mass_high = massArray[i+1]

                dataIP  = dataIPallm[where( (dataIPallm.mstar GE mass_low) AND (dataIPallm.mstar LE mass_high) )]
                dataAll = dataAllallm[where( (dataAllallm.mstar GE mass_low) AND (dataAllallm.mstar LE mass_high) )]

                msuffix   = 'M' + strtrim(i+1,1)
                masslabel = 'M* = [' + decimal(mass_low,2) + ', ' + decimal(mass_high,2) + ']; Med z = ' + $
                        decimal(median(dataIP[where(dataIP.SFQ EQ 1)].zprimus),2) + ' (SF) ' + $
                        decimal(median(dataIP[where(dataIP.SFQ EQ 0)].zprimus),2) + ' (Q)'
                legend = [legend, masslabel]
;		legend = [legend, 'M* = [' + decimal(massArray[i],2) + ', ' + decimal(massArray[i+1],2) + ']']

		mrange = 'M' + strtrim(i+1,1)
		data = MRDFITS('~/conformity/results/match_IP_sample/massBins/latefrac_' + mrange + '_targ_weight_matchedIPsample_dR1Mpc_BSE.fits', 1)

		Rmin		= data.rmin
		Rmax		= data.rmax
		total_IPSF	= data.n_tot_IPSF
		total_IPQ	= data.n_tot_IPQ
		late_IPSF	= data.n_late_IPSF
		late_IPQ	= data.n_late_IPQ

		frac_IPSF	= float(late_IPSF)/total_IPSF
		frac_IPQ	= float(late_IPQ)/total_IPQ

		norm_diff	= (frac_IPSF - frac_IPQ)/((frac_IPSF + frac_IPQ)/2.)

		OPLOT, Rmin, norm_diff, LINESTYLE=0, COLOR=colors[i]
;		PRINT, Rmin, Rmax
	ENDFOR
	
	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, '0.2 < z < 1.0', ALIGNMENT=1.0
	XYOUTS, xmin+0.05*xr, ymin+0.90*yr, 'Matched SF and Q IP Samples', ALIGNMENT=0.0
	XYOUTS, xmin+0.05*xr, ymin+0.85*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=0.0
	LEGEND, legend, LINESTYLE=[0,0,0], COLOR=colors, BOX=0, /BOTTOM, /RIGHT
;	LEGEND, ['Unweighted totals (default)', 'Weighted by TARG_WEIGHT', 'Annulus median of all IPs'], COLOR=colors, LINESTYLE=[0,0,0], BOX=0, /TOP, /RIGHT

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
