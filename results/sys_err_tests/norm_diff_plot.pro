;HAVE LATEFRAC DATA SAVED

;conservative_mass_cutoff (+0.0dex)
;variable_mass+0.5dex
;variable_mass+0.8dex
;variable_mass+1.0dex

PRO norm_diff_plot, outputFormat
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/norm_diff_plot', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=10, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	!P.CHARSIZE=1.2
	
	xmin = 0
	xmax = 15
	xr = xmax - xmin

	ymin = 0
	ymax = 0.1
	yr = ymax - ymin

	PLOT, findgen(10), findgen(10), xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle=textoidl('Normalized % Difference'), $
		title='Normalized % diff. between LT- and ET-IP late-type fraction', /NODATA
	
	var_mass_paths 	  = []
	targ_weight_paths = []
	median_paths 	  = []
	FOREACH dm,['0.5','0.8','1.0'] DO BEGIN
		var_mass_paths    = [var_mass_paths, '~/results/variable_mass+' + dm + 'dex/latefrac_data_hist/latefrac_0.2_1.0_zerodSFQ_IP_dz2.0.fits']
		targ_weight_paths = [targ_weight_paths, '~/results/variable_mass+' + dm + 'dex/latefrac_data_hist/latefrac_0.2_1.0_targ_weight_zerodSFQ_IP_dz2.0.fits']
		median_paths 	  = [median_paths, '~/results/variable_mass+' + dm + 'dex/latefrac_data_hist/latefrac_0.2_1.0_median_zerodSFQ_IP_dz2.0.fits']
	ENDFOREACH

	dataPaths    = ['~/results/conservative_mass_cutoff/latefrac_data_hist/latefrac_0.2_1.0_zerodSFQ_IP_dz2.0_dm0.0.fits', var_mass_paths, $
			'~/results/conservative_mass_cutoff/latefrac_data_hist/latefrac_0.2_1.0_targ_weight_zerodSFQ_IP_dz2.0_dm0.0.fits', targ_weight_paths, $
			'~/results/conservative_mass_cutoff/latefrac_data_hist/latefrac_0.2_1.0_median_zerodSFQ_IP_dz2.0_dm0.0.fits', median_paths]

	massDiffs = ['0.0', '0.5', '0.8', '1.0']
	linestyles = [1,2,3,0]
	colors = [cgColor('Magenta'), cgColor('Blue'), cgColor('Green')]
;	colors = [cgColor('Magenta'), cgColor('Purple'), cgColor('Blue'), cgColor('Green')]

	FOR i=0,7 DO BEGIN
		data = mrdfits(dataPaths[i],1)

		Rmax		= data.rmax
		total_IPSF	= data.n_tot_IPSF
		total_IPQ	= data.n_tot_IPQ
		late_IPSF	= data.n_late_IPSF
		late_IPQ	= data.n_late_IPQ

		frac_IPSF	= float(late_IPSF)/total_IPSF
		frac_IPQ	= float(late_IPQ)/total_IPQ

		norm_diff	= (frac_IPSF - frac_IPQ)/((frac_IPSF + frac_IPQ)/2.)

		IF (i LE 3) THEN BEGIN
			color = colors[0]
			linestyle = linestyles[i]
		ENDIF ELSE BEGIN
			color = colors[1]
			linestyle = linestyles[i-4]
		ENDELSE

		OPLOT, Rmax, norm_diff, LINESTYLE=linestyle, COLOR=color
	ENDFOR

	FOR i=8,11 DO BEGIN
		data = mrdfits(dataPaths[i],1)

		dR = data[1].rmax - data[0].rmax
		Rmax		= data.rmax
		frac_IPSF	= data.frac_IPSF
		frac_IPQ	= data.frac_IPQ

		norm_diff	= (frac_IPSF - frac_IPQ)/((frac_IPSF + frac_IPQ)/2.)

		OPLOT, Rmax-dR, norm_diff, LINESTYLE=linestyles[i-8], COLOR=colors[2]
	ENDFOR
	
	XYOUTS, xmin+0.95*xr, ymin+0.75*yr, 'SF/Q IP mass comp. diff.', ALIGNMENT=1.0
	LEGEND, massDiffs + ' dex', LINESTYLE=linestyles, POSITION=[xmin+0.67*xr, ymin+0.75*yr], BOX=0;, ALIGNMENT=1.0; /CENTER, /RIGHT
	LEGEND, ['Unweighted totals (default)', 'Weighted by TARG_WEIGHT', 'Annulus median of all IPs'], COLOR=colors, LINESTYLE=[0,0,0], BOX=0, /TOP, /RIGHT

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'

;LATEFRAC DATA GENERATED WHEN PROGRAM IS RUN (CHANGE)

;conservative_mass_median
;variable_mass+0.5dex_median
;variable_mass+0.8dex_median
;variable_mass+1.0dex_median

;conservative_mass_targ_weight 	* TO DO *
;variable_mass+0.5dex_targ_weight
;variable_mass+0.8dex_targ_weight
;variable_mass+1.0dex_targ_weight

END
