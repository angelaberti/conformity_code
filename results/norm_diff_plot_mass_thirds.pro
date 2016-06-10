; uses histograms and delta_z <= dz_coeff*0.005*(1+z)

PRO norm_diff_plot_mass_thirds, outputFormat
	input_dm = 1.0
	string_dm = strtrim(strcompress(string(input_dm, format='(f20.1)')),1)

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/norm_diff_mass_thirds_' + string_dm + 'dex', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	!P.CHARSIZE=1.1

	dz_coeff = 2.
	zmin = 0.2
	zmax = 1.0

	n_annuli 	= 60
	dRproj		= 0.25

	Rproj_array = float(dRproj)*findgen(n_annuli+1)

	xmin = 0
        xmax = n_annuli*dRproj
        xr = xmax-xmin

        ymin = -0.05
        ymax = 0.20
        yr = ymax-ymin

;	!P.MULTI = 0

	colors = ['Green', 'Cyan', 'Blue', 'Magenta', 'Red']

	PLOT, findgen(10), findgen(10), xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Normalized % diff. in late-type frac between LT/ET IPs', /NODATA, $
		title = 'Signal strength IP mass thirds ' + textoidl('(var. IP mass diff. 1.0 dex; \Deltaz=2.0)')

	medSF = ['10.35', '10.75', '11.08']
	medQ  =	['10.39', '10.77', '11.10']

	LEGEND, textoidl('Median M_{IP}: ') + medSF + ' (SF) ' + medQ + ' (Q)', LINESTYLE=[0,0,0], color=cgColor(colors[0:2]), BOX=0, /TOP, /RIGHT
;	LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], BOX=0, /BOTTOM, /LEFT
;	XYOUTS, xmin+0.05*xr, ymin+0.05*yr, '0.2 < z < 1.0', ALIGNMENT=0.0
;	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, IPtypeLabel + ' IPs', ALIGNMENT=1.0

	norm_diff_T1 = get_norm_diff(1, string_dm)
	norm_diff_T2 = get_norm_diff(2, string_dm)
	norm_diff_T3 = get_norm_diff(3, string_dm)

	FOR i=0,2 DO BEGIN
		OPLOT, Rproj_array+0.5*dRproj, get_norm_diff(i+1, string_dm), color=cgColor(colors[i])
	ENDFOR
;	OPLOT, Rproj_array, norm_diff_T2, color=cgColor(colors[1])
;	OPLOT, Rproj_array, norm_diff_T3, color=cgColor(colors[2])

	SAVE, Rproj_array, norm_diff_T1, norm_diff_T2, norm_diff_T3, FILE='norm_diff_mass_thirds.sav'

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

FUNCTION get_norm_diff, third, string_dm
	latefracdata 	= '~/conformity/results/mass_bins/latefrac_data_hist/variable_mass+' + string_dm + 'dex/latefrac_0.2_1.0_dataIP_T' + strtrim(strcompress(third),1) + '.fits'
	data		= mrdfits(latefracdata, 1)

	fracIPSF = data.n_late_IPSF/data.n_tot_IPSF
	fracIPQ  = data.n_late_IPQ/data.n_tot_IPQ
		
	norm_diff = (fracIPSF-fracIPQ)/((fracIPSF+fracIPQ)/2)

	RETURN, norm_diff	
END
