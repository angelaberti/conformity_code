; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)
; FD = field-dependent SFR vs. mass line

PRO latefracplot_targ_weight_FD, input_dm, outputFormat;, zmin, zmax;, dz_coeff, printEvery

        string_dm = strtrim(strcompress(string(input_dm, format='(f20.1)')),1)

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/latefracplot_targ_weight_' + string_dm + 'dex_zall_FD', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

        IPdataPath = '~/field-dependent_SFR_vs_mass/'
        zerodInputFile = 'zerodSFQ_IP_dz' + strtrim(strcompress(string(dz_coeff, format='(f20.1)')),1) + '_' + string_dm + 'dex.fits'

;	PRINT, 'Input data: ' + IPdataPath + zerodInputFile

        data = mrdfits('~/field-dependent_SFR_vs_mass/latefrac_' + strtrim(strcompress(string(zmin, format='(f20.1)')), 1) + '_' $
				+ strtrim(strcompress(string(zmax, format='(f20.1)')) ,1) + '_targ_weight_' + zerodInputFile,1)
	Rmax 		= 15.
	dRproj		= 0.25
	n_annuli 	= Ceil(Rmax/dRproj)
	printEvery	= 100

;	data = zerodInput[where((zerodInput.zprimus GE zmin) and (zerodInput.zprimus LE zmax))]
;	data = data[where(data.targ_weight GE 1.)]

	!p.multi=0
	!p.charsize=1.25

	Rproj_array = float(dRproj)*(findgen(n_annuli+1)+0.5)

	xmin = 0
	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = 0.35
	ymax = 1
	yr = ymax-ymin

;	get_targ_weight, data, dz_coeff, 1, neigh_all_grand_total, neigh_late_grand_total, n_annuli, dRproj, printEvery

	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPSF = data.n_late_IPSF
	n_tot_IPQ  = data.n_tot_IPQ
	n_late_IPQ = data.n_late_IPQ

	PLOT, Rproj_array, float(n_late_IPSF)/n_tot_IPSF, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', $
	  title = 'M13 Mass Limit + ' + string_dm + ' (0.0) dex for SF (Q) IP ' + textoidl('(\Deltaz=2.0)'), LINESTYLE=0
	LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], BOX=0, /BOTTOM, /LEFT
	OPLOT, Rproj_array, float(n_late_IPQ)/n_tot_IPQ, LINESTYLE=2

	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, '0.2 < z < 1.0', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'

END
