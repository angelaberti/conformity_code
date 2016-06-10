; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)
; PNZ = permute neighbor redshifts

PRO latefracplot_variable_targ_weight_PNZ, outputFormat;, zmin, zmax;, dz_coeff, printEvery

	string_dm = '1.0'

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/latefracplot_targ_weight_1.0dex_zall_PNZ', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

        IPdataPath = '~/conformity/results/variable_mass+1.0dex/IP_data/'
        zerodInputFile = 'zerodSFQ_IP_dz' + strtrim(strcompress(string(dz_coeff, format='(f20.1)')),1) + '.fits'

	PRINT, 'Input data: ', zerodInputFile

        zerodInput = mrdfits(IPdataPath + zerodInputFile, 1)
	dataIP = zerodInput[where(zerodInput.IP EQ 1)]

	Rmax 		= 15.
	dRproj		= 0.25
	n_annuli 	= Ceil(Rmax/dRproj)
	printEvery	= 100

	!p.multi=0
	!p.charsize=1.1

	Rproj_array = float(dRproj)*(findgen(n_annuli+1)+0.5)

	xmin = 0
	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = 0.5
	ymax = 1
	yr = ymax-ymin

;	get_targ_weight, data, dz_coeff, 1, neigh_all_grand_total, neigh_late_grand_total, n_annuli, dRproj, printEvery

	data = mrdfits('~/conformity/results/permute_neighbor_z/latefrac_' + strtrim(strcompress(string(zmin, format='(f20.1)')), 1) + '_' + strtrim(strcompress(string(zmax, format='(f20.1)')) ,1) + '_targ_weight_' + zerodInputFile, 1)

	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPSF = data.n_late_IPSF
	n_tot_IPQ  = data.n_tot_IPQ
	n_late_IPQ = data.n_late_IPQ

	PLOT, findgen(10), findgen(10), xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', $
	  title = 'M13 Mass Limit + 1.0 (0.0) dex for SF (Q) IPs ' + textoidl('(\Deltaz=2.0)'), /NODATA
	LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], BOX=0, /BOTTOM, /LEFT
	OPLOT, Rproj_array, n_late_IPSF/n_tot_IPSF, LINESTYLE=0
	OPLOT, Rproj_array, n_late_IPQ/n_tot_IPQ, LINESTYLE=2

        XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'SF IP: ' + bigint(n_elements(dataIP[where((dataIP.SFQ EQ 1) AND (dataIP.targ_weight GE 1.))])) + $
                ', Med M* ' + decimal(median(dataIP[where(dataIP.SFQ EQ 1)].mstar),2) + $
                ', Med z ' + decimal(median(dataIP[where(dataIP.SFQ EQ 1)].zprimus),2), ALIGNMENT=1.0
        XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Q IP: ' + bigint(n_elements(dataIP[where((dataIP.SFQ EQ 0) AND (dataIP.targ_weight GE 1.))])) + $
                ', Med M* ' + decimal(median(dataIP[where(dataIP.SFQ EQ 0)].mstar),2) + $
                ', Med z ' + decimal(median(dataIP[where(dataIP.SFQ EQ 0)].zprimus),2), ALIGNMENT=1.0

	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, '0.2 < z < 1.0', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.85*yr, 'NEIGHBOR REDSHIFTS PERMUTED WITHIN FIELD', ALIGNMENT=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
