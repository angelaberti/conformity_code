; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

PRO latefracplot_matchedIP_zbinsEqualThirds_samePlot, outputFormat;, zmin, zmax;, dz_coeff, printEvery
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/latefracplot_matchedIP_zbinsEqualThirds', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

        dataIPallz  = MRDFITS('~/conformity/results/match_IP_sample/matchedIPsample.fits', 1)
	dataAllallz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

	Rmax 		= 5.
	dRproj		= 0.25
	n_annuli 	= Ceil(Rmax/dRproj)
	printEvery	= 100

        dataAllallz = dataAllallz[where(dataAllallz.targ_weight GE 1.)]
        dataIPallz  = dataIPallz[where(dataIPallz.targ_weight GE 1.)]

	zArray = [0.2, 0.502986, 0.699939, 1.0]
	PRINT, zArray

	legend = []
;	FOR i=0,2 DO legend = [legend, 'z = [' + strtrim(string(zArray[i], format='(f20.2)'),1) + ', ' + strtrim(string(zArray[i+1], format='(f20.2)'),1) + ']']
        FOR i=0,2 DO BEGIN
                dataIP = dataIPallz[where((dataIPallz.zprimus GE zArray[i]) AND (dataIPallz.zprimus LE zArray[i+1]))]
                legend = [legend, 'z = [' + decimal(zArray[i], 2) + ', ' + decimal(zArray[i+1], 2) + ']; Med M*, z: ' + $
                decimal(median(dataIP[where(dataIP.SFQ EQ 1)].mstar), 2) + ', ' + decimal(median(dataIP[where(dataIP.SFQ EQ 1)].zprimus), 2) + ' (SF) ' + $
		decimal(median(dataIP[where(dataIP.SFQ EQ 0)].mstar), 2) + ', ' + decimal(median(dataIP[where(dataIP.SFQ EQ 0)].zprimus), 2) + ' (Q)']
        ENDFOR

	!p.multi=0
	!p.charsize=1.1

	Rproj_array = float(dRproj)*(findgen(n_annuli+1)+0.5)

	xmin = 0
	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = 0.5
	ymax = 1
	yr = ymax-ymin

;	colors = ['Magenta', 'Blue', 'Cyan', 'Green']
	colors = ['Magenta', 'Blue', 'Green']

	FOR i=0,2 DO BEGIN
		z_low   = zArray[i]
		z_high  = zArray[i+1]

;		zsuffix = strtrim(string(0.2+0.2*i,format='(f20.1)'),1) + '_' + strtrim(string(0.2+0.2*(i+1),format='(f20.1)'),1)
		zsuffix = 'T' + strtrim(i+1,1)
;		zlabel  = strtrim(string(0.2+0.2*i,format='(f20.1)'),1) + ' < z < ' + strtrim(string(0.2+0.2*(i+1),format='(f20.1)'),1)

		lfdata = MRDFITS('~/conformity/results/match_IP_sample/latefrac_' + zsuffix + '_targ_weight_matchedIPsample_BSE.fits',1)
		
		n_tot_IPSF  = lfdata.n_tot_IPSF
		n_late_IPSF = lfdata.n_late_IPSF
		n_tot_IPQ  = lfdata.n_tot_IPQ
		n_late_IPQ = lfdata.n_late_IPQ
		errors_IPSF = lfdata.errors_IPSF
		errors_IPQ  = lfdata.errors_IPQ

;		dataIP  = dataIPallz[where( (dataIPallz.zprimus GE z_low) AND (dataIPallz.zprimus LE z_high) )]
;		dataAll = dataAllallz[where( (dataAllallz.zprimus GE z_low) AND (dataAllallz.zprimus LE z_high) )]

		IF i EQ 0 THEN BEGIN
			PLOT, findgen(10), findgen(10), xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', $
			  title = 'Matched SF and Q IP Samples ' + textoidl('(\Deltaz=2.0)'), /NODATA
			XYOUTS, xmin+0.95*xr, ymin+0.9375*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
			XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0
			LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], BOX=0, /TOP, /LEFT
			LEGEND, legend, LINESTYLE=[0,0,0], color=cgColor(colors), BOX=0, /BOTTOM, /LEFT
		ENDIF

		OPLOT, Rproj_array, n_late_IPSF/float(n_tot_IPSF), LINESTYLE=0, color=cgColor(colors[i])
                IF (i EQ 0) OR (i EQ 2) THEN $
                ERRPLOT, Rproj_array, n_late_IPSF/float(n_tot_IPSF)+errors_IPSF, n_late_IPSF/float(n_tot_IPSF)-errors_IPSF, color=cgColor(colors[i])

	        OPLOT, Rproj_array, n_late_IPQ/float(n_tot_IPQ), LINESTYLE=2, color=cgColor(colors[i])
                IF (i EQ 0) OR (i EQ 2) THEN $
                ERRPLOT, Rproj_array, n_late_IPQ/float(n_tot_IPQ)+errors_IPQ, n_late_IPQ/float(n_tot_IPQ)-errors_IPQ, color=cgColor(colors[i])

;		XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'SF IP: ' + bigint(n_elements(dataIP[where((dataIP.SFQ EQ 1) AND (dataIP.targ_weight GE 1.))])), ALIGNMENT=1.0
;		XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Q IP: ' + bigint(n_elements(dataIP[where((dataIP.SFQ EQ 0) AND (dataIP.targ_weight GE 1.))])), ALIGNMENT=1.0
;		XYOUTS, xmin+0.95*xr, ymin+0.05*yr, zlabel, ALIGNMENT=1.0
	ENDFOR

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
