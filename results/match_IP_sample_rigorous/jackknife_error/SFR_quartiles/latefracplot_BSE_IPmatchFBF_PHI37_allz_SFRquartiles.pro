PRO latefracplot_BSE_IPmatchFBF_PHI37_allz_SFRquartiles, outputFormat;, zmin, zmax;, dz_coeff, printEvery
	PHI = 3.7
	stringPHI = strtrim(string(PHI, format='(f20.1)'),2)

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	Rmax 		= 5.
	dRproj		= 1.
	n_annuli 	= Ceil(Rmax/dRproj)
	printEvery	= 100

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 0.50 THEN string_dR = '_dR500kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

	!p.multi=0
	!p.charsize=1.25

	Rproj_array = float(dRproj)*(findgen(n_annuli+1))	
;	PRINT, Rproj_array

	xmin = 0
	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = 0.71
	ymax = 0.88
	yr = ymax-ymin

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/latex/figures/latefracplot_BSE_SFRquartiles', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

  colors = ['red','orange','cyan','blue']
  FOR i=0,3 DO BEGIN
	data = MRDFITS('~/results/match_IP_sample_rigorous/jackknife_error/SFR_quartiles/latefrac_BSE_quart'+STRTRIM(i+1,2)+'.fits', 1)

	n_tot_IP  = data.n_tot_IP
	n_late_IP = data.n_late_IP
	errors_IP = data.errors_IP
	frac_IP   = n_late_IP/n_tot_IP

    IF (i EQ 0) THEN $
	PLOT, Rproj_array+0.5*dRproj, float(n_late_IP)/n_tot_IP, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', /NODATA

	OPLOT, Rproj_array+0.5*dRproj, frac_IP, LINESTYLE=0, color=cgColor(colors[i])
	ERRPLOT, Rproj_array+0.5*dRproj, frac_IP-errors_IP, frac_IP+errors_IP, color=cgColor(colors[i])

  ENDFOR
	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
