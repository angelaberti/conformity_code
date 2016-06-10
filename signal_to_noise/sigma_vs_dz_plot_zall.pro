PRO sigma_vs_dz_plot_zall, outputFormat ;, dm, zmin, zmax, show_legend, inputFiles
	zmin = 0.2
	zmax = 1.0
	dm = 0.5
	show_legend = 1

	IF (string(outputFormat) eq 'ps') THEN BEGIN
		PS_OPEN, '../figures/sigma_vs_dz_zall', THICK=5, /ENCAP
		DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
		SET_PLOT, 'X'
	ENDELSE

	!p.multi = 0
	!p.charsize = 1.5

	xmin	= 0.5
	xmax	= 3.0
	xr 	= xmax - xmin

	ymin	= 0.
	ymax	= 10
	yr	= ymax - ymin

	dz_coeffArray = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

	sigmaRange000_075 = [] ; sigmaRange(inputFile,0,2)
	sigmaRange075_200 = [] ; sigmaRange(inputFile,3,7)
	sigmaRange075_300 = [] ; sigmaRange(inputFile,3,11)
	sigmaRange075_400 = [] ; sigmaRange(inputFile,3,15)

; build up arrays of sigma values for each dz_coeff
	FOREACH element, dz_coeffArray DO BEGIN
		dz = element
;		data = mrdfits(element, 1)
		inputFile = '../results/default_parameters/latefrac_data_hist/latefrac_' + strtrim(string(zmin, format='(f20.1)'),1) + '_' $
					+ strtrim(string(zmax, format='(f20.1)'),1) $
					+ '_zerodSFQ_IP_dz' + strtrim(strcompress(string(dz, format='(f20.1)')),1) + '_dm' $
					+ strtrim(string(dm, format='(f20.1)'),1) + '.fits'
		print, inputFile
		data = mrdfits(inputFile, 1)

		frac_IPSF = data.n_late_IPSF/data.n_tot_IPSF
		frac_IPQ  = data.n_late_IPQ/data.n_tot_IPQ

		dfrac_IPSF = poissonError(data.n_late_IPSF, data.n_tot_IPSF)
		dfrac_IPQ  = poissonError(data.n_late_IPQ, data.n_tot_IPQ)

		sigmaRange000_075 = [sigmaRange000_075, sigmaRange(inputFile,0,2)]
		sigmaRange075_200 = [sigmaRange075_200, sigmaRange(inputFile,3,7)]
		sigmaRange075_300 = [sigmaRange075_300, sigmaRange(inputFile,3,11)]
		sigmaRange075_400 = [sigmaRange075_400, sigmaRange(inputFile,3,15)]
	ENDFOREACH

                nx=1.
                ny=1.

                xm=0.1
                ym=0.15

                dx = (1 - 1.5*xm)/nx
                dy = (1 - 2.0*ym)/ny

                FOR i=0,nx*ny-1 DO BEGIN
                        IF (i le (nx-1)) THEN BEGIN
                                xtitle = textoidl("\Deltaz / (0.005(1+z_{IP}))")
                                xtickformat = ''
                        ENDIF ELSE BEGIN
                                xtitle = ''
                                xtickformat = '(A1)'
                        ENDELSE
                        IF (i mod nx eq 0) THEN BEGIN
                                ytitle = textoidl("\sigma over range of projected radii")
                                ytickformat = ''
                        ENDIF ELSE BEGIN
                                ytitle = ''
                                ytickformat = '(A1)'
                        ENDELSE

                        pos = [ xm + dx*float(i mod nx), ym + dy*float(floor(float(i)/nx)), $
                                xm + dx*(1. + float(i mod nx)), ym + dy*(1. + float(floor(float(i)/nx))) ]
                        PRINT, pos

	plot, dz_coeffArray, sigmaRange000_075, xrange=[xmin,xmax], yrange=[ymin,ymax], linestyle=0, position=pos, $
;	  xtitle=textoidl("\Deltaz / (0.005(1+z_{IP}))"), ytitle=textoidl("\sigma over range of projected radii"), 
	  xtitle=xtitle, ytitle=ytitle, title='S/N vs. Cylinder Depth (default=1.0)'
	oplot, dz_coeffArray, sigmaRange075_200, linestyle=1
	oplot, dz_coeffArray, sigmaRange075_300, linestyle=2
	oplot, dz_coeffArray, sigmaRange075_400, linestyle=3

	IF show_legend eq 1 THEN BEGIN
		LEGEND, ['R < 0.75 Mpc', '0.75 to 2 Mpc', '0.75 to 3 Mpc', '0.75 to 4 Mpc'], $
			linestyle=[0,1,2,3], box=0, /BOTTOM, /LEFT
	ENDIF

        xyouts, xmin+0.95*xr, ymin+0.05*yr, strtrim(string(zmin, format='(f20.1)'),1) + " < z < " + strtrim(string(zmax, format='(f20.1)'),1), $
	  ALIGNMENT=1.0

	ENDFOR

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
