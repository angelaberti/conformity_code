PRO SN_vs_dzCoeff, outputFormat ;, dm, zmin, zmax, show_legend, inputFiles
;	zmin = 0.2
;	zmax = 1.0
;	dm = 0.5
;	show_legend = 1

	IF (STRING(outputFormat) EQ 'ps') THEN BEGIN
		PS_OPEN, '~/figures/matchedIPsampleFBF/SN_vs_dzCoeff', THICK=3, /ENCAP
		DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
		SET_PLOT, 'X'
	ENDELSE

	!p.multi = [0,2,2]
	!p.charsize = 1.1

        nx=2.
        ny=2.

        xm=0.08
        ym=0.1

	dx = (1 - 1.5*xm)/nx
	dy = (1 - 2.0*ym)/ny

	xmin	= 0.25
	xmax	= 3.25
	xr 	= xmax - xmin

	ymin	= -0.5
	ymax	= 6
	yr	= ymax - ymin

	dz_coeffArray = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

  FOR i=0,nx*ny-1 DO BEGIN
	zmin = 0.2*FLOAT(i+1)
	zmax = 0.2*FLOAT(i+2)

	sigmaArray0_1 = []
	sigmaArray1_2 = []
	sigmaArray1_3 = []
	sigmaArray3_5 = []

	IF (i EQ 1) THEN (zrange = 'allz') ELSE $
	IF (i EQ 0) THEN (zrange = 'T3') ELSE $
	IF (i EQ 2) THEN (zrange = 'T1') ELSE $
	IF (i EQ 3) THEN (zrange = 'T2')
	;(zrange = 'T' + STRTRIM(STRING(i,format='(i20)'),2))
	PRINT, i, ' z range: ' + zrange

	zlabels = ['high z', 'all z', 'low z', 'mid z']
; build up arrays of sigma values for each dz
	FOREACH dz,dz_coeffArray DO BEGIN
		string_dzCoeff = STRTRIM(STRING(dz, format='(f20.1)'),2)

		inputFile = '~/results/match_IP_sample_rigorous/jackknife_error/cylinderDepthTest/normsig_' + zrange + $
			'_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_dzCoeff' + string_dzCoeff + '.fits'
		data = MRDFITS(inputFile, 1)

		frac_IPSF	= data.n_late_IPSF/data.n_tot_IPSF
		frac_IPQ	= data.n_late_IPQ/data.n_tot_IPQ

		normsig		= data.normsig		
		normsig_errors	= data.normsig_errors

		sigmaArray0_1 = [sigmaArray0_1, getSigmaRange_JKE(data, 0, 0)]
		sigmaArray1_2 = [sigmaArray1_2, getSigmaRange_JKE(data, 1, 1)]
		sigmaArray1_3 = [sigmaArray1_3, getSigmaRange_JKE(data, 1, 2)]
		sigmaArray3_5 = [sigmaArray3_5, getSigmaRange_JKE(data, 3, 4)]

;		sigmaArray0_1 = [sigmaArray0_1, getSigmaRange_JKE(data, 0, 0)]
;		sigmaArray1_2 = [sigmaArray1_2, getSigmaRange_JKE(data, 1, 1)]
;		sigmaArray1_3 = [sigmaArray1_3, getSigmaRange_JKE(data, 1, 2)]
;		sigmaArray3_5 = [sigmaArray3_5, getSigmaRange_JKE(data, 1, 4)]
	ENDFOREACH

        IF (i LE (nx-1)) THEN BEGIN
		xtitle = textoidl('\Deltaz / (0.005(1+z_{IP}))')
                xtickformat = ''
        ENDIF ELSE BEGIN
                xtitle = ''
                xtickformat = '(A1)'
        ENDELSE
        IF (i mod nx EQ 0) THEN BEGIN
                ytitle = textoidl('S/N over R_{proj} (Mpc)')
                ytickformat = ''
        ENDIF ELSE BEGIN
                ytitle = ''
                ytickformat = '(A1)'
        ENDELSE
	
	IF i EQ 3 THEN XYOUTS, xmax, 1.025*ymax,'S/N vs Cylinder Depth', ALIGNMENT=0.5, charsize=1.5

        pos = [ xm + dx*FLOAT(i mod nx), ym + dy*FLOAT(floor(FLOAT(i)/nx)), $
                xm + dx*(1. + FLOAT(i mod nx)), ym + dy*(1. + FLOAT(floor(FLOAT(i)/nx))) ]
        PRINT, pos

	PLOT, dz_coeffArray, sigmaArray0_1, xrange=[xmin,xmax], yrange=[ymin,ymax], LINESTYLE=0, position=pos, $
;	  xtitle=textoidl('\Deltaz / (0.005(1+z_{IP}))'), ytitle=textoidl('\sigma over range of projected radii'), 
	  xtickformat=xtickformat, ytickformat=ytickformat, xtitle=xtitle, ytitle=ytitle
	OPLOT, dz_coeffArray, sigmaArray1_2, LINESTYLE=1
	OPLOT, dz_coeffArray, sigmaArray1_3, LINESTYLE=2
	OPLOT, dz_coeffArray, sigmaArray3_5, LINESTYLE=3

	IF i EQ nx*ny-1 THEN BEGIN
		LEGEND, [textoidl('R_{proj} < 1 Mpc'), '1-2 Mpc', '1-3 Mpc', '3-5 Mpc'], $
			LINESTYLE=[0,1,2,3], BOX=0, /TOP, /RIGHT
	ENDIF

	XYOUTS, xmin+0.05*xr, ymin+0.925*yr, zlabels[i], ALIGNMENT=0.0

  ENDFOR
	IF (STRING(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
