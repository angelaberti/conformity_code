; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; latefrac:	calculates late-type fractions
; BSE:		uses bootstrap error method

PRO latefracplot_BSE_allz, outputFormat;, zmin, zmax;, dz_coeff, printEvery
	!P.FONT=0
	SFcolor = cgcolor('blue')
	Qcolor  = cgcolor('red')

  IPdatafiles = [ $
; MATCHED MASS ONLY
	'~/conformity/results/match_IP_sample_rigorous/unmatchedExample/matchMassOnly/matchedIPsampleMassOnlyFBF_PHI3.7.fits', $
; MATCHED Z AND MASS
	'~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', $
; NO MATCHING
	'~/conformity/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits', $
;	'~/conformity/results/match_IP_sample_rigorous/unmatchedExample/noMatching/unmatchedSampleFBF_PHI3.7.fits', $
; MATCHED Z ONLY
	'~/conformity/results/match_IP_sample_rigorous/unmatchedExample/matchRedshiftOnly/matchedIPsampleRedshiftOnlyFBF_PHI3.7.fits' $
  ]

  datafiles = [ $
; MATCHED MASS ONLY
	'~/conformity/results/match_IP_sample_rigorous/unmatchedExample/latefrac_matchMassOnly_allz_dR1Mpc_BSE.fits', $
; MATCHED Z AND MASS
	'~/conformity/results/match_IP_sample_rigorous/latefrac_allz_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits', $
; NO MATCHING
	'~/conformity/results/match_IP_sample_rigorous/unmatchedExample/latefrac_conservative_allz_dR1Mpc_BSE.fits', $
;	'~/conformity/results/conservative_mass_cutoff/latefrac_data_hist/latefrac_0.2_1.0_targ_weight_zerodSFQ_IP_dz2.0_dm0.0.fits', $
; MATCHED Z ONLY
	'~/conformity/results/match_IP_sample_rigorous/unmatchedExample/latefrac_matchRedshiftOnly_allz_dR1Mpc_BSE.fits' $
  ]

  titles = [$;'Default Parameters', $
	textoidl('(c) Matched M_{*} Only'), $
	textoidl('(d) Matched M_{*} & Redshift'), $
	'(a) No Matching', $
	'(b) Matched Redshift Only']

;  tags = ['default', 'singleMass', 'var1.0dex']
;  tags = ['default', 'singleMass', 'conservative', 'conservative0.3dex', 'var1.0dex']
;  tags = ['matchMassOnly', '', 'noMatching', 'matchRedshiftOnly']
  tags = ['matchMassOnly', '', 'conservative', 'matchRedshiftOnly']

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	Rmax 		= 15.
	dRproj		= 1.
	n_annuli 	= CEIL(Rmax/dRproj)
	printEvery	= 100

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

	ERASE
	!P.MULTI=4
	!P.CHARSIZE=2.25
	charsz=1.1

	Rproj_array = FINDGEN(Rmax)+1
	Rplot_array = FINDGEN(Rmax) + 0.5*dRproj
;	Rproj_array = [0.0, 0.5, 1.0+FINDGEN(Rmax)]
;	Rplot_array = [0.25, 0.75, 1.5+FINDGEN(Rmax-1)]

	PRINT, Rproj_array
	PRINT, Rplot_array

	xmin = 0
	xmax = 15
;	xmax = n_annuli*dRproj
	xr = xmax-xmin

;	ymin = 0.6
;	ymax = 0.89
	ymin = 0.72
	ymax = 0.81
	yr = ymax-ymin

	z_low = 0.2
	z_high = 1.0

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/latex/figures/unmatchedIPsampleCompare', THICK=5, /ENCAP
		DEVICE, /INCH, XS=8, YS=8, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	dataAll_allz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)
; eliminate data with targ_weight < 1
        dataAll_allz = dataAll_allz[where(dataAll_allz.targ_weight GE 1.)]

	nx=2.
	ny=2.

	xm=0.10
	ym=0.10

	dx = (1 - 1.5*xm)/nx
	dy = (1 - 2.0*ym)/ny

	dm = 0.005

  FOR i=0,nx*ny-1 DO BEGIN
        IF (i LE (nx-1)) THEN BEGIN
;		xtitle = 'Projected Radius (Mpc)'
		xtitle = textoidl('R_{proj} (Mpc)')
                xtickformat = ''
        ENDIF ELSE BEGIN
                xtitle = ''
                xtickformat = '(A1)'
        ENDELSE
        IF (i mod nx EQ 0) THEN BEGIN
;		ytitle = 'Late-type Fraction'
		ytitle = textoidl('f_{late}')
                ytickformat = ''
        ENDIF ELSE BEGIN
                ytitle = ''
                ytickformat = '(A1)'
        ENDELSE

        pos = [ xm + dx*FLOAT(i mod nx) + dm, ym + dy*FLOAT(FLOOR(FLOAT(i)/nx)) + dm, $
                xm + dx*(1. + FLOAT(i mod nx)) - dm, ym + dy*(1. + FLOAT(FLOOR(FLOAT(i)/nx))) - dm ]
;	PRINT, pos

;  FOR i=0,n_elements(IPdatafiles)-1 DO BEGIN
	PRINT, ''
	PRINT, IPdatafiles[i]
	IPstats, IPdatafiles[i]

	SFcolor = cgcolor('blue')
	Qcolor  = cgcolor('red')

;	IF (i EQ 5) THEN dataIP_allz = MRDFITS('~/conformity/results/match_IP_sample_rigorous/latefrac_allz_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits', 1) ELSE $
;	IF (i EQ 5) THEN dataIP_allz = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1) ELSE $
	  dataIP_allz  = MRDFITS(IPdatafiles[i], 1)

	; eliminate data with targ_weight < 1
        dataIP_allz  = dataIP_allz[where((dataIP_allz.targ_weight GE 1.) AND (dataIP_allz.IP EQ 1))]

	dataAll = dataAll_allz
	dataIP = dataIP_allz

	zlabel = 'z = [' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']'

	data = MRDFITS(datafiles[i], 1)

	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPSF = data.n_late_IPSF
	errors_IPSF = data.errors_IPSF
	frac_IPSF   = n_late_IPSF/n_tot_IPSF

	n_tot_IPQ  = data.n_tot_IPQ
	n_late_IPQ = data.n_late_IPQ
	errors_IPQ = data.errors_IPQ
	frac_IPQ   = n_late_IPQ/n_tot_IPQ

;	PLOT, Rplot_array, frac_IPSF, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', POSITION=pos, $
	PLOT, Rplot_array, frac_IPSF, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle=xtitle, ytitle=ytitle, xtickformat=xtickformat, ytickformat=ytickformat, POSITION=pos, $
	  /NODATA
;	  title = 'Field-by-Field Matched SF & Q IP Samples ' + textoidl('(\Deltaz=2.0)'), /NODATA
	IF (i EQ 0) THEN XYOUTS, xmin+0.05*xr, ymin+0.90*yr, textoidl('(c) Matched M_{*} Only'), ALIGNMENT=0.0, CHARSIZE=charsz
	IF (i EQ 1) THEN XYOUTS, xmin+0.05*xr, ymin+0.90*yr, textoidl('(d) Matched M_{*} & Redshift'), ALIGNMENT=0.0, CHARSIZE=charsz
	IF (i GT 1) THEN XYOUTS, xmin+0.05*xr, ymin+0.90*yr, titles[i], ALIGNMENT=0.0, CHARSIZE=charsz

	OPLOT, Rplot_array, frac_IPSF, LINESTYLE=0, COLOR=SFcolor
	ERRPLOT, Rplot_array, frac_IPSF-errors_IPSF, frac_IPSF+errors_IPSF, COLOR=SFcolor

	OPLOT, Rplot_array, frac_IPQ, LINESTYLE=2, COLOR=Qcolor
	ERRPLOT, Rplot_array, frac_IPQ-errors_IPQ, frac_IPQ+errors_IPQ, COLOR=Qcolor

;	XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'SF IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 1)])), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Q IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 0)])), ALIGNMENT=1.0

	zStatsSF = textoidl('z_{med,SF} = ') +	decimal(MEDIAN(dataIP[WHERE(dataIP.SFQ EQ 1)].zprimus),2)
	zStatsQ  = textoidl('z_{med,Q} = ') +	decimal(MEDIAN(dataIP[WHERE(dataIP.SFQ EQ 0)].zprimus),2)
	mStatsSF = textoidl('M_{*med,SF} = ') +	decimal(MEDIAN(dataIP[WHERE(dataIP.SFQ EQ 1)].mstar),2)
	mStatsQ  = textoidl('M_{*med,Q} = ') + 	decimal(MEDIAN(dataIP[WHERE(dataIP.SFQ EQ 0)].mstar),2)

	XYOUTS, xmin+0.40*xr, ymin+0.15*yr, zStatsSF, ALIGNMENT=1.0, CHARSIZE=charsz
	XYOUTS, xmin+0.40*xr, ymin+0.075*yr, zStatsQ, ALIGNMENT=1.0, CHARSIZE=charsz
	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, mStatsSF, ALIGNMENT=1.0, CHARSIZE=charsz
	XYOUTS, xmin+0.95*xr, ymin+0.075*yr, mStatsQ, ALIGNMENT=1.0, CHARSIZE=charsz

	IF (i EQ 2) THEN $
;		XYOUTS, xmin+0.15*xr, ymin+0.825*yr, zlabel, ALIGNMENT=0.0, CHARSIZE=charsz
		XYOUTS, xmin+0.95*xr, ymin+0.90*yr, zlabel, ALIGNMENT=1.0, CHARSIZE=charsz
	IF (i EQ 3) THEN $
;		LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], COLOR=[SFcolor, Qcolor], CHARSIZE=0.9*charsz, BOX=0, POSITION=[4,0.795]
		LEGEND, ['', '', 'SF IP', 'Q IP'], LINESTYLE=[1,1,0,2], COLOR=[cgColor('white'), cgColor('white'), SFcolor, Qcolor], CHARSIZE=charsz, BOX=0, /TOP, /RIGHT
  ENDFOR

	IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
