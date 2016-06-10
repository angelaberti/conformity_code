; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; latefrac:	calculates late-type fractions
; BSE:		uses bootstrap error method

PRO IPhist_allz, outputFormat;, zmin, zmax;, dz_coeff, printEvery
	!P.FONT=0
	SFcolor = cgcolor('blue')
	Qcolor  = cgcolor('red')

  IPdatafiles = [ $
	'~/conformity/results/match_IP_sample_rigorous/unmatchedExample/matchMassOnly/matchedIPsampleMassOnlyFBF_PHI3.7.fits', $
	'~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', $
	'~/conformity/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits', $
	'~/conformity/results/single_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_singleMass.fits'] ; matched z only

  titles = [$;'Default Parameters', $
	'(c) Matched ' + textoidl('M_{\ast}') + ' Only', $
	'(d) Matched ' + textoidl('M_{\ast}') + ' & Redshift', $
	'(a) No Matching', $
	'(b) Matched Redshift Only']

;  tags = ['default', 'singleMass', 'var1.0dex']
;  tags = ['default', 'singleMass', 'conservative', 'conservative0.3dex', 'var1.0dex']
  tags = ['matchMassOnly', '', 'conservative', 'singleMass']

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	Rmax 		= 5.
	dRproj		= 1.
	n_annuli 	= CEIL(Rmax/dRproj)
	printEvery	= 100

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

	ERASE
;	!P.MULTI=[4,2,2]
	!P.MULTI=4
	!P.CHARSIZE=2.5
	charsz=1.25

	Rproj_array = [0.0, 0.5, 1.0+FINDGEN(Rmax)]	
	Rplot_array = [0.25, 0.75, 1.5+FINDGEN(Rmax)]	
;	Rproj_array = FLOAT(dRproj)*(findgen(n_annuli+1))	
;	PRINT, Rproj_array

	xmin = 8
	xmax = 12
;	xmin = 0.2
;	xmax = 1.0
	xr = xmax-xmin

	ymin = 0
	ymax = 9000
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

  FOR i=0,nx*ny-1 DO BEGIN
        IF (i LE (nx-1)) THEN BEGIN
		xtitle = 'Stellar Mass'
;		xtitle = 'z
                xtickformat = ''
        ENDIF ELSE BEGIN
                xtitle = ''
                xtickformat = '(A1)'
        ENDELSE
        IF (i mod nx EQ 0) THEN BEGIN
                ytitle = 'Number'
                ytickformat = ''
        ENDIF ELSE BEGIN
                ytitle = ''
                ytickformat = '(A1)'
        ENDELSE

        pos = [ xm + dx*FLOAT(i mod nx), ym + dy*FLOAT(FLOOR(FLOAT(i)/nx)), $
                xm + dx*(1. + FLOAT(i mod nx)), ym + dy*(1. + FLOAT(FLOOR(FLOAT(i)/nx))) ]
;	PRINT, pos

;  FOR i=0,n_elements(IPdatafiles)-1 DO BEGIN
	PRINT, ''
	PRINT, IPdatafiles[i]
	IPstats, IPdatafiles[i]

	SFcolor = cgcolor('blue')
	Qcolor  = cgcolor('red')

	dataIP_allz  = MRDFITS(IPdatafiles[i], 1)

	; eliminate data with targ_weight < 1
        dataIP_allz  = dataIP_allz[where((dataIP_allz.targ_weight GE 1.) AND (dataIP_allz.IP EQ 1))]

	dataAll = dataAll_allz
	dataAllSF = dataAll[where(dataAll.SFQ EQ 1)]
	dataAllQ  = dataAll[where(dataAll.SFQ EQ 0)]

	dataIP = dataIP_allz
	dataIPSF = dataIP[where(dataIP.SFQ EQ 1)]
	dataIPQ  = dataIP[where(dataIP.SFQ EQ 0)]

	zlabel = 'z=[' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']'

	IF (i EQ 1) THEN data = MRDFITS('~/conformity/results/match_IP_sample_rigorous/latefrac_allz_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE_halo12.fits', 1) ELSE $
	  data = MRDFITS('~/conformity/results/match_IP_sample_rigorous/unmatchedExample/latefrac_' + tags[i] + '_allz_dR1Mpc_BSE.fits', 1)

	b=0.25

	PLOTHIST, dataAll.mstar, bin=b, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle=xtitle, ytitle=ytitle, xtickformat=xtickformat, ytickformat=ytickformat, POSITION=pos, $
;	PLOTHIST, dataAll.mstar, bin=b, xtitle=xtitle, ytitle=ytitle, xtickformat=xtickformat, ytickformat=ytickformat, POSITION=pos, $
	  /NODATA
	XYOUTS, xmin+0.05*xr, ymin+0.90*yr, titles[i], ALIGNMENT=0.0, CHARSIZE=charsz

	PLOTHIST, dataAllSF.mstar, bin=b, COLOR=SFcolor, LINESTYLE=2, /OVERPLOT
	PLOTHIST, dataAllQ.mstar, bin=b, COLOR=Qcolor, LINESTYLE=2, /OVERPLOT
	PLOTHIST, dataIPSF.mstar, bin=b, COLOR=SFcolor, /OVERPLOT
	PLOTHIST, dataIPQ.mstar, bin=b, COLOR=Qcolor, /OVERPLOT

	IF (i EQ 2) THEN BEGIN
;		XYOUTS, xmin+0.15*xr, ymin+0.825*yr, zlabel, ALIGNMENT=0.0, CHARSIZE=charsz
		XYOUTS, xmin+0.05*xr, ymin+0.05*yr, zlabel, ALIGNMENT=0.0, CHARSIZE=charsz
;		LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], COLOR=[SFcolor, Qcolor], CHARSIZE=charsz, BOX=0, /CENTER, /RIGHT
;		LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], COLOR=[SFcolor, Qcolor], CHARSIZE=charsz, BOX=0, POS=[3,0.85]
	ENDIF
  ENDFOR

	IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
