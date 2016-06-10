; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; latefrac:	calculates late-type fractions
; BSE:		uses bootstrap error method

PRO test, outputFormat;, zmin, zmax;, dz_coeff, printEvery
	!P.FONT=0
	SFcolor = cgcolor('blue')
	Qcolor  = cgcolor('red')

  IPdatafiles = [ $
	'~/conformity/results/match_IP_sample_rigorous/unmatchedExample/noMatching/unmatchedSampleFBF_PHI3.7.fits', $
	'~/conformity/results/match_IP_sample_rigorous/unmatchedExample/matchMassOnly/matchedIPsampleMassOnlyFBF_PHI3.7.fits', $
	'~/conformity/results/match_IP_sample_rigorous/unmatchedExample/matchRedshiftOnly/matchedIPsampleRedshiftOnlyFBF_PHI3.7.fits', $
	'~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits']
;	'~/conformity/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits', $

  titles = [$;'Default Parameters', $
	'(a) No Matching', $
	'(b) Matched ' + textoidl('M_{\ast}') + ' Only', $
	'(c) Matched Redshift Only', $
	'(d) Matched ' + textoidl('M_{\ast}') + ' & Redshift']

  tags = ['noMatching', 'matchMassOnly', 'matchRedshiftOnly', '']

	dz_coeff = 2.0
	zmin	= 0.2
	zmax	= 1.0

	dRproj	= 1.

	IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
	IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

	ERASE
;	!P.MULTI=12
;	!P.MULTI=[12,3,4]
	!P.CHARSIZE=1
	charsz=1.25

	z_low = 0.2
	z_high = 1.0

	NX=3
	NY=4

	IF (string(outputFormat) eq 'ps') THEN BEGIN
		axisColor = cgColor('Black')
		PS_OPEN, '~/latex/figures/unmatchedIPsampleCompare_wHist', THICK=5, /ENCAP
		DEVICE, /INCH, XS=6, YS=8, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
		axisColor = cgColor('White')
		SET_PLOT, 'X'
		THICK=1
	ENDELSE

	ERASE
;	!P.MULTI=NX*NY

; OUTER MARGINS
        XMouter=0.05
        YMouter=0.10
; INNER MARGINS
        XMinner = 0.10
        YMinner = 0.02

        DX = (1 - 1.5*XMouter)/NX
        DY = (1 - 2.0*YMouter)/NY

  POSITION_ARRAY = []
  FOR i=0,NX*NY-1 DO BEGIN
        POS = [ XMouter + DX*FLOAT(i mod NX) + XMinner/2., $
                YMouter + DY*FLOAT(FLOOR(FLOAT(i)/NX)) + YMinner/2., $
                XMouter + DX*(1. + FLOAT(i mod NX)) - XMinner/2., $
                YMouter + DY*(1. + FLOAT(FLOOR(FLOAT(i)/NX))) - YMinner/2. ]
	PRINT, i, POS
	POSITION_ARRAY = [POSITION_ARRAY, POS]
  ENDFOR

;  FOR i=0,NX*NY-1 DO BEGIN
;	IF (i LE (NX-1)) THEN BEGIN
		xtitle = 'xtitle'
		xtickformat = ''
;	ENDIF ELSE BEGIN
;		xtitle = ''
;		xtickformat = '(A1)'
;	ENDELSE
;	IF (i mod NX EQ 0) THEN BEGIN
		ytitle = 'ytitle'
		ytickformat = ''
;	ENDIF ELSE BEGIN
;		ytitle = ''
;		ytickformat = '(A1)'
;	ENDELSE

 FOR i=0,n_elements(IPdatafiles)-1 DO BEGIN

	dataAll_allz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)
; eliminate data with targ_weight < 1
	dataAll_allz = dataAll_allz[where(dataAll_allz.targ_weight GE 1.)]

	PRINT, ''
	PRINT, IPdatafiles[i]
	IPstats, IPdatafiles[i]

	SFcolor = cgcolor('blue')
	Qcolor  = cgcolor('red')

	dataIP_allz  = MRDFITS(IPdatafiles[i], 1)

; eliminate data with targ_weight < 1
	dataIP_allz  = dataIP_allz[where((dataIP_allz.targ_weight GE 1.) AND (dataIP_allz.IP EQ 1))]

	dataAll = dataAll_allz
	dataIP	= dataIP_allz

	zlabel = 'z=[' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']'

	IF (i EQ 3) THEN datafile = '~/conformity/results/match_IP_sample_rigorous/unmatchedExample/halo12test/latefrac_allz_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE_halo12.fits' ELSE $
	  datafile = '~/conformity/results/match_IP_sample_rigorous/unmatchedExample/halo12test/latefrac_' + tags[i] + '_allz_dR1Mpc_BSE.fits'
	data = MRDFITS(datafile, 1)
	PRINT, datafile

	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPSF = data.n_late_IPSF
	errors_IPSF = data.errors_IPSF
	frac_IPSF   = n_late_IPSF/n_tot_IPSF

	n_tot_IPQ  = data.n_tot_IPQ
	n_late_IPQ = data.n_late_IPQ
	errors_IPQ = data.errors_IPQ
	frac_IPQ   = n_late_IPQ/n_tot_IPQ

; 1st PLOT
	Rmax = 9.

	Rproj_array = [0.0, 0.5, 1.0+FINDGEN(Rmax)]
	Rplot_array = [0.25, 0.75, 1.5+FINDGEN(Rmax-1)]

	PRINT, Rproj_array
	PRINT, Rplot_array

	xmin = 0
	xmax = 9
	xr = xmax-xmin

	ymin = 0.71
	ymax = 0.88
	yr = ymax-ymin

	PRINT, POSITION_ARRAY[3*i]

	PLOT, Rplot_array, frac_IPSF, xrange=[xmin,xmax], yrange=[ymin,ymax], $
		xtitle=xtitle, ytitle=ytitle, xtickformat=xtickformat, ytickformat=ytickformat;, $
;		POSITION=POSITION_ARRAY[3*i], /NODATA
;	XYOUTS, xmin+0.05*xr, ymin+0.90*yr, titles[i], ALIGNMENT=0.0, CHARSIZE=charsz

	OPLOT, Rplot_array, frac_IPSF, LINESTYLE=0, COLOR=SFcolor
	ERRPLOT, Rplot_array, frac_IPSF-errors_IPSF, frac_IPSF+errors_IPSF, COLOR=SFcolor

	OPLOT, Rplot_array, frac_IPQ, LINESTYLE=2, COLOR=Qcolor
	ERRPLOT, Rplot_array, frac_IPQ-errors_IPQ, frac_IPQ+errors_IPQ, COLOR=Qcolor

; 2nd PLOT: MASS HISTOGRAMS
	dataAllSF = dataAll[where(dataAll.SFQ EQ 1)]
	dataAllQ  = dataAll[where(dataAll.SFQ EQ 0)]

	dataIPSF = dataIP[where(dataIP.SFQ EQ 1)]
	dataIPQ  = dataIP[where(dataIP.SFQ EQ 0)]

	xmin = 9.1
	xmax = 11.6
	xr = xmax-xmin
	ymin = 0
	ymax = 4000
	yr = ymax-ymin
	b = 0.25

	PLOTHIST, dataIP.mstar, bin=b, xrange=[xmin,xmax], yrange=[ymin,ymax], color=axisColor, $
		xtitle='Stellar Mass', ytitle='Number', POSITION=POSITION_ARRAY[3*i+1]
;	XYOUTS, xmin+0.05*xr, ymin+0.90*yr, titles[i], ALIGNMENT=0.0, CHARSIZE=charsz

;	PLOTHIST, dataAllSF.mstar, bin=b, COLOR=SFcolor, LINESTYLE=2, /OVERPLOT
;	PLOTHIST, dataAllQ.mstar, bin=b, COLOR=Qcolor, LINESTYLE=2, /OVERPLOT
	PLOTHIST, dataIPSF[where((dataIPSF.mstar GE 9.1) AND (dataIPSF.mstar LE 11.6))].mstar, bin=b, COLOR=SFcolor, /OVERPLOT
	PLOTHIST, dataIPQ[where((dataIPQ.mstar GE 9.1) AND (dataIPQ.mstar LE 11.6))].mstar, bin=b, COLOR=Qcolor, /OVERPLOT

; 3rd PLOT: REDSHIFT HISTROGRAMS
	xmin = 0.175
	xmax = 1.025
	xr = xmax-xmin
	ymin = 0
	ymax = 1200
	yr = ymax-ymin
	b = 0.05

	PLOTHIST, dataIP.zprimus, bin=b, xrange=[xmin,xmax], yrange=[ymin,ymax], color=axisColor, $
		xtitle='z', ytitle='Number', POSITION=POSITION_ARRAY[3*i+2]

;	PLOTHIST, dataAllSF.zprimus, bin=b, COLOR=SFcolor, LINESTYLE=2, /OVERPLOT
;	PLOTHIST, dataAllQ.zprimus, bin=b, COLOR=Qcolor, LINESTYLE=2, /OVERPLOT
	PLOTHIST, dataIPSF.zprimus, bin=b, COLOR=SFcolor, /OVERPLOT
	PLOTHIST, dataIPQ.zprimus, bin=b, COLOR=Qcolor, /OVERPLOT

;	IF (i EQ 2) THEN $
;		XYOUTS, xmin+0.15*xr, ymin+0.825*yr, zlabel, ALIGNMENT=0.0, CHARSIZE=charsz
;		XYOUTS, xmin+0.95*xr, ymin+0.05*yr, zlabel, ALIGNMENT=1.0, CHARSIZE=charsz
;	IF (i EQ 2) THEN $
;		LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], COLOR=[SFcolor, Qcolor], CHARSIZE=charsz, BOX=0, /CENTER, /RIGHT
;		LEGEND, ['', 'Late-type IP', 'Early-type IP'], LINESTYLE=[1,0,2], COLOR=[cgColor('white'), SFcolor, Qcolor], CHARSIZE=charsz, BOX=0, /TOP, /RIGHT
  ENDFOR

	IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
