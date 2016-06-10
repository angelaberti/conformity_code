; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; latefrac:	calculates late-type fractions
; BSE:		uses bootstrap error method

PRO latefracplot_BSE_allz_wHist, outputFormat
	!P.FONT=0
	SFcolor = cgcolor('blue')
	Qcolor  = cgcolor('red')

  IPdatafiles = [ $
; MATCHED Z AND MASS
	'~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', $
; MATCHED MASS ONLY
	'~/conformity/results/match_IP_sample_rigorous/unmatchedExample/matchMassOnly/matchedIPsampleMassOnlyFBF_PHI3.7.fits', $
; MATCHED Z ONLY
	'~/conformity/results/match_IP_sample_rigorous/unmatchedExample/matchRedshiftOnly/matchedIPsampleRedshiftOnlyFBF_PHI3.7.fits', $
; NO MATCHING
;	'~/conformity/results/match_IP_sample_rigorous/unmatchedExample/noMatching/unmatchedSampleFBF_PHI3.7.fits', $
	'~/conformity/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits' $
  ]

  datafiles = [ $
; MATCHED Z AND MASS
	'~/conformity/results/match_IP_sample_rigorous/latefrac_allz_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits', $
; MATCHED MASS ONLY
	'~/conformity/results/match_IP_sample_rigorous/unmatchedExample/latefrac_matchMassOnly_allz_dR1Mpc_BSE.fits', $
; MATCHED Z ONLY
	'~/conformity/results/match_IP_sample_rigorous/unmatchedExample/latefrac_matchRedshiftOnly_allz_dR1Mpc_BSE.fits', $
; NO MATCHING
;	'~/conformity/results/conservative_mass_cutoff/latefrac_data_hist/latefrac_0.2_1.0_targ_weight_zerodSFQ_IP_dz2.0_dm0.0.fits', $
	'~/conformity/results/match_IP_sample_rigorous/unmatchedExample/latefrac_conservative_allz_dR1Mpc_BSE.fits' $
  ]

  titles = [$;'Default Parameters', $
	'Matched ' + textoidl('M_{\ast}') + ' & Redshift', $
	'Matched ' + textoidl('M_{\ast}') + ' Only', $
	'Matched Redshift Only', $
	'No Matching']

	dz_coeff = 2.0

	Rmax 		= 15.
	dRproj		= 1.
	n_annuli 	= CEIL(Rmax/dRproj)
	printEvery	= 100

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

	ERASE
	!P.MULTI=12
	!P.CHARSIZE=1.75
	charsz=0.75

	Rproj_array = FINDGEN(Rmax)+1
	Rplot_array = FINDGEN(Rmax) + 0.5*dRproj
;	Rproj_array = [0.0, 0.5, 1.0+FINDGEN(Rmax)]
;	Rplot_array = [0.25, 0.75, 1.5+FINDGEN(Rmax-1)]

;	PRINT, Rproj_array
;	PRINT, Rplot_array

	xmin = 0
	xmax = 15
;	xmax = n_annuli*dRproj
	xr = xmax-xmin

;	ymin = 0.6
;	ymax = 0.89
	ymin = 0.71
	ymax = 0.88
	yr = ymax-ymin

	z_low = 0.2
	z_high = 1.0

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/latex/figures/unmatchedIPsampleCompare_wHist', THICK=5, /ENCAP
		DEVICE, /INCH, XS=6, YS=8, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	nx=3.
	ny=4.

	xm=0.10
	ym=0.10

	dx = (1 - 1.5*xm)/nx
	dy = (1 - 2.0*ym)/ny

	dm = 0.005

	posArray = []
  FOR i=0,nx*ny-1 DO BEGIN
        pos = [ xm + dx*FLOAT(i mod nx) + dm, ym + dy*FLOAT(FLOOR(FLOAT(i)/nx)) + dm, $
                xm + dx*(1. + FLOAT(i mod nx)) - dm, ym + dy*(1. + FLOAT(FLOOR(FLOAT(i)/nx))) - dm ]
	posArray = [[posArray],[pos]]
  ENDFOR
;	PRINT, posArray

	dataAll_allz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)
        dataAll_allz = dataAll_allz[where(dataAll_allz.targ_weight GE 1.)]

  FOR i=0,n_elements(IPdatafiles)-1 DO BEGIN
	PRINT, IPdatafiles[i]
	IPstats, IPdatafiles[i]
	dataIP_allz = MRDFITS(IPdatafiles[i],1)
	dataIP_allz = dataIP_allz[where((dataIP_allz.targ_weight GE 1.) AND (dataIP_allz.IP EQ 1))]

	dataAll = dataAll_allz
	dataIP = dataIP_allz

	dataIPSF = dataIP[WHERE(dataIP.SFQ EQ 1)]
	dataIPQ  = dataIP[WHERE(dataIP.SFQ EQ 0)]

	zlabel = 'z=[0.20, 1.00]'
;	zlabel = 'z=[' + decimal(zmin,2) + ', ' + decimal(zmax,2) + ']'

	data = MRDFITS(datafiles[i], 1)

	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPSF = data.n_late_IPSF
	errors_IPSF = data.errors_IPSF
	frac_IPSF   = n_late_IPSF/n_tot_IPSF

	n_tot_IPQ  = data.n_tot_IPQ
	n_late_IPQ = data.n_late_IPQ
	errors_IPQ = data.errors_IPQ
	frac_IPQ   = n_late_IPQ/n_tot_IPQ

  IF (i EQ 0) THEN BEGIN
	xtitle_left	= 'Projected Radius (Mpc)'
	xtitle_mid	= 'z'
	xtitle_right	= 'Stellar Mass'
	xtickformat	= ''
  ENDIF ELSE BEGIN
	xtitle_left	= ''
	xtitle_mid	= ''
	xtitle_right	= ''
	xtickformat	= '(A1)' 
  ENDELSE
	ytitle		= 'Late-type Fraction'
	ytickformat	= ''

;	PLOT, Rplot_array, frac_IPSF, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', POSITION=posArray[*,nx*i], /NODATA
	PLOT, Rplot_array, frac_IPSF, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle=xtitle_left, ytitle=ytitle, xtickformat=xtickformat, ytickformat=ytickformat, POSITION=posArray[*,nx*i], $
	  /NODATA
	XYOUTS, xmin+0.05*xr, ymin+0.90*yr, titles[i MOD nx], ALIGNMENT=0.0, CHARSIZE=charsz

	IF (i EQ (ny-1)) THEN $
;		XYOUTS, xmin+0.15*xr, ymin+0.825*yr, zlabel, ALIGNMENT=0.0, CHARSIZE=charsz
		XYOUTS, xmin+0.95*xr, ymin+0.05*yr, zlabel, ALIGNMENT=1.0, CHARSIZE=charsz
	IF (i EQ (ny-1)) THEN $
;		LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], COLOR=[cgColor('blue'), cgColor('red')], CHARSIZE=charsz, BOX=0, /CENTER, /RIGHT
		LEGEND, ['', 'Late-type IP', 'Early-type IP'], LINESTYLE=[1,0,2], COLOR=[cgColor('white'), cgColor('blue'), cgColor('red')], CHARSIZE=charsz, BOX=0, /TOP, /RIGHT

	OPLOT, Rplot_array, frac_IPSF, LINESTYLE=0, COLOR=cgColor('blue')
	ERRPLOT, Rplot_array, frac_IPSF-errors_IPSF, frac_IPSF+errors_IPSF, COLOR=cgColor('blue')

	OPLOT, Rplot_array, frac_IPQ, LINESTYLE=2, COLOR=cgColor('red')
	ERRPLOT, Rplot_array, frac_IPQ-errors_IPQ, frac_IPQ+errors_IPQ, COLOR=cgColor('red')

;; Z HISTOGRAM
	bin = 0.05
	zmin = 0.2
	zmax = 1.0
;	xrange_z = [zmin-bin/2,zmax+bin/2]
	xrange_z = [zmin,zmax]
	yrange_z = [0, 1400]

	zbins = zmin+bin*FINDGEN(CEIL((zmax-zmin)/bin))
	PRINT, 'z bins: ', zbins

;	ydataSF = HISTOGRAM(dataIPSF.zprimus, bin=bin)
;	ydataQ  = HISTOGRAM(dataIPQ.zprimus, bin=bin)

	ydataSF = []
	ydataQ	= []
	zPoints = []
    FOR j=0,CEIL((zmax-zmin)/bin) DO BEGIN
	z_low  = zmin+j*bin
	z_high = zmin+(j+1)*bin 	
;	z_low  = MIN(dataIP.zprimus)+j*bin
;	z_high = MIN(dataIP.zprimus)+(j+1)*bin 	
	zPoints = [zPoints, z_low]
	ydataSF = [ydataSF, n_elements(dataIPSF[WHERE((dataIPSF.zprimus GE z_low) AND (dataIPSF.zprimus LT z_high), /NULL)]) ]
	ydataQ  = [ydataQ, n_elements(dataIPQ[WHERE((dataIPQ.zprimus GE z_low) AND (dataIPQ.zprimus LT z_high), /NULL)]) ]
    ENDFOR
	PRINT, TRANSPOSE([[zPoints],[ydataSF],[ydataQ]])

	PLOT, FINDGEN(10), xrange=xrange_z, yrange=yrange_z, $
	  POSITION=posArray[*,nx*i+1], xtitle=xtitle_mid, xtickformat=xtickformat, ytitle='', ytickformat='(A1)', /NODATA
	OPLOT, zPoints+0.5*bin, ydataSF, PSYM=10, COLOR=cgColor('blue')
	OPLOT, zPoints+0.5*bin, ydataQ, PSYM=10, COLOR=cgColor('red')

;	PLOTHIST, dataIPSF.zprimus, bin=bin, COLOR=cgColor('blue'), xrange=xrange_z, yrange=yrange_mass, $
;	  POSITION=posArray[*,nx*i+1], xtitle=xtitle_mid, xtickformat=xtickformat, ytitle='', ytickformat='(A1)'
;	PLOTHIST, dataIPQ.zprimus, bin=bin, COLOR=cgColor('red'), /OVERPLOT

;; MASS HISTOGRAM
	bin = 0.2
	Mmin = 9.1
	Mmax = 11.6
	xrange_mass = [Mmin,Mmax]
	yrange_mass = [0, 3000]

	massbins = Mmin+bin*FINDGEN(CEIL((Mmax-Mmin)/bin))
	PRINT, 'Mass bins: ', massbins

;    IF (i EQ 0) OR (i EQ 1) THEN BEGIN
;	dataIPSF = dataIPSF[WHERE((dataIPSF.mstar GE Mmin) AND (dataIPSF.mstar LE Mmax))]
;	dataIPQ  = dataIPQ[WHERE((dataIPQ.mstar GE Mmin) AND (dataIPQ.mstar LE Mmax))]
;    ENDIF

;  PRINT, i, MINMAX(dataIPSF.mstar), MINMAX(dataIPQ.mstar)

;	ydataSF = HISTOGRAM(dataIPSF.mstar, bin=bin)
;	ydataQ  = HISTOGRAM(dataIPQ.mstar, bin=bin)

	ydataSF = []
	ydataQ	= []
	massPoints = []
    FOR j=0,CEIL((Mmax-Mmin)/bin) DO BEGIN
	m_low  = Mmin+j*bin
	m_high = Mmin+(j+1)*bin 	
;	m_low  = MIN(dataIP.mstar)+j*bin
;	m_high = MIN(dataIP.mstar)+(j+1)*bin 	
	massPoints = [massPoints, m_low]
	ydataSF = [ydataSF, n_elements(dataIPSF[WHERE((dataIPSF.mstar GE m_low) AND (dataIPSF.mstar LT m_high), /NULL)]) ]
	ydataQ  = [ydataQ, n_elements(dataIPQ[WHERE((dataIPQ.mstar GE m_low) AND (dataIPQ.mstar LT m_high), /NULL)]) ]
    ENDFOR
	PLOT, FINDGEN(10), xrange=xrange_mass, yrange=yrange_mass, $
	  POSITION=posArray[*,nx*i+2], xtitle=xtitle_right, xtickformat=xtickformat, ytitle='', ytickformat='(A1)', /NODATA
	OPLOT, massPoints + 0.5*bin, ydataSF, PSYM=10, COLOR=cgColor('blue')
	OPLOT, massPoints + 0.5*bin, ydataQ, PSYM=10, COLOR=cgColor('red')

;	PLOTHIST, dataIPSF.mstar, bin=bin, COLOR=cgColor('blue'), xrange=xrange_mass, yrange=yrange_mass, $
;	  POSITION=posArray[*,nx*i+2], xtitle=xtitle_right, xtickformat=xtickformat, ytitle='', ytickformat='(A1)'
;	PLOTHIST, dataIPQ.mstar, bin=bin, COLOR=cgColor('red'), /OVERPLOT

  ENDFOR

	IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
