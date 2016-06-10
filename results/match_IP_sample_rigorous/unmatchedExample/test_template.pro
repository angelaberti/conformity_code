PRO test_template, outputFormat
	!P.FONT=0

	SFcolor = cgcolor('blue')
	Qcolor  = cgcolor('red')

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	Rmax 		= 5.
	dRproj		= 1.
	n_annuli 	= CEIL(Rmax/dRproj)
	printEvery	= 100

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

	NX=1
	NY=3

	ERASE
	!P.MULTI=NX*NY
	!P.CHARSIZE=1.75
	charsz=1.4

	Rproj_array = FLOAT(dRproj)*(findgen(n_annuli+1))	
;	PRINT, Rproj_array

	xmin = 0
	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = 0.71
	ymax = 0.88
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

	dataAll_allz = MRDFITS('~/results/zerodSFQ_all_cart.fits', 1)
	; eliminate data with targ_weight < 1
        dataAll_allz = dataAll_allz[where(dataAll_allz.targ_weight GE 1.)]

; OUTER MARGINS
	XMouter=0.10
	YMouter=0.10
; INNER MARGINS
	XMinner = 0.02
	YMinner = 0.02

	DX = (1 - 1.5*XMouter)/NX
	DY = (1 - 2.0*YMouter)/NY

  FOR i=0,NX*NY-1 DO BEGIN

; EQUAL SIZE GRID
	POS = [ XMouter + DX*FLOAT(i mod NX) + XMinner/2., $
		YMouter + DY*FLOAT(FLOOR(FLOAT(i)/NX)) + YMinner/2., $
		XMouter + DX*(1. + FLOAT(i mod NX)) - XMinner/2., $
		YMouter + DY*(1. + FLOAT(FLOOR(FLOAT(i)/NX))) - YMinner/2. ]

; UNEQUAL HEIGHT UPPER AND LOWER PANELS (2 ROWS)
;	S = 0.67
;  IF (i LE NX-1) THEN BEGIN
;	POS = [ XMouter + DX*FLOAT(i mod NX) + XMinner/2., $
;		YMouter + DY*FLOAT(FLOOR(FLOAT(i)/NX)) + YMinner/2., $
;		XMouter + DX*(1. + FLOAT(i mod NX)) - XMinner/2., $
;		S*(YMouter + DY*(1. + FLOAT(FLOOR(FLOAT(i)/NX)))) - YMinner/2. ]
;  ENDIF ELSE BEGIN
;	POS = [ XMouter + DX*FLOAT(i mod NX) + XMinner/2., $
;		S*(YMouter + DY*FLOAT(FLOOR(FLOAT(i)/NX))) + YMinner/2., $
;		XMouter + DX*(1. + FLOAT(i mod NX)) - XMinner/2., $
;		YMouter + DY*(1. + FLOAT(FLOOR(FLOAT(i)/NX))) - YMinner/2. ]
;  ENDELSE
	PRINT, i, POS
	
        IF (i LE (NX-1)) THEN BEGIN
                xtitle = 'xtitle'
                xtickformat = ''
        ENDIF ELSE BEGIN
                xtitle = ''
                xtickformat = '(A1)'
        ENDELSE
        IF (i mod NX EQ 0) THEN BEGIN
                ytitle = 'ytitle'
                ytickformat = ''
        ENDIF ELSE BEGIN
                ytitle = ''
                ytickformat = '(A1)'
        ENDELSE
	
;	PLOT, Rproj_array+0.5*dRproj, frac_IPSF, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle=xtitle, ytitle=ytitle, xtickformat=xtickformat, ytickformat=ytickformat, POSITION=POS, $
	PLOT, FINDGEN(10), FINDGEN(10), xtitle=xtitle, ytitle=ytitle, xtickformat=xtickformat, ytickformat=ytickformat, POSITION=POS, /NODATA
	XYOUTS, 5, 5, i+1, ALIGNMENT=0.0, CHARSIZE=charsz
  ENDFOR

	IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
