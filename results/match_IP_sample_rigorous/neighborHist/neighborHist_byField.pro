PRO neighborHist_byField, binwidth, outputFormat
	PHI = 3.7
	stringPHI = strtrim(string(PHI, format='(f20.1)'),2)

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	; read in data files
	dataIPallz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
	dataAllallz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

	Rmax 		= 15.
	dRproj		= 0.25
	n_annuli 	= CEIL(Rmax/dRproj)
	printEvery	= 100

	IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
	IF dRproj EQ 0.50 THEN string_dR = '_dR500kpc'
	IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

	; eliminate data with targ_weight < 1
	dataAllallz = dataAllallz[WHERE(dataAllallz.targ_weight GE 1.)]
	dataIPallz  = dataIPallz[WHERE(dataIPallz.targ_weight GE 1.)]

	dataIPSF = dataIPallz[WHERE(dataIPallz.SFQ EQ 1)]
	dataIPQ	 = dataIPallz[WHERE(dataIPallz.SFQ EQ 0)]
	dataSF	 = dataAllallz[WHERE(dataAllallz.SFQ EQ 1)]
	dataQ	 = dataAllallz[WHERE(dataAllallz.SFQ EQ 0)]

;	!P.MULTI=[10,2,5]
;	!P.MULTI=[5,1,5]
	ERASE
	!P.MULTI=10
	!P.CHARSIZE=1.75
	charsz=1

	Rproj_array = float(dRproj)*(findgen(n_annuli+1))	
;	PRINT, Rproj_array

	xmin = 8
	xmax = 11.5
	xr = xmax-xmin

	ymin = 0.5
	ymax = 1
	yr = ymax-ymin

	posArray = grid_array(2,5)

	!P.FONT=0
	IF (string(outputFormat) EQ 'ps') THEN BEGIN
		PS_OPEN, '~/latex/figures/neighborHist_byField', THICK=5, /ENCAP
		DEVICE, /INCH, XS=10, YS=10
	ENDIF ELSE BEGIN
		SET_PLOT, 'X'
		THICK=1
	ENDELSE

	zerodSFQ_cart = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits',1)
	fields = zerodSFQ_cart[uniq(zerodSFQ_cart.field, sort(zerodSFQ_cart.field))].field

	ERASE
;	binwidth = 0.25

	dataIP  = dataIPallz
	dataAll = dataAllallz

  FOR f=0,n_elements(fields)-1 DO BEGIN
	dataIPfield = dataIPallz[WHERE(dataIPallz.field EQ fields[f])]
	dataAllfield = dataAllallz[WHERE(dataAllallz.field EQ fields[f])]

	stringField = STRTRIM(fields[f],2)
	stringFieldOut = ['CDFS', 'XMM-CFHTLS', 'COSMOS', 'ES1', 'XMM-SXDS']

	h = HIST_POINTS(dataSF[WHERE(dataSF.field EQ fields[f])].mstar, binwidth)
    IF (f EQ 0) THEN $
	PLOT, h[*,0], h[*,1], /NODATA, xrange=[xmin,xmax], yrange=[0,1200], xtitle='Stellar Mass', ytitle='Number', POSITION=posArray[*,2*f] $
    ELSE $
	PLOT, h[*,0], h[*,1], /NODATA, xrange=[xmin,xmax], yrange=[0,1200], xtickformat='(A1)', ytitle='Number', POSITION=posArray[*,2*f]

	OPLOT, h[*,0], h[*,1], COLOR=cgColor('blue'), LINESTYLE=2
	h = HIST_POINTS(dataQ[WHERE(dataQ.field EQ fields[f])].mstar, binwidth)
	OPLOT, h[*,0], h[*,1], PSYM=0, COLOR=cgColor('red'), LINESTYLE=2
	h = HIST_POINTS(dataIPSF[WHERE(dataIPSF.field EQ fields[f])].mstar, binwidth)
	OPLOT, h[*,0], h[*,1], PSYM=0, COLOR=cgColor('blue')
	h = HIST_POINTS(dataIPQ[WHERE(dataIPQ.field EQ fields[f])].mstar, binwidth)
	OPLOT, h[*,0], h[*,1], PSYM=0, COLOR=cgColor('red')

    IF (f EQ 4) THEN $
	LEGEND, ['SF IP', 'Q IP', 'all SF', 'all Q'], COLOR=cgColor(['blue','red','blue','red']), LINESTYLE=[0,0,2,2], /TOP, /LEFT, BOX=0, $
		NUMBER=2, PSPACING=3, CHARSIZE=0.8*charsz
  
;	get_hist_data, dataIPfield, dataAllfield, dz_coeff, 1, neighborMassesIPSF
;	SAVE, neighborMassesIPSF, FILENAME='neighborMassesIPSF_'+stringField+'.sav'
	RESTORE, 'neighborMassesIPSF_'+stringField+'.sav'
	histPointsIPSF	= HIST_POINTS(neighborMassesIPSF, binwidth)
;	histPointsIPSF	= HIST_POINTS(neighborMassesIPSF[UNIQ(neighborMassesIPSF,SORT(neighborMassesIPSF))], binwidth)
	x_IPSF = histPointsIPSF[*,0]
	y_IPSF = histPointsIPSF[*,1]

;	get_hist_data, dataIPfield, dataAllfield, dz_coeff, 0, neighborMassesIPQ
;	SAVE, neighborMassesIPQ, FILENAME='neighborMassesIPQ_'+stringField+'.sav'
	RESTORE, 'neighborMassesIPQ_'+stringField+'.sav'
	histPointsIPQ	= HIST_POINTS(neighborMassesIPQ, binwidth)
;	histPointsIPQ	= HIST_POINTS(neighborMassesIPQ[UNIQ(neighborMassesIPQ,SORT(neighborMassesIPQ))], binwidth)
	x_IPQ = histPointsIPQ[*,0]
	y_IPQ = histPointsIPQ[*,1]

;	PLOT, FINDGEN(10), xrange=[MIN([MIN(x_IPSF),MIN(x_IPQ)])-binwidth/2., MAX([MAX(x_IPSF),MAX(x_IPQ)])+binwidth/2.], $

    IF (f EQ 4) THEN BEGIN
	title = 'Neighbors Within 0 < R < 1 Mpc'
    ENDIF ELSE BEGIN
	title=''
    ENDELSE

    IF (f EQ 0) THEN BEGIN
	xtitle = 'Stellar Mass'
	xtickformat = ''
    ENDIF ELSE BEGIN
	xtitle=''
	xtickformat='(A1)'
    ENDELSE
	PLOT, FINDGEN(10), xrange=[xmin,xmax], $
		yrange=[MIN([MIN(y_IPSF),MIN(y_IPQ)]), 800], POSITION=posArray[*,2*f+1], $
		xtitle=xtitle, /NODATA, xminor=4, title=title, xtickformat=xtickformat, YSTYLE=4
	AXIS, YAXIS=0, YTICKFORMAT='(A1)'
	AXIS, YAXIS=1
    IF (f EQ 4) THEN $
	LEGEND, ['neigh of SF IP', 'neigh of Q IP'], COLOR=cgColor(['blue','red']), LINESTYLE=[0,0], /TOP, /RIGHT, BOX=0, $
		NUMBER=2, PSPACING=3, CHARSIZE=0.8*charsz

;    ENDIF
;	OPLOT, histPointsIPSF[*,0], histPointsIPSF[*,1], LINESTYLE=f, COLOR=cgColor('blue')
;	OPLOT, histPointsIPQ[*,0], histPointsIPQ[*,1], LINESTYLE=f, COLOR=cgColor('red')
;	OPLOT, histPointsIPSF[*,0], histPointsIPSF[*,1], LINESTYLE=f, COLOR=cgColor('BLU'+STRING(f+3))
;	OPLOT, histPointsIPQ[*,0], histPointsIPQ[*,1], LINESTYLE=f, COLOR=cgColor('RED'+STRING(f+3))
	OPLOT, histPointsIPSF[*,0], histPointsIPSF[*,1], LINESTYLE=0, COLOR=cgColor('blue')
	OPLOT, histPointsIPQ[*,0], histPointsIPQ[*,1], LINESTYLE=0, COLOR=cgColor('red')
	XYOUTS, xmin+0.05*xr, ymin+0.80*800, stringFieldOut[f], ALIGNMENT=0.0, CHARSIZE=charsz

  ENDFOR

;	LEGEND, ['SF IP neighbors', 'Q IP neighbors'], COLOR=cgColor(['blue','red']), LINESTYLE=[0,0], /TOP, /LEFT, BOX=0, $
;		NUMBER=2, PSPACING=3

	IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

PRO get_hist_data, dataIP, dataAll, dz_coeff, IPSFstatus, neighborMasses
	data_isIP = dataIP
	data = dataAll
	
        IP = data_isIP[WHERE(data_isIP.SFQ EQ IPSFstatus)]


	IF (IPSFstatus EQ 1) THEN (IPtype = 'Star-forming') ELSE (IPtype = 'Quiescent')
;	IF (IPSFstatus EQ 1) THEN (IP = IP[UNIQ(IP.objname)])
	PRINT, IPtype + ' IP galaxies: ' + strtrim(n_elements(IP),1)

	neighborMasses = []
	FOR i=0,(n_elements(IP)-1) DO BEGIN
;	FOR i=0,199 DO BEGIN
                currentIP = IP[i]

                IPneighbors = data[ WHERE( $
			(data.field EQ currentIP.field) AND $
			(data.objname NE currentIP.objname) AND $
			(ABS(data.zprimus - currentIP.zprimus) LE dz_coeff*0.005*(1.+currentIP.zprimus)) AND $
			(SQRT((data.xprop-currentIP.xprop)^2 + (data.yprop-currentIP.yprop)^2) LE 1.), /NULL ) ]
		
		IF (IPneighbors NE !NULL) THEN (neighborMasses = [neighborMasses, IPneighbors.mstar])
	IF (i MOD 100 EQ 0) THEN PRINT, i

        ENDFOR
END

FUNCTION HIST_POINTS, inputData, binwidth

;	PRINT, MINMAX(inputData)

	xpts = FLOOR(100*MIN(inputData))/100.-(FLOOR(100*MIN(inputData)) MOD (100*binwidth))/100. $
		+ binwidth*INDGEN( CEIL( ( MAX(inputData)-MIN(inputData) )/binwidth ) )
	ypts = HISTOGRAM(inputData, bin=binwidth)

;	PLOT, xpts+0.5*binwidth, HISTOGRAM(inputData, bin=binwidth), PSYM=10, xrange=MINMAX(inputData), yrange=[0,MAX(HISTOGRAM(inputData,bin=binwidth))]
;	PLOT, xpts, HISTOGRAM(inputData, bin=binwidth), PSYM=10, xrange=[MIN(inputData)-binwidth,MAX(inputData)+binwidth], yrange=1.05*[0,MAX(HISTOGRAM(inputData,bin=binwidth))]

	points = []
	FOR i=0,n_elements(xpts)-2 DO BEGIN
		points = [[points],$
		[xpts[i],ypts[i]],$
		[xpts[i+1],ypts[i]]]
	ENDFOR

	histPoints = [[xpts[0],0],$
		[points],$
		[xpts[-1],ypts[-1]],$
		[xpts[-1]+binwidth,ypts[-1]],$
		[xpts[-1]+binwidth,0]]

;	PRINT, TRANSPOSE(histPoints)

;	OPLOT, points[0,*], points[1,*], LINESTYLE=0, COLOR=COLOR

	RETURN, TRANSPOSE(histPoints)
END
