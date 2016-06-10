PRO neighborHist_SFtype, binwidth, outputFormat
	PHI = 3.7
	stringPHI = strtrim(string(PHI, format='(f20.1)'),2)

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	; read in data files
	dataIPallz  = MRDFITS('~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
	dataAllallz = MRDFITS('~/results/zerodSFQ_all_cart.fits', 1)

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

	!P.MULTI=[2,1,2]
;	!P.MULTI=0
	!P.CHARSIZE=1.2
	!P.FONT=0

	Rproj_array = float(dRproj)*(findgen(n_annuli+1))	
;	PRINT, Rproj_array

	xmin = 8
	xmax = 11.5
	xr = xmax-xmin

	ymin = 0
	ymax = 2800
	yr = ymax-ymin

	IF (string(outputFormat) EQ 'ps') THEN BEGIN
		PS_OPEN, '~/latex/figures/neighborHist_SFtype', THICK=5, /ENCAP
		DEVICE, /INCH, XS=6, YS=6, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
		SET_PLOT, 'X'
		THICK=1
	ENDELSE

	ERASE

	dataIP  = dataIPallz
	dataAll = dataAllallz

;	get_hist_data, dataIP, dataAll, dz_coeff, 1, 1, neighborMassesIPSF_SF
;	SAVE, neighborMassesIPSF_SF, FILENAME='neighborMassesIPSF_SF.sav'
	RESTORE, 'neighborMassesIPSF_SF.sav'
	histPointsIPSF_SF	= HIST_POINTS(neighborMassesIPSF_SF, binwidth)
	x_IPSF_SF = histPointsIPSF_SF[*,0]
	y_IPSF_SF = histPointsIPSF_SF[*,1]

;	get_hist_data, dataIP, dataAll, dz_coeff, 1, 0, neighborMassesIPSF_Q
;	SAVE, neighborMassesIPSF_Q, FILENAME='neighborMassesIPSF_Q.sav'
	RESTORE, 'neighborMassesIPSF_Q.sav'
	histPointsIPSF_Q	= HIST_POINTS(neighborMassesIPSF_Q, binwidth)
	x_IPSF_Q = histPointsIPSF_Q[*,0]
	y_IPSF_Q = histPointsIPSF_Q[*,1]

;	get_hist_data, dataIP, dataAll, dz_coeff, 0, 1, neighborMassesIPQ_SF
;	SAVE, neighborMassesIPQ_SF, FILENAME='neighborMassesIPQ_SF.sav'
	RESTORE, 'neighborMassesIPQ_SF.sav'
	histPointsIPQ_SF	= HIST_POINTS(neighborMassesIPQ_SF, binwidth)
	x_IPQ_SF = histPointsIPQ_SF[*,0]
	y_IPQ_SF = histPointsIPQ_SF[*,1]

;	get_hist_data, dataIP, dataAll, dz_coeff, 0, 0, neighborMassesIPQ_Q
;	SAVE, neighborMassesIPQ_Q, FILENAME='neighborMassesIPQ_Q.sav'
	RESTORE, 'neighborMassesIPQ_Q.sav'
	histPointsIPQ_Q		= HIST_POINTS(neighborMassesIPQ_Q, binwidth)
	x_IPQ_Q = histPointsIPQ_Q[*,0]
	y_IPQ_Q = histPointsIPQ_Q[*,1]

	posArray = grid_array(1,2)

;	PLOT, FINDGEN(10), xrange=[MIN([MIN(x_IPSF),MIN(x_IPQ)])-binwidth/2., MAX([MAX(x_IPSF),MAX(x_IPQ)])+binwidth/2.], $
	PLOT, FINDGEN(10), xrange=[xmin,xmax], $
		yrange=[0, 2800], POSITION=posArray[*,1], xtickformat='(A1)', $
		ytitle=textoidl('N_{neighbors}'), /NODATA, xminor=4 ;, title = 'Neighbors Within 0 < R < 1 Mpc'
;		xtitle = textoidl('log ( M_{*}/M' + sunsymbol() + ' )'), ytitle=textoidl('N_{neighbors}'), /NODATA, xminor=4 ;, title = 'Neighbors Within 0 < R < 1 Mpc'
	XYOUTS, xmin+0.95*xr, ymin+0.85*yr, 'Star-forming IPs', ALIGNMENT=1.0
	OPLOT, histPointsIPSF_SF[*,0], histPointsIPSF_SF[*,1], LINESTYLE=0, COLOR=cgColor('blue')
	OPLOT, histPointsIPSF_Q[*,0], histPointsIPSF_Q[*,1], LINESTYLE=2, COLOR=cgColor('red')
	LEGEND, [textoidl('SF neighbors; M_{*med} = ') + decimal(MEDIAN(neighborMassesIPSF_SF),3), $
		textoidl(' Q neighbors; M_{*med} = ') + decimal(MEDIAN(neighborMassesIPSF_Q),3)], COLOR=cgColor(['blue', 'red']), LINESTYLE=[0,2], /TOP, /LEFT, BOX=0, $
		NUMBER=2, PSPACING=3, CHARSIZE=1.1

	PLOT, FINDGEN(10), xrange=[xmin,xmax], $
		yrange=[0, 2800], POSITION=posArray[*,0], $
		xtitle = textoidl('log ( M_{*}/M' + sunsymbol() + ' )'), ytitle=textoidl('N_{neighbors}'), /NODATA, xminor=4 
	XYOUTS, xmin+0.95*xr, ymin+0.85*yr, 'Quiescent IPs', ALIGNMENT=1.0
	OPLOT, histPointsIPQ_SF[*,0], histPointsIPQ_SF[*,1], LINESTYLE=0, COLOR=cgColor('blue')
	OPLOT, histPointsIPQ_Q[*,0], histPointsIPQ_Q[*,1], LINESTYLE=2, COLOR=cgColor('red')
;	LEGEND, ['SF IP; all neigh', 'Q IP; all neigh', 'SF IP; unique neigh', 'Q IP; unique neigh'], COLOR=cgColor(['BLU4','RED4', 'blue', 'red']), LINESTYLE=[4,3,0,2], /TOP, /LEFT, BOX=0, $
	LEGEND, [textoidl('SF neighbors; M_{*med} = ') + decimal(MEDIAN(neighborMassesIPQ_SF),3), $
		textoidl(' Q neighbors; M_{*med} = ') + decimal(MEDIAN(neighborMassesIPQ_Q),3)], COLOR=cgColor(['blue', 'red']), LINESTYLE=[0,2], /TOP, /LEFT, BOX=0, $
		NUMBER=2, PSPACING=3, CHARSIZE=1.1

;	PRINT, 'Median SF IP neighbor mass: ', MEDIAN(neighborMassesIPSF)
;	PRINT, 'Median Q IP neighbor mass: ', MEDIAN(neighborMassesIPQ)

	IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

PRO get_hist_data, dataIP, dataAll, dz_coeff, IPSFstatus, neighSFstatus, neighborMasses
	data_isIP = dataIP
	data = dataAll
	
        IP = data_isIP[WHERE(data_isIP.SFQ EQ IPSFstatus)]

	IF (IPSFstatus EQ 1) THEN (IPtype = 'Star-forming') ELSE (IPtype = 'Quiescent')
;	IF (IPSFstatus EQ 1) THEN (IP = IP[UNIQ(IP.objname)])
	PRINT, IPtype + ' IP galaxies: ' + strtrim(n_elements(IP),1)

	neighborMasses = []
	FOR i=0,(n_elements(IP)-1) DO BEGIN
;	FOR i=0,99 DO BEGIN
                currentIP = IP[i]

                IPneighbors = data[ WHERE( $
			(data.field EQ currentIP.field) AND (data.SFQ EQ neighSFstatus) AND $
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
