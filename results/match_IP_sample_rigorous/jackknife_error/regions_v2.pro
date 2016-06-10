PRO regions_v2, outputFormat
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsampleFBF/jackknife_regions', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=8
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

; read in data files
	IPdataAllFields  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)
	allDataAllFields = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

; eliminate data with targ_weight < 1
	allDataAllFields = allDataAllFields[where(allDataAllFields.targ_weight GE 1.)]
	IPdataAllFields = IPdataAllFields[where(IPdataAllFields.targ_weight GE 1.)]

	!p.multi=[0,2,2]
;	!p.multi=[2,2,1]
;	!p.multi=1
        !p.charsize=0.8
;	!X.MARGIN=[8,4]
;	!Y.MARGIN=[5,5]

	xmin = 0
	xmax = 2.5
	xr = xmax-xmin

	ymin = 0
	ymax = 2.3
	yr = ymax-ymin

        theta = 2*!PI*FINDGEN(17)/16.
        W = 0.1
        USERSYM, W*COS(theta), W*SIN(theta)

;	fields = ['cdfs', 'cfhtls_xmm', 'cosmos', 'es1', 'xmm']
	fields = ['cdfs', 'cfhtls_xmm', 'cosmos', 'es1']

	ERASE
	FOR j=0,n_elements(fields)-1 DO	plotField, allDataAllFields, IPdataAllFields, fields[j]

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

PRO plotField, allData, IPdata, field

	regionColors = cgcolor(['cyan', 'green', 'yellow', 'orange'])

;	All galaxies within given field
	allDataf = allData[where(strtrim(allData.field,2) EQ field)]
;	All IP within given field
	IPdataf  = IPdata[where(strtrim(IPdata.field,2) EQ field)]

	xmin = MIN(allDataf.RA)
	xmax = MAX(allDataf.RA)
;	xmax = xmin + 2.5
	xrange = xmax-xmin
	ymin = MIN(allDataf.DEC)
	ymax = MAX(allDataf.DEC)
;	ymax = ymin + 2.3
	yrange = ymax - ymin

	range = MAX([xrange,yrange])
	xmax = xmin + range
	ymax = ymin + range

	linecol = cgColor('red')

	PLOT, allDataf.RA, allDataf.DEC, title=field, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='RA', ytitle='DEC', /NODATA

	randPath = '~/edge_effects_test/random/'
	IF (field NE 'cfhtls_xmm') THEN randFieldData = MRDFITS(randPath + field + '_random.fits',1) ELSE $
		randFieldData = MRDFITS(randPath + 'xmm_random.fits',1)
;	OPLOT, randFieldData.RA, randfieldData.DEC, psym=8, color=cgcolor('green')
;	OPLOT, allDataf.RA, allDataf.DEC, psym=8, color=cgcolor('green')
	orderedRA  = allDataf[SORT(allDataf.RA)].RA
	orderedDEC = allDataf[SORT(allDataf.DEC)].DEC

  IF (field EQ 'cdfs') THEN BEGIN
	RA1 	= 0.9 + MIN(allDataf.RA)
	RA2 	= MEDIAN(allDataf[where(allDataf.RA GE RA1)].RA) 

;	OPLOT, [RA1,RA1], [ymin,ymax], LINESTYLE=0, THICK=2, color=linecol
;	OPLOT, [RA2,RA2], [ymin,ymax], LINESTYLE=0, THICK=2, color=linecol

	rand1 = randFieldData[WHERE(randFieldData.RA LT RA1)]
	OPLOT, rand1.RA, rand1.DEC, psym=3, color=regionColors[0]
	rand2 = randFieldData[WHERE( (randFieldData.RA GE RA1) AND (randFieldData.RA LE RA2) )]
	OPLOT, rand2.RA, rand2.DEC, psym=3, color=regionColors[1]
	rand3 = randFieldData[WHERE(randFieldData.RA GT RA2)]
	OPLOT, rand3.RA, rand3.DEC, psym=3, color=regionColors[2]

	numR01 = n_elements(allDataf[where(allDataf.RA LT RA1)])
	numR02 = n_elements(allDataf[where( (allDataf.RA GE RA1) AND (allDataf.RA LE RA2) )])
	numR03 = n_elements(allDataf[where(allDataf.RA GT RA2)])

	PRINT, field, numR01, numR02, numR03
	PRINT, 'RA: ', RA1, RA2
  ENDIF

  IF (field EQ 'cfhtls_xmm') THEN BEGIN
	DEC1 	= -5.0
	DEC2	= -4.347
	RA1 	= MEDIAN(allDataf[where(allDataf.DEC GE DEC1)].RA) 

;	OPLOT, [xmin,xmax], [DEC1,DEC1], LINESTYLE=0, THICK=2, color=linecol
;	OPLOT, [xmin,xmax], [DEC2,DEC2], LINESTYLE=0, THICK=2, color=linecol

	XMM = allData[WHERE(allData.field EQ 'xmm       ')]

	xmmDEC1 = MAX(XMM.DEC)
	xmmDEC2 = MAX(XMM[WHERE(XMM.RA LT 34.2)].DEC)
	xmmDEC3 = MIN(XMM[WHERE(XMM.RA GT 34.8)].DEC)
	xmmDEC4 = MIN(XMM.DEC)

	xmmRA1  = MIN(XMM.RA)
	xmmRA2a = MIN(XMM[WHERE(XMM.DEC GT -4.7)].RA)
	xmmRA2b = MIN(XMM[WHERE(XMM.DEC LT -5.3)].RA)
	xmmRA3a = MAX(XMM[WHERE(XMM.DEC GT -4.7)].RA)
	xmmRA3b = MAX(XMM[WHERE(XMM.DEC LT -5.3)].RA)
	xmmRA4  = MAX(XMM.RA)

	randXMM = randFieldData

	randXMM1 = randXMM[WHERE( (randXMM.RA GE xmmRA2a) AND (randXMM.RA LE xmmRA3a) AND (randXMM.DEC LE xmmDEC1) AND (randXMM.DEC GE xmmDEC2) )]
	randXMM2 = randXMM[WHERE( (randXMM.RA GE xmmRA1) AND (randXMM.RA LE xmmRA4) AND (randXMM.DEC LE xmmDEC2) AND (randXMM.DEC GE xmmDEC3) )]
	randXMM3 = randXMM[WHERE( (randXMM.RA GE xmmRA2b) AND (randXMM.RA LE xmmRA3b) AND (randXMM.DEC LE xmmDEC3) AND (randXMM.DEC GE xmmDEC4) )]

	rand1 = randFieldData[WHERE(randFieldData.DEC LT DEC1)]
	OPLOT, rand1.RA, rand1.DEC, psym=3, color=regionColors[0]
	rand2 = randFieldData[WHERE( (randFieldData.DEC GE DEC1) AND (randFieldData.DEC LE DEC2) )]
	OPLOT, rand2.RA, rand2.DEC, psym=3, color=regionColors[1]
	rand3 = randFieldData[WHERE(randFieldData.DEC GT DEC2)]
	OPLOT, rand3.RA, rand3.DEC, psym=3, color=regionColors[2]

	OPLOT, randXMM1.RA, randXMM1.DEC, psym=3, color=regionColors[3]
	OPLOT, randXMM2.RA, randXMM2.DEC, psym=3, color=regionColors[3]
	OPLOT, randXMM3.RA, randXMM3.DEC, psym=3, color=regionColors[3]
;	OPLOT, allData[WHERE(allData.field EQ 'xmm       ')].RA, allData[WHERE(allData.field EQ 'xmm       ')].DEC, psym=8, color=regionColors[3]
	OPLOT, IPdata[WHERE(IPdata.field EQ 'xmm       ')].RA, IPdata[WHERE(IPdata.field EQ 'xmm       ')].DEC, psym=8

	numR04 = n_elements(allDataf[where(allDataf.DEC LT DEC1)])
	numR05 = n_elements(allDataf[where( (allDataf.DEC GE DEC1) AND (allDataf.DEC LE DEC2) )])
	numR06 = n_elements(allDataf[where(allDataf.DEC GT DEC2)])

	PRINT, field, numR04, numR05, numR06
	PRINT, 'DEC: ', DEC1, DEC2
  ENDIF

  IF (field EQ 'cosmos') THEN BEGIN
	DEC1 	= MEDIAN(allDataf.DEC)

;	OPLOT, [xmin,xmax], [DEC1,DEC1], LINESTYLE=0, THICK=2, color=linecol

	rand1 = randFieldData[WHERE(randFieldData.DEC LE DEC1)]
	OPLOT, rand1.RA, rand1.DEC, psym=3, color=regionColors[0]
	rand2 = randFieldData[WHERE(randFieldData.DEC GT DEC1)]
	OPLOT, rand2.RA, rand2.DEC, psym=3, color=regionColors[1]

	numR07 = n_elements(allDataf[where(allDataf.DEC LE DEC1)])
	numR08 = n_elements(allDataf[where(allDataf.DEC GT DEC1)])

	PRINT, field, numR07, numR08
	PRINT, 'DEC: ', DEC1
  ENDIF

  IF (field EQ 'es1') THEN BEGIN
	PRINT, field, n_elements(allDataf)
	OPLOT, randFieldData.RA, randFieldData.DEC, psym=3, color=regionColors[0]
  ENDIF

  IF (field EQ 'XMM') THEN BEGIN
	PRINT, field, n_elements(allDataf)
  ENDIF

	OPLOT, IPdataf.RA, IPdataf.DEC, psym=8

END
