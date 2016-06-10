PRO regions, outputFormat
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsampleFBF/jackknife_regions', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=9, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

; read in data files
	IPdataAllFields  = MRDFITS('~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)
	allDataAllFields = MRDFITS('~/results/zerodSFQ_all_cart.fits', 1)

; eliminate data with targ_weight < 1
	allDataAllFields = allDataAllFields[where(allDataAllFields.targ_weight GE 1.)]
	IPdataAllFields = IPdataAllFields[where(IPdataAllFields.targ_weight GE 1.)]

	!p.multi=[0,3,2]
;	!p.multi=[2,2,1]
;	!p.multi=1
        !p.charsize=1.5
;	!X.MARGIN=[8,4]
;	!Y.MARGIN=[5,5]

	xmin = 0
	xmax = 2.5
	xr = xmax-xmin

	ymin = 0
	ymax = 2.3
	yr = ymax-ymin

        theta = 2*!PI*FINDGEN(17)/16.
        W = 0.25
        USERSYM, W*COS(theta), W*SIN(theta)

	fields = ['cdfs', 'cfhtls_xmm', 'cosmos', 'es1', 'xmm']
;	fields = ['cfhtls_xmm', 'cosmos', 'xmm']

;	allData = allDataAllFields
;	IPdata = IPdataAllFields

	ERASE
	FOR j=0,n_elements(fields)-1 DO	plotField, allDataAllFields, IPdataAllFields, fields[j]
;	FOR j=0,n_elements(fields)-1 DO	plotField, allData, IPdata, fields[j]

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

PRO plotField, allData, IPdata, field

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
	OPLOT, allDataf.RA, allDataf.DEC, psym=8, color=cgcolor('green')
	OPLOT, IPdataf.RA, IPdataf.DEC, psym=8
	orderedRA  = allDataf[SORT(allDataf.RA)].RA
	orderedDEC = allDataf[SORT(allDataf.DEC)].DEC

  IF (field EQ 'cdfs') THEN BEGIN
	RA1 	= 0.9 + MIN(allDataf.RA)
	RA2 	= MEDIAN(allDataf[where(allDataf.RA GE RA1)].RA) 

	OPLOT, [RA1,RA1], [ymin,ymax], LINESTYLE=0, THICK=2, color=linecol
	OPLOT, [RA2,RA2], [ymin,ymax], LINESTYLE=0, THICK=2, color=linecol

	numR01 = n_elements(allDataf[where(allDataf.RA LT RA1)])
	numR02 = n_elements(allDataf[where( (allDataf.RA GE RA1) AND (allDataf.RA LE RA2) )])
	numR03 = n_elements(allDataf[where(allDataf.RA GT RA2)])

;	IPdataf[where(IPdataf.RA LT RA1)], allDataf)

;	IPdata = IPdataf[where( (IPdataf.RA GE RA1) AND (IPdataf.RA LE RA2) )]
;	IPdata = IPdataf[where(IPdataf.RA GT RA2)]

	PRINT, field, numR01, numR02, numR03
	PRINT, 'RA: ', RA1, RA2
  ENDIF

  IF (field EQ 'cfhtls_xmm') THEN BEGIN
	DEC1 	= -5.0
	DEC2	= -4.347
	RA1 	= MEDIAN(allDataf[where(allDataf.DEC GE DEC1)].RA) 

	OPLOT, [xmin,xmax], [DEC1,DEC1], LINESTYLE=0, THICK=2, color=linecol
	OPLOT, [xmin,xmax], [DEC2,DEC2], LINESTYLE=0, THICK=2, color=linecol

	numR04 = n_elements(allDataf[where(allDataf.DEC LT DEC1)])
	numR05 = n_elements(allDataf[where( (allDataf.DEC GE DEC1) AND (allDataf.DEC LE DEC2) )])
	numR06 = n_elements(allDataf[where(allDataf.DEC GT DEC2)])

	PRINT, field, numR04, numR05, numR06
	PRINT, 'DEC: ', DEC1, DEC2
  ENDIF

  IF (field EQ 'cosmos') THEN BEGIN
	DEC1 	= MEDIAN(allDataf.DEC)

	OPLOT, [xmin,xmax], [DEC1,DEC1], LINESTYLE=0, THICK=2, color=linecol

	numR07 = n_elements(allDataf[where(allDataf.DEC LE DEC1)])
	numR08 = n_elements(allDataf[where(allDataf.DEC GT DEC1)])

	PRINT, field, numR07, numR08
	PRINT, 'DEC: ', DEC1
  ENDIF

  IF (field EQ 'es1') THEN BEGIN
	condR09 = 'strtrim(IPdata.field,2) EQ field'
	PRINT, field, n_elements(allDataf)
  ENDIF

  IF (field EQ 'xmm') THEN BEGIN
	condR10 = 'strtrim(IPdata.field,2) EQ field'
	PRINT, field, n_elements(allDataf)
  ENDIF

END
