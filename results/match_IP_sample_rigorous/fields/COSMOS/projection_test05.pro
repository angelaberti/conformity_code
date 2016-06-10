PRO projection_test05, outputFormat
	dz_coeff = 2.0

	zmin = 0.68
	zmax = 1.

	; read in data files
	dataIP  = MRDFITS('~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)
	data	= MRDFITS('~/results/zerodSFQ_all_cart.fits', 1)

	ERASE
	!p.multi=[6,3,2]
;	!p.multi=[4,2,2]
;	!p.multi=1
	!p.charsize=3

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsampleFBF/COSMOS_projection_test', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=8
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	theta = 2*!PI*FINDGEN(65)/64.

	redRange1 = [0.65, 0.76]
	redRange2 = [0.82, 1.]
	nbins 	  = 20
	addPoints = 1

	partition, data, redRange1, 'XPROP', 'YPROP', nbins, addPoints
	partition, data, redRange1, 'ZPRIMUS', 'YPROP', nbins, addPoints
	partition, data, redRange1, 'ZPRIMUS', 'XPROP', nbins, addPoints
	partition, data, redRange2, 'XPROP', 'YPROP', nbins, addPoints
	partition, data, redRange2, 'ZPRIMUS', 'YPROP', nbins, addPoints
	partition, data, redRange2, 'ZPRIMUS', 'XPROP', nbins, addPoints

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END


PRO partition, data, redRange, xfield, yfield, nbins, addPoints
	dataAll = data

	theta = 2*!PI*FINDGEN(33)/32.
	W = 0.25
	USERSYM, W*COS(theta), W*SIN(theta)

	SFcolor = cgColor('blue')
	Qcolor  = cgColor('red')

	dataAll	= dataAll[WHERE( (dataAll.field EQ 'cosmos    ') AND (dataAll.zprimus GE redRange[0]) AND (dataAll.zprimus LE redRange[1]) )]	

	tnames = TAG_NAMES(dataAll)
	dataXX = dataAll.(WHERE(STRCMP(tnames, xfield) EQ 1))
	dataYY = dataAll.(WHERE(STRCMP(tnames, yfield) EQ 1))

	xmin = MIN(dataXX)
	xmax = MAX(dataXX)
	xr   = xmax-xmin
	ymin = MIN(dataYY)
	ymax = MAX(dataYY)
	yr   = ymax-ymin

	dataSF = dataAll[WHERE(dataAll.SFQ EQ 1)]
	dataQ  = dataAll[WHERE(dataAll.SFQ EQ 0)]

	dataXX_SF = dataSF.(WHERE(STRCMP(tnames, xfield) EQ 1))
	dataYY_SF = dataSF.(WHERE(STRCMP(tnames, yfield) EQ 1))
	dataXX_Q  = dataQ.(WHERE(STRCMP(tnames, xfield) EQ 1))
	dataYY_Q  = dataQ.(WHERE(STRCMP(tnames, yfield) EQ 1))

	PLOT, [1], [1], /NODATA, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle=xfield, ytitle=yfield

	neighSFall_dx = dataXX_SF - MEAN(dataXX_SF)
	neighQall_dx  = dataXX_Q - MEAN(dataXX_Q)
	neighSFall_dy = dataYY_SF - MEAN(dataYY_SF)
	neighQall_dy  = dataYY_Q - MEAN(dataYY_Q)

  IF (xfield EQ 'ZPRIMUS') THEN BEGIN
	neighSFall_dx = dataXX_SF
	neighQall_dx  = dataXX_Q
  ENDIF
  IF (yfield EQ 'ZPRIMUS') THEN BEGIN
	neighSFall_dy = dataYY_SF
	neighQall_dy  = dataYY_Q
  ENDIF

	neighSFall = TRANSPOSE([[neighSFall_dx], [neighSFall_dy]])
	neighQall  = TRANSPOSE([[neighQall_dx], [neighQall_dy]])
	
	xdivs = xmin+FINDGEN(nbins+1)*(xr/nbins)
	ydivs = ymin+FINDGEN(nbins+1)*(yr/nbins)
	
	numArraySF = []
	numArrayQ  = []
	boxTotalArray = []
	FOR x_index=0,n_elements(xdivs)-2 DO BEGIN
		in_x_range_SF = neighSFall_dy[WHERE( (neighSFall_dx GE xdivs[x_index]) AND (neighSFall_dx LT xdivs[x_index+1]) )]
		in_x_range_Q  = neighQall_dy[WHERE( (neighQall_dx GE xdivs[x_index]) AND (neighQall_dx LT xdivs[x_index+1]) )]

		numColSF = []
		numColQ  = []
		boxTotalCol = []
		FOR y_index=0,n_elements(ydivs)-2 DO BEGIN
			in_box_SF = n_elements(in_x_range_SF[WHERE( (in_x_range_SF GE ydivs[y_index]) AND (in_x_range_SF LT ydivs[y_index+1]) )])
			in_box_Q  = n_elements(in_x_range_Q[WHERE( (in_x_range_Q GE ydivs[y_index]) AND (in_x_range_Q LT ydivs[y_index+1]) )])
			boxTotal = in_box_SF + in_box_Q
			
			numColSF = [numColSF, in_box_SF]
			numColQ  = [numColQ, in_box_Q]
			boxTotalCol = [boxTotalCol, boxTotal]
		ENDFOR
		numArraySF = [[numArraySF], [numColSF]]
		numArrayQ  = [[numArrayQ], [numColQ]]
		boxTotalArray = [[boxTotalArray], [boxTotalCol]]
	ENDFOR

	boxWeightArray = boxTotalArray-MIN(boxTotalArray)
	boxWeightArray = boxTotalArray/FLOAT(MAX(boxTotalArray))
	
	boxFracArraySF = numArraySF/FLOAT(boxTotalArray)
	boxFracArrayQ  = numArrayQ/FLOAT(boxTotalArray)

;	blueArray = CEIL(255*boxWeightArray*boxFracArraySF)
;	redArray  = CEIL(255*boxWeightArray*boxFracArrayQ)
	blueArray = 25+CEIL(255*boxWeightArray*boxFracArraySF/1.1)
	redArray  = 25+CEIL(255*boxWeightArray*boxFracArrayQ/1.1)

	blueArray = TRANSPOSE(blueArray)
	redArray  = TRANSPOSE(redArray)
	
	FOR i=0,n_elements(xdivs)-2 DO BEGIN
		FOR j=0,n_elements(ydivs)-2 DO BEGIN
			rect = [ [xdivs[i],ydivs[j]], [xdivs[i+1],ydivs[j]], [xdivs[i+1],ydivs[j+1]], [xdivs[i],ydivs[j+1]] ]
;			PRINT, rect[*,0]
			POLYFILL, rect, COLOR=color24([redArray[i,j],0,blueArray[i,j]])
		ENDFOR
	ENDFOR

;	FOR i=0,n_elements(xdivs)-2 DO OPLOT, [xdivs[i],xdivs[i]], [ymin,ymax], LINESTYLE=1, THICK=1 
;	FOR i=0,n_elements(ydivs)-2 DO OPLOT, [xmin,xmax], [ydivs[i],ydivs[i]], LINESTYLE=1, THICK=1 

  IF (addPoints EQ 1) THEN BEGIN
	OPLOT, neighSFall_dx, neighSFall_dy, COLOR=SFcolor, psym=8
	OPLOT, neighQall_dx, neighQall_dy, COLOR=Qcolor, psym=8
  ENDIF

	XYOUTS, xmin+0.05*xr, ymin+0.05*yr, 'z=['+decimal(redRange[0],2)+', '+decimal(redRange[1],2)+']', CHARSIZE=1.5
END
