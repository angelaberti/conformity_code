PRO projection_test04, outputFormat
	dz_coeff = 2.0

	zmin = 0.68
	zmax = 1.

	; read in data files
	dataIP  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)
	dataAll = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

	dataAll	  = dataAll[where( (dataAll.field EQ 'cosmos    ') AND (dataAll.zprimus GE zmin) AND (dataAll.zprimus LE zmax) )]
	dataAllSF = dataAll[where(dataAll.SFQ EQ 1)]
	dataAllQ  = dataAll[where(dataAll.SFQ EQ 0)]

	dataIP	  = dataIP[where( (dataIP.field EQ 'cosmos    ') AND (dataIP.zprimus GE zmin) AND (dataIP.zprimus LE zmax) )]
	dataIPSF = dataIP[where(dataIP.SFQ EQ 1)]
	dataIPQ  = dataIP[where(dataIP.SFQ EQ 0)]

	SFcolor = cgColor('blue')
	Qcolor  = cgColor('red')

	ERASE
	!p.multi=[2,2,1]
	!p.charsize=1.25

	xmin = -10
	xmax = -1.*xmin
	xr = xmax-xmin

	ymin = -10
	ymax = -1.*ymin
	yr = ymax-ymin
	
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsampleFBF/COSMOS_projection_test', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=8
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	theta = 2*!PI*FINDGEN(65)/64.

	PLOT, [1], [1], xrange=[xmin,xmax], yrange=[ymin,ymax], title='Q IP', /NODATA
;	partition, dataIPQ, dataAll, dz_coeff, SFcolor, Qcolor, xmin, xmax, ymin, ymax
	plotNeigh, dataIPQ, dataAll, dz_coeff, SFcolor, Qcolor;, neighSFdistsAll, neighQdistsAll, dz_coeff
;	FOR i=2*xmin,2*xmax-1 DO OPLOT, [i/2.,i/2.], [ymin,ymax], LINESTYLE=0, THICK=3;, COLOR=cgcolor('black')
;	FOR i=2*ymin,2*ymax-1 DO OPLOT, [xmin,xmax], [i/2.,i/2.], LINESTYLE=0, THICK=3;, COLOR=cgcolor('black')
	FOR i=1,5 DO OPLOT, 2.*i*COS(theta), 2.*i*SIN(theta), LINESTYLE=0, THICK=1
	FOR i=1,5 DO OPLOT, (2.*i-1)*COS(theta), (2.*i-1)*SIN(theta), LINESTYLE=2, THICK=1

	PLOT, [1], [1], xrange=[xmin,xmax], yrange=[ymin,ymax], title='SF IP', /NODATA
;	partition, dataIPSF, dataAll, dz_coeff, SFcolor, Qcolor, xmin, xmax, ymin, ymax
	plotNeigh, dataIPSF, dataAll, dz_coeff, SFcolor, Qcolor;, neighSFdistsAll, neighQdistsAll, dz_coeff
	FOR i=1,5 DO OPLOT, 2.*i*COS(theta), 2.*i*SIN(theta), LINESTYLE=0, THICK=1
	FOR i=1,5 DO OPLOT, (2.*i-1)*COS(theta), (2.*i-1)*SIN(theta), LINESTYLE=2, THICK=1

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

PRO getDists, IPdata, dataAll, neighSFdistsAll, neighQdistsAll, dz_coeff
	neighSFdistsAll = []
	neighQdistsAll  = []

	FOR i=0,n_elements(IPdata)-1 DO BEGIN
		currentIP = IPdata[i]
	
		; keep galaxies within dz=dz+coeff*0.005*(1+z) of IP candidate (and in same field)
		neigh = dataAll[ where( (dataAll.objname NE currentIP.objname) AND $
                        (ABS(dataAll.zprimus - currentIP.zprimus) LE dz_coeff*0.005*(1.+currentIP.zprimus)) ) ]

		neighSF = neigh[where(neigh.SFQ EQ 1)]
		neighQ  = neigh[where(neigh.SFQ EQ 0)]

		neighSFdists = [SQRT((neighSF.xprop-currentIP.xprop)^2 + (neighSF.yprop-currentIP.yprop)^2)]
		neighQdists = [SQRT((neighQ.xprop-currentIP.xprop)^2 + (neighQ.yprop-currentIP.yprop)^2)]

		neighSFdistsAll = [neighSFdistsAll, neighSFdists]
		neighQdistsAll  = [neighQdistsAll, neighQdists]

;		OPLOT, neighSF.xprop-currentIP.xprop, neighSF.yprop-currentIP.yprop, COLOR=SFcolor, psym=3
;		OPLOT, neighQ.xprop-currentIP.xprop, neighQ.yprop-currentIP.yprop, COLOR=Qcolor, psym=3
	ENDFOR
END


PRO partition, IPdata, dataAll, dz_coeff, SFcolor, Qcolor, xmin, xmax, ymin, ymax
	theta = 2*!PI*FINDGEN(33)/32.
	W = 0.25
	USERSYM, W*COS(theta), W*SIN(theta)

	neighSFall_dx = []
	neighSFall_dy = []
	neighQall_dx  = []
	neighQall_dy  = []

	FOR i=0,n_elements(IPdata)-1 DO BEGIN
		currentIP = IPdata[i]
	
		; keep galaxies within dz=dz+coeff*0.005*(1+z) of IP candidate (and in same field)
		neigh = dataAll[ where( (dataAll.objname NE currentIP.objname) AND $
                        (ABS(dataAll.zprimus - currentIP.zprimus) LE dz_coeff*0.005*(1.+currentIP.zprimus)) ) ]

		neighSF = neigh[where(neigh.SFQ EQ 1)]
		neighQ  = neigh[where(neigh.SFQ EQ 0)]

		neighSFall_dx = [neighSFall_dx, neighSF.xprop-currentIP.xprop]
		neighSFall_dy = [neighSFall_dy, neighSF.yprop-currentIP.yprop]
		neighQall_dx  = [neighQall_dx, neighQ.xprop-currentIP.xprop]
		neighQall_dy  = [neighQall_dy, neighQ.yprop-currentIP.yprop]
	ENDFOR

	neighSFall = TRANSPOSE([[neighSFall_dx], [neighSFall_dy]])
	neighQall  = TRANSPOSE([[neighQall_dx], [neighQall_dy]])

	xdivs = xmin+FINDGEN(4*(xmax-xmin)+1)/4
	ydivs = ymin+FINDGEN(4*(ymax-ymin)+1)/4

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

;	blueArray = 50+CEIL(255*boxWeightArray*boxFracArraySF/1.25)
;	redArray  = 50+CEIL(255*boxWeightArray*boxFracArrayQ/1.25)
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

;	OPLOT, neighSFall_dx, neighSFall_dy, COLOR=SFcolor, psym=8
;	OPLOT, neighQall_dx, neighQall_dy, COLOR=Qcolor, psym=8

END


PRO plotNeigh, IPdata, dataAll, dz_coeff, SFcolor, Qcolor;, neighSFdistsAll, neighQdistsAll, dz_coeff
;	neighSFdistsAll = []
;	neighQdistsAll  = []

	theta = 2*!PI*FINDGEN(33)/32.
	W = 0.25
	USERSYM, W*COS(theta), W*SIN(theta)

	FOR i=0,n_elements(IPdata)-1 DO BEGIN
		currentIP = IPdata[i]
	
		; keep galaxies within dz=dz+coeff*0.005*(1+z) of IP candidate (and in same field)
		neigh = dataAll[ where( (dataAll.objname NE currentIP.objname) AND $
                        (ABS(dataAll.zprimus - currentIP.zprimus) LE dz_coeff*0.005*(1.+currentIP.zprimus)) ) ]

		neighSF = neigh[where(neigh.SFQ EQ 1)]
		neighQ  = neigh[where(neigh.SFQ EQ 0)]

		OPLOT, neighSF.xprop-currentIP.xprop, neighSF.yprop-currentIP.yprop, COLOR=SFcolor, psym=8
		OPLOT, neighQ.xprop-currentIP.xprop, neighQ.yprop-currentIP.yprop, COLOR=Qcolor, psym=8

;		OPLOT, ABS(neighQ.xprop-currentIP.xprop), ABS(neighQ.yprop-currentIP.yprop), COLOR=Qcolor, psym=8
;		OPLOT, ABS(neighSF.xprop-currentIP.xprop), ABS(neighSF.yprop-currentIP.yprop), COLOR=SFcolor, psym=8
	ENDFOR
END
