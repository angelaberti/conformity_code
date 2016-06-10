PRO projection_test02, outputFormat
	dz_coeff = 2.0

	; read in data files
	dataIP  = MRDFITS('~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)
	dataAll = MRDFITS('~/results/zerodSFQ_all_cart.fits', 1)

	SFcolor = cgColor('blue')
	Qcolor  = cgColor('red')

	ERASE
	!p.multi=[6,3,2]
	!p.charsize=1.25

;	zmin = 0.65
;	zmax = 0.75
	zmin = .35
	zmax = .4

	dataAll	  = dataAll[where( (dataAll.field EQ 'cosmos    ') AND (dataAll.zprimus GE zmin) AND (dataAll.zprimus LE zmax) )]
	dataAllSF = dataAll[where(dataAll.SFQ EQ 1)]
	dataAllQ  = dataAll[where(dataAll.SFQ EQ 0)]

	dataIP	 = dataIP[where( (dataIP.field EQ 'cosmos    ') AND (dataIP.zprimus GE zmin) AND (dataIP.zprimus LE zmax) )]
	dataIPSF = dataIP[where(dataIP.SFQ EQ 1)]
	dataIPQ  = dataIP[where(dataIP.SFQ EQ 0)]

	xmin = MIN(dataAll.xprop)
	xmax = MAX(dataAll.xprop)
	xr = xmax-xmin

	ymin = MIN(dataAll.yprop)
	ymax = MAX(dataAll.yprop)
	yr = ymax-ymin
	
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsampleFBF/COSMOS_projection_test', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=8
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	theta = 2*!PI*FINDGEN(17)/16.
	W = 0.25
	USERSYM, W*COS(theta), W*SIN(theta)

	PLOT, [1], [1], xrange=[xmin,xmax], yrange=[ymin,ymax], /NODATA
	OPLOT, dataAllSF.xprop, dataAllSF.yprop, PSYM=8, color=SFcolor
	OPLOT, dataAllQ.xprop, dataAllQ.yprop, PSYM=8, color=Qcolor

	PLOT, [1], [1], xrange=[zmin,zmax], yrange=[ymin,ymax], /NODATA
	OPLOT, dataAllSF.zprimus, dataAllSF.yprop, PSYM=8, color=SFcolor
	OPLOT, dataAllQ.zprimus, dataAllQ.yprop, PSYM=8, color=Qcolor

	PLOT, [1], [1], xrange=[zmin,zmax], yrange=[xmin,xmax], /NODATA
	OPLOT, dataAllSF.zprimus, dataAllSF.xprop, PSYM=8, color=SFcolor
	OPLOT, dataAllQ.zprimus, dataAllQ.xprop, PSYM=8, color=Qcolor

;	OPLOT, dataIPSF.xprop, dataIPSF.yprop, PSYM=4, color=cgcolor('cyan')
;	OPLOT, dataIPQ.xprop, dataIPQ.yprop, PSYM=4, color=cgcolor('magenta')

;	PLOT, [1], [1], xrange=[xmin,xmax], yrange=[ymin,ymax], title='SF IP', /NODATA
;	plotNeigh, dataIPSF, dataAll, dz_coeff, SFcolor, Qcolor;, neighSFdistsAll, neighQdistsAll, dz_coeff
;	FOR i=1,8 DO OPLOT, 2.*i*COS(theta), 2.*i*SIN(theta), LINESTYLE=0, THICK=3
;	FOR i=1,8 DO OPLOT, (2.*i-1)*COS(theta), (2.*i-1)*SIN(theta), LINESTYLE=2, THICK=2

;	PLOT, [1], [1], xrange=[xmin,xmax], yrange=[ymin,ymax], title='Q IP', /NODATA
;	plotNeigh, dataIPQ, dataAll, dz_coeff, SFcolor, Qcolor;, neighSFdistsAll, neighQdistsAll, dz_coeff
;	FOR i=1,8 DO OPLOT, 2.*i*COS(theta), 2.*i*SIN(theta), LINESTYLE=0, THICK=3
;	FOR i=1,8 DO OPLOT, (2.*i-1)*COS(theta), (2.*i-1)*SIN(theta), LINESTYLE=2, THICK=2

;	getDists, dataIPSF, dataAll, neighSFdistsAll, neighQdistsAll, dz_coeff
;	PLOTHIST, neighSFdistsAll, bin=1, color=SFcolor, xrange=[0,30];, yrange=[0,22000]
;	PLOTHIST, neighQdistsAll, bin=1, color=Qcolor, /OVERPLOT
;	SFdists = HISTOGRAM(neighSFdistsAll, binsize=1)
;	Qdists  = HISTOGRAM(neighQdistsAll, binsize=1)

;	OPLOT, 0.5+FINDGEN(30), FLOAT(SFdists)/(SFdists+Qdists), LINESTYLE=0, color=SFcolor

;	getDists, dataIPQ, dataAll, neighSFdistsAll, neighQdistsAll, dz_coeff
;	SFdists = HISTOGRAM(neighSFdistsAll, binsize=1)
;	Qdists  = HISTOGRAM(neighQdistsAll, binsize=1)
;	OPLOT, 0.5+FINDGEN(30), FLOAT(SFdists)/(SFdists+Qdists), LINESTYLE=0, color=Qcolor

;	getDists, dataIPSF, dataAll, neighSFdistsAll, neighQdistsAll, dz_coeff
;	PLOTHIST, neighSFdistsAll, bin=1, color=SFcolor, xrange=[0,30], yrange=[0,22000]
;	PLOTHIST, neighQdistsAll, bin=1, color=Qcolor, /OVERPLOT

;	getDists, dataIPQ, dataAll, neighSFdistsAll, neighQdistsAll, dz_coeff
;	PLOTHIST, neighSFdistsAll, bin=1, color=SFcolor, LINESTYLE=2, /OVERPLOT
;	PLOTHIST, neighQdistsAll, bin=1, color=Qcolor, LINESTYLE=2, /OVERPLOT

;	XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'SF IP: ' + bigint(n_elements(UNIQ(dataIP[where(dataIP.SFQ EQ 1)].objname))), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Q IP: ' + bigint(n_elements(dataIP[where(dataIP.SFQ EQ 0)])), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, zlabel, ALIGNMENT=1.0

;	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0
;	XYOUTS, xmin+0.05*xr, ymin+0.925*yr, 'COSMOS only', ALIGNMENT=0.0

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
