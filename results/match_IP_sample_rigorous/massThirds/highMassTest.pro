PRO highMassTest, outputFormat
; read in data files
        dataIPallFields  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)
	dataAllallFields = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

; eliminate data with targ_weight < 1
        dataIPallFields  = dataIPallFields[where(dataIPallFields.targ_weight GE 1.)]
        dataAllallFields  = dataAllallFields[where(dataAllallFields.targ_weight GE 1.)]

; divide IP population into 3 mass bins with equal numbers of IPs
	lower_index = CEIL(n_elements(dataIPallFields)/3.)
	upper_index = 2*FLOOR(n_elements(dataIPallFields)/3.)

	orderedMasses = dataIPallFields[SORT(dataIPallFields.mstar)].mstar

	lowerMass = orderedMasses(lower_index)
	upperMass = orderedMasses(upper_index)

	massArray = [min(dataIPallFields.mstar), lowerMass, upperMass, max(dataIPallFields.mstar)]
	PRINT, massArray

; isolate highest mass third
;	dataAllallFields = dataAllallFields[where(dataAllallFields.mstar GE upperMass)]
	dataIPallFields  = dataIPallFields[where(dataIPallFields.mstar GE upperMass)]

        dataIPallFieldsSF = dataIPallFields[where(dataIPallFields.SFQ EQ 1)]
        dataIPallFieldsQ  = dataIPallFields[where(dataIPallFields.SFQ EQ 0)]
	PRINT, n_elements(dataIPallFieldsSF), n_elements(dataIPallFieldsQ)

        SFcolor = cgColor('blue')
        Qcolor  = cgColor('red')

	!p.multi=[15,5,3]
        !p.charsize=2

;	xmin = MIN(dataAll.xprop)
;	xmax = MAX(dataAll.xprop)
	xmin = -5
	xmax = 5
	xr = xmax-xmin

;	ymin = MIN(dataAll.yprop)
;	ymax = MAX(dataAll.yprop)
	ymin = -5
	ymax = 5
	yr = ymax-ymin

	IF (string(outputFormat) eq 'ps') THEN BEGIN
                PS_OPEN, '~/figures/matchedIPsampleFBF/highMassTest', THICK=3, /ENCAP
                DEVICE, /INCH, XS=8, YS=8
        ENDIF ELSE BEGIN
                SET_PLOT, 'X'
                THICK=1
        ENDELSE

        theta = 2*!PI*FINDGEN(17)/16.
        W = 0.25
        USERSYM, W*COS(theta), W*SIN(theta)

	fields = ['cdfs', 'cfhtls_xmm', 'cosmos', 'es1', 'xmm']
;	fields = ['cdfs', 'cfhtls_xmm', 'cosmos', 'xmm']
;	fields = ['cfhtls_xmm', 'cosmos', 'xmm']

	ERASE
	FOR j=0,n_elements(fields)-1 DO BEGIN
		getDists, fields[j], dataIPallFields, dataAllallFields, neighSFdistsAll, neighQdistsAll, 2.0, SFcolor, Qcolor
	ENDFOR
	FOR j=0,n_elements(fields)-1 DO BEGIN
		PLOT, [1], [1], xrange=[xmin,xmax], yrange=[ymin,ymax], title='Late-type IP', xtitle='Relative Mpc', ytitle='Relative Mpc', /NODATA
		plotNeigh, fields[j], dataIPallFieldsSF, dataAllallFields, 2.0, SFcolor, Qcolor;, neighSFdistsAll, neighQdistsAll, dz_coeff
	ENDFOR
	FOR j=0,n_elements(fields)-1 DO BEGIN
		PLOT, [1], [1], xrange=[xmin,xmax], yrange=[ymin,ymax], title='Early-type IP', xtitle='Relative Mpc', ytitle='Relative Mpc', /NODATA
		plotNeigh, fields[j], dataIPallFieldsQ, dataAllallFields, 2.0, SFcolor, Qcolor;, neighSFdistsAll, neighQdistsAll, dz_coeff
	ENDFOR

;	OPLOT, dataAllSF.xprop, dataAllSF.yprop, PSYM=8, color=SFcolor
;	OPLOT, dataAllQ.xprop, dataAllQ.yprop, PSYM=8, color=Qcolor
;	OPLOT, dataIPSF.xprop, dataIPSF.yprop, PSYM=4, color=cgcolor('cyan')
;	OPLOT, dataIPQ.xprop, dataIPQ.yprop, PSYM=4, color=cgcolor('magenta')
END

PRO getDists, field, IPdata, dataAll, neighSFdistsAll, neighQdistsAll, dz_coeff, SFcolor, Qcolor
	IPSFneighDistsAll   = []
	IPSFneighSFdistsAll = []
	IPSFneighQdistsAll  = []
	IPQneighDistsAll   = []
	IPQneighSFdistsAll = []
	IPQneighQdistsAll  = []

	fieldIPSF = IPdata[where( (strtrim(IPdata.field,2) EQ field) AND (IPdata.SFQ EQ 1) )]
	FOR i=0,n_elements(fieldIPSF)-1 DO BEGIN
		currentIP = fieldIPSF[i]

                ; keep galaxies within dz=dz+coeff*0.005*(1+z) of IP candidate (and in same field)
                neigh = dataAll[ where( (dataAll.objname NE currentIP.objname) AND (strtrim(dataAll.field,2) EQ field) AND $
			(ABS(dataAll.zprimus - currentIP.zprimus) LE dz_coeff*0.005*(1.+currentIP.zprimus)) ) ]
                neighDists = [SQRT((neigh.xprop-currentIP.xprop)^2 + (neigh.yprop-currentIP.yprop)^2)]

                neighSF = neigh[where(neigh.SFQ EQ 1)]
                neighQ  = neigh[where(neigh.SFQ EQ 0)]
                neighSFdists = [SQRT((neighSF.xprop-currentIP.xprop)^2 + (neighSF.yprop-currentIP.yprop)^2)]
                neighQdists = [SQRT((neighQ.xprop-currentIP.xprop)^2 + (neighQ.yprop-currentIP.yprop)^2)]

		IPSFneighDistsAll   = [IPSFneighDistsAll, neighDists]
		IPSFneighSFdistsAll = [IPSFneighSFdistsAll, neighSFdists]
		IPSFneighQdistsAll  = [IPSFneighQdistsAll, neighQdists]
        ENDFOR

	fieldIPQ  = IPdata[where( (strtrim(IPdata.field,2) EQ field) AND (IPdata.SFQ EQ 0) )]
	FOR i=0,n_elements(fieldIPQ)-1 DO BEGIN
		currentIP = fieldIPQ[i]

                ; keep galaxies within dz=dz+coeff*0.005*(1+z) of IP candidate (and in same field)
                neigh = dataAll[ where( (dataAll.objname NE currentIP.objname) AND (strtrim(dataAll.field,2) EQ field) AND $
			(ABS(dataAll.zprimus - currentIP.zprimus) LE dz_coeff*0.005*(1.+currentIP.zprimus)) ) ]
                neighDists = [SQRT((neigh.xprop-currentIP.xprop)^2 + (neigh.yprop-currentIP.yprop)^2)]

                neighSF = neigh[where(neigh.SFQ EQ 1)]
                neighQ  = neigh[where(neigh.SFQ EQ 0)]

                neighSFdists = [SQRT((neighSF.xprop-currentIP.xprop)^2 + (neighSF.yprop-currentIP.yprop)^2)]
                neighQdists = [SQRT((neighQ.xprop-currentIP.xprop)^2 + (neighQ.yprop-currentIP.yprop)^2)]

		IPQneighDistsAll   = [IPQneighDistsAll, neighDists]
		IPQneighSFdistsAll = [IPQneighSFdistsAll, neighSFdists]
		IPQneighQdistsAll  = [IPQneighQdistsAll, neighQdists]
        ENDFOR

;	IPSFtotals 	= HISTOGRAM(IPSFneighSFdistsAll, bin=1) + HISTOGRAM(IPSFneighQdistsAll, bin=1)
	IPSFtotals 	= HISTOGRAM(IPSFneighDistsAll, bin=1)	
	IPSFlatefrac	= HISTOGRAM(IPSFneighSFdistsAll, bin=1)/FLOAT(IPSFtotals)
;	IPQtotals 	= HISTOGRAM(IPQneighSFdistsAll, bin=1) + HISTOGRAM(IPQneighQdistsAll, bin=1)
	IPQtotals 	= HISTOGRAM(IPQneighDistsAll, bin=1)	
	IPQlatefrac	= HISTOGRAM(IPQneighSFdistsAll, bin=1)/FLOAT(IPQtotals)

	PRINT, field, IPSFtotals[0:4]
	PRINT, field, IPSFlatefrac[0:4]
	PRINT, field, IPQtotals[0:4]
	PRINT, field, IPQlatefrac[0:4]
	PRINT, ''

	PLOT, 0.5+FINDGEN(5), IPSFlatefrac, xrange=[0,5], yrange=[0.25,0.5], LINESTYLE=0, xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction'
	XYOUTS, 2.5, 0.48, strtrim(field,2), ALIGNMENT=0.5, charsize=1.5
	IF (strtrim(field,2) EQ 'cdfs') THEN LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[2,0], BOX=0, /BOTTOM, /LEFT, charsize=1.25
	OPLOT, 0.5+FINDGEN(5), IPQlatefrac, LINESTYLE=2

;	PLOTHIST, IPSFneighSFdistsAll, bin=1.0, xrange=[0,5], color=SFcolor, LINESTYLE=0
;	PLOTHIST, IPSFneighQdistsAll, bin=1.0, color=Qcolor, LINESTYLE=0, /OVERPLOT
;	PLOTHIST, IPQneighSFdistsAll, bin=1.0, color=SFcolor, LINESTYLE=2, /OVERPLOT
;	PLOTHIST, IPQneighQdistsAll, bin=1.0, color=Qcolor, LINESTYLE=2, /OVERPLOT
END

PRO plotNeigh, field, IPdata, dataAll, dz_coeff, SFcolor, Qcolor;, neighSFdistsAll, neighQdistsAll, dz_coeff
;	IPdata  = IPdata[where(strtrim(IPdata.field,2) EQ field, /NULL)]
;	dataAll = dataAll[where(strtrim(dataAll.field,2) EQ field, /NULL)]

;	PRINT, strtrim(field,2), n_elements(IPdata), n_elements(dataAll)

        theta = 2*!PI*FINDGEN(33)/32.
        W = 0.25
        USERSYM, W*COS(theta), W*SIN(theta)

	fieldIP = IPdata[where(strtrim(IPdata.field,2) EQ field)]
        FOR i=0,n_elements(fieldIP)-1 DO BEGIN
                currentIP = fieldIP[i]

                ; keep galaxies within dz=dz+coeff*0.005*(1+z) of IP candidate (and in same field)
                neigh = dataAll[ where( (dataAll.objname NE currentIP.objname) AND (strtrim(dataAll.field,2) EQ field) AND $
                        (ABS(dataAll.zprimus - currentIP.zprimus) LE dz_coeff*0.005*(1.+currentIP.zprimus)) ) ]

                neighSF = neigh[where(neigh.SFQ EQ 1)]
                neighQ  = neigh[where(neigh.SFQ EQ 0)]
		
;		IF (i EQ 0) THEN PLOT, neighSF.xprop-currentIP.xprop, neighSF.yprop-currentIP.yprop, COLOR=SFcolor, psym=8 ELSE $
                OPLOT, [neighSF.xprop]-currentIP.xprop, [neighSF.yprop]-currentIP.yprop, COLOR=SFcolor, psym=8
                OPLOT, [neighQ.xprop]-currentIP.xprop, [neighQ.yprop]-currentIP.yprop, COLOR=Qcolor, psym=8

		FOR j=1,5 DO OPLOT, j*COS(theta), j*SIN(theta), LINESTYLE=2, THICK=2

;               OPLOT, ABS(neighQ.xprop-currentIP.xprop), ABS(neighQ.yprop-currentIP.yprop), COLOR=Qcolor, psym=8
;               OPLOT, ABS(neighSF.xprop-currentIP.xprop), ABS(neighSF.yprop-currentIP.yprop), COLOR=SFcolor, psym=8
        ENDFOR
END

