; takes given IP sample (conservative, variable_mass+1.0dex, etc.) and selects SF and Q IP samples $
; with matching mean mstar and redshift

; currently using M13 mass limit

PRO matchSamplesAllFields, outputFormat
	seed = 1
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsamples_M13highMassCut', /ENCAP, THICK=3
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	dataAll = mrdfits('~/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits',1)
	data = dataAll[where(dataAll.IP eq 1)]
		
	!P.MULTI = 0
	!P.CHARSIZE = 1.1

;	!x.margin=[6,2.5] ; left, right
;	!y.margin=[3.5,2.5] ; bottom, top

	xmin = 0.2
	xmax = 1.0
	xr = xmax-xmin

	ymin = 9.1
	ymax = 11.6
	yr = ymax-ymin

;	IF outputFormat eq 'ps' THEN axisColor = cgColor('Black') ELSE axisColor = cgColor('Black')
	Qcolor  = cgColor('Red')
	SFcolor = cgColor('Blue')

	mass_binwidth 	= 0.15
        z_binwidth	= 0.05
	n_levels 	= 10	
	
	xtitle = 'Redshift'	
	ytitle = textoidl('Log (M / M_{\odot})')

	zmin = 0.2
	zmax = 1.0
	zArray = zmin + z_binwidth*findgen((zmax-zmin)/z_binwidth + 1)

	PLOT, findgen(10), findgen(10), psym=3, xrange=[xmin, xmax], yrange=[ymin, ymax], xtitle=xtitle, ytitle=ytitle, title='Matched SF and Q IP Samples', /NODATA

	dataSF	= data[where(data.SFQ EQ 1)]
	dataQ	= data[where(data.SFQ EQ 0)]

;	OPLOT, dataSF.zprimus, dataSF.mstar, PSYM=3, COLOR=SFcolor
;	OPLOT, dataQ.zprimus, dataQ.mstar, PSYM=3, COLOR=Qcolor

	cont = getContours(data, 1, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	CONTOUR, cont, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=2, COLOR=SFcolor, /OVERPLOT
	cont = getContours(data, 0, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	CONTOUR, cont, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=2, COLOR=Qcolor, /OVERPLOT

        mass    = [10.8,        10.9,   11.0,   11.1,   11.2,   11.3,   11.4]
        mass_z2 = [10.8,        10.9,   11.0,   11.1,           11.3]
        phi_z2  = [-2.894,	-3.180, -3.380, -3.370,         -4.600]
        phi_z3  = [-2.865,	-2.905, -3.113, -3.430, -3.670, -4.040, -4.150]
        phi_z4  = [-2.762,	-2.978, -3.117, -3.347, -3.540, -3.830, -4.230]
        phi_z5  = [-2.185,	-2.980, -3.139, -3.368, -3.539, -3.930, -4.080]
        phi_z65 = [-2.744,	-2.788, -3.000, -3.148, -3.286, -3.671, -4.040]
        phi_z8  = [-2.908,	-3.011, -3.113, -3.279, -3.405, -3.630, -3.920]

        phiArray = [[phi_z3], [phi_z4], [phi_z5], [phi_z65], [phi_z8]]

        targPhi = -3.9
        mArray	= [INTERPOL(mass_z2, phi_z2, targPhi)]
	zarray	= [0.25, 0.35, 0.45, 0.55, 0.725, 0.9]

        FOR i=0,n_elements(zArray)-2 DO mArray = [mArray, INTERPOL(mass, phiArray[*,i], targPhi)]

        massFunc = INTERPOL(mArray, zArray, data.zprimus)

;	highMass = [11.0, 11.0, 11.0, 11.1, 11.2, 11.4]

;	dataQcut  = dataQ[where(dataQ.mstar LE highcut(dataQ.zprimus))]
	dataQcut  = dataQ[where(dataQ.mstar LE INTERPOL(mArray, zarray, dataQ.zprimus))]
;	dataSFcut = dataSF[where(dataSF.mstar LE highcut(dataSF.zprimus))]
	dataSFcut = dataSF[where(dataSF.mstar LE INTERPOL(mArray, zarray, dataSF.zprimus))]

	OPLOT, dataQcut.zprimus, dataQcut.mstar, PSYM=3, COLOR=Qcolor
	OPLOT, dataSFcut.zprimus, dataSFcut.mstar, PSYM=3, COLOR=SFcolor

	contQ = getContours(dataQcut, 0, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
;	normContQ = contQ/float(total(contQ))
;	FOR i=0,n_elements(normContQ[*,0])-1 DO PRINT, strtrim(string(normContQ[*,n_elements(normContQ[*,0])-1-i],format='(f20.7)'),1)
;	CONTOUR, contQ, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=3, COLOR=Qcolor, /OVERPLOT

	contSF = getContours(dataSFcut, 1, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
;	normContSF = contSF/float(total(contSF))
;	FOR i=0,n_elements(contSF[*,0])-1 DO PRINT, ypts[i], strtrim(string(contSF[*,i],format='(f20.7)'),1)
;	FOR i=0,n_elements(normContSF[*,0])-1 DO PRINT, ypts[i], strtrim(string(normContSF[*,n_elements(normContSF[*,0])-1-i],format='(f20.7)'),1)
;	CONTOUR, contSF, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=3, COLOR=SFcolor, /OVERPLOT

;	[ ROW, COLUMN ]
;	ROW: 0.2 + 0.05 to 1.0
;	COL: 9.1 + 0.15 to 11.5
	n_rows = ceil((1.0-0.2)/z_binwidth)
	c_cols = ceil((11.5-9.1)/mass_binwidth)

	QsampleAll	= []
	SFselectAll	= []
	index = 0
	FOR r=0,n_elements(contQ[0,*])-1 DO BEGIN
		min_m = 9.1 + mass_binwidth*r
		max_m = 9.1 + mass_binwidth*(r+1)

		FOR c=0,n_elements(contQ[*,0])-1 DO BEGIN
			index += 1
			min_z = 0.2 + z_binwidth*c
			max_z = 0.2 + z_binwidth*(c+1)

			Qsample  = dataQcut[WHERE( (dataQcut.zprimus GE min_z) AND (dataQcut.zprimus LE max_z) AND $
				(dataQcut.mstar GE min_m) AND (dataQcut.mstar LE max_m), /NULL )]
			SFsample = dataSFcut[WHERE( (dataSFcut.zprimus GE min_z) AND (dataSFcut.zprimus LE max_z) AND $
				(dataSFcut.mstar GE min_m) AND (dataSFcut.mstar LE max_m), /NULL )]
			targetNum = n_elements(Qsample)

			IF (targetNum NE 0) THEN BEGIN
				SFselect	= SFsample[ROUND(n_elements(SFsample)*RANDOMU(seed,targetNum))]
				SFselectAll	= [SFselectAll, SFselect]
				QsampleAll	= [QsampleAll, Qsample]
			ENDIF ELSE BEGIN
				SFselect = !NULL
			ENDELSE

			IF (index mod 10 EQ 0) THEN PRINT, min_z, ' ', min_m, ' targetNum: ', strcompress(targetNum), '  SFselect: ', strcompress(n_elements(SFselect)), ' of ', strcompress(n_elements(SFsample))
		ENDFOR
	ENDFOR

	PRINT, n_elements(QsampleAll), n_elements(SFselectAll)
	PRINT, minmax(QsampleAll.mstar), median(QsampleAll.mstar)
	PRINT, minmax(SFselectAll.mstar), median(SFselectAll.mstar)
	PRINT, minmax(QsampleAll.zprimus), median(QsampleAll.zprimus)
	PRINT, minmax(SFselectAll.zprimus), median(SFselectAll.zprimus)

	XYOUTS, xmin+0.95*xr, ymin+0.17*yr, 'median Q redshift: ' + twoDecimal(median(QsampleAll.zprimus)), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.13*yr, 'median SF redshift: ' + twoDecimal(median(SFselectAll.zprimus)), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.09*yr, 'median Q IP mass: ' + twoDecimal(median(QsampleAll.mstar)), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, 'median SF IP mass: ' + twoDecimal(median(SFselectAll.mstar)), ALIGNMENT=1.0

	; contours for Q sample within matching criteria
;	contQ = getContours(QsampleAll, 0, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
;	CONTOUR, contQ, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=2, /OVERPLOT

	; contours for matched SF sample
	contSF = getContours(SFselectAll, 1, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	CONTOUR, contSF, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=5, COLOR=cgColor('Black'), /OVERPLOT

	zz = 0.2 + 0.01*FINDGEN(80)
;	OPLOT, zz, highcut(zz), color=cgcolor('White'), LINESTYLE=2, THICK=3
	OPLOT, data[SORT(data.zprimus)].zprimus, massFunc[SORT(massFunc)], color=cgcolor('Green'), LINESTYLE=2, THICK=3
	
;	OPLOT, QsampleAll.zprimus, QsampleAll.mstar, PSYM=3, COLOR=cgColor('Orange')
;	OPLOT, SFselectAll.zprimus, SFselectAll.mstar, PSYM=3, COLOR=cgColor('Cyan')

	allIP = [QsampleAll, SFselectAll]

	MWRFITS, allIP, '~/results/match_IP_sample_rigorous/matchedIPsample_M13highMassCut.fits', /CREATE

;	FOR i=0,n_elements(xpts)-1 DO OPLOT, [xpts[i],xpts[i]], [min(ypts),max(ypts)], LINESTYLE=1, THICK=2
;	FOR j=0,n_elements(ypts)-1 DO OPLOT, [min(xpts),max(xpts)], [ypts[j],ypts[j]], LINESTYLE=1, THICK=2

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'	
END

FUNCTION getContours, data, SFstatus, n_levels, levels, xbinwidth, ybinwidth, xpts, ypts
;	xbinwidth = 0.05
;	ybinwidth = 0.2

	Qcolor  = cgColor('Red')
	SFcolor = cgColor('Blue')
	IF SFstatus EQ 1 THEN SFQcolor=SFcolor ELSE SFQcolor=Qcolor

	dataSFQ = data[where(data.SFQ eq SFstatus)]
	mass	= dataSFQ.mstar
	z	= dataSFQ.zprimus

	mass_width 	= ybinwidth
	z_width	= xbinwidth
		
;	n_massbins = ceil((max(mass)-min(mass))/mass_width)
	n_massbins = ceil((11.5-9.1)/mass_width)
        n_zbins    = floor((1.0-0.2)/z_width)

	hist2d = []
	FOR s_index=0,n_zbins DO BEGIN
		z_low  = 0.2 + z_width*float(s_index)
        	z_high = 0.2 + z_width*float(s_index+1)

		in_z_range = dataSFQ[where((dataSFQ.zprimus GE z_low) AND (dataSFQ.zprimus LE z_high), /NULL)]

                z_row = []
                FOR m_index=0,n_massbins DO BEGIN
;			m_low   = min(mass) + mass_width*float(m_index)
;			m_high  = min(mass) + mass_width*float(m_index+1)
			m_low   = 9.1 + mass_width*float(m_index)
			m_high  = 9.1 + mass_width*float(m_index+1)

                        IF in_z_range NE !NULL THEN $
                        in_unit = n_elements(in_z_range[where((in_z_range.mstar GE m_low) AND (in_z_range.mstar LE m_high), /NULL)]) $
                        ELSE in_unit = 0

                        z_row = [z_row, in_unit]
                ENDFOR

                hist2d = [[hist2d], [z_row]]
        ENDFOR

;	masspts	= min(mass) + mass_width*(findgen(n_massbins+1))
	masspts	= 9.1 + mass_width*(findgen(n_massbins+1))
	zpts	= 0.2 + z_width*(findgen(n_zbins+1))
	
;	PRINT, 'mass bins: ', n_massbins
;	PRINT, 'mass pts:  ', n_elements(masspts)
;	PRINT, 'z bins:    ', n_zbins
;	PRINT, 'z pts:     ', n_elements(zpts)

	xpts = zpts
	ypts = masspts	

	levels = float(max(hist2d))*findgen(n_levels)/(n_levels-1)

;	PRINT, minmax(hist2d)
;	PRINT, 'Levels: ', levels

	RETURN, TRANSPOSE(hist2d)
END

FUNCTION highcut, zz
	z0 	= 0.9
	a 	= -1.2
	b 	= 0
	c	= 11.4
;	PRINT, twoDecimal(a) + '*(z - ' + twoDecimal(z0) + ')^2 + ' + twoDecimal(b) + '*z + ' + twoDecimal(c)
	RETURN, a*(zz - z0)^2 + b*zz + c
END

FUNCTION twoDecimal, input
	RETURN, strtrim(string(input, format='(f20.2)'),1)
END
