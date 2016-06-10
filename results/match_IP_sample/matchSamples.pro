; takes given IP sample (conservative, variable_mass+1.0dex, etc.) and selects SF and Q IP samples $
; with matching mean mstar and redshift

; currently using M13 mass limit

PRO matchSamples, outputFormat
	seed = 1
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsamples', /ENCAP, THICK=3
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	dataAll = mrdfits('~/conformity/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits',1)
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

;	IF outputFormat eq 'ps' THEN axisColor = cgColor('Black') ELSE axisColor = cgColor('White')
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

	cont = getContours(data, 1, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	CONTOUR, cont, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=2, COLOR=SFcolor, /OVERPLOT
	cont = getContours(data, 0, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	CONTOUR, cont, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=2, COLOR=Qcolor, /OVERPLOT

	x = findgen(11)/10
	OPLOT, x, highcut(x), LINESTYLE=2, THICK=3
	FOREACH xx,x DO PRINT, xx, highcut(xx)

	dataSF	= data[where(data.SFQ EQ 1)]
	dataQ	= data[where(data.SFQ EQ 0)]
	dataQcut  = dataQ[where(dataQ.mstar LE highcut(dataQ.zprimus))]
	dataSFcut = dataSF[where(dataSF.mstar LE highcut(dataSF.zprimus))]
;	OPLOT, dataQ.zprimus, dataQ.mstar, PSYM=3, COLOR=Qcolor
;	OPLOT, dataSF.zprimus, dataSF.mstar, PSYM=3, COLOR=SFcolor

	contQ = getContours(dataQcut, 0, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	normContQ = contQ/float(total(contQ))
;	FOR i=0,n_elements(normContQ[*,0])-1 DO PRINT, strtrim(string(normContQ[*,n_elements(normContQ[*,0])-1-i],format='(f20.7)'),1)
;	CONTOUR, contQ, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=3, COLOR=Qcolor, /OVERPLOT

	contSF = getContours(dataSFcut, 1, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	normContSF = contSF/float(total(contSF))
;	PRINT, '             ', strtrim(string(xpts,format='(f20.7)'),2)
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

;			IF (index mod 10 EQ 0) THEN PRINT, min_z, ' ', min_m, ' targetNum: ', strcompress(targetNum), '  SFselect: ', strcompress(n_elements(SFselect)), ' of ', strcompress(n_elements(SFsample))
		ENDFOR
	ENDFOR

	PRINT, n_elements(QsampleAll), n_elements(SFselectAll)
	PRINT, minmax(QsampleAll.mstar), median(QsampleAll.mstar)
	PRINT, minmax(SFselectAll.mstar), median(SFselectAll.mstar)
	PRINT, minmax(QsampleAll.zprimus), median(QsampleAll.zprimus)
	PRINT, minmax(SFselectAll.zprimus), median(SFselectAll.zprimus)

	XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'median Q redshift: ' + twoDecimal(median(QsampleAll.zprimus)), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'median SF redshift: ' + twoDecimal(median(SFselectAll.zprimus)), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'median Q IP mass: ' + twoDecimal(median(QsampleAll.mstar)), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, 'median SF IP mass: ' + twoDecimal(median(SFselectAll.mstar)), ALIGNMENT=1.0

	contQ = getContours(QsampleAll, 0, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	CONTOUR, contQ, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=2, /OVERPLOT

	contSF = getContours(SFselectAll, 1, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	CONTOUR, contSF, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=5, COLOR=cgColor('Black'), /OVERPLOT

;	OPLOT, QsampleAll.zprimus, QsampleAll.mstar, PSYM=3, COLOR=cgColor('Orange')
;	OPLOT, SFselectAll.zprimus, SFselectAll.mstar, PSYM=3, COLOR=cgColor('Cyan')

;	allIP = [QsampleAll, SFselectAll]

;	MWRFITS, allIP, 'matchedIPsample.fits', /CREATE

;	FOR i=0,n_elements(xpts)-1 DO OPLOT, [xpts[i],xpts[i]], [min(ypts),max(ypts)], LINESTYLE=1, THICK=2
;	FOR j=0,n_elements(ypts)-1 DO OPLOT, [min(xpts),max(xpts)], [ypts[j],ypts[j]], LINESTYLE=1, THICK=2

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'	
END

FUNCTION getContours, data, SFstatus, n_levels, levels, xbinwidth, ybinwidth, xpts, ypts
	Qcolor  = cgColor('Red')
	SFcolor = cgColor('Blue')
	IF SFstatus EQ 1 THEN SFQcolor=SFcolor ELSE SFQcolor=Qcolor

	dataSFQ = data[where(data.SFQ eq SFstatus)]
	mass	= dataSFQ.mstar
	z	= dataSFQ.zprimus

	mass_binwidth 	= ybinwidth
	z_binwidth	= xbinwidth
		
;	n_massbins = ceil((max(mass)-min(mass))/mass_binwidth)
	n_massbins = ceil((11.5-9.1)/mass_binwidth)
        n_zbins    = floor((1.0-0.2)/z_binwidth)

	hist2d = []
	FOR s_index=0,n_zbins DO BEGIN
		z_low  = 0.2 + z_binwidth*float(s_index)
        	z_high = 0.2 + z_binwidth*float(s_index+1)

		in_z_range = dataSFQ[where((dataSFQ.zprimus GE z_low) AND (dataSFQ.zprimus LE z_high), /NULL)]

                z_row = []
                FOR m_index=0,n_massbins DO BEGIN
;			m_low   = min(mass) + mass_binwidth*float(m_index)
;			m_high  = min(mass) + mass_binwidth*float(m_index+1)
			m_low   = 9.1 + mass_binwidth*float(m_index)
			m_high  = 9.1 + mass_binwidth*float(m_index+1)

                        IF in_z_range NE !NULL THEN $
                        in_unit = n_elements(in_z_range[where((in_z_range.mstar GE m_low) AND (in_z_range.mstar LE m_high), /NULL)]) $
                        ELSE in_unit = 0

                        z_row = [z_row, in_unit]
                ENDFOR

                hist2d = [[hist2d], [z_row]]
        ENDFOR

;	masspts	= min(mass) + mass_binwidth*(findgen(n_massbins+1))
	masspts	= 9.1 + mass_binwidth*(findgen(n_massbins+1))
	zpts	= 0.2 + z_binwidth*(findgen(n_zbins+1))
	
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
