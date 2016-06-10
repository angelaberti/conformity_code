; FBF = field-by-field (median redshift and stellar mass matching of SF and Q IP populations done within each field)

; takes given IP sample (conservative, variable_mass+1.0dex, etc.) and selects SF and Q IP samples $
; with matching mean mstar and redshift

; currently using M13 mass limit

PRO matchSamples_FBF, outputFormat
	seed = 1
	dataAll = mrdfits('~/conformity/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits',1)
	dataIPallFields = dataAll[where(dataAll.IP eq 1)]

	fields = dataAll[uniq(dataAll.field, sort(dataAll.field))].field
		
	QsampleAllFields =  []
	SFselectAllFields = []

  FOR f=0,n_elements(fields)-1 DO BEGIN
	data = dataIPallFields[WHERE(dataIPallFields.field EQ fields[f])]
	PRINT, 'FIELD: ' + strcompress(fields[f])
	PRINT, 'IPs: ' + strcompress(n_elements(data))

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsamples_' + strcompress(strtrim(fields[f],2)), /ENCAP, THICK=3
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	!P.MULTI = 0
	!P.CHARSIZE = 1.1

;	!x.margin=[6,2.5] ; left, right
;	!y.margin=[3.5,2.5] ; bottom, top

	xmin = 0.2
	xmax = 1.0
;	xmin = 0.5
;	xmax = 0.55

	xr = xmax-xmin

	ymin = 9.1
	ymax = 11.6
;	ymin = 11.2
;	ymax = 11.35
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
	; PLOT GRID
;	FOR i=0,ceil((zmax-zmin)/z_binwidth) DO OPLOT, [zmin+i*z_binwidth, zmin+i*z_binwidth], [ymin,ymax], linestyle=2, thick=1	
;	FOR j=0,ceil((ymax-ymin)/mass_binwidth) DO OPLOT, [zmin,zmax], [ymin+j*mass_binwidth, ymin+j*mass_binwidth], linestyle=2, thick=1	

	cont = getContours(data, 1, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	CONTOUR, cont, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=2, COLOR=SFcolor, /OVERPLOT
	cont = getContours(data, 0, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	CONTOUR, cont, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=2, COLOR=Qcolor, /OVERPLOT

	x = findgen(11)/10
	OPLOT, x, highcut(x), LINESTYLE=2, THICK=2
	FOREACH xx,x DO PRINT, xx, highcut(xx)

	dataSF	= data[where(data.SFQ EQ 1)]
	dataQ	= data[where(data.SFQ EQ 0)]

;	OPLOT, dataSF.zprimus, dataSF.mstar, psym=3, color=cgColor('Cyan')
;	OPLOT, dataQ.zprimus, dataQ.mstar, psym=3, color=cgColor('Red')

	dataQcut  = dataQ[where(dataQ.mstar LE highcut(dataQ.zprimus))]
	dataSFcut = dataSF[where(dataSF.mstar LE highcut(dataSF.zprimus))]
	PRINT, 'num elements dataSFcut: ', n_elements(dataSFcut)

	contQ = getContours(dataQcut, 0, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	normContQ = contQ/float(total(contQ))

	contSF = getContours(dataSFcut, 1, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	normContSF = contSF/float(total(contSF))

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

			IF (n_elements(SFsample) EQ 0) OR (targetNum EQ 0) THEN BEGIN
				SFselect = !NULL
				Qsample  = !NULL
			ENDIF ELSE BEGIN
				SFselect = SFsample[ROUND(n_elements(SFsample)*RANDOMU(seed,targetNum))]
			ENDELSE

			SFselectAll	= [SFselectAll, SFselect]
			QsampleAll	= [QsampleAll, Qsample]

;			IF (index mod 10 EQ 0) THEN PRINT, min_z, ' ', min_m, ' targetNum: ', strcompress(targetNum), '  SFselect: ', strcompress(n_elements(SFselect)), ' of ', strcompress(n_elements(SFsample))
		ENDFOR
	ENDFOR

	PRINT, fields[f], n_elements(QsampleAll), n_elements(SFselectAll)
;	PRINT, minmax(QsampleAll.mstar), median(QsampleAll.mstar)
;	PRINT, minmax(SFselectAll.mstar), median(SFselectAll.mstar)
;	PRINT, minmax(QsampleAll.zprimus), median(QsampleAll.zprimus)
;	PRINT, minmax(SFselectAll.zprimus), median(SFselectAll.zprimus)

	XYOUTS, xmin+0.95*xr, ymin+0.17*yr, 'median Q redshift: ' + twoDecimal(median(QsampleAll.zprimus)), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.13*yr, 'median SF redshift: ' + twoDecimal(median(SFselectAll.zprimus)), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.09*yr, 'median Q IP mass: ' + twoDecimal(median(QsampleAll.mstar)), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, 'median SF IP mass: ' + twoDecimal(median(SFselectAll.mstar)), ALIGNMENT=1.0

	QsampleAllFields =  [QsampleAllFields, QsampleAll]
	SFselectAllFields = [SFselectAllFields, SFselectAll]

	contQ = getContours(QsampleAll, 0, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	CONTOUR, contQ, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=2, /OVERPLOT

	contSF = getContours(SFselectAll, 1, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	CONTOUR, contSF, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=5, COLOR=cgColor('Black'), /OVERPLOT

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
  ENDFOR

	contQ = getContours(QsampleAllFields, 0, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	CONTOUR, contQ, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=2, /OVERPLOT

	contSF = getContours(SFselectAllFields, 1, 10, levels, z_binwidth, mass_binwidth, xpts, ypts)
	CONTOUR, contSF, xpts, ypts, LEVELS=levels, C_LINESTYLE=0, C_THICK=5, COLOR=cgColor('Black'), /OVERPLOT

	allIP = [QsampleAllFields, SFselectAllFields]

	MWRFITS, allIP, 'matchedIPsample_FBF.fits', /CREATE
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
