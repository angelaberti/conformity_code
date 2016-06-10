PRO IPhist_latefrac_vs_z, outputFormat
	!P.FONT=0

        IF (string(outputFormat) EQ 'ps') THEN BEGIN
		PS_OPEN, '~/latex/figures/IPhist_latefrac_vs_z', THICK=5, /ENCAP
		DEVICE, /INCH, XS=8, YS=8, SET_FONT='Palatino-Roman'
        ENDIF ELSE BEGIN
                SET_PLOT, 'X'
        ENDELSE

	bin = 0.05

	ERASE
	!P.MULTI=2
	!P.CHARSIZE=1.75
	posArray = grid_array_deluxe_unequal(1,2,0.6,0.05,0)

	zmin = 0.2
	zmax = 1.0

	ymin = 0.0
	ymax = 1400
	yr = ymax-ymin	

	zbins = zmin+bin*FINDGEN(CEIL((zmax-zmin)/bin))

	dataAll = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)
	dataIP  = MRDFITS('~/conformity/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits',1)

	dataAll   = dataAll[WHERE(dataAll.targ_weight GE 1.)]
	dataAllSF = dataAll[WHERE(dataAll.SFQ EQ 1)]
	dataAllQ  = dataAll[WHERE(dataAll.SFQ EQ 0)]

	dataIP 	 = dataIP[WHERE((dataIP.targ_weight GE 1.) AND (dataIP.IP EQ 1))]
	dataIPSF = dataIP[WHERE(dataIP.SFQ EQ 1)]
	dataIPQ  = dataIP[WHERE(dataIP.SFQ EQ 0)]

	PLOT, INDGEN(10), xrange=[zmin-bin/2,zmax+bin/2], yrange=[ymin,ymax], PSYM=10, xtickformat='(A1)', ytitle=textoidl('N_{gal}'), POSITION=posArray[*,1], /NODATA
;	LEGEND, ['Late-type IP','Early-type IP'], COLOR=[cgColor('blue'),cgColor('red')], BOX=0, /TOP, /RIGHT, LINESTYLE=[0,0], NUMBER=2, PSPACING=3
	plot_zhist, dataIPSF, 0, cgColor('blue'), bin, zbins, 45
;	plot_zhist, dataIPQ, 0, cgColor('white'), bin, zbins, 45
	plot_zhist, dataIPQ, 3, cgColor('red'), bin, zbins, 135
;	XYOUTS, 0.6, 400, 'Early-type IP', ALIGNMENT=0.5, FILL_BACKGROUND=1, FILL_COLOR=cgColor('white')
	LEGEND, ['SF IP', 'Q IP'], LINESTYLE=[0,3], COLOR=cGcolor(['blue','red']), /TOP, /RIGHT, BOX=0, NUMBER=2, PSPACING=3
;	LEGEND, ['Late-type IP'], /MIDDLE, /CENTER
	
	bluePoints=[[0.85,ymin+0.90*yr],$
		[0.85,ymin+0.95*yr],$
		[1.0,ymin+0.95*yr],$
		[1.0,ymin+0.90*yr]]
	redPoints=[[0.85,ymin+0.8*yr],$
		[0.85,ymin+0.85*yr],$
		[1.0,ymin+0.85*yr],$
		[1.0,ymin+0.80*yr]]

	BPT=TRANSPOSE(bluePoints)
	RPT=TRANSPOSE(redPoints)
;	POLYFILL, bluePoints, COLOR=cgColor('blue'), ORIENTATION=45
;	OPLOT, [BPT[*,0],BPT[0,0]], [BPT[*,1],BPT[0,1]], COLOR=cgColor('blue'), PSYM=0
;	POLYFILL, redPoints, COLOR=cgColor('red'), ORIENTATION=135
;	OPLOT, [RPT[*,0],RPT[0,0]], [RPT[*,1],RPT[0,1]], COLOR=cgColor('red'), PSYM=0

;	XYOUTS, 0.82, ymin+0.9*yr, 'SF IP', ALIGNMENT=1.0
;	XYOUTS, 0.82, ymin+0.8*yr, 'Q IP', ALIGNMENT=1.0

	SFarray = HISTOGRAM(dataIPSF.zprimus, bin=bin)
	Qarray  = HISTOGRAM(dataIPQ.zprimus, bin=bin)
;	plot_zhist, dataIPSF, 0, cgColor('blue'), bin, zbins, 45
	frac = FLOAT(SFarray)/(SFarray + Qarray)

;	PRINT, TRANSPOSE([[SFarray],[Qarray],[frac]])

	PLOT, [zmin,zbins+bin/2,zmax], [frac[0],frac,frac[-1]], xrange=[zmin-bin/2,zmax+bin/2], yrange=[0.55,0.85], LINESTYLE=0, xtitle='Redshift', ytitle=textoidl('SF fraction'), POSITION=posArray[*,0], /NODATA, yticks=3, yminor=2
;	PLOT, [zmin,zbins+bin/2,zmax], [frac[0],frac,frac[-1]], xrange=[zmin-bin/2,zmax+bin/2], yrange=[0.5,1], LINESTYLE=0, xtitle='Redshift', ytitle='Late-type Fraction', POSITION=posArray[*,0], COLOR=cgColor('gray')
	POLYFILL, [zmin,zmin,zbins+bin/2,zmax,zmax], [0.56,frac[0],frac,frac[-1],0.56], COLOR=cgColor('gray')
	
        IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
       	SET_PLOT, 'X'
END

PRO plot_zhist, data, LS, COLOR, bin, zbins, angle
	zmin = 0.2
	zmax = 1.0

;	histSF = HISTOGRAM(dataIPSF.zprimus, bin=bin)
;	histQ  = HISTOGRAM(dataIPQ.zprimus, bin=bin)

	histData = HISTOGRAM(data.zprimus, bin=bin)
;	OPLOT, zbins+bin/2, histData, COLOR=COLOR, PSYM=10

;	PRINT, TRANSPOSE([[zbins],[histData]])

	points = []
	FOR i=0,n_elements(zbins)-2 DO BEGIN
		points = [[points],$
			[zbins[i],histData[i]],$
			[zbins[i+1],histData[i]]]
	ENDFOR
	
	points = [[zbins[0],10],$
		[points],$
		[zbins[-1],histData[-1]],$
		[zbins[-1]+bin,histData[-1]],$
		[zbins[-1]+bin,10]]

	OPLOT, points[0,*], points[1,*], LINESTYLE=LS, COLOR=COLOR
;	POLYFILL, points[0,*], points[1,*], COLOR=COLOR, /LINE_FILL, ORIENTATION=angle
END
