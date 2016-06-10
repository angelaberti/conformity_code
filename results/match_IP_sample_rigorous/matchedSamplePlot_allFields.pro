PRO matchedSamplePlot_allFields, outputFormat
	!P.FONT=0

	targPhi = -3.7
	stringPHI = strtrim(string(-1.*targPhi, format='(f20.1)'), 2)

	seed = 1
	dataAll = mrdfits('~/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits',1)
        data = dataAll[where(dataAll.IP eq 1)]

        fields = dataAll[uniq(dataAll.field, sort(dataAll.field))].field

	A = FINDGEN(17) * (!PI*2/16.)
	W = 0.5
	USERSYM, W*COS(A), W*SIN(A), /FILL
	PSYM=3
	selectSym=3

        mass    = [10.8,        10.9,   11.0,   11.1,   11.2,   11.3,   11.4]
        mass_z2 = [10.8,        10.9,   11.0,   11.1,           11.3]
        phi_z2  = [-2.894,	-3.180, -3.380, -3.370,         -4.600]
        phi_z3  = [-2.865,      -2.905, -3.113, -3.430, -3.670, -4.040, -4.150]
        phi_z4  = [-2.762,      -2.978, -3.117, -3.347, -3.540, -3.830, -4.230]
        phi_z5  = [-2.185,	-2.980, -3.139, -3.368, -3.539, -3.930, -4.080]
        phi_z65 = [-2.744,	-2.788, -3.000, -3.148, -3.286, -3.671, -4.040]
        phi_z8  = [-2.908,	-3.011, -3.113, -3.279, -3.405, -3.630, -3.920]

        phiArray = [[phi_z3], [phi_z4], [phi_z5], [phi_z65], [phi_z8]]

	mArray  = [INTERPOL(mass_z2, phi_z2, targPhi)]
	interp_zArray = [0.25, 0.35, 0.45, 0.55, 0.725, 0.9]

        FOR i=0,n_elements(interp_zArray)-2 DO mArray = [mArray, INTERPOL(mass, phiArray[*,i], targPhi)]

;	PRINT, 'mArray: ', mArray

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/latex/figures/matchedSamplePlot_allFields', /ENCAP, THICK=5
	        DEVICE, /INCH, XS=6, YS=5, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	!P.MULTI = 0
	!P.CHARSIZE = 1.5

;	!x.margin=[6,2.5] ; left, right
;	!y.margin=[3.5,2.5] ; bottom, top

	xmin = 0.2
	xmax = 1.0
	xr = xmax-xmin

	yminPlot = 8.7
	ymaxPlot = 11.7
	yrPlot = ymaxPlot-yminPlot

	IF outputFormat eq 'ps' THEN SFselectColor = cgColor('Black') ELSE SFselectColor = cgColor('White')
	Qcolor  = cgColor('Red')
	SFcolor = cgColor('Blue')
	textColor = cgColor('Black')
;	textColor = cgColor('White')

	mass_binwidth	= 0.2
        z_binwidth	= 0.05
	n_levels	= 10

	xtitle = 'Redshift'	
;	ytitle = textoidl('log (M_{stellar} / M_{\odot})')
	ytitle = textoidl('log ( M_{*}/M') + sunsymbol() + ' )'

	zmin = 0.2
	zmax = 1.0
	zArray = zmin + z_binwidth*findgen((zmax-zmin)/z_binwidth + 1)

	; SET UP PLOT AREA
	PLOT, data.zprimus, data.mstar, xrange=[xmin, xmax], yrange=[yminPlot, ymaxPlot], xtitle=xtitle, ytitle=textoidl('log (M_{*}/M_{\odot} )'), COLOR=textColor, /NODATA
;	title='Matched SF and Q IP Samples', COLOR=textColor, /NODATA

;	LEGEND, ['Late-type IPs','Early-type IPs'], LINESTYLE=[0,2], COLOR=[cgColor('blue'),cgColor('red')], BOX=0, /BOTTOM, /RIGHT
	LEGEND, ['SF IP','Q IP'], LINESTYLE=[0,2], COLOR=[cgColor('blue'),cgColor('red')], BOX=0, /BOTTOM, /RIGHT, THICK=5, NUMBER=2, PSPACING=3, CHARSIZE=1.25, textcolors=cgColor('Black')
;	LEGEND, ['SF IP','Q IP'], LINESTYLE=[0,2], COLOR=[cgColor('blue'),cgColor('red')], BOX=0, POS=[0.52,8.65], THICK=5, NUMBER=2, PSPACING=3

	IPselect = MRDFITS('~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)
	massArray_IPselect = FLOOR(10.*MIN(IPselect.mstar))/10. + mass_binwidth*FINDGEN((CEIL(10.*MAX(IPselect.mstar))/10. - FLOOR(10.*MIN(IPselect.mstar))/10.)/mass_binwidth + 1)	
;	PRINT, 'IP selection mass array: ', massArray_IPselect

	SFselect = IPselect[where(IPselect.SFQ EQ 1)]
	Qselect  = IPselect[where(IPselect.SFQ EQ 0)]

	dataSF = data[where(data.SFQ EQ 1)]
	dataQ  = data[where(data.SFQ EQ 0)]

	; SCATTER PLOT OF ALL PRIMUS SOURCES WITH TARG_WEIGHT >= 1
	allPRIMUS = MRDFITS('~/results/zerodSFQ_all_cart.fits', 1)
	allPRIMUS = allPRIMUS[WHERE(allPRIMUS.TARG_WEIGHT GE 1)]
	allPRIMUS_SF = allPRIMUS[WHERE(allPRIMUS.SFQ EQ 1)]
	allPRIMUS_Q  = allPRIMUS[WHERE(allPRIMUS.SFQ EQ 0)]
;	OPLOT, allPRIMUS_SF.zprimus, allPRIMUS_SF.mstar, psym=psym, COLOR=cgColor('cyan')
;	OPLOT, allPRIMUS_Q.zprimus, allPRIMUS_Q.mstar, psym=psym, COLOR=cgColor('orange')

	; SCATTER PLOT OF ALL IP IN ALL FIELDS
;	OPLOT, dataSF.zprimus, dataSF.mstar, psym=psym, COLOR=SFcolor
;	OPLOT, dataQ.zprimus, dataQ.mstar, psym=psym, COLOR=Qcolor

;	levels = 60.*INDGEN(n_levels)
;	levels = [5,levels[1:*]]

;	FOREACH elem,massArray_IPselect DO OPLOT, [zmin,zmax], [elem,elem], LINESTYLE=0
;	FOREACH elem,zArray DO OPLOT, [elem,elem], [yminPlot,ymaxPlot], LINESTYLE=0

;	cont = getContours(IPselect, z_binwidth, mass_binwidth, xpts, ypts)
	z_binwidth = 0.05
	mass_binwidth = 0.2
	cont = getContours(IPselect, z_binwidth, mass_binwidth, xpts, ypts)
;;	CONTOUR, MIN_CURVE_SURF(cont, XGRID=[0.2, z_binwidth], YGRID=[MIN(IPselect.mstar), mass_binwidth], GS=[z_binwidth,mass_binwidth]), $
;;		xpts+0.25*z_binwidth, ypts+0.0*mass_binwidth, $
;;		/FILL, LEVELS=levels, /OVERPLOT
	levels = LEVELS(cont, 6, 0.05)
	PRINT, 'All levels: ', levels
	CONTOUR, MIN_CURVE_SURF(cont, XGRID=[0.2, z_binwidth], YGRID=[MIN(IPselect.mstar), mass_binwidth], GS=[z_binwidth,mass_binwidth]), $
		xpts+0.25*z_binwidth, ypts+0.25*mass_binwidth, LEVELS=levels, /OVERPLOT, /FILL

;	contSF = getContours(allPRIMUS_SF, z_binwidth, mass_binwidth, xpts, ypts)
	contSF = getContours(dataSF, z_binwidth, mass_binwidth, xpts, ypts)
;	levels=[20, 0.8*MAX(contSF)*(1+INDGEN(5))/6.]

	levels = LEVELS(contSF, 6, 0.1)
	PRINT, 'SF levels: ', levels
;	CONTOUR, MIN_CURVE_SURF(contSF, XGRID=[0.2, z_binwidth], YGRID=[MIN(allPRIMUS_SF.mstar), mass_binwidth], GS=[z_binwidth,mass_binwidth]), xpts+0.25*z_binwidth, ypts+0.25*mass_binwidth, C_THICK=5, COLOR=SFcolor, /OVERPLOT, LEVELS=levels
	CONTOUR, MIN_CURVE_SURF(contSF, XGRID=[0.2, z_binwidth], YGRID=[MIN(dataSF.mstar)-mass_binwidth, mass_binwidth], GS=[z_binwidth,mass_binwidth]), xpts+0.25*z_binwidth, ypts+mass_binwidth, C_THICK=5, COLOR=SFcolor, /OVERPLOT, LEVELS=levels

;	FOR i=0,n_elements(xpts)-1 DO BEGIN
;	  FOR j=0,n_elements(ypts)-1 DO BEGIN
;	    IF (contSF[i,j] GT 0) THEN XYOUTS, xpts[i]+0.5*z_binwidth, ypts[j]+0.33*mass_binwidth, STRTRIM(contSF[i,j],2), ALIGNMENT=0.5, CHARSIZE=1.5
;	  ENDFOR
;	ENDFOR

;	contQ = getContours(allPRIMUS_Q, z_binwidth, mass_binwidth, xpts, ypts)
	contQ = getContours(dataQ, z_binwidth, mass_binwidth, xpts, ypts)
;	levels=[20, 0.8*MAX(contQ)*(1+INDGEN(5))/6.]
	levels = LEVELS(contQ, 6, 0.05)
	PRINT, 'Q levels:  ', levels
;	CONTOUR, MIN_CURVE_SURF(contQ, XGRID=[0.2, z_binwidth], YGRID=[MIN(allPRIMUS_Q.mstar), mass_binwidth], GS=[z_binwidth,mass_binwidth]), xpts+0.25*z_binwidth, ypts+0.25*mass_binwidth, C_THICK=5, COLOR=Qcolor, C_LINESTYLE=2, /OVERPLOT, LEVELS=levels
	CONTOUR, MIN_CURVE_SURF(contQ, XGRID=[0.2, z_binwidth], YGRID=[MIN(dataQ.mstar), mass_binwidth], GS=[z_binwidth,mass_binwidth]), xpts+0.25*z_binwidth, ypts+0.0*mass_binwidth, C_THICK=5, COLOR=Qcolor, C_LINESTYLE=2, /OVERPLOT, LEVELS=levels

	; HIGHLIGHT SELECTED SOURCES
;	OPLOT, IPselect.zprimus, IPselect.mstar, psym=8, COLOR=cgcolor('blue')

	; DETERMINE HIGH MASS CUTOFF FOR SF IP SOURCES
        massFunc = INTERPOL(mArray, interp_zArray, data.zprimus)

	; CUT SF AND Q IP SAMPLES BASED ON (SF) HIGH MASS CUTOFF (ABOVE)
	dataQcut  = dataQ[where(dataQ.mstar LE INTERPOL(mArray, interp_zArray, dataQ.zprimus))]
	dataSFcut = dataSF[where(dataSF.mstar LE INTERPOL(mArray, interp_zArray, dataSF.zprimus))]

	; SCATTER PLOT OF CUT SOURCES (BELOW HIGH MASS CUTOFF)
;	OPLOT, dataQcut.zprimus, dataQcut.mstar, psym=8, COLOR=Qcolor
;	OPLOT, dataSFcut.zprimus, dataSFcut.mstar, psym=8, COLOR=SFcolor

;	XYOUTS, xmin+0.05*xr, yminPlot+0.95*yrPlot, STRTRIM(fields[f],2), ALIGNMENT=0.0

;	XYOUTS, xmin+0.95*xr, yminPlot+0.13*yrPlot, 'Unique IPs: ' + strtrim(n_elements(UNIQ(SFselect.objname)),2) + ' (SF) ' + strtrim(n_elements(Qselect),2) + ' (Q)', ALIGNMENT=1.0, COLOR=textColor
;	XYOUTS, xmin+0.95*xr, yminPlot+0.09*yrPlot, 'Median redshift: ' + decimal(median(SFselect.zprimus),2) + ' (SF) ' + decimal(median(Qselect.zprimus),2) + ' (Q)', ALIGNMENT=1.0, COLOR=textColor
;	XYOUTS, xmin+0.95*xr, yminPlot+0.05*yrPlot, 'Median mass: ' + decimal(median(SFselect.mstar),2) + ' (SF) ' + decimal(median(Qselect.mstar),2) + ' (Q)', ALIGNMENT=1.0, COLOR=textColor

;	OPLOT, data[SORT(data.zprimus)].zprimus, massFunc[SORT(massFunc)], LINESTYLE=2, THICK=5, color=textColor

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

FUNCTION twoDecimal, input
	RETURN, strtrim(string(input, format='(f20.2)'),1)
END

FUNCTION levels, grid_data, n_levels, minPercent

	minlev = CEIL(minPercent*MAX(grid_data))

	width = CEIL((MAX(grid_data)-minlev)/FLOAT(n_levels))
	levels = [minlev+MIN(grid_data)+width*FINDGEN(n_levels)]
;	PRINT, levels

	RETURN, levels
END
