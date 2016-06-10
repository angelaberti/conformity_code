PRO ipplot, zerodInputFile, f, X, Y, IPtype, SFtype

	ymin = -40
	ymax =  40
	yr = ymax-ymin

        A = FINDGEN(17) * (!PI*2/16.)
        USERSYM, 0.3*COS(A), 0.3*SIN(A), /FILL

	zerodInput = mrdfits(zerodInputFile, 1)

	fields	= zerodInput[uniq(zerodInput.field, sort(zerodInput.field))].field
	n_fields = n_elements(fields)

	data	= zerodInput[where(zerodInput.field eq f)]

	dataCut = data[where((data.IP eq IPtype) AND (data.SFQ eq SFtype))]
;	dataIP	= data[where(data.IP eq 1)]
;	dataNIP	= data[where(data.IP eq 0)]

;	IF (IPtype eq 1) THEN data=dataIP ELSE $
;	IF (IPtype eq 0) THEN data=dataNIP

	tnames  = TAG_NAMES(zerodInput)
	tindex  = where(STRCMP(tnames, X, /FOLD_CASE) eq 1)
	data_x  = dataCut.(tindex)
	tindex  = where(STRCMP(tnames, Y, /FOLD_CASE) eq 1)
	data_y  = dataCut.(tindex)

	IF (X eq 'dz') THEN BEGIN
		xmin=500
		xmax=2400
		dist = " along line-of-sight"
	ENDIF ELSE BEGIN
		xmin=-40
		xmax=40
		dist = " (projected)"
	ENDELSE
	xr = xmax-xmin

	IF IPtype eq 1 THEN IPstatus = "IP" ELSE IPstatus = "non-IP"
	IF SFtype eq 1 THEN SFQstatus = "Star-forming (late-type)" ELSE SFQstatus = "Quiescent (early-type)"

	plot, data_x, data_y, psym=8, xrange=[xmin,xmax], yrange=[ymin,ymax], $
		xtitle=textoidl('Mpc h^{-1}'), ytitle=textoidl('Mpc h^{-1}')
	xyouts, xmin+0.05*xr, ymin+0.90*yr, "Field: " + f
	xyouts, xmin+0.05*xr, ymin+0.85*yr, SFQstatus + " " + IPstatus + dist
;	xyouts, xmin+0.05*xr, ymin+0.80*yr, SFQstatus
END
