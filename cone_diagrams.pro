PRO cone_diagrams, outputFormat
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/latex/figures/cone_diagrams', /ENCAP, THICK=5
	        DEVICE, /INCH, XS=10, YS=10, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

;	!x.margin=[8,2.5]; left, right
;	!y.margin=[4,2] ; bottom, top

	data = mrdfits('~/results/zerodSFQ_all_cart.fits',1)

;	fields = data[uniq(data.field, sort(data.field))].field
	fields = ['xmm       ','cfhtls_xmm','es1       ','cosmos    ','cdfs      ']
;	fieldnames = ['CDFS', 'XMM-CFHTLS', 'COSMOS', 'ES1', 'XMM-SXDS']
	fieldnames = ['XMM-SXDS', 'XMM-CFHTLS', 'ES1', 'COSMOS', 'CDFS']

	!P.FONT = 0
	!P.MULTI = [5,1,5]
	!P.CHARSIZE = 3.0
	
	ERASE

  FOR f=0,n_elements(fields)-1 DO BEGIN
	dataf = data[WHERE(data.field EQ fields[f])]
	PRINT, 'FIELD: ' + strcompress(fields[f])

        nx=1.
        ny=5.

        xm=0.12
        ym=0.12

        dx = (1 - 1.5*xm)/nx
        dy = (1 - 2.0*ym)/ny

	pos = [ xm + dx*float(f mod nx), ym + dy*float(floor(float(f)/nx)), $
		xm + dx*(1. + float(f mod nx)), ym + dy*(1. + float(floor(float(f)/nx))) ]

	xmin = dcomovinglos(0.175, /Mpc)
	xmax = dcomovinglos(1.025, /Mpc)
	xr = xmax-xmin

;	ymin = -1.1*MAX(ABS(dataf.yprop))
;	ymax = 1.1*MAX(ABS(dataf.yprop))
	ymin = -48
	ymax = 48
	yr = ymax-ymin

	A = FINDGEN(33)*(!PI*2/32.)
	W = 0.15
	USERSYM, W*COS(A), W*SIN(A), /FILL
	sym=8

	IF f EQ n_elements(fields)-1 THEN title='' ELSE title=''
	IF f EQ 2 THEN ytitle = 'Relative Distance (Mpc)' ELSE ytitle=''
	IF f EQ 0 THEN BEGIN
		xtitle = 'Physical Distance (Mpc)'
		xtickformat = ''
	ENDIF ELSE BEGIN
		xtitle=''
		xtickformat = '(A1)'
	ENDELSE

	datafQ  = dataf[where(dataf.SFQ EQ 0)]
	datafSF = dataf[where(dataf.SFQ EQ 1)]

	XTICKN=[0.2, 0.4, 0.6, 0.8, 1.0]

	PLOT, dcomovinglos(datafSF.zprimus, /Mpc), datafSF.yprop, position=pos, xrange=[xmin, xmax], yrange=[ymin, ymax], xtitle=xtitle, ytitle=ytitle, title=title, xtickformat=xtickformat, /NODATA
	IF f EQ 4 THEN AXIS, XAXIS=1, XTICKS=4, XTICKV=dcomovinglos(XTICKN, /Mpc), XTICKN=decimal(XTICKN,2), xtitle = 'Redshift'
	OPLOT, dcomovinglos(datafSF.zprimus, /Mpc), datafSF.yprop, psym=sym, color=cgcolor('Blue');, symsize=.1, thick=1
	OPLOT, dcomovinglos(datafQ.zprimus, /Mpc), datafQ.yprop, psym=sym, color=cgcolor('Red');, symsize=.1, thick=1
	XYOUTS, xmin+0.025*xr, ymin+0.80*yr, strcompress(strtrim(fieldnames[f],1)), ALIGNMENT=0.0, size=1.5
  ENDFOR
	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
