; dx, dy, dz in zerodSFQ_*.fits are in Mpc/h

PRO latefrac_VisualTest_v3, outputFormat;, z, dz, plotEvery
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/visual_test', THICK=2, /ENCAP
	        DEVICE, /INCH, XS=6, YS=9
		axisColor = cgColor('Black')
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		axisColor = cgColor('White')
	ENDELSE

	dz_coeff = 2.0
	plotEvery = 400

;	!p.multi=[24,4,6]
	!p.multi=24
;	!x.margin=[3,4] ; left, right
;	!y.margin=[3.5,1] ; bottom, top

	data = mrdfits('~/results/variable_mass+1.0dex/IP_data/zerodSFQ_IP_dz2.0.fits', 1)
	dataIP  = data[where(data.IP eq 1)]	

;	IP_index = 1

	halfRange = 4. ; Mpc

	!p.charsize = 2
	xyoutscharsize = 0.8

	textColor = axisColor

	IPcolorSF = cgColor('Blue')
	IPcolorQ  = cgColor('Red')
	IPsym 	  = 4

	nonIPcolorSF = cgColor('Blue')
	nonIPcolorQ = cgColor('Red')
	nonIPsym   = 1

	dashedLineColor = cgColor('Green')

        A = FINDGEN(33)*(!PI*2/32.)

	NIPplot = 0
	
	ERASE

	total_neighbors = 0

	nx=4.
        ny=6.

        xm=0.08
        ym=0.08

        dx = (1 - 1.5*xm)/nx
        dy = (1 - 2.0*ym)/ny

	PRINT, 'n_elements(dataIP): ', n_elements(dataIP)

;	FOR IP_index=0,n_elements(dataIP)-1 DO BEGIN
;	FOR IP_index=0,19 DO BEGIN
;	  IF (IP_index mod plotEvery eq 0) THEN BEGIN

	FOR i=0,nx*ny-1 DO BEGIN

		IP_index = plotEvery*(i+1)
	
                IF (i le (nx-1)) THEN BEGIN
                        xtitle = 'Mpc'
                        xtickformat = ''
                ENDIF ELSE BEGIN
                        xtitle = ''
                        xtickformat = '(A1)'
                ENDELSE
                IF (i mod nx eq 0) THEN BEGIN
                        ytitle = 'Mpc'
                        ytickformat = ''
                ENDIF ELSE BEGIN
                        ytitle = ''
                        ytickformat = '(A1)'
                ENDELSE

                pos = [ xm + dx*float(i mod nx), ym + dy*float(floor(float(i)/nx)), $
                        xm + dx*(1. + float(i mod nx)), ym + dy*(1. + float(floor(float(i)/nx))) ]
		PRINT, pos

		currentIP = dataIP[IP_index]

		dataLocal = data[where( (data.field eq currentIP.field) AND ( ABS(currentIP.zprimus - data.zprimus) le dz_coeff*0.005*(1.+currentIP.zprimus) ) AND $
			( SQRT( (currentIP.xprop-data.xprop)^2 + (currentIP.yprop-data.yprop)^2 ) le 4. ) AND $
			(currentIP.objname NE data.objname) )]

		zdiff = ABS(dataLocal.zprimus - currentIP.zprimus)
		PRINT, 'Delta z: ' + strtrim(dz_coeff*0.005*(1. + currentIP.zprimus),1) + '    Max z diff of neighbors: ' + strtrim(MAX(zdiff),1)

		dataHyperLocal = dataLocal[where(SQRT( (currentIP.xprop-dataLocal.xprop)^2 + (currentIP.yprop-dataLocal.yprop)^2 ) le 0.5, /NULL) ] 
		IF dataHyperLocal NE !NULL THEN BEGIN
			PRINT, n_elements(dataHyperLocal)
			PRINT, 'IP mass: ' + strtrim(currentIP.mstar)
			FOREACH element,dataHyperLocal DO PRINT, element.mstar
			PRINT, ''
		ENDIF

		Nneighbors = strcompress(n_elements(dataLocal))
		total_neighbors += n_elements(dataLocal)

		NIPplot += 1

		xoff = currentIP.xprop
		yoff = currentIP.yprop

		PLOT, [currentIP.xprop]-xoff, [currentIP.yprop]-yoff, psym=IPsym, color=axisColor, $
			xrange=[currentIP.xprop-xoff-halfRange, currentIP.xprop-xoff+halfRange], $
                        yrange=[currentIP.yprop-yoff-halfRange, currentIP.yprop-yoff+halfRange], $
			xtitle=xtitle, ytitle=ytitle, xtickformat=xtickformat, ytickformat=ytickformat, position=pos, /NODATA
		OPLOT, [currentIP.xprop], [currentIP.yprop], psym=IPsym, color=IPcolor

		FOR j=1,4 DO BEGIN
			OPLOT, currentIP.xprop-xoff+float(j)*COS(A), currentIP.yprop-yoff+float(j)*SIN(A), linestyle=4, color=dashedLineColor
		ENDFOR

		OPLOT, [(dataLocal[where(dataLocal.SFQ eq 1)].xprop)]-xoff, [(dataLocal[where(dataLocal.SFQ eq 1)].yprop)]-yoff, psym=nonIPsym, color=nonIPcolorSF
		OPLOT, [(dataLocal[where(dataLocal.SFQ eq 0)].xprop)]-xoff, [(dataLocal[where(dataLocal.SFQ eq 0)].yprop)]-yoff, psym=nonIPsym, color=nonIPcolorQ

;               XYOUTS, currentIP.xprop-xoff+halfRange, currentIP.yprop-yoff-halfRange, textoidl("IP M_{\odot} ") + strtrim(string(currentIP.mstar, format='(f20.2)'),1), color=IPcolor, ALIGNMENT=1.0, charsize=xyoutscharsize
                XYOUTS, currentIP.xprop-xoff-0.9*halfRange, currentIP.yprop-yoff+0.85*halfRange, strtrim(currentIP.field,1), color=IPcolor, charsize=xyoutscharsize, ALIGNMENT=0.0
;		XYOUTS, currentIP.xprop-xoff-0.9*halfRange, currentIP.yprop-yoff-0.90*halfRange, textoidl('N_{neigh}: ') + strtrim(Nneighbors,1), color=IPcolor, charsize=xyoutscharsize, ALIGNMENT=0
                XYOUTS, currentIP.xprop-xoff-0.9*halfRange, currentIP.yprop-yoff-0.90*halfRange, strtrim(Nneighbors,1), color=IPcolor, charsize=xyoutscharsize, ALIGNMENT=0
;	  ENDIF
	ENDFOR
	
	PRINT, 'Total neighbors: ', total_neighbors
;	PRINT, 'Total IP: ', n_elements(dataIP)
	PRINT, 'Total IP: ', NIPplot
	PRINT, 'Average neighbors per IP: ', float(total_neighbors)/NIPplot

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
