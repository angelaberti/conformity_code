PRO visTest, IPtype
	dataAll	= mrdfits('../results/default_parameters/IP_data/zerodSFQ_IP_dz1.0_dm0.5.fits',1)

	ERASE

;	!p.multi=[0,2,2]
;	!p.position=aspect(1.0)
	!p.multi=0
;	!x.margin=[15,2.5] ; left, right
;	!y.margin=[3.5,3] ; bottom, top

	showNum = 1000
;	IPtype = 1

	FOR i=1,4 DO BEGIN
		plot_coadd, dataAll, 0.2*float(i), 0.2*float(i+1), IPtype, showNum, i
	ENDFOR
END

PRO plot_coadd, dataAll, zmin, zmax, IPtype, showNum, plot_index

	Rmax = 0.5 ; Mpc
	dz_coeff = 1.0

	data	= dataAll[where( (dataAll.zprimus ge zmin) AND (dataAll.zprimus le zmax) )]
	dataIP	= data[where( (data.IP eq 1) AND (data.SFQ eq IPtype) )]

	dw = 0.475
	dh = 0.475

	pos = [ 0.1 + dw*float(1-(plot_index mod 2)), $
		0.1 + dh*float(1-((plot_index+1)/2 mod 2)), $
		0.475 + dw*float(1-(plot_index mod 2)), $
		0.475 + dh*float(1-((plot_index+1)/2 mod 2))]

	PLOT, [0], [0], xrange=[-1.*Rmax,Rmax], yrange=[-1.*Rmax,Rmax], position=pos, xtitle='Mpc', ytitle='Mpc', /NODATA
	XYOUTS, -0.9*Rmax, 0.85*Rmax, strtrim(string(zmin,format='(f20.1)'),1) + " < z < " + strtrim(string(zmax,format='(f20.1)'),1)

	PRINT, "IP type: ", strcompress(IPtype)
	IF IPtype eq 1 THEN IPlabel='SF IPs: ' ELSE IPlabel='Q IPs'
	XYOUTS, 0.9*Rmax, 0.85*Rmax, IPlabel + strcompress(n_elements(dataIP)), ALIGNMENT=1.0

	FOR i=0,n_elements(dataIP)-1 DO BEGIN
;	FOR i=0,showNum-1 DO BEGIN
		currentIP = dataIP[i]
		xOffset = currentIP.xprop
		yOffset = currentIP.yprop

		neighborCriteria = where(( data.field eq currentIP.field ) AND $
					 ( ABS(data.zprimus - currentIP.zprimus) le dz_coeff*0.005*(1. + currentIP.zprimus) ) AND $
					 ( SQRT((currentIP.xprop-data.xprop)^2 + (currentIP.yprop-data.yprop)^2) le Rmax ), /NULL)

		neighbors = data[neighborCriteria]
		SFneighbors = neighbors[where(neighbors.SFQ eq 1, /NULL)]
		Qneighbors  = neighbors[where(neighbors.SFQ eq 0, /NULL)]
;		PRINT, i, " total neighbors:", n_elements(neighbors), " SF:", n_elements(SFneighbors), " Q:", strcompress(n_elements(Qneighbors))

		IF n_elements(SFneighbors) ne 0 THEN BEGIN
			OPLOT, [SFneighbors.xprop-xOffset], [SFneighbors.yprop-yOffset], psym=3, color=cgColor('Cyan')
		ENDIF
		IF n_elements(Qneighbors) ne 0 THEN BEGIN
			OPLOT, [Qneighbors.xprop-xOffset], [Qneighbors.yprop-yOffset], psym=3, color=cgColor('Orange')
		ENDIF
	ENDFOR

	A = FINDGEN(33)*2*!PI/33.
	xcirc = 0.25*COS(A)
	ycirc = 0.25*SIN(A)
	OPLOT, xcirc, ycirc, linestyle=2, thick=3
END
