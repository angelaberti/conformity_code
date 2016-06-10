; dx, dy, dz in zerodSFQ_*.fits are in Mpc/h

PRO latefrac_VisualTest_coadd_v3, outputFormat
  IF string(outputFormat) eq 'ps' THEN BEGIN
	SET_PLOT, 'ps'
	DEVICE, file="latefrac_VisualTest_coadd.ps", /landscape ; $
;		+ strtrim(string(z,format='(f20.2)'),1) $
;		+ "_" + strtrim(string(dz,format='(f20.2)'),1) + ".ps", /landscape
  ENDIF ELSE BEGIN
	SET_PLOT, 'X'
  ENDELSE

;	!p.multi=[0,2,2]
	!p.multi=0
;	!x.margin=[15,2.5] ; left, right
;	!y.margin=[3.5,3] ; bottom, top

	field = 'all'
	plotEvery = 1

	SFQ = 0
	maxIP = 9

	dataAll = mrdfits('../results/default_parameters/IP_data/zerodSFQ_IP_dz1.0_dm0.5.fits', 1)

	FOR i=1,1 DO BEGIN
		plot_coadd, dataAll, 0.60+0.0*float(i), 0.40, field, SFQ, plotEvery, maxIP
	ENDFOR

  IF string(outputFormat) eq 'ps' THEN BEGIN
	DEVICE, /CLOSE
  ENDIF
	SET_PLOT, 'X'
END
PRO plot_coadd, dataAll, z, dz, field, SFQ, plotEvery, maxIP

;	PRINT, "Field: ", field
	
	data_zcut = dataAll[where( (dataAll.zprimus ge z-dz) AND $
				   (dataAll.zprimus le z+dz) )]
  	IF (field eq "all") THEN (data = data_zcut) ELSE (data = data_zcut[where(data_zcut.field eq field)])

	dataIP  = data[where((data.IP eq 1) AND (data.SFQ eq SFQ))]	

;	PRINT, "IP galaxies in redshift range: ", n_elements(dataIP)

	H0 	 = 100.*redh100()

	halfRange = 0.5 ; Mpc

	!p.charsize = 1.1
	xyoutscharsize = 1.0

	axisColor = cgColor('White')
	textColor = cgColor('White')
;	axisColor = cgColor('Black')
;	textColor = cgColor('Black')

;	IPcolorSF = cgColor('Cyan')
;	IPcolorQ  = cgColor('Red')
	IPcolorSF = axisColor
	IPcolorQ  = axisColor
	IPsym 	  = 4

	nonIPcolorSF = cgColor('Cyan')
	nonIPcolorQ  = cgColor('Red')
	nonIPsym     = 3

	dashedLineColor = axisColor

	plot_index	= 0
	n_rings		= 4
	ring_width	= 0.25 ; Mpc

;	PRINT, "ring width: ", ring_width

	ring_totals_early = 0.0*findgen(float(n_rings))
	ring_totals_late  = 0.0*findgen(float(n_rings))

	IF SFQ eq 1 THEN SFQstatus="Late" ELSE SFQstatus="Early"

	FOR IP_index=0,n_elements(dataIP)-1 DO BEGIN
;	FOR IP_index=0,maxIP-1 DO BEGIN
		; only plot certain fraction of IP neighbors
		IF (IP_index mod plotEvery eq 0) THEN BEGIN
			currentIP = dataIP[IP_index]

			rings = ring_width*findgen(float(n_rings)+1.)
;			PRINT, "rings:", rings

			distances = ((currentIP.xprop-data.xprop)^2 + (currentIP.yprop-data.yprop)^2)^0.5
			dataLocal = data[where( (data.field eq currentIP.field) AND $
						(ABS(data.zprimus - currentIP.zprimus) le 0.005*(1.+currentIP.zprimus)) AND $
						(distances gt 0.) AND (distances le halfRange) )]
		
			FOR s=0,n_elements(rings)-2 DO BEGIN
				data_ring 	= dataLocal[where( (distances ge rings[s]) AND (distances le rings[s+1]) )]
				
				subtot_early	= n_elements(data_ring[where(data_ring.SFQ eq 0)])
				subtot_late 	= n_elements(data_ring[where(data_ring.SFQ eq 1)])
				
				ring_totals_early[s] += subtot_early
				ring_totals_late[s]  += subtot_late
			ENDFOR			

		      ; SET UP AXES BASED ON FIRST IP
			IF (IP_index eq 0) THEN BEGIN
			      ; SHIFT POSITION OF FIRST IP TO ORIGIN
				xOffsetMaster = dataIP[0].xprop
				yOffsetMaster = dataIP[0].yprop

				xmin	= dataIP[0].xprop - xOffsetMaster - halfRange
				xmax	= dataIP[0].xprop - xOffsetMaster + halfRange
				xr	= xmax-xmin

	                        ymin2	= dataIP[0].yprop - yOffsetMaster - halfRange
	                        ymax2	= dataIP[0].yprop - yOffsetMaster + halfRange
	                        yr2	= ymax2-ymin2

			      ; DRAW AXES
	                        PLOT, [currentIP.xprop - xOffsetMaster], [currentIP.yprop - yOffsetMaster], psym=IPsym, color=axisColor, $
	                                xrange=[xmin, xmax], yrange=[ymin2, ymax2], $
					xtitle="Relative Mpc", ytitle="Relative Mpc", title = SFQstatus + "-type IPs"
			      ; PLOT FIRST IP
	                        OPLOT, [currentIP.xprop - xOffsetMaster], [currentIP.yprop - yOffsetMaster], psym=IPsym, color=IPcolor
	
			      ; OVERPLOT SF AND Q NON-IP GALAXIES OF FIRST IP
				OPLOT, [(dataLocal[where(dataLocal.SFQ eq 1)].xprop - xOffsetMaster)], [(dataLocal[where(dataLocal.SFQ eq 1)].yprop - yOffsetMaster)], psym=nonIPsym, color=cgColor('Red')
				OPLOT, [(dataLocal[where(dataLocal.SFQ eq 0)].xprop - xOffsetMaster)], [(dataLocal[where(dataLocal.SFQ eq 0)].yprop - yOffsetMaster)], psym=nonIPsym, color=cgColor('Cyan')

				plot_index += 1
			ENDIF

		      ; OVERPLOT SUBSEQUENT IP AND NON-IP GALAXIES RELATIVE TO FIRST
			xOffset = currentIP.xprop
			yOffset = currentIP.yprop

		      ; PLOT IP
;			OPLOT, [currentIP.xprop - xOffset], [currentIP.yprop - yOffset], psym=IPsym, color=IPcolor
		
		      ; PLOT SF AND Q NON-IP GALAXIES
                        OPLOT, [(dataLocal[where(dataLocal.SFQ eq 1)].xprop) - xOffset], [(dataLocal[where(dataLocal.SFQ eq 1)].yprop) - yOffset], psym=nonIPsym, color=nonIPcolorSF
                        OPLOT, [(dataLocal[where(dataLocal.SFQ eq 0)].xprop) - xOffset], [(dataLocal[where(dataLocal.SFQ eq 0)].yprop) - yOffset], psym=nonIPsym, color=nonIPcolorQ

			plot_index += 1
		ENDIF
	ENDFOR

	A = 2.*!PI*FINDGEN(33)/32.
	xcirc = 0.25*COS(A)
	ycirc = 0.25*SIN(A)
	OPLOT, xcirc, ycirc, linestyle=1, color=axisColor

	XYOUTS, xmin+0.02*xr, ymin2+0.95*yr2, strtrim(plot_index,1) + "/" + strtrim(n_elements(dataIP),1) + " IPs shown"
	XYOUTS, xmin+0.75*xr, ymin2+0.95*yr2, strtrim(string(z-dz,format='(f20.1)'),1) + " < z < " + strtrim(string(z+dz,format='(f20.1)'),1)

	FOR s=0,n_rings-1 DO BEGIN
		PRINT, "Ring: ", strtrim(s+1,1), "  Total early: ", ring_totals_early[s], "  Total late: ", ring_totals_late[s]
	ENDFOR

	RETURN
END
