; dx, dy, dz in zerodSFQ_*.fits are in Mpc/h

PRO latefrac_VisualTest_coadd_v2
;	CD, '/IP_data'

;	SET_PLOT, 'ps'
;	DEVICE, file="latefrac_VisualTest_coadd.ps", /landscape ; $
;		+ strtrim(string(z,format='(f20.2)'),1) $
;		+ "_" + strtrim(string(dz,format='(f20.2)'),1) + ".ps", /landscape

	!p.multi=[0,2,2]
;	!p.multi=0
	!x.margin=[15,2.5] ; left, right
	!y.margin=[3.5,3] ; bottom, top

	vel = 1500.
	field = 'all'
	plotEvery = 1

	dataAll = mrdfits('../results/default_parameters/IP_data/zerodSFQ_IP_dv1500kms_dm0.5.fits', 1)
;	dataAll = mrdfits('IP_data/zerodSFQ_IP_dv5000kms_dm0.5.fits', 1)

	plot_coadd, dataAll, 0.30, 0.10, vel, field, 1, plotEvery
	plot_coadd, dataAll, 0.70, 0.10, vel, field, 1, plotEvery
	plot_coadd, dataAll, 0.30, 0.10, vel, field, 0, plotEvery
	plot_coadd, dataAll, 0.70, 0.10, vel, field, 0, plotEvery

;	DEVICE, /CLOSE
;	SET_PLOT, 'X'
END

PRO plot_coadd, dataAll, z, dz, delta_z, field, SFQ, plotEvery

;	PRINT, "Field: ", field
	
	data_zcut = dataAll[where( (dataAll.zprimus ge z-dz) AND $
				   (dataAll.zprimus le z+dz) )]
  	IF (field eq "all") THEN (data = data_zcut) ELSE (data = data_zcut[where(data_zcut.field eq field)])

	dataIP  = data[where((data.IP eq 1) AND (data.SFQ eq SFQ))]	

;	PRINT, "IP galaxies in redshift range: ", n_elements(dataIP)

;	IP_index = 1
;	delta_z  = 1500. ; km/s
	delta_z = float(delta_z)
	H0 	 = 100.*redh100()

	halfRange = 4. ; Mpc

	!p.charsize = 1.1
	xyoutscharsize = 1.0

	axisColor = cgColor('White')
	textColor = cgColor('White')
;	axisColor = cgColor('Black')
;	textColor = cgColor('Black')

;	IPcolorSF = cgColor('Blue')
;	IPcolorQ  = cgColor('Red')
	IPcolorSF = axisColor
	IPcolorQ  = axisColor
	IPsym 	  = 4

	nonIPcolorSF = cgColor('Blue')
	nonIPcolorQ = cgColor('Red')
	nonIPsym   = 3

	dashedLineColor = axisColor

        A = FINDGEN(33)*(!PI*2/32.)

	plot_index = 0
	n_slices = 10
	vslice_width = 2.*delta_z/float(n_slices)

;	PRINT, "velocity slice width: ", vslice_width

	slice_totals_early = 0.0*findgen(float(n_slices))
	slice_totals_late  = 0.0*findgen(float(n_slices))

;	FOR IP_index=0,n_elements(dataIP)-1 DO BEGIN
	FOR IP_index=0,99 DO BEGIN
		IF (IP_index mod plotEvery eq 0) THEN BEGIN
			currentIP = dataIP[IP_index]

			vslices = (H0*currentIP.dz/redh100() - delta_z) + vslice_width*findgen(float(n_slices)+1.)
			PRINT, vslices			

			dataLocal = data[where( (data.field eq currentIP.field) AND $
						(ABS(H0*data.dz/redh100() - H0*currentIP.dz/redh100()) le delta_z) AND $
						( (((currentIP.dx-data.dx)/redh100())^2 + ((currentIP.dy-data.dy)/redh100())^2)^0.5 le halfRange ) AND $
						(currentIP.objname ne data.objname) )]
			FOR s=0,n_elements(vslices)-2 DO BEGIN
				data_vslice = dataLocal[where( (H0*dataLocal.dz/redh100() ge vslices[s]) AND (H0*dataLocal.dz/redh100() le vslices[s+1]) )]
				
				subtot_early = n_elements(data_vslice[where(data_vslice.SFQ eq 0)])
				subtot_late  = n_elements(data_vslice[where(data_vslice.SFQ eq 1)])
				
				slice_totals_early[s] += subtot_early
				slice_totals_late[s]  += subtot_late

;				PRINT, s, " Total galaxies in slice: ",  strtrim(n_elements(data_vslice),1), vslices[s], vslices[s+1], subtot_early, slice_totals_early[s], subtot_late, slice_totals_late[s]
			ENDFOR			

			IF currentIP.SFQ eq 1 THEN IPcolor = IPcolorSF ELSE IPcolor = IPcolorQ

                        ymin2 = H0*dataIP[0].dz/redh100()-delta_z
                        ymax2 = H0*dataIP[0].dz/redh100()+delta_z
                        yr2 = ymax2-ymin2

;                        xyouts, currentIP.dx/redh100(), currentIP.dy/redh100()-0.9*halfRange, textoidl("     IP M_{\odot} ") + strtrim(string(currentIP.mstar, format='(f20.2)'),1), color=IPcolor
;                        xyouts, currentIP.dx/redh100()-0.9*halfRange, currentIP.dy/redh100()+0.8*halfRange, "z=" + strtrim(string(currentIP.zprimus, format='(f20.3)'),1), color=axisColor
;                        xyouts, currentIP.dx/redh100()-0.9*halfRange, currentIP.dy/redh100()+0.6*halfRange, "Field: " + strtrim(currentIP.field,1), color=axisColor, charsize=xyoutscharsize

		      ; SET UP AXES BASED ON FIRST IP
			IF (IP_index eq 0) THEN BEGIN
			      ; DRAW AXES
				IF SFQ eq 1 THEN SFQstatus="late" ELSE SFQstatus="early"

				xmin = currentIP.dx/redh100()-halfRange
				xmax = currentIP.dx/redh100()+halfRange
				xr = xmax-xmin

	                        PLOT, [currentIP.dx/redh100()], [H0*currentIP.dz/redh100()], psym=IPsym, color=axisColor, $
	                                xrange=[currentIP.dx/redh100()-halfRange, currentIP.dx/redh100()+halfRange], $
	                                yrange=[ymin2, ymax2], xtitle="Relative Mpc", ytitle=textoidl("km s^{-1}"), $
					title=  "z="+strtrim(string(z,format='(f20.2)'),1) + textoidl("\pm") + strtrim(string(dz,format='(f20.2)'),1) + $
						"   IP: "+SFQstatus+ $
						"   Field: "+string(field), /nodata
			      ; PLOT FIRST IP
	                        OPLOT, [currentIP.dx/redh100()], [H0*currentIP.dz/redh100()], psym=IPsym, color=IPcolor
	
			      ; OVERPLOT SF AND Q NON-IP GALAXIES OF FIRST IP
;				OPLOT, [(dataLocal[where(dataLocal.SFQ eq 1)].dx)/redh100()], [(H0*dataLocal[where(dataLocal.SFQ eq 1)].dz)/redh100()], psym=nonIPsym, color=nonIPcolorSF
;				OPLOT, [(dataLocal[where(dataLocal.SFQ eq 0)].dx)/redh100()], [(H0*dataLocal[where(dataLocal.SFQ eq 0)].dz)/redh100()], psym=nonIPsym, color=nonIPcolorQ

		      ; OVERPLOT SUBSEQUENT IP AND NON-IP GALAXIES RELATIVE TO FIRST
			ENDIF ELSE BEGIN
				xOffset = (dataIP[IP_index].dx-dataIP[0].dx)/redh100()
				yOffset = H0*(dataIP[IP_index].dz-dataIP[0].dz)/redh100()

			      ; PLOT IP
;				OPLOT, [currentIP.dx/redh100()-xOffset], [H0*currentIP.dz/redh100()-yOffset], psym=IPsym, color=IPcolor
		
			      ; PLOT SF AND Q NON-IP GALAXIES
	                        OPLOT, [(dataLocal[where(dataLocal.SFQ eq 1)].dx)/redh100()-xOffset], [(H0*dataLocal[where(dataLocal.SFQ eq 1)].dz)/redh100()-yOffset], psym=nonIPsym, color=nonIPcolorSF
	                        OPLOT, [(dataLocal[where(dataLocal.SFQ eq 0)].dx)/redh100()-xOffset], [(H0*dataLocal[where(dataLocal.SFQ eq 0)].dz)/redh100()-yOffset], psym=nonIPsym, color=nonIPcolorQ
			ENDELSE

			plot_index += 1
		ENDIF

;                IF (dataLocalInner ne !null) THEN BEGIN
;                        i = 0
;                        FOREACH element, dataLocalInner DO BEGIN
;                                xyouts, currentIP.dx/redh100(), ymin2+0.10*yr2*(i+1), $
;                                  textoidl("non-IP M_{\odot} ") + strtrim(string(element.mstar, format='(f20.2)'),1), charsize=xyoutscharsize
;                                i += 1
;                        ENDFOREACH
;                ENDIF

;		ENDIF
	ENDFOR

        OPLOT, [dataIP[0].dx/redh100()-halfRange, dataIP[0].dx/redh100()+halfRange], [ymin2+(delta_z-1500.), ymin2+(delta_z-1500.)], linestyle=3, color=dashedLineColor
        OPLOT, [dataIP[0].dx/redh100()-halfRange, dataIP[0].dx/redh100()+halfRange], [ymax2-(delta_z-1500.), ymax2-(delta_z-1500.)], linestyle=3, color=dashedLineColor

        OPLOT, [dataIP[0].dx/redh100()-0.5, dataIP[0].dx/redh100()-0.5], [ymin2, ymax2], linestyle=3, color=dashedLineColor
        OPLOT, [dataIP[0].dx/redh100()+0.5, dataIP[0].dx/redh100()+0.5], [ymin2, ymax2], linestyle=3, color=dashedLineColor

	XYOUTS, xmin-0.30*xr, ymin2-0.05*yr2, textoidl("\Deltav=") + strtrim(string(delta_z, format='(i20)'),1) + textoidl(" km s^{-1}")
	XYOUTS, xmin-0.30*xr, ymin2-0.15*yr2, strtrim(plot_index,1) + "/" + strtrim(n_elements(dataIP),1)
	XYOUTS, xmin-0.30*xr, ymin2-0.20*yr2, "IPs shown"

;	PRINT, strtrim(plot_index,1) + " of " + strtrim(n_elements(dataIP),1) + " IP galaxies within z=" $
;		+ strtrim(string(z,format='(f20.2)'),1) + "+/-" + strtrim(string(dz,format='(f20.3)'),1) + " plotted."
	
	PRINT, "z = "+strtrim(string(z,format='(f20.2)'),1) + "+/-" + strtrim(string(dz,format='(f20.2)'),1) + "    IP type: " + SFQstatus, "  Cylinder depth (km/s): ", delta_z

	FOR s=0,n_slices-1 DO BEGIN
		PRINT, "Slice: ", strtrim(s+1,1), "  Total early: ", slice_totals_early[s], "  Total late: ", slice_totals_late[s]
	ENDFOR

	RETURN
END
