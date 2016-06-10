PRO latefrac_VisualTest

	SET_PLOT, 'ps'
	DEVICE, file="latefrac_VisualTest.ps", /landscape

	!p.multi=[0,3,2]

	!x.margin=[7,2.5] ; left, right
	!y.margin=[3.5,1] ; bottom, top

	data = mrdfits("zerodSFQ_IP_dv1500kms_dm0.5.fits", 1)
	dataIP = data[where(data.IP eq 1)]	

	IP_index = 1
	delta_z  = 1500. ; km/s
	HO 	 = 100.*redh100()

	halfRange = 4. ; Mpc

	xyoutscharsize = 1

;	axisColor = cgColor('White')
;	textColor = cgColor('White')
	axisColor = cgColor('Black')
	textColor = cgColor('Black')


	IPcolorSF = cgColor('Blue')
	IPcolorQ  = cgColor('Red')
	IPsym 	  = 4

	nonIPcolorSF = cgColor('Blue')
	nonIPcolorQ = cgColor('Red')
	nonIPsym   = 1

        A = FINDGEN(33)*(!PI*2/32.)
	
	FOR IP_index=0,499 DO BEGIN
		currentIP = dataIP[IP_index]
		dataLocal = data[where( (data.field eq currentIP.field) AND $
					(ABS(HO*data.dz/redh100() - HO*currentIP.dz/redh100()) le delta_z) AND $
					( (((currentIP.dx-data.dx)/redh100())^2 + ((currentIP.dy-data.dy)/redh100())^2)^0.5 le halfRange ) AND $
					(currentIP.objname ne data.objname) )]
			dataLocalInner = dataLocal[where(( ((currentIP.dx-dataLocal.dx)/redh100())^2 + ((currentIP.dy-dataLocal.dy)/redh100())^2 )^0.5 le 0.5, /null)] 

		IF n_elements(dataLocalInner) gt 0 THEN BEGIN
			print, n_elements(dataLocal), currentIP.dx/redh100(), currentIP.dy/redh100(), dataLocal[0].dx/redh100(), dataLocal[0].dy/redh100()

			IF currentIP.SFQ eq 1 THEN IPcolor = IPcolorSF ELSE IPcolor = IPcolorQ

			PLOT, [currentIP.dx/redh100()], [currentIP.dy/redh100()], psym=IPsym, color=axisColor, $
				xrange=[currentIP.dx/redh100()-halfRange, currentIP.dx/redh100()+halfRange], $
				yrange=[currentIP.dy/redh100()-halfRange, currentIP.dy/redh100()+halfRange], $
				xtitle="Mpc", ytitle="Mpc", /nodata
			OPLOT, [currentIP.dx/redh100()], [currentIP.dy/redh100()], psym=IPsym, color=IPcolor

			OPLOT, currentIP.dx/redh100()+0.5*COS(A), currentIP.dy/redh100()+0.5*SIN(A), linestyle=4, color=cgColor('Green')

			OPLOT, [(dataLocal[where(dataLocal.SFQ eq 1)].dx)/redh100()], [(dataLocal[where(dataLocal.SFQ eq 1)].dy)/redh100()], psym=nonIPsym, color=nonIPcolorSF
			OPLOT, [(dataLocal[where(dataLocal.SFQ eq 0)].dx)/redh100()], [(dataLocal[where(dataLocal.SFQ eq 0)].dy)/redh100()], psym=nonIPsym, color=nonIPcolorQ
			
			IF dataLocalInner ne !null THEN BEGIN
				xyouts, currentIP.dx/redh100()-0.9*halfRange, currentIP.dy/redh100()-0.9*halfRange, "IP mass: " + strtrim(string(currentIP.mstar, format='(f20.2)'),1), color=IPcolor, charsize=xyoutscharsize
			
				i = 0
				FOREACH element, dataLocalInner DO BEGIN
					xyouts, currentIP.dx/redh100()-0.9*halfRange, currentIP.dy/redh100()-0.9*halfRange+(i+1)*0.50, strtrim(string(element.mstar, format='(f20.2)'),1), charsize=xyoutscharsize
					i += 1
				ENDFOREACH
			ENDIF

;			BREAK
		ENDIF	
	ENDFOR
	
	DEVICE, /CLOSE
	SET_PLOT, 'X'

END

