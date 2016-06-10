; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

PRO latefrac_default_params, dz_coeff, zmin, zmax ; , printEvery
	path = 'IP_data/'
	zerodInputFile = 'zerodSFQ_IP_dz' + strtrim(strcompress(string(dz_coeff, format='(f20.1)')),1) + '_dm0.5.fits'

	zerodInputRaw = mrdfits(path+zerodInputFile, 1)

	n_annuli 	= 16
	dRproj		= 0.25
	printEvery	= 1000

	zerodInput = zerodInputRaw[where((zerodInputRaw.zprimus ge zmin) and (zerodInputRaw.zprimus le zmax))]
	PRINT, "total galaxies: ", n_elements(zerodInputRaw), $
	  " zmin:", zmin, $
	  " zmax:", zmax, $
	  " galaxies within z range:", n_elements(zerodInput)

	Rproj_array = float(dRproj)*findgen(n_annuli+1)

	fields	= zerodInput[uniq(zerodInput.field, sort(zerodInput.field))].field
	n_fields = n_elements(fields)

	zerodInput_isIP = zerodInput[where(zerodInput.IP eq 1)]

	n_tot_IPSF 	= 0.0*findgen(n_annuli)
	n_late_IPSF 	= 0.0*findgen(n_annuli)
	n_tot_IPQ 	= 0.0*findgen(n_annuli)
	n_late_IPQ 	= 0.0*findgen(n_annuli)

	outputIndex	= 0L

; LOOP THROUGH FIELDS
	FOR f=0,(n_fields-1) DO BEGIN 
		currentField	= fields[f]
		currentField_IP	= zerodInput_isIP[ where( zerodInput_isIP.field eq currentField, /null ) ]
		n_IP		= n_elements(currentField_IP)
		currentField_all= zerodInput[ where( zerodInput.field eq currentField, /null ) ]
		n_all		= n_elements(currentField_all)
		PRINT, "Field: ", currentField, $
		       "total galaxies in field:", n_all, $
		       "IP galaxies in field: ", n_IP

; LOOP THROUGH ALL IP GALAXIES IN FIELD
		FOR i=0,(n_IP-1) DO BEGIN 
			currentIP 	= currentField_IP[i]
			c		= 299792.458 ; c in km/s
;			vcurr		= (100.0*redh100())*(currentIP.dz/redh100()/(1+currentIP.zprimus)) 	; radial velocity (km/s)
;					; (H0 in km/s/Mpc )*(physical distance in Mpc (NOT Mpc/h))
			IF outputIndex mod printEvery eq 0 THEN BEGIN
				PRINT, outputIndex
;				print, "Index:" + string(strcompress(outputIndex)) + ",  vcurr:" + string(strcompress(vcurr)) + ",  SFQ:" + string(strcompress(currentIP.SFQ))
			ENDIF
; REDSHIFT SPACE CUT
;			nonIP_dzcut	= currentField_all[ where( (currentIP.objname ne currentField_all.objname) $
;			  , /null ) ]
;			  AND ( ABS( vcurr - (100.0*redh100())*(currentField_all.dz/redh100()/(1+currentField_all.zprimus)) ) lt float(deltav) ), /null ) ]
			
			nonIP_dzcut	= currentField_all[ where( ABS(currentField_all.zprimus - currentIP.zprimus) le dz_coeff*0.005*(1.+currentIP.zprimus), /null ) ]

			n_dzcut		= n_elements(nonIP_dzcut)
;			print, "n_dzcut:" + string(strcompress(n_dzcut))

; COUNTS TOTAL # OF NON-IP AND # OF LATE-TYPE NON-IP IN EACH PROJECTED RADIAL ANNULUS
			FOR j=0,(n_annuli-1) DO BEGIN
				IF (nonIP_dzcut ne !null) THEN BEGIN
					annulus_all = nonIP_dzcut[ where( $
					  ( ( (currentIP.xprop-nonIP_dzcut.xprop)^2 + (currentIP.yprop-nonIP_dzcut.yprop)^2 )^0.5 ge Rproj_array[j] ) AND $
					  ( ( (currentIP.xprop-nonIP_dzcut.xprop)^2 + (currentIP.yprop-nonIP_dzcut.yprop)^2 )^0.5 lt Rproj_array[j+1] ), /null ) ]
				ENDIF ELSE BEGIN
					annulus_all = []
				ENDELSE
				IF ( annulus_all ne !null ) THEN ( annulus_late = annulus_all[where( annulus_all.SFQ eq 1, /null )] ) ELSE ( annulus_late=[] )
	
				IF (currentIP.SFQ eq 1) THEN BEGIN
					n_tot_IPSF[j] 	+= n_elements(annulus_all)
					n_late_IPSF[j] 	+= n_elements(annulus_late)
					IF outputIndex mod printEvery eq 0 THEN BEGIN
						print, Rproj_array[j], " ", Rproj_array[j+1], " ", n_elements(annulus_all), " ", n_elements(annulus_late), $
						  " ", n_tot_IPSF[j], " ", n_late_IPSF[j]
					ENDIF
				ENDIF ELSE BEGIN
					n_tot_IPQ[j] 	+= n_elements(annulus_all)
					n_late_IPQ[j] 	+= n_elements(annulus_late)
					IF outputIndex mod printEvery eq 0 THEN BEGIN
						print, Rproj_array[j], " ", Rproj_array[j+1], " ", n_elements(annulus_all), " ", n_elements(annulus_late), $
						  " ", n_tot_IPQ[j], " ", n_late_IPQ[j]
					ENDIF
				ENDELSE
			ENDFOR
			
			outputIndex += 1
		ENDFOR
	ENDFOR

	templateRow 	= create_struct('Rmin', 0.0, 'Rmax', 0.0, 'n_tot_IPSF', 0.0, 'n_late_IPSF', 0.0, 'n_tot_IPQ', 0.0, 'n_late_IPQ', 0.0)
	outputStruct 	= replicate(templateRow, n_annuli)

	FOR i=0,n_annuli-1 DO BEGIN
		newRow	= create_struct('Rmin', Rproj_array[i], 'Rmax', Rproj_array[i+1], $
		  'n_tot_IPSF', n_tot_IPSF[i], 'n_late_IPSF', n_late_IPSF[i], 'n_tot_IPQ', n_tot_IPQ[i], 'n_late_IPQ', n_late_IPQ[i])
		print, newRow
		outputStruct[i] = newRow
	ENDFOR

	mwrfits, outputStruct, "latefrac_data/latefrac_" + strtrim(strcompress(string(zmin, format='(f20.1)')), 1) + "_" $
					   + strtrim(strcompress(string(zmax, format='(f20.1)')) ,1) + "_" $
					   + zerodInputFile, /create
END
