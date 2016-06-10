; uses histograms and delta_z <= dz_coeff*0.005*(1+z)

PRO latefrac_hist_byField, field; zmin, zmax, dz_coeff
	dz_coeff = 2.
	zmin = 0.2
	zmax = 1.0
	path = 'IP_data/'
	zerodInputFile = 'zerodSFQ_IP_dz' + strtrim(strcompress(string(dz_coeff, format='(f20.1)')),1) + '_dm0.0.fits'

	PRINT, "Input data: ", zerodInputFile

	zerodInput = mrdfits(path+zerodInputFile, 1)

	n_annuli 	= 16
	dRproj		= 0.25
	printEvery	= 100

	data = zerodInput[where((zerodInput.zprimus ge zmin) AND $
				(zerodInput.zprimus le zmax) AND $
				(strtrim(zerodInput.field) eq field) )]

	PRINT, "Field: ", field, $ 
	  " total galaxies: ", n_elements(zerodInput), $
	  " zmin:", zmin, $
	  " zmax:", zmax, $
	  " galaxies within z range:", n_elements(data)

	Rproj_array = float(dRproj)*findgen(n_annuli+1)

	data_isIP = data[where(data.IP eq 1)]

	n_tot_IPSF 	= 0.0*findgen(n_annuli)
	n_late_IPSF 	= 0.0*findgen(n_annuli)
	n_tot_IPQ 	= 0.0*findgen(n_annuli)
	n_late_IPQ 	= 0.0*findgen(n_annuli)

	outputIndex	= 0L

	c = 299792.458 ; c in km/s

; LOOP THROUGH ALL STAR-FORMING IP GALAXIES
	IPSF 	= data_isIP[where(data_isIP.SFQ eq 1)]
	PRINT, 'SF IP galaxies:', n_elements(IPSF)

	FOR i=0,(n_elements(IPSF)-1) DO BEGIN 
		currentIP 	= IPSF[i]

		; keep galaxies within dz=dz+coeff*0.005*(1+z) of IP candidate (and in same field) 
		dz_slice_all 	= data[ where( (data.field eq currentIP.field) AND $
			(ABS(data.zprimus - currentIP.zprimus) le dz_coeff*0.005*(1.+currentIP.zprimus)) AND $
			(data.objname ne currentIP.objname) ) ]
		dz_slice_all_dists = ((currentIP.xprop - dz_slice_all.xprop)^2 + $
				   (currentIP.yprop - dz_slice_all.yprop)^2)^(0.5)
		; add stats for current IP to n_tot_IPSF
		dR_array_tot 	= histogram(dz_slice_all_dists, bin=dRproj, min=0., nbins=n_annuli)
		n_tot_IPSF 	= n_tot_IPSF + dR_array_tot[0:15]
		
		dz_slice_late	= dz_slice_all[where(dz_slice_all.SFQ eq 1)]
		dz_slice_late_dists = ((currentIP.xprop - dz_slice_late.xprop)^2 + $
				   (currentIP.yprop - dz_slice_late.yprop)^2)^(0.5)
		dR_array_late	= histogram(dz_slice_late_dists, bin=dRproj, min=0., nbins=n_annuli)
		n_late_IPSF	= n_late_IPSF + dR_array_late[0:15]

		IF outputIndex mod printEvery eq 0 THEN PRINT, outputIndex

		outputIndex += 1
	ENDFOR

; LOOP THROUGH ALL QUIESCENT GALAXIES
	IPQ 	= data_isIP[where(data_isIP.SFQ eq 0)]
	PRINT, 'Q IP galaxies:', n_elements(IPQ)

	FOR i=0,(n_elements(IPQ)-1) DO BEGIN 
		currentIP 	= IPQ[i]

		; keep galaxies within dz=dz+coeff*0.005*(1+z) of IP candidate (and in same field) 
		dz_slice_all 	= data[ where( (data.field eq currentIP.field) AND $
			(ABS(data.zprimus - currentIP.zprimus) le dz_coeff*0.005*(1.+currentIP.zprimus)) AND $
			(data.objname ne currentIP.objname) ) ]
		dz_slice_all_dists = ((currentIP.xprop - dz_slice_all.xprop)^2 + $
				   (currentIP.yprop - dz_slice_all.yprop)^2)^(0.5)
		; add stats for current IP to n_tot_IPQ
		dR_array_tot 	= histogram(dz_slice_all_dists, bin=dRproj, min=0., nbins=n_annuli)
		n_tot_IPQ 	= n_tot_IPQ + dR_array_tot[0:15]
		
		dz_slice_late	= dz_slice_all[where(dz_slice_all.SFQ eq 1)]
		dz_slice_late_dists = ((currentIP.xprop - dz_slice_late.xprop)^2 + $
				   (currentIP.yprop - dz_slice_late.yprop)^2)^(0.5)
		dR_array_late	= histogram(dz_slice_late_dists, bin=dRproj, min=0., nbins=n_annuli)
		n_late_IPQ	= n_late_IPQ + dR_array_late[0:15]

		IF outputIndex mod printEvery eq 0 THEN PRINT, outputIndex

		outputIndex += 1
	ENDFOR

	PRINT, Rproj_array
	PRINT, 'n_tot_IPSF:', n_tot_IPSF
	PRINT, 'n_late_IPSF:', n_late_IPSF
	PRINT, 'n_tot_IPQ:', n_tot_IPQ
	PRINT, 'n_late_IPQ:', n_late_IPQ

        templateRow     = create_struct('Rmin', 0.0, 'Rmax', 0.0, 'n_tot_IPSF', 0.0, 'n_late_IPSF', 0.0, 'n_tot_IPQ', 0.0, 'n_late_IPQ', 0.0)
        outputStruct    = replicate(templateRow, n_annuli)

	FOR i=0,n_annuli-1 DO BEGIN
		newRow	= create_struct('Rmin', Rproj_array[i], 'Rmax', Rproj_array[i+1], $
		  'n_tot_IPSF', n_tot_IPSF[i], 'n_late_IPSF', n_late_IPSF[i], 'n_tot_IPQ', n_tot_IPQ[i], 'n_late_IPQ', n_late_IPQ[i])
		print, newRow
		outputStruct[i] = newRow
	ENDFOR

	mwrfits, outputStruct, "latefrac_data_hist_byField/latefrac_" + strtrim(strcompress(string(zmin, format='(f20.1)')), 1) + "_" $
					   + strtrim(strcompress(string(zmax, format='(f20.1)')) ,1) + "_" $
					   + field + "_" + zerodInputFile, /create
END
