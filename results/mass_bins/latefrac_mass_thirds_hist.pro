; uses histograms and delta_z <= dz_coeff*0.005*(1+z)

PRO latefrac_mass_thirds_hist, input_dm, third;, zmin, zmax;, dz_coeff
;	input_dm = 1.0

	dz_coeff = 2.
	zmin = 0.2
	zmax = 1.0

	string_third = strtrim(strcompress(string(third, format='(i20)')),1)
	string_dm = strtrim(strcompress(string(input_dm, format='(f20.1)')),1)

	allDataPath	= '~/results/variable_mass+' + string_dm + 'dex/IP_data/'
	zerodInputFile	= 'zerodSFQ_IP_dz2.0.fits'

	zerodInput = mrdfits(allDataPath + zerodInputFile, 1)

	PRINT, 'All data: ' + allDataPath + zerodInputFile

	n_annuli 	= 60
	dRproj		= 0.25
	printEvery	= 1000

	data = zerodInput[where((zerodInput.zprimus ge zmin) and (zerodInput.zprimus le zmax))]

	PRINT, 'total galaxies: ', n_elements(zerodInput), $
	  ' zmin:', zmin, $
	  ' zmax:', zmax, $
	  ' galaxies within z range:', n_elements(data)

	Rproj_array = float(dRproj)*findgen(n_annuli+1)

	data_isIP = mrdfits('~/results/mass_bins/IP_data/variable_mass+' + string_dm + 'dex/dataIP_T' + string_third + '.fits', 1)

	n_tot_IPSF 	= 0.0*findgen(n_annuli)
	n_late_IPSF 	= 0.0*findgen(n_annuli)
	n_tot_IPQ 	= 0.0*findgen(n_annuli)
	n_late_IPQ 	= 0.0*findgen(n_annuli)

	outputIndex	= 0L

	c = 299792.458 ; c in km/s

; LOOP THROUGH ALL STAR-FORMING IP GALAXIES
	IPSF 	= data_isIP[where(data_isIP.SFQ eq 1)]
	PRINT, 'SF IP galaxies:', n_elements(IPSF)

	subtotal_neighbors_IPSF	= 0L*indgen(n_annuli)
	FOR i=0,(n_elements(IPSF)-1) DO BEGIN 
		currentIP 	= IPSF[i]

		; keep galaxies within dz=dz+coeff*0.005*(1+z) of IP candidate (and in same field) 
		dz_slice_all 	= data[ where( (data.field eq currentIP.field) AND (data.objname ne currentIP.objname) AND $
			(ABS(data.zprimus - currentIP.zprimus) le dz_coeff*0.005*(1.+currentIP.zprimus)) ) ]
		dz_slice_all_dists = SQRT( (currentIP.xprop - dz_slice_all.xprop)^2 + (currentIP.yprop - dz_slice_all.yprop)^2 )
	
		FOR n=0,n_annuli-1 DO BEGIN
			currentMinR = dRproj*float(n)
			currentMaxR = dRproj*float(n+1)
			annulus_dists = dz_slice_all_dists[where( (dz_slice_all_dists GE currentMinR) AND (dz_slice_all_dists LE currentMaxR), /NULL )]
;			subtotal_annuli_dists = dz_slice_all_dists[where( (dz_slice_all_dists LE currentMaxR), /NULL )]
;			PRINT, currentMinR, ' ', currentMaxR, ' ', n_elements(annulus_dists), ' ', n_elements(subtotal_annuli_dists)
			subtotal_neighbors_IPSF[n] += n_elements(dz_slice_all_dists[where( (dz_slice_all_dists LE currentMaxR), /NULL) ])
;			subtotal_neighbors_IPSF[n] += n_elements(annulus_dists)
		ENDFOR
;		PRINT, strcompress(outputIndex), ') ', subtotal_neighbors_IPSF

		; add stats for current IP to n_tot_IPSF
		dR_array_tot 	= histogram(dz_slice_all_dists, bin=dRproj, min=0., nbins=n_annuli)
		n_tot_IPSF 	= n_tot_IPSF + dR_array_tot[0:n_annuli-1]
		
		dz_slice_late	= dz_slice_all[where(dz_slice_all.SFQ eq 1)]
		dz_slice_late_dists = ((currentIP.xprop - dz_slice_late.xprop)^2 + $
				   (currentIP.yprop - dz_slice_late.yprop)^2)^(0.5)
		dR_array_late	= histogram(dz_slice_late_dists, bin=dRproj, min=0., nbins=n_annuli)
		n_late_IPSF	= n_late_IPSF + dR_array_late[0:n_annuli-1]

		IF outputIndex mod printEvery eq 0 THEN PRINT, outputIndex
		outputIndex += 1
	ENDFOR

; LOOP THROUGH ALL QUIESCENT GALAXIES
	IPQ 	= data_isIP[where(data_isIP.SFQ eq 0)]
	PRINT, 'Q IP galaxies:', n_elements(IPQ)

	subtotal_neighbors_IPQ 	= 0L*indgen(n_annuli)
	FOR i=0,(n_elements(IPQ)-1) DO BEGIN 
		currentIP 	= IPQ[i]

		; keep galaxies within dz=dz+coeff*0.005*(1+z) of IP candidate (and in same field) 
		dz_slice_all 	= data[ where( (data.field eq currentIP.field) AND (data.objname ne currentIP.objname) AND $
			(ABS(data.zprimus - currentIP.zprimus) le dz_coeff*0.005*(1.+currentIP.zprimus)) ) ]
		dz_slice_all_dists = SQRT( (currentIP.xprop - dz_slice_all.xprop)^2 + (currentIP.yprop - dz_slice_all.yprop)^2 )
	
		FOR n=0,n_annuli-1 DO BEGIN
			currentMinR = dRproj*float(n)
			currentMaxR = dRproj*float(n+1)
			annulus_dists = dz_slice_all_dists[where( (dz_slice_all_dists GE currentMinR) AND (dz_slice_all_dists LE currentMaxR), /NULL )]
;			subtotal_annuli_dists = dz_slice_all_dists[where( (dz_slice_all_dists LE currentMaxR), /NULL )]
;			PRINT, currentMinR, ' ', currentMaxR, ' ', n_elements(annulus_dists), ' ', n_elements(subtotal_annuli_dists)
			subtotal_neighbors_IPQ[n] += n_elements(dz_slice_all_dists[where( (dz_slice_all_dists LE currentMaxR), /NULL) ])
;			subtotal_neighbors_IPQ[n] += n_elements(annulus_dists)
		ENDFOR
;		PRINT, strcompress(outputIndex), ') ', subtotal_neighbors_IPQ

		; add stats for current IP to n_tot_IPQ
		dR_array_tot 	= histogram(dz_slice_all_dists, bin=dRproj, min=0., nbins=n_annuli)
		n_tot_IPQ 	= n_tot_IPQ + dR_array_tot[0:n_annuli-1]
		
		dz_slice_late	= dz_slice_all[where(dz_slice_all.SFQ eq 1)]
		dz_slice_late_dists = ((currentIP.xprop - dz_slice_late.xprop)^2 + $
				   (currentIP.yprop - dz_slice_late.yprop)^2)^(0.5)
		dR_array_late	= histogram(dz_slice_late_dists, bin=dRproj, min=0., nbins=n_annuli)
		n_late_IPQ	= n_late_IPQ + dR_array_late[0:n_annuli-1]

		IF outputIndex mod printEvery eq 0 THEN PRINT, outputIndex
		outputIndex += 1
	ENDFOR

;	PRINT, Rproj_array
;	PRINT, 'n_tot_IPSF:', n_tot_IPSF
;	PRINT, 'n_late_IPSF:', n_late_IPSF
;	PRINT, 'n_tot_IPQ:', n_tot_IPQ
;	PRINT, 'n_late_IPQ:', n_late_IPQ

;	templateRow     = create_struct('Rmin', 0.0, 'Rmax', 0.0, 'n_tot_IPSF', 0.0, 'n_late_IPSF', 0.0, 'n_tot_IPQ', 0.0, 'n_late_IPQ', 0.0)
;	outputStruct    = replicate(templateRow, n_annuli)

	outputStruct = mrd_struct( ['Rmin','Rmax','n_tot_IPSF','n_late_IPSF','subtotal_neighbors_IPSF','n_tot_IPQ','n_late_IPQ','subtotal_neighbors_IPQ'], ['0.0', '0.0', '0.0', '0.0', '0L', '0.0', '0.0', '0L'], n_annuli )
;		[Rproj_array[i], Rproj_array[i+1], n_tot_IPSF[i], n_late_IPSF[i], subtotal_neighbors_IPSF[i], n_tot_IPQ[i], n_late_IPQ[i], subtotal_neighbors_IPQ[i]] )

	FOR i=0,n_annuli-1 DO BEGIN
		outputStruct[i].Rmin = Rproj_array[i]
		outputStruct[i].Rmax = Rproj_array[i+1]
	ENDFOR

	outputStruct.n_tot_IPSF		= n_tot_IPSF
	outputStruct.n_late_IPSF	= n_late_IPSF
	outputStruct.subtotal_neighbors_IPSF = subtotal_neighbors_IPSF
	outputStruct.n_tot_IPQ		= n_tot_IPQ
	outputStruct.n_late_IPQ		= n_late_IPQ
	outputStruct.subtotal_neighbors_IPQ = subtotal_neighbors_IPQ

	FOREACH element,outputStruct DO PRINT, element
;	newRow = create_struct( ['Rmin','Rmax','n_tot_IPSF','n_late_IPSF','subtotal_neighbors_IPSF','n_tot_IPQ','n_late_IPQ','subtotal_neighbors_IPQ'], $
;					[Rproj_array[i], Rproj_array[i+1], n_tot_IPSF[i], n_late_IPSF[i], subtotal_neighbors_IPSF[i], n_tot_IPQ[i], n_late_IPQ[i], subtotal_neighbors_IPQ[i]] )
;	print, newRow
;	outputStruct[i] = newRow

	mwrfits, outputStruct, '~/results/mass_bins/latefrac_data_hist/variable_mass+' + string_dm + 'dex/latefrac_' $
				+ strtrim(strcompress(string(zmin, format='(f20.1)')), 1) + '_' $
				+ strtrim(strcompress(string(zmax, format='(f20.1)')) ,1) + '_dataIP_T' + string_third + '.fits', /create
END
