; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)
; PNZ = permute neighbor redshifts

PRO latefrac_conservative_targ_weight_PNZ, outputFormat;, zmin, zmax;, dz_coeff, printEvery

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/latefracplot_targ_weight_conservative_mass_zall_PNZ', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

        IPdataPath = '~/conformity/results/conservative_mass_cutoff/IP_data/'
        zerodInputFile = 'zerodSFQ_IP_dz' + strtrim(strcompress(string(dz_coeff, format='(f20.1)')),1) + '_dm0.0.fits'

	PRINT, 'Input data: ', zerodInputFile

        zerodInput = mrdfits(IPdataPath + zerodInputFile, 1)

	Rmax 		= 15.
	dRproj		= 0.25
	n_annuli 	= Ceil(Rmax/dRproj)
	printEvery	= 100

	data = zerodInput[where((zerodInput.zprimus GE zmin) and (zerodInput.zprimus LE zmax) AND (zerodInput.targ_weight GE 1.))]

	fields = data[uniq(data.field, sort(data.field))].field

	tempData = MRD_STRUCT('zshuffle', '0.0', n_elements(data))

        combine_structs, data, tempData, dataPlus
	data = dataPlus

        allData = []
        FOR f=0,n_elements(fields)-1 DO BEGIN
                fieldData = data[where(data.field EQ fields[f])]
                fieldData.zshuffle = SHUFFLE(fieldData.zprimus)
                allData = [allData, fieldData]
        ENDFOR

	data = allData

	PLOT, data.zprimus, data.zshuffle, psym=3

	!p.multi=0
	!p.charsize=1.1

	Rproj_array = float(dRproj)*(findgen(n_annuli+1)+0.5)

	xmin = 0
	xmax = n_annuli*dRproj
	xr = xmax-xmin

	ymin = 0.35
	ymax = 1
	yr = ymax-ymin

	get_targ_weight, data, dz_coeff, 1, neigh_all_grand_total, neigh_late_grand_total, n_annuli, dRproj, printEvery

	n_tot_IPSF  = neigh_all_grand_total
	n_late_IPSF = neigh_late_grand_total

	PLOT, Rproj_array, float(neigh_late_grand_total)/neigh_all_grand_total, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Late-type Fraction', $
	  title = 'M13 Mass Limit ' + textoidl('(\Deltaz=2.0)'), /NODATA
	LEGEND, ['Late-type IP', 'Early-type IP'], LINESTYLE=[0,2], BOX=0, /BOTTOM, /LEFT
	OPLOT, Rproj_array, float(neigh_late_grand_total)/neigh_all_grand_total, LINESTYLE=0

	get_targ_weight, data, dz_coeff, 0, neigh_all_grand_total, neigh_late_grand_total, n_annuli, dRproj, printEvery

	n_tot_IPQ  = neigh_all_grand_total
	n_late_IPQ = neigh_late_grand_total

	OPLOT, Rproj_array, float(neigh_late_grand_total)/neigh_all_grand_total, LINESTYLE=2

	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, '0.2 < z < 1.0', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.85*yr, 'NEIGHBOR REDSHIFTS PERMUTED WITHIN FIELD', ALIGNMENT=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'

        templateRow     = create_struct('Rmin', 0.0, 'Rmax', 0.0, 'n_tot_IPSF', 0.0, 'n_late_IPSF', 0.0, 'n_tot_IPQ', 0.0, 'n_late_IPQ', 0.0)
        outputStruct    = replicate(templateRow, n_annuli)

        FOR i=0,n_annuli-1 DO BEGIN
                newRow  = create_struct('Rmin', Rproj_array[i], 'Rmax', Rproj_array[i+1], $
                  'n_tot_IPSF', n_tot_IPSF[i], 'n_late_IPSF', n_late_IPSF[i], 'n_tot_IPQ', n_tot_IPQ[i], 'n_late_IPQ', n_late_IPQ[i])
                print, newRow
                outputStruct[i] = newRow
        ENDFOR

        mwrfits, outputStruct, '~/conformity/results/permute_neighbor_z/latefrac_' + strtrim(strcompress(string(zmin, format='(f20.1)')), 1) + '_' $
				+ strtrim(strcompress(string(zmax, format='(f20.1)')) ,1) + '_targ_weight_' + zerodInputFile, /create
END

PRO get_targ_weight, data, dz_coeff, IPSFstatus, neigh_all_grand_total, neigh_late_grand_total, n_annuli, dRproj, printEvery
	data_isIP = data[where(data.IP EQ 1)]

	allIP_total_array = [] ; list of all (weighted) neighbors around all IPs by annulus
	allIP_late_array  = [] ; list of late-type (weighted) neighbors around all IPs by annulus

        IP = data_isIP[where(data_isIP.SFQ EQ IPSFstatus)]
	IF (IPSFstatus EQ 1) THEN (IPtype = 'Star-forming') ELSE (IPtype = 'Quiescent')
	PRINT, IPtype + ' IP galaxies: ' + strtrim(n_elements(IP),1)

	outputIndex = 0L

	FOR i=0,(n_elements(IP)-1) DO BEGIN
;	FOR i=0,19 DO BEGIN
                currentIP = IP[i]

                ; keep galaxies within dz=dz+coeff*0.005*(1+z) of IP candidate (and in same field)
                all_dz_neigh = data[ where( (data.field EQ currentIP.field) AND (data.objname NE currentIP.objname) AND $
                        (ABS(data.zshuffle - currentIP.zprimus) le dz_coeff*0.005*(1.+currentIP.zprimus)) ) ]
		all_dz_neigh_dists = [SQRT((all_dz_neigh.xprop-currentIP.xprop)^2 + (all_dz_neigh.yprop-currentIP.yprop)^2)]

		currentIP_total_array = [] ; weight total of all neighbors by annulus for current IP
		currentIP_late_array  = [] ; weight total of late-type neighbors by annulus for current IP

		; GET LATE FRACTION FOR EACH ANNULUS
                FOR n=0,n_annuli-1 DO BEGIN
			currentRmin = dRproj*float(n)
			currentRmax = dRproj*float(n+1)
			
			annulus_neigh = all_dz_neigh[where((all_dz_neigh_dists GE currentRmin) AND (all_dz_neigh_dists LE currentRmax), /NULL)]
			IF (annulus_neigh NE !NULL) THEN $
				annulus_neigh_SF = annulus_neigh[where(annulus_neigh.SFQ EQ 1, /NULL)] ELSE $
				annulus_neigh_SF = []

			IF (annulus_neigh NE !NULL) THEN $
				currentIP_annulus_total = TOTAL(annulus_neigh.targ_weight) ELSE $
				currentIP_annulus_total = 0
			IF (annulus_neigh_SF NE !NULL) THEN $
				currentIP_annulus_late = TOTAL(annulus_neigh_SF.targ_weight) ELSE $
				currentIP_annulus_late = 0
			
			currentIP_total_array = [currentIP_total_array, currentIP_annulus_total]
			currentIP_late_array  = [currentIP_late_array, currentIP_annulus_late]
		ENDFOR

		allIP_total_array = [[allIP_total_array], [currentIP_total_array]]
		allIP_late_array  = [[allIP_late_array], [currentIP_late_array]]

                IF outputIndex mod printEvery EQ 0 THEN PRINT, outputIndex
		outputIndex += 1
        ENDFOR

	neigh_all_grand_total  = []
	neigh_late_grand_total = []

	; ARRAY[ COL, ROW ]
	FOR n=0,n_annuli-1 DO BEGIN
		neigh_all_grand_total_annulus  = TOTAL(TRANSPOSE([allIP_total_array[n,*]]))
		neigh_late_grand_total_annulus = TOTAL(TRANSPOSE([allIP_late_array[n,*]]))
		neigh_all_grand_total 	= [neigh_all_grand_total, neigh_all_grand_total_annulus]
		neigh_late_grand_total	= [neigh_late_grand_total, neigh_late_grand_total_annulus]
	ENDFOR

	PRINT, 'all neighbor weighted totals: ', neigh_all_grand_total
	PRINT, 'late neighbor weighted totals: ', neigh_late_grand_total
END
