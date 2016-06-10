PRO normsig_SFRquarts_JKE_IPmatchFBF_PHI37_allz, outputFormat
	PHI = 3.7
	stringPHI = strtrim(string(PHI, format='(f20.1)'),2)

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

; IP samples are highest and lowest SFR/M* quartiles, with matched redshift and z distributions
;	IPdataAllz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
;	IPdataAllz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/SFR_quartiles/matchedIPsampleFBF.fits', 1)
	IPdata_low  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/SFR_quartiles/matchedIPsampleFBF_PHI3.7_lowQuart.fits', 1)
	IPdata_high = MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/SFR_quartiles/matchedIPsampleFBF_PHI3.7_highQuart.fits', 1)
	allDataAllz = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)

	Rmax 		= 15.
	dRproj		= 1.
	n_annuli 	= CEIL(Rmax/dRproj)
	printEvery	= 100

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

; eliminate data with targ_weight < 1
        allDataAllz	= allDataAllz[WHERE(allDataAllz.targ_weight GE 1.)]
	IPdata_high	= IPdata_high[WHERE(IPdata_high.targ_weight GE 1.)]
	IPdata_low	= IPdata_low[WHERE(IPdata_low.targ_weight GE 1.)]
;       dd 		= IPdataAllz[WHERE(IPdataAllz.targ_weight GE 1.)]
;       dd_unique	= dd[UNIQ(dd.objname, SORT(dd.objname))]
;	PRINT, 'All IP: ', n_elements(dd)
;	PRINT, 'Unique IP: ', n_elements(dd_unique)

; USE ONLY UNIQUE IPs TO DETERMINE QUARTILES
;	SFRperM_med     = MEDIAN(dd_unique.SFR-dd_unique.mstar)
;	dd_lowHalf	= dd_unique[WHERE((dd_unique.SFR-dd_unique.mstar) LE SFRperM_med)]
;	dd_highHalf     = dd_unique[WHERE((dd_unique.SFR-dd_unique.mstar) GE SFRperM_med)]

; USE ALL IPs TO DETERMINE QUARTILES
;	SFRperM_med     = MEDIAN(dd.SFR-dd.mstar)
;	dd_lowHalf	= dd[WHERE((dd.SFR-dd.mstar) LE SFRperM_med)]
;	dd_highHalf     = dd[WHERE((dd.SFR-dd.mstar) GE SFRperM_med)]

;	SFRperM_low     = MEDIAN(dd_lowHalf.SFR-dd_lowHalf.mstar)
;	SFRperM_high    = MEDIAN(dd_highHalf.SFR-dd_highHalf.mstar)

;	PRINT, SFRperM_low, SFRperM_med, SFRperM_high
; LOWEST & HIGHEST QUARTILES
;	IPdataAllz	= dd[WHERE( ((dd.SFR-dd.mstar) LE SFRperM_low) OR ((dd.SFR-dd.mstar) GE SFRperM_high) )]
; MIDDLE QUARTILES
;	IPdataAllz	= dd[WHERE( ((dd.SFR-dd.mstar) GT SFRperM_low) AND ((dd.SFR-dd.mstar) LT SFRperM_high) )]

;	PRINT, n_elements(IPdataAllz[WHERE(IPdataAllz.SFQ EQ 1)]), n_elements(IPdataAllz[WHERE(IPdataAllz.SFQ EQ 0)])

;	PLOT, IPdataAllz.mstar, IPdataAllz.SFR, psym=3

	!p.multi=0
	!p.charsize=1.25

	Rproj_array = FLOAT(dRproj)*(findgen(n_annuli+1))	
;	PRINT, 'Rproj (Mpc): ', Rproj_array

	xmin = 0
;	xmax = n_annuli*dRproj
	xmax = 5
	xr = xmax-xmin

	ymin = -0.2
	ymax = 0.2
	yr = ymax-ymin

	z_low = 0.2
	z_high = 1.0

	zsuffix = 'allz'
	zlabel  = strtrim(string(z_low,format='(f20.2)'),1) + ' < z < ' + strtrim(string(z_high,format='(f20.2)'),1)

	IF (string(outputFormat) eq 'ps') THEN BEGIN
;	        PS_OPEN, '~/figures/matchedIPsampleFBF/normsigplot_JKE_IPmatchFBF_PHI37_allz', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

;	dataIP  = IPdataAllz
	dataAll = allDataAllz

; ***** High Quartile IPs ***** 
	get_weights_JKE, IPdata_high, dataAll, dz_coeff, allTotals_allRegions_IPhigh, lateTotals_allRegions_IPhigh, n_annuli, dRproj, printEvery

; ***** Low Quartile IPs *****
	get_weights_JKE, IPdata_low, dataAll, dz_coeff, allTotals_allRegions_IPlow, lateTotals_allRegions_IPlow, n_annuli, dRproj, printEvery

	get_normsig_with_error, allTotals_allRegions_IPhigh, lateTotals_allRegions_IPhigh, totalArray_IPhigh, lateArray_IPhigh, $
		allTotals_allRegions_IPlow, lateTotals_allRegions_IPlow, totalArray_IPlow, lateArray_IPlow, normsigArray, errorArray, n_annuli

	IPhigh_totalArray	= totalArray_IPhigh
	IPhigh_lateArray	= lateArray_IPhigh
	IPhigh_latefracArray	= IPhigh_lateArray/IPhigh_totalArray

	IPlow_totalArray	= totalArray_IPlow
	IPlow_lateArray		= lateArray_IPlow
	IPlow_latefracArray	= IPlow_lateArray/IPlow_totalArray

	PLOT, Rproj_array+0.5*dRproj, normsigArray, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Normalized Conformity Signal', $
	  title = 'Field-by-Field Matched SF & Q IP Samples ' + textoidl('(\Deltaz=2.0)'), /NODATA

	OPLOT, Rproj_array+0.5*dRproj, normsigArray;, LINESTYLE=0
	OPLOT, [xmin,xmax], [0,0], LINESTYLE=2
	ERRPLOT, Rproj_array+0.5*dRproj, normsigArray-errorArray, normsigArray+errorArray

;	XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'SF IP: ' + bigint(n_elements(dataIP[WHERE(dataIP.SFQ EQ 1)])), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Q IP: ' + bigint(n_elements(dataIP[WHERE(dataIP.SFQ EQ 0)])), ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
;	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, zlabel, ALIGNMENT=1.0

;	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'

	templateRow = create_struct('Rmin', 0.0, 'Rmax', 0.0, $
		'n_tot_IPhigh', 0.0, 'n_late_IPhigh', 0.0, $
		'n_tot_IPlow', 0.0, 'n_late_IPlow', 0.0, $
		'normsig', 0.0, 'normsig_errors', 0.0)
	outputStruct    = replicate(templateRow, n_annuli)

	FOR j=0,n_annuli-1 DO BEGIN
		newRow  = create_struct('Rmin', Rproj_array[j], 'Rmax', Rproj_array[j+1], $
			'n_tot_IPhigh', totalArray_IPhigh[j], 'n_late_IPhigh', lateArray_IPhigh[j], $
			'n_tot_IPlow', totalArray_IPlow[j], 'n_late_IPlow', lateArray_IPlow[j], $
			'normsig', normsigArray[j], 'normsig_errors', errorArray[j])
		outputStruct[j] = newRow
	ENDFOR

;	MWRFITS, outputStruct, '~/conformity/results/match_IP_sample_rigorous/jackknife_error/SFR_quartiles/normsig_matchedIPsampleFBF_uniqIPlowuarts.fits', /CREATE
	MWRFITS, outputStruct, '~/conformity/results/match_IP_sample_rigorous/jackknife_error/SFR_quartiles/normsig_matchedIPsampleFBF_highLowQuarts.fits', /CREATE
END


PRO get_weights_JKE, dataIP, dataAll, dz_coeff, allTotals_allRegions, lateTotals_allRegions, n_annuli, dRproj, printEvery
; DEFINE REGIONS AND THEIR IP SAMPLES
	R01_IP = dataIP[WHERE( (STRTRIM(dataIP.field,2) EQ 'cdfs') AND (dataIP.RA LT 52.71101) )]
	R02_IP = dataIP[WHERE( (STRTRIM(dataIP.field,2) EQ 'cdfs') AND (dataIP.RA GE 52.71101) AND (dataIP.RA LE 53.47923) )]
	R03_IP = dataIP[WHERE( (STRTRIM(dataIP.field,2) EQ 'cdfs') AND (dataIP.RA GT 53.47923) )]
	R04_IP = dataIP[WHERE( (STRTRIM(dataIP.field,2) EQ 'cfhtls_xmm') AND (dataIP.DEC LT -5.0) )]
	R05_IP = dataIP[WHERE( (STRTRIM(dataIP.field,2) EQ 'cfhtls_xmm') AND (dataIP.DEC GE -5.0) AND (dataIP.DEC LE -4.3470) )]
	R06_IP = dataIP[WHERE( (STRTRIM(dataIP.field,2) EQ 'cfhtls_xmm') AND (dataIP.DEC GT -4.3470) )]
	R07_IP = dataIP[WHERE( (STRTRIM(dataIP.field,2) EQ 'cosmos') AND (dataIP.DEC LE 2.3545399) )]
	R08_IP = dataIP[WHERE( (STRTRIM(dataIP.field,2) EQ 'cosmos') AND (dataIP.DEC GT 2.3545399) )]
	R09_IP = dataIP[WHERE( STRTRIM(dataIP.field,2) EQ 'es1' )]
	R10_IP = dataIP[WHERE( STRTRIM(dataIP.field,2) EQ 'xmm' )]

	regions = LIST(R01_IP, R02_IP, R03_IP, R04_IP, R05_IP, R06_IP, R07_IP, R08_IP, R09_IP, R10_IP)
;	regions = LIST(R07_IP, R08_IP, R09_IP, R10_IP)

	lateTotals_allRegions = []
	allTotals_allRegions  = []

  FOR r=0,n_elements(regions)-1 DO BEGIN
	regIPall = regions[r]
	regionIP = regIPall

;	IF (IPquartile EQ 'high') THEN (IPtype = 'Highest SFR/M* quartile') ELSE (IPtype = 'Lowest SFR/M* quartile')
	PRINT, 'IP galaxies in region ' + strtrim(r,2) + ': ' + strtrim(n_elements(regionIP),1)

	region_totalArray = [] ; list of all (weighted) neighbors around all region IPs by annulus
	region_lateArray  = [] ; list of late-type (weighted) neighbors around all region IPs by annulus
	outputIndex = 0L

	FOR i=0,(n_elements(regionIP)-1) DO BEGIN
;	FOR i=0,19 DO BEGIN
		currentIP = regionIP[i]

	; keep galaxies within dz=dz+coeff*0.005*(1+z) of IP (and in same field)
		all_dz_neigh = dataAll[ WHERE( (dataAll.field EQ currentIP.field) AND (dataAll.objname NE currentIP.objname) AND $
			(ABS(dataAll.zprimus - currentIP.zprimus) LE dz_coeff*0.005*(1.+currentIP.zprimus)) ) ]
		all_dz_neigh_dists = [SQRT((all_dz_neigh.xprop-currentIP.xprop)^2 + (all_dz_neigh.yprop-currentIP.yprop)^2)]

		currentIP_totalArray = [] ; weight total of all neighbors by annulus for current IP
		currentIP_lateArray  = [] ; weight total of late-type neighbors by annulus for current IP

		; GET LATE FRACTION FOR EACH ANNULUS
		FOR n=0,n_annuli-1 DO BEGIN
			currentRmin = dRproj*FLOAT(n)
			currentRmax = dRproj*FLOAT(n+1)
			
			annulus_neigh = all_dz_neigh[WHERE((all_dz_neigh_dists GE currentRmin) AND (all_dz_neigh_dists LE currentRmax), /NULL)]
			IF (annulus_neigh NE !NULL) THEN $
				annulus_neighSF = annulus_neigh[WHERE(annulus_neigh.SFQ EQ 1, /NULL)] ELSE $
				annulus_neighSF = []
			IF (annulus_neigh NE !NULL) THEN $
				currentIP_annulusTotal = TOTAL(annulus_neigh.targ_weight) ELSE $
				currentIP_annulusTotal = 0
			IF (annulus_neighSF NE !NULL) THEN $
				currentIP_annulusLate = TOTAL(annulus_neighSF.targ_weight) ELSE $
				currentIP_annulusLate = 0
			
			currentIP_totalArray = [currentIP_totalArray, currentIP_annulusTotal]
			currentIP_lateArray  = [currentIP_lateArray, currentIP_annulusLate]
		ENDFOR

		region_totalArray = [[region_totalArray], [currentIP_totalArray]]
		region_lateArray  = [[region_lateArray], [currentIP_lateArray]]

;		IF outputIndex mod printEvery EQ 0 THEN PRINT, outputIndex
		outputIndex += 1
        ENDFOR
;	PRINT, 'region_totalArray: ', region_totalArray

	region_allTotal  = []
	region_lateTotal = []

	; ARRAY[ COL, ROW ]
	FOR n=0,n_annuli-1 DO BEGIN
		region_allTotal_annulus  = TOTAL(TRANSPOSE([region_totalArray[n,*]]))
		region_lateTotal_annulus = TOTAL(TRANSPOSE([region_lateArray[n,*]]))
		region_allTotal		 = [region_allTotal, region_allTotal_annulus]
		region_lateTotal	 = [region_lateTotal, region_lateTotal_annulus]
	ENDFOR

	allTotals_allRegions 	= [[allTotals_allRegions], [region_allTotal]]
	lateTotals_allRegions	= [[lateTotals_allRegions], [region_lateTotal]]
  ENDFOR
END


FUNCTION getAnnuliWeightTotalArray, regionData
; returns total weight by annulus for specified regions (regionData)
; regionData can be lists of all neighbor or late-type neighbor weights (for either SF or Q IP) by annulus, one list per region

	annuliWeightTotalArray = []

	n_annuli = n_elements(regionData[*,0])	

	FOR n=1,n_annuli DO BEGIN
		annulusWeightTotal = 0
		FOR r=0,n_elements(regionData[0,*])-1 DO BEGIN
			regionList = regionData[*,r]	; row r, all columns
			annulusWeightTotal = (annulusWeightTotal + regionList[n-1])
		ENDFOR
		annuliWeightTotalArray = [annuliWeightTotalArray, annulusWeightTotal]
	ENDFOR

	RETURN, annuliWeightTotalArray
END


PRO get_normsig_with_error, allTotals_allRegions_IPhigh, lateTotals_allRegions_IPhigh, totalArray_IPhigh, lateArray_IPhigh, $
	allTotals_allRegions_IPlow, lateTotals_allRegions_IPlow, totalArray_IPlow, lateArray_IPlow, normsigArray, errorArray, n_annuli
	n_regions	 = n_elements(allTotals_allRegions_IPhigh[0,*])
	allRegionIndices = INDGEN(n_regions)

; SF IP LATE-TYPE FRACTION WITH ENTIRE SAMPLE (ALL REGIONS)
	totalArray_IPhigh 	= getAnnuliWeightTotalArray(allTotals_allRegions_IPhigh)
	lateArray_IPhigh		= getAnnuliWeightTotalArray(lateTotals_allRegions_IPhigh)
	latefracArray_IPhigh	= lateArray_IPhigh/totalArray_IPhigh
; Q IP LATE-TYPE FRACTION WITH ENTIRE SAMPLE (ALL REGIONS)
	totalArray_IPlow		= getAnnuliWeightTotalArray(allTotals_allRegions_IPlow)
	lateArray_IPlow		= getAnnuliWeightTotalArray(lateTotals_allRegions_IPlow)
	latefracArray_IPlow	= lateArray_IPlow/totalArray_IPlow
; NORM SIGNAL WITH ENTIRE SAMPLE (ALL REGIONS)
	normsigArray = (latefracArray_IPhigh - latefracArray_IPlow)/((latefracArray_IPhigh + latefracArray_IPlow)/2.)

; GET ARRAY OF LATE-TYPE FRACTIONS FOR ALL SUBREGIONS
	latefracArray_IPhigh_subRegionArray = []
	latefracArray_IPlow_subRegionArray  = []
	normsigArray_subRegionArray	  = []
	FOR ex=0,n_regions-1 DO BEGIN
		PRINT, 'Excluded region: ', ex
		subRegionIndices = allRegionIndices[WHERE(allRegionIndices NE ex)]
		
		totalNeigh_IPhigh_subRegionArray = []
		totalNeigh_IPlow_subRegionArray  = []
		lateNeigh_IPhigh_subRegionArray  = []
		lateNeigh_IPlow_subRegionArray   = []

		FOREACH index,subRegionIndices DO BEGIN
			; CREATE ARRAY OF SUBREGION DATA FOR **ALL** ANNULUS NEIGHBORS
			totalNeigh_IPhigh_subRegionArray = [[totalNeigh_IPhigh_subRegionArray], [allTotals_allRegions_IPhigh[*,index]]]
			totalNeigh_IPlow_subRegionArray  = [[totalNeigh_IPlow_subRegionArray], [allTotals_allRegions_IPlow[*,index]]]
			; CREATE ARRAY OF SUBREGION DATA FOR **LATE** ANNULUS NEIGHBORS
			lateNeigh_IPhigh_subRegionArray  = [[lateNeigh_IPhigh_subRegionArray], [lateTotals_allRegions_IPhigh[*,index]]]
			lateNeigh_IPlow_subRegionArray   = [[lateNeigh_IPlow_subRegionArray], [lateTotals_allRegions_IPlow[*,index]]]
		ENDFOREACH

		; GET **ALL** ANNULUS NEIGHBOR WEIGHT TOTALS FOR GIVEN SUBREGION
		totalNeigh_IPhigh_subRegion = getAnnuliWeightTotalArray(totalNeigh_IPhigh_subRegionArray)
		totalNeigh_IPlow_subRegion  = getAnnuliWeightTotalArray(totalNeigh_IPlow_subRegionArray)

		; GET **LATE** ANNULUS NEIGHBOR WEIGHT TOTALS FOR GIVEN SUBREGION
		lateNeigh_IPhigh_subRegion  = getAnnuliWeightTotalArray(lateNeigh_IPhigh_subRegionArray)
		lateNeigh_IPlow_subRegion   = getAnnuliWeightTotalArray(lateNeigh_IPlow_subRegionArray)

		PRINT, 'totalNeigh_IPlow_subRegion: ', totalNeigh_IPlow_subRegion
		PRINT, 'lateNeigh_IPlow_subRegion: ', lateNeigh_IPlow_subRegion
		PRINT, 'totalNeigh_IPhigh_subRegion: ', totalNeigh_IPhigh_subRegion
		PRINT, 'lateNeigh_IPhigh_subRegion: ', lateNeigh_IPhigh_subRegion
;		PRINT, ''

		; COMPUTE SUBREGION LATE-TYPE FRACTION AND ADD TO ARRAY OF ALL SUBREGION LATE-TYPE FRACTIONS
		latefracArray_IPhigh_subRegion	= lateNeigh_IPhigh_subRegion/totalNeigh_IPhigh_subRegion
		latefracArray_IPlow_subRegion 	= lateNeigh_IPlow_subRegion/totalNeigh_IPlow_subRegion

		latefracArray_IPhigh_subRegionArray = [[latefracArray_IPhigh_subRegionArray], [latefracArray_IPhigh_subRegion]]
		latefracArray_IPlow_subRegionArray  = [[latefracArray_IPlow_subRegionArray], [latefracArray_IPlow_subRegion]]

		; COMPUTE SUBREGION NORM SIGNAL AND ADD TO ARRAY OF ALL SUBREGION NORM SIGNALS
		normsigArray_subRegion		= (latefracArray_IPhigh_subRegion - latefracArray_IPlow_subRegion)/((latefracArray_IPhigh_subRegion + latefracArray_IPlow_subRegion)/2.)
		normsigArray_subRegionArray	= [[normsigArray_subRegionArray], [normsigArray_subRegion]]
	ENDFOR

; COMPUTE ERROR OF NORMALIZED SIGNAL (ALL SUBREGIONS)
	diffArray = []
	FOR r=0,n_regions-1 DO diffArray = [[diffArray], [normsigArray_subRegionArray[*,r] - normsigArray]]
	errorArray = []
	FOR c=0,n_annuli-1 DO errorArray = [errorArray, SQRT( ((n_regions-1.)/n_regions)*TOTAL( diffArray[c,*]^2 ) )]

	PRINT, '((n_regions-1)^0.5)*(STDDEV(normalized signals of all subregions)) for each annulus'
	FOR c=0,n_annuli-1 DO PRINT, SQRT(n_regions-1)*STDDEV(normsigArray_subRegionArray[c,*])

	PRINT, 'errorArray: ', errorArray
END
