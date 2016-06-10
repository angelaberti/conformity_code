; uses my original (non-histogram) method and delta_z <= dz_coeff*0.005*(1+z)

; normsig:	calculates normalized conformity signal (and error)
; JKE:		uses jackknife error method
; IPmatchFBF:	SF and Q IP populations have been selected to have same mass and redshift distribution WITHIN EACH FIELD
; PHI37:	IP matching procedure uses high mass cutoff threshold based on Moustakas13 stellar mass function results
; zbinsEqualThirds: late-type fractions are computed for three redshift ranges with equal # IP

PRO normsig_JKE_IPmatchFBFnoCosmos_PHI37_zbinsEqualThirds, outputFormat
	PHI = 3.7
	stringPHI = strtrim(string(PHI, format='(f20.1)'),2)

	dz_coeff = 2.0
	zmin = 0.2
	zmax = 1.0

	; read in data files
	IPdataAllz  = MRDFITS('~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)
	allDataAllz = MRDFITS('~/results/zerodSFQ_all_cart.fits', 1)

	Rmax 		= 15.
	dRproj		= 1.
	n_annuli 	= CEIL(Rmax/dRproj)
	printEvery	= 100

        IF dRproj EQ 0.25 THEN string_dR = '_dR250kpc'
        IF dRproj EQ 1.00 THEN string_dR = '_dR1Mpc'

; eliminate data with targ_weight < 1
        allDataAllz = allDataAllz[WHERE( (allDataAllz.targ_weight GE 1.) AND (allDataAllz.field NE 'cosmos    ') )]
        IPdataAllz  = IPdataAllz[WHERE( (IPdataAllz.targ_weight GE 1.) AND (IPdataAllz.field NE 'cosmos    ') )]

; divide IP population into 3 redshift bins with equal numbers of IPs
	orderedRedshifts = IPdataAllz[SORT(IPdataAllz.zprimus)].zprimus

	lower_index = CEIL(n_elements(IPdataAllz)/3.)
	upper_index = 2*FLOOR(n_elements(IPdataAllz)/3.)

	lower_z = orderedRedshifts[lower_index]
	upper_z = orderedRedshifts[upper_index]

	zArray = [0.2, lower_z, upper_z, 1.0]
	PRINT, stringPHI, ' ', zArray

	!p.multi=0
	!p.charsize=1.25

	Rproj_array = FLOAT(dRproj)*(findgen(n_annuli+1))	
;	PRINT, Rproj_array

	xmin = 0
;	xmax = n_annuli*dRproj
	xmax = 5
	xr = xmax-xmin

	ymin = -0.2
	ymax = 0.2
	yr = ymax-ymin

  FOR i=0,2 DO BEGIN
	z_low   = zArray[i]
;	z_low = 0.2
	z_high  = zArray[i+1]
;	z_high = 1.0

	zsuffix = 'T' + strtrim(i+1,1)
;	zsuffix = 'allz'
	zlabel = 'z = [' + decimal(z_low,2) + ', ' + decimal(z_high,2) + ']'
;	zlabel  = strtrim(string(z_low,format='(f20.2)'),1) + ' < z < ' + strtrim(string(z_high,format='(f20.2)'),1)

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsampleFBF/normsigplot_JKE_IPmatchFBFnoCosmos_PHI37_zbinsEqualThirds', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
		THICK=1
	ENDELSE

	dataIP  = IPdataAllz[WHERE( (IPdataAllz.zprimus GE z_low) AND (IPdataAllz.zprimus LE z_high) )]
	dataAll = allDataAllz[WHERE( (allDataAllz.zprimus GE z_low) AND (allDataAllz.zprimus LE z_high) )]
;	dataIP  = IPdataAllz
;	dataAll = allDataAllz

; ***** SF IP ***** 
	get_weights_JKE, dataIP, dataAll, dz_coeff, 1, allTotals_allRegions_IPSF, lateTotals_allRegions_IPSF, n_annuli, dRproj, printEvery

; ***** Q IP *****
	get_weights_JKE, dataIP, dataAll, dz_coeff, 0, allTotals_allRegions_IPQ, lateTotals_allRegions_IPQ, n_annuli, dRproj, printEvery

	get_normsig_with_error, allTotals_allRegions_IPSF, lateTotals_allRegions_IPSF, totalArray_IPSF, lateArray_IPSF, $
		allTotals_allRegions_IPQ, lateTotals_allRegions_IPQ, totalArray_IPQ, lateArray_IPQ, normsigArray, errorArray, n_annuli

	IPSF_totalArray		= totalArray_IPSF
	IPSF_lateArray		= lateArray_IPSF
	IPSF_latefracArray	= IPSF_lateArray/IPSF_totalArray

	IPQ_totalArray		= totalArray_IPQ
	IPQ_lateArray		= lateArray_IPQ
	IPQ_latefracArray	= IPQ_lateArray/IPQ_totalArray

	PLOT, Rproj_array+0.5*dRproj, normsigArray, xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Projected Radius (Mpc)', ytitle='Normalized Conformity Signal', $
	  title = 'Field-by-Field Matched SF & Q IP Samples ' + textoidl('(\Deltaz=2.0)'), /NODATA

	OPLOT, Rproj_array+0.5*dRproj, normsigArray;, LINESTYLE=0
	OPLOT, [xmin,xmax], [0,0], LINESTYLE=2
	ERRPLOT, Rproj_array+0.5*dRproj, normsigArray-errorArray, normsigArray+errorArray

	XYOUTS, xmin+0.95*xr, ymin+0.20*yr, 'SF IP: ' + bigint(n_elements(dataIP[WHERE(dataIP.SFQ EQ 1)])), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.15*yr, 'Q IP: ' + bigint(n_elements(dataIP[WHERE(dataIP.SFQ EQ 0)])), ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.10*yr, 'Binwidth: ' + strtrim(string(dRproj,format='(f20.2)'),1) + ' Mpc', ALIGNMENT=1.0
	XYOUTS, xmin+0.95*xr, ymin+0.05*yr, zlabel, ALIGNMENT=1.0

	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, 'Weighted by TARG_WEIGHT', ALIGNMENT=1.0

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'

	templateRow = create_struct('Rmin', 0.0, 'Rmax', 0.0, $
		'n_tot_IPSF', 0.0, 'n_late_IPSF', 0.0, $
		'n_tot_IPQ', 0.0, 'n_late_IPQ', 0.0, $
		'normsig', 0.0, 'normsig_errors', 0.0)
	outputStruct    = replicate(templateRow, n_annuli)

	FOR j=0,n_annuli-1 DO BEGIN
		newRow  = create_struct('Rmin', Rproj_array[j], 'Rmax', Rproj_array[j+1], $
			'n_tot_IPSF', totalArray_IPSF[j], 'n_late_IPSF', lateArray_IPSF[j], $
			'n_tot_IPQ', totalArray_IPQ[j], 'n_late_IPQ', lateArray_IPQ[j], $
			'normsig', normsigArray[j], 'normsig_errors', errorArray[j])
		outputStruct[j] = newRow
	ENDFOR

	MWRFITS, outputStruct, '~/results/match_IP_sample_rigorous/jackknife_error/noCosmos/zThirds/normsig_' + zsuffix + '_targ_weight_IPmatchFBFnoCosmos_PHI' + stringPHI + string_dR + '_JKE.fits', /CREATE
  ENDFOR
END


PRO get_weights_JKE, dataIP, dataAll, dz_coeff, IPSFstatus, allTotals_allRegions, lateTotals_allRegions, n_annuli, dRproj, printEvery
; DEFINE REGIONS AND THEIR IP SAMPLES
	R01_IP = dataIP[WHERE( (STRTRIM(dataIP.field,2) EQ 'cdfs') AND (dataIP.RA LT 52.71101) )]
	R02_IP = dataIP[WHERE( (STRTRIM(dataIP.field,2) EQ 'cdfs') AND (dataIP.RA GE 52.71101) AND (dataIP.RA LE 53.47923) )]
	R03_IP = dataIP[WHERE( (STRTRIM(dataIP.field,2) EQ 'cdfs') AND (dataIP.RA GT 53.47923) )]
	R04_IP = dataIP[WHERE( (STRTRIM(dataIP.field,2) EQ 'cfhtls_xmm') AND (dataIP.DEC LT -5.0) )]
	R05_IP = dataIP[WHERE( (STRTRIM(dataIP.field,2) EQ 'cfhtls_xmm') AND (dataIP.DEC GE -5.0) AND (dataIP.DEC LE -4.3470) )]
	R06_IP = dataIP[WHERE( (STRTRIM(dataIP.field,2) EQ 'cfhtls_xmm') AND (dataIP.DEC GT -4.3470) )]
;	R07_IP = dataIP[WHERE( (STRTRIM(dataIP.field,2) EQ 'cosmos') AND (dataIP.DEC LE 2.3545399) )]
;	R08_IP = dataIP[WHERE( (STRTRIM(dataIP.field,2) EQ 'cosmos') AND (dataIP.DEC GT 2.3545399) )]
	R09_IP = dataIP[WHERE( STRTRIM(dataIP.field,2) EQ 'es1' )]
	R10_IP = dataIP[WHERE( STRTRIM(dataIP.field,2) EQ 'xmm' )]

;	regions = LIST(R01_IP, R02_IP, R03_IP, R04_IP, R05_IP, R06_IP, R07_IP, R08_IP, R09_IP, R10_IP)
	regions = LIST(R01_IP, R02_IP, R03_IP, R04_IP, R05_IP, R06_IP, R09_IP, R10_IP)

	lateTotals_allRegions = []
	allTotals_allRegions  = []

  FOR r=0,n_elements(regions)-1 DO BEGIN
	regIPall = regions[r]
	regionIP = regIPall[WHERE(regIPall.SFQ EQ IPSFstatus)]

	IF (IPSFstatus EQ 1) THEN (IPtype = 'Star-forming') ELSE (IPtype = 'Quiescent')
	PRINT, IPtype + ' IP galaxies in region ' + strtrim(r,2) + ': ' + strtrim(n_elements(regionIP),1)

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
;	PRINT, region_totalArray

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


PRO get_normsig_with_error, allTotals_allRegions_IPSF, lateTotals_allRegions_IPSF, totalArray_IPSF, lateArray_IPSF, $
	allTotals_allRegions_IPQ, lateTotals_allRegions_IPQ, totalArray_IPQ, lateArray_IPQ, normsigArray, errorArray, n_annuli
	n_regions	 = n_elements(allTotals_allRegions_IPSF[0,*])
	allRegionIndices = INDGEN(n_regions)

; SF IP LATE-TYPE FRACTION WITH ENTIRE SAMPLE (ALL REGIONS)
	totalArray_IPSF 	= getAnnuliWeightTotalArray(allTotals_allRegions_IPSF)
	lateArray_IPSF		= getAnnuliWeightTotalArray(lateTotals_allRegions_IPSF)
	latefracArray_IPSF	= lateArray_IPSF/totalArray_IPSF
; Q IP LATE-TYPE FRACTION WITH ENTIRE SAMPLE (ALL REGIONS)
	totalArray_IPQ		= getAnnuliWeightTotalArray(allTotals_allRegions_IPQ)
	lateArray_IPQ		= getAnnuliWeightTotalArray(lateTotals_allRegions_IPQ)
	latefracArray_IPQ	= lateArray_IPQ/totalArray_IPQ
; NORM SIGNAL WITH ENTIRE SAMPLE (ALL REGIONS)
	normsigArray = (latefracArray_IPSF - latefracArray_IPQ)/((latefracArray_IPSF + latefracArray_IPQ)/2.)

; GET ARRAY OF LATE-TYPE FRACTIONS FOR ALL SUBREGIONS
	latefracArray_IPSF_subRegionArray = []
	latefracArray_IPQ_subRegionArray  = []
	normsigArray_subRegionArray	  = []
	FOR ex=0,n_regions-1 DO BEGIN
		PRINT, 'Excluded region: ', ex
		subRegionIndices = allRegionIndices[WHERE(allRegionIndices NE ex)]
		
		totalNeigh_IPSF_subRegionArray = []
		totalNeigh_IPQ_subRegionArray  = []
		lateNeigh_IPSF_subRegionArray  = []
		lateNeigh_IPQ_subRegionArray   = []
		FOREACH index,subRegionIndices DO BEGIN
			; CREATE ARRAY OF SUBREGION DATA FOR **ALL** ANNULUS NEIGHBORS
			totalNeigh_IPSF_subRegionArray = [[totalNeigh_IPSF_subRegionArray], [allTotals_allRegions_IPSF[*,index]]]
			totalNeigh_IPQ_subRegionArray  = [[totalNeigh_IPQ_subRegionArray], [allTotals_allRegions_IPQ[*,index]]]
			; CREATE ARRAY OF SUBREGION DATA FOR **LATE** ANNULUS NEIGHBORS
			lateNeigh_IPSF_subRegionArray  = [[lateNeigh_IPSF_subRegionArray], [lateTotals_allRegions_IPSF[*,index]]]
			lateNeigh_IPQ_subRegionArray   = [[lateNeigh_IPQ_subRegionArray], [lateTotals_allRegions_IPQ[*,index]]]
		ENDFOREACH

		; GET **ALL** ANNULUS NEIGHBOR WEIGHT TOTALS FOR GIVEN SUBREGION
		totalNeigh_IPSF_subRegion = getAnnuliWeightTotalArray(totalNeigh_IPSF_subRegionArray)
		totalNeigh_IPQ_subRegion  = getAnnuliWeightTotalArray(totalNeigh_IPQ_subRegionArray)
		; GET **LATE** ANNULUS NEIGHBOR WEIGHT TOTALS FOR GIVEN SUBREGION
		lateNeigh_IPSF_subRegion  = getAnnuliWeightTotalArray(lateNeigh_IPSF_subRegionArray)
		lateNeigh_IPQ_subRegion   = getAnnuliWeightTotalArray(lateNeigh_IPQ_subRegionArray)

		PRINT, totalNeigh_IPQ_subRegion
		PRINT, lateNeigh_IPQ_subRegion
		PRINT, ''

		; COMPUTE SUBREGION LATE-TYPE FRACTION AND ADD TO ARRAY OF ALL SUBREGION LATE-TYPE FRACTIONS
		latefracArray_IPSF_subRegion	= lateNeigh_IPSF_subRegion/totalNeigh_IPSF_subRegion
		latefracArray_IPQ_subRegion 	= lateNeigh_IPQ_subRegion/totalNeigh_IPQ_subRegion

		latefracArray_IPSF_subRegionArray = [[latefracArray_IPSF_subRegionArray], [latefracArray_IPSF_subRegion]]
		latefracArray_IPQ_subRegionArray  = [[latefracArray_IPQ_subRegionArray], [latefracArray_IPQ_subRegion]]

		; COMPUTE SUBREGION NORM SIGNAL AND ADD TO ARRAY OF ALL SUBREGION NORM SIGNALS
		normsigArray_subRegion		= (latefracArray_IPSF_subRegion - latefracArray_IPQ_subRegion)/((latefracArray_IPSF_subRegion + latefracArray_IPQ_subRegion)/2.)
		normsigArray_subRegionArray	= [[normsigArray_subRegionArray], [normsigArray_subRegion]]
	ENDFOR

; COMPUTE ERROR OF NORMALIZED SIGNAL (ALL SUBREGIONS)
	diffArray = []
	FOR r=0,n_regions-1 DO diffArray = [[diffArray], [normsigArray_subRegionArray[*,r] - normsigArray]]
	errorArray = []
	FOR c=0,n_annuli-1 DO errorArray = [errorArray, SQRT( ((n_regions-1.)/n_regions)*TOTAL( diffArray[c,*]^2 ) )]

	PRINT, '((n_regions-1)^0.5)*(STDDEV(normalized signals of all subregions)) for each annulus'
	FOR c=0,n_annuli-1 DO PRINT, SQRT(n_regions-1)*STDDEV(normsigArray_subRegionArray[c,*])

	PRINT, 'errorArray'
	PRINT, errorArray
END
