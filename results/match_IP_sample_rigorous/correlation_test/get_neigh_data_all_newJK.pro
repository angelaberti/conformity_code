PRO get_neigh_data_all_newJK;, minIPmstar
	minIPmstar = 10.7
	m = 2

	dR = 0.1 ; Mpc
	dz_coeff = 2.0
	targPhi = -3.7
	stringPHI = strtrim(string(-1.*targPhi, format='(f20.1)'), 2)

	dataAll = MRDFITS('~/conformity/results/conservative_mass_cutoff/allAboveMassCompLim-0.0.fits',1)
	dataAll = dataAll[WHERE(dataAll.targ_weight GE 1.)]

	dataIP 	= MRDFITS('~/conformity/results/match_IP_sample_rigorous/correlation_test/allAboveMassCompLim+0.5_IP.fits',1)
	dataIP  = dataIP[WHERE((dataIP.IP EQ 1) AND (dataIP.targ_weight GE 1.))]

	fields = dataAll[uniq(dataAll.field, sort(dataAll.field))].field

	es1Lim		= [0.3, 0.4, 0.5]
	cdfsLim		= [0.3, 0.4, 0.5]
	cosmosLim	= [0.5, 0.65, 0.8]
	xmmLim		= [0.5, 0.65, 0.8]
	cfhtls_xmmLim	= [0.65, 0.8, 1.0]

	maxIPmstar = minIPmstar + 0.3

	dataIP_allFields = []
  FOR f=0,n_elements(fields)-1 DO BEGIN
	IF (STRTRIM(fields[f],2) EQ 'es1') OR (STRTRIM(fields[f],2) EQ 'cdfs') THEN (zlims = es1Lim)
	IF (STRTRIM(fields[f],2) EQ 'cosmos') OR (STRTRIM(fields[f],2) EQ 'xmm') THEN (zlims = cosmosLim)
	IF (STRTRIM(fields[f],2) EQ 'cfhtls_xmm') THEN (zlims = cfhtls_xmmLim)

	dataIPfield = dataIP[WHERE((dataIP.mstar GE minIPmstar) AND (dataIP.mstar LE maxIPmstar) AND $
		(dataIP.field EQ fields[f]) AND (dataIP.zprimus LE zlims[m]))]

	PRINT, fields[f], zlims[m], n_elements(dataIPfield)

	dataIP_allFields = [dataIP_allFields, dataIPfield]
	PRINT, 'IP: ', n_elements(dataIP_allFields)
  ENDFOR
	dataIP = dataIP_allFields

	mbin = 0.1
	mArray = minIPmstar-0.5+FINDGEN((maxIPmstar-minIPmstar+0.5)/mbin+1)*mbin

	medSSFRarray = []
	FOR i=0,n_elements(mArray)-2 DO BEGIN
		mlow	= mArray[i]
		mhigh	= mArray[i+1]

		subdata = dataAll[WHERE((dataAll.mstar GE mlow) AND (dataAll.mstar LT mhigh))]
		medSSFRarray = [medSSFRArray, MEDIAN(subdata.SFR-subdata.mstar)]
	ENDFOR
;	PLOT, HISTOGRAM(dataIP.SFR-dataIP.mstar, bin=0.1)

;	PRINT, 'medSSFRarray: ', medSSFRarray
;	PRINT, n_elements(mArray), n_elements(medSSFRarray)

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

	templateRow  = CREATE_STRUCT('JK_region', 0, 'zIP', 0.0, 'IP_type', 0, 'IP_weight', 0.0, $
		'IP_mstar', 0.0, 'IP_SFR', 0.0, $
		'n_neighbors', 0L, 'neigh_all_weight_total', 0.0, 'n_neigh_late', 0L, 'neigh_late_weight_total', 0.0)
	neighData    = REPLICATE(templateRow, n_elements(dataIP))
;	neighData    = REPLICATE(templateRow, 1000)

  count = 0
  FOR r=0,n_elements(regions)-1 DO BEGIN
	regionIP = regions[r]
	PRINT, 'IP in region ', r, ': ', n_elements(regionIP) 

	FOR i=0,n_elements(regionIP)-1 DO BEGIN
;	FOR i=0,99 DO BEGIN
	  currentIP = regionIP[i]

	  neigh = dataAll[ WHERE( (dataAll.field EQ currentIP.field) AND (dataAll.objname NE currentIP.objname) AND $
		(ABS(dataAll.zprimus - currentIP.zprimus) LE dz_coeff*0.005*(1.+currentIP.zprimus)) AND $
		((currentIP.mstar - dataAll.mstar) GT 0.0) AND ((currentIP.mstar - dataAll.mstar) LE 0.5) ) ]
	  neigh_dists	= [SQRT((neigh.xprop-currentIP.xprop)^2 + (neigh.yprop-currentIP.yprop)^2)]

	  neighInRange	= neigh[WHERE((neigh_dists GE 0.3) AND (neigh_dists LE 4.), /NULL)]
	
	IF (neighInRange EQ !NULL) THEN BEGIN
	  newRow  = CREATE_STRUCT('JK_region', r, 'zIP', currentIP.zprimus, 'IP_type', currentIP.SFQ, $
		'IP_weight', currentIP.targ_weight, 'IP_mstar', currentIP.mstar, 'IP_SFR', currentIP.SFR, $
		'n_neighbors', 0L, 'neigh_all_weight_total', 0.0, $
		'n_neigh_late', 0L, 'neigh_late_weight_total', 0.0)
	ENDIF ELSE BEGIN
	  newRow  = CREATE_STRUCT('JK_region', r, 'zIP', currentIP.zprimus, 'IP_type', currentIP.SFQ, $
		'IP_weight', currentIP.targ_weight, 'IP_mstar', currentIP.mstar, 'IP_SFR', currentIP.SFR, $
		'n_neighbors', n_elements(neighInRange), 'neigh_all_weight_total', TOTAL(neighInRange.targ_weight), $
		'n_neigh_late', n_elements(neighInRange[WHERE(neighInRange.SFQ EQ 1, /NULL)]), $
		'neigh_late_weight_total', TOTAL(neighInRange[WHERE(neighInRange.SFQ EQ 1)].targ_weight))
	ENDELSE
;	  PRINT, count, newRow
	  neighData[count] = newRow

	  IF (count mod 100 EQ 0) THEN PRINT, count
	  count += 1
	ENDFOR

  ENDFOR	
	IPmassRange = decimal(minIPmstar,1) + '_' + decimal(maxIPmstar,1)
;	PRINT, IPmassRange
	MWRFITS, neighData, '~/conformity/results/match_IP_sample_rigorous/correlation_test/M13/neighborDataAll_' + IPmassRange + '.fits', /CREATE
END
