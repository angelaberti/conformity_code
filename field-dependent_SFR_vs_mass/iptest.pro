; define 'velocity' difference with dz=0.005(1+z)
; input_dm is amount above Moustakas13 mass limit for SF IPs

PRO iptest, input_dm ;, zerodInputFile, dz_coeff, printEvery ; finds isolated primaries in a given field
	string_dm = strtrim(strcompress(string(input_dm, format='(f20.1)')), 1)

	zerodInputRaw = mrdfits('~/field-dependent_SFR_vs_mass/zerodSFQ_all_cart.fits', 1)
	dz_coeff = 2.
	printEvery = 100
	
	zerodInput = zerodInputRaw[where( (zerodInputRaw.zprimus GE 0.2) and (zerodInputRaw.zprimus LE 1.0) )]

	fields	 = zerodInput[uniq(zerodInput.field, sort(zerodInput.field))].field
	n_fields = n_elements(fields)

	templateRow = create_struct(zerodInput[0], 'IP', 0) ; default IP status is 0 (not IP)
	zerodOutput = replicate(templateRow, n_elements(zerodInput))

;	print, 'zerodInputRaw:', n_elements(zerodInputRaw)
;	print, 'zerodInput (0.2 <= z <= 1.0):', n_elements(zerodInput)
;	print, 'zerodOutput:', n_elements(zerodOutput)
	
	outputIndex = 0L

; make initial mass cut based on field, z, and star-forming status to get sample to test for IP status
;	IPcands    = []
;	nonIPcands = []
	IPcands    = mrdfits('~/field-dependent_SFR_vs_mass/allAboveMassCompLim_var' + string_dm + 'dex.fits', 1)
	nonIPcands = mrdfits('~/field-dependent_SFR_vs_mass/allBelowMassCompLim_var' + string_dm + 'dex.fits', 1)

;	FOR i=0,n_elements(zerodInput)-1 DO BEGIN
;		IF (zerodInput[i].SFQ eq 1) THEN (dm = input_dm) ELSE (dm = 0.0)
;		massLimit = mlimit(massStruct, fields, zerodInput[i]) + dm ; dm is tolerance in dex beyond Moustakas limit
;		IF (zerodInput[i].mstar GE massLimit) $
;		  THEN (IPcands	   = [IPcands,    zerodInput[i]]) $
;		  ELSE (nonIPcands = [nonIPcands, zerodInput[i]])
;		IF (i mod printEvery eq 0) and (zerodInput[i].mstar GE massLimit) THEN BEGIN
;			PRINT, zerodInput[i].field, $
;			       '  redshift:', zerodInput[i].zprimus, $
;			       '  stellar mass:', zerodInput[i].mstar, $
;			       '  SF status:', zerodInput[i].SFQ, $
;			       '  dm:', dm, $
;			       '  mass limit:', massLimit, $
;			       '  IP candidates:', n_elements(IPcands)
;		ENDIF	
;	END

;	mwrfits, IPcands, '~/field-dependent_SFR_vs_mass/allAboveMassCompLim_var' + string_dm + 'dex.fits', /create
;	mwrfits, nonIPcands, '~/field-dependent_SFR_vs_mass/allBelowMassCompLim_var' + string_dm + 'dex.fits', /create

	PRINT, 'IPcands:', n_elements(IPcands)
	PRINT, 'nonIPcands:', n_elements(nonIPcands)

; ADD IP STATUS OF NON-IP CANDIDATES TO OUTPUT ARRAY
	FOR i=0,(n_elements(nonIPcands)-1) DO BEGIN 
		currentGal 	= nonIPcands[i]
		newRow 		= create_struct(currentGal, 'IP', 0)
		zerodOutput[outputIndex] = newRow
		outputIndex += 1
	ENDFOR

; LOOP THROUGH FIELDS
	FOR f=0,(n_fields-1) DO BEGIN
		currentField		= fields[f]
;		PRINT, strtrim(currentField)
		currentFieldGals 	= zerodInput[ where(zerodInput.field eq currentField, /null) ]
		n_gal			= n_elements(currentFieldGals)
		currentFieldIPcands 	= IPcands[ where(IPcands.field eq currentField, /null) ]

; LOOP THROUGH ALL IP CANDIDATES IN FIELD
		FOR i=0,(n_elements(currentFieldIPcands)-1) DO BEGIN 
			currentIPcand 	= currentFieldIPcands[i]
			currentIPcandStatus = 1							; 1=True, 0=False

; TEST IP STATUS BY MAKING CUTS TO SET OF OTHER GALAXIES IN FIELD
		      ; RADIAL VELOCITY CUT
			testGals_dzcut		= currentFieldGals[ where( (currentIPcand.objname ne currentFieldGals.objname) AND $
			  (ABS(currentFieldGals.zprimus - currentIPcand.zprimus) LE dz_coeff*0.005*(1.+currentIPcand.zprimus)), /null ) ]
;			  (ABS( vcurrIPcand-(100.0*redh100())*(currentFieldGals.dz/redh100()/(1 + currentFieldGals.zprimus)) ) lt float(dv)), /null ) ]
			n_dzcut		= n_elements(testGals_dzcut)
	
		      ; PROJECTED SEPARATION CUT
			IF testGals_dzcut ne !null THEN BEGIN
				testGals_pscut	= testGals_dzcut[ where( $
				  ( (currentIPcand.xprop-testGals_dzcut.xprop)^2 $
				  + (currentIPcand.yprop-testGals_dzcut.yprop)^2 )^0.5 lt 0.5, /null ) ] ; 500 kpc (not kpc/h)
				n_pscut 	= n_elements(testGals_pscut)
			ENDIF ELSE BEGIN
				testGals_pscut	= !null
				n_pscut 	= 0
			ENDELSE

		      ; MASS CUT
			IF testGals_pscut ne !null THEN BEGIN
				testGals_mcut	= testGals_pscut[ where(testGals_pscut.mstar gt (currentIPcand.mstar-alog10(2.0)), /null) ]
				n_mcut		= n_elements(testGals_mcut)
			ENDIF ELSE BEGIN
				testGals_mcut	= !null
				n_mcut		= 0
			ENDELSE

			IF (n_mcut eq 0) THEN (currentIPcandStatus = 1) ELSE (currentIPcandStatus = 0)

			IF outputIndex mod printEvery eq 0 THEN BEGIN
				PRINT, ''
				PRINT, 'Index:' + strcompress(string(outputIndex))
				PRINT, 'Galaxies within dz =' + strtrim(strcompress(string(dz_coeff, format='(f20.1)'))) + '*0.005(1+z) =' + strcompress(string(dz_coeff*0.005*(1.+currentIPcand.zprimus))) + ':' + strcompress(string(n_dzcut))
				PRINT, '...and within 500 kpc:' + strcompress(string(n_pscut))
				PRINT, '...and above mass threshold of IP candidate:' + strcompress(string(n_mcut))
				PRINT, 'IP status:' + strcompress(string(currentIPcandStatus))
			ENDIF

			newRow = create_struct(currentIPcand, 'IP', currentIPcandStatus)
			zerodOutput[outputIndex] = newRow
			outputIndex += 1
		ENDFOR
	ENDFOR

	PRINT, ''
	PRINT, 'Total galaxies:   ' + strcompress(string(n_elements(zerodInput)))
	PRINT, 'IP candidates:    ' + strcompress(string(n_elements(IPcands)))
	PRINT, 'Total IP galaxies:' + strcompress(string(n_elements(zerodOutput[where(zerodOutput.IP eq 1)])))

;	IF outputIndex ne n_elements(zerodInput_zcut) THEN message, 'outputIndex mismatch'

	mwrfits, zerodOutput, '~/field-dependent_SFR_vs_mass/zerodSFQ_IP_dz' + strtrim(strcompress(string(dz_coeff, format='(f20.1)')),1) + '_' + string_dm + 'dex.fits', /create
END
