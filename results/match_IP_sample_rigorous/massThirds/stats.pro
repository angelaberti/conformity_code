PRO stats
	zmin = 0.2
	zmax = 1.0

	; read in data files
        dataIPallz  = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)

	; eliminate data with targ_weight < 1
        dataIPallz  = dataIPallz[where(dataIPallz.targ_weight GE 1.)]

; EQUAL IP
	; divide IP population into 3 mass bins with equal numbers of IPs
;	lower_index = ceil(n_elements(dataIPallz)/3.)
;	upper_index = 2*floor(n_elements(dataIPallz)/3.)

;	orderedMasses = dataIPallz[SORT(dataIPallz.mstar)].mstar

;	lowerMass = orderedMasses(lower_index)
;	upperMass = orderedMasses(upper_index)

;	massArray = [min(dataIPallz.mstar), lowerMass, upperMass, max(dataIPallz.mstar)]


; EQUAL RANGE

;	minMass = MIN(dataIPallz.mstar)
;	maxMass = MAX(dataIPallz.mstar)
;	massRange = maxMass - minMass
;	lowerMass = minMass + massRange/3.
;	upperMass = maxMass - massRange/3.

;	massArray = [minMass, lowerMass, upperMass, maxMass]
;	PRINT, massArray

; HALF DEX

	massArray = [10., 10.5, 11., 11.5]


;	PRINT, 'Lower index: ', lower_index, ' Upper index: ', upper_index
;	PRINT, 'Lower mass: ', lowerMass, ' Upper mass: ', upperMass
;	PRINT, massArray

  FOR i=0,2 DO BEGIN
	mass_low   = massArray[i]
	mass_high  = massArray[i+1]

;	massSuffix = 'M' + strtrim(i+1,1) + 'equalIP'
;	masslabel  = textoidl('M_{*} = [') + strtrim(string(mass_low,format='(f20.2)'),1) + ', ' + strtrim(string(mass_high,format='(f20.2)'),1) + ']'

	dataIP  = dataIPallz[where( (dataIPallz.mstar GE mass_low) AND (dataIPallz.mstar LE mass_high) )]
	PRINT, n_elements(dataIP[where(dataIP.SFQ EQ 1)]), n_elements(dataIP[where(dataIP.SFQ EQ 0)])
  ENDFOR
END
