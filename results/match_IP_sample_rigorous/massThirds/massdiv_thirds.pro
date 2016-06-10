; T1 lowest masses, T3 highest masses

PRO massdiv_thirds
;	data	= mrdfits('~/results/variable_mass+' + string_dm + 'dex/IP_data/zerodSFQ_IP_dz2.0.fits',1)
	dataIP	= MRDFITS('~/results/match_IP_sampe_rigorous/matchedIPsampleFBF_M13highMassCut.fits', 1)

	; initialize
	dataIP_T1 = []
	dataIP_T2 = []
	dataIP_T3 = []

	lower_index = ceil(n_elements(dataIP)/3.)
	upper_index = 2*floor(n_elements(dataIP)/3.)

	orderedMasses = dataAbove[SORT(dataIP.mstar)].mstar
;	PRINT, orderedMasses

	lowerMass = orderedMasses(lower_index)
	upperMass = orderedMasses(upper_index)

	PRINT, 'Lower index: ', lower_index, ' Upper index: ', upper_index
	PRINT, 'Lower mass: ', lowerMass, ' Upper mass: ', upperMass

	dataIP_T1 = dataIP[where(dataIP.mstar LE lowerMass)]
	dataIP_T2 = dataIP[where( (dataIP.mstar GE lowerMass) AND (dataIP.mstar LE upperMass))]
	dataIP_T3 = dataIP[where(dataIP.mstar GE upperMass)]

	PRINT, n_elements(dataIP_T1)
	PRINT, n_elements(dataIP_T2)
	PRINT, n_elements(dataIP_T3)

	MWRFITS, dataIP_T1, '~/results/match_IP_sample_rigorous/mass_bins/IP_data/matchedIPsampleFBF_M13highMassCut_T1.fits', /CREATE
	MWRFITS, dataIP_T2, '~/results/match_IP_sample_rigorous/mass_bins/IP_data/matchedIPsampleFBF_M13highMassCut_T3.fits', /CREATE
	MWRFITS, dataIP_T3, '~/results/match_IP_sample_rigorous/mass_bins/IP_data/matchedIPsampleFBF_M13highMassCut_T3.fits', /CREATE
END
