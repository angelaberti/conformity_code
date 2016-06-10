PRO sigmaRange_BSE_singleField
	int250kpc = [ [0,3], [4,7], [4,11], [12,19] ]
	int500kpc = [ [0,1], [2,3], [2,5], [6,9] ]
	int1Mpc   = [ [0,0], [1,1], [1,2], [3,4] ]

	intervals = int1Mpc

	allResults = []
	datapath = '~/conformity/results/match_IP_sample_rigorous/fields/'

	fields = ['cdfs', 'cfhtls_xmm', 'cosmos', 'es1', 'xmm']

	FOR j=0,4 DO BEGIN
		data = MRDFITS(datapath + 'latefrac_allz_' + fields[j] + '_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits',1)
		fieldResults = fields[j]
		FOR i=0,3 DO BEGIN
			fieldResults = [fieldResults, getSigmaRange_BSE(data, intervals[0,i], intervals[1,i])]
		ENDFOR	
		allResults = [ [allResults], [fieldResults] ]
	ENDFOR

	FOR i=0,n_elements(fields)-1 DO PRINT, allResults[*,i] 
END

FUNCTION getSigmaRange_BSE, data, lower_index, upper_index
	n_late_IPSF = data.n_late_IPSF
	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPQ  = data.n_late_IPQ
	n_tot_IPQ   = data.n_tot_IPQ
 
;	frac_IPSF = n_late_IPSF/n_tot_IPSF
;	frac_IPQ  = n_late_IPQ/n_tot_IPQ

	dfrac_IPSF = data.errors_IPSF
	dfrac_IPQ  = data.errors_IPQ

	; wmean = weighted mean
	wmean_frac_IPSF = TOTAL(n_late_IPSF[lower_index:upper_index])/TOTAL(n_tot_IPSF[lower_index:upper_index])
	dfrac_IPSFmean  = MEAN(dfrac_IPSF[lower_index:upper_index])/(FLOAT(n_elements(dfrac_IPSF[lower_index:upper_index])))^0.5

	wmean_frac_IPQ  = TOTAL(n_late_IPQ[lower_index:upper_index])/TOTAL(n_tot_IPQ[lower_index:upper_index])
	dfrac_IPQmean   = MEAN(dfrac_IPQ[lower_index:upper_index])/(FLOAT(n_elements(dfrac_IPQ[lower_index:upper_index])))^0.5

	signalOverRange = (wmean_frac_IPSF - wmean_frac_IPQ)/((wmean_frac_IPSF + wmean_frac_IPQ)/2)
	sigmaOverRange	= ABS(signalOverRange)/((dfrac_IPSFmean^2 + dfrac_IPQmean^2)^(0.5))

	dR = data[0].Rmax - data[0].Rmin
	R_low  = data[lower_index].Rmin
	R_high = data[upper_index].Rmax

	result = decimal(R_low, 2) + '  ' + decimal(R_high, 2) + '  ' + decimal(signalOverRange, 4) + '  ' + decimal(sigmaOverRange, 4)

	RETURN, result
END

FUNCTION decimal, input, places
	format = '(f20.' + strtrim(places, 1) + ')'
	RETURN, strtrim(string(input, format=format),1)
END
