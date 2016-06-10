PRO sigmaRange_BSE_halo12, dataSet
; 0	0.0-0.5
; 1	0.5-1.0
; 2	1.0-2.0
; 3	2.0-3.0
; 4	3.0-4.0
; 5	4.0-5.0

; want 0.0-0.5, 0.5-3.0, 3.0-5.0

	intervals = [ [0,0], $
		      [1,3], $
		      [4,5] ]

	datapath = '~/results/match_IP_sample_rigorous/unmatchedExample/halo12test/'

	datafiles = [ $
		'latefrac_allz_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE_halo12.fits', $
		'latefrac_noMatching_allz_dR1Mpc_BSE.fits' $
	]

  FOREACH datafile,datafiles DO BEGIN
	data = MRDFITS(datapath + datafile, 1)

	PRINT, datafile
	PRINT, getSigmaRange_BSE(data, intervals[0,0], intervals[1,0])
	PRINT, getSigmaRange_BSE(data, intervals[0,1], intervals[1,1])
	PRINT, getSigmaRange_BSE(data, intervals[0,2], intervals[1,2])
  ENDFOREACH
END

FUNCTION getSigmaRange_BSE, data, lower_index, upper_index
	n_late_IPSF = data.n_late_IPSF
	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPQ  = data.n_late_IPQ
	n_tot_IPQ   = data.n_tot_IPQ
 
	dfrac_IPSF = data.errors_IPSF
	dfrac_IPQ  = data.errors_IPQ

	; wmean = weighted mean
	wmean_frac_IPSF = TOTAL(n_late_IPSF[lower_index:upper_index])/TOTAL(n_tot_IPSF[lower_index:upper_index])
	dfrac_IPSFmean  = MEAN(dfrac_IPSF[lower_index:upper_index])/(FLOAT(n_elements(dfrac_IPSF[lower_index:upper_index])))^0.5

	wmean_frac_IPQ  = TOTAL(n_late_IPQ[lower_index:upper_index])/TOTAL(n_tot_IPQ[lower_index:upper_index])
	dfrac_IPQmean   = MEAN(dfrac_IPQ[lower_index:upper_index])/(FLOAT(n_elements(dfrac_IPQ[lower_index:upper_index])))^0.5

	signalOverRange = (wmean_frac_IPSF - wmean_frac_IPQ)/((wmean_frac_IPSF + wmean_frac_IPQ)/2)
	errorOverRange	= SQRT(dfrac_IPSFmean^2 + dfrac_IPQmean^2)
;	sigmaOverRange	= ABS(signalOverRange)/((dfrac_IPSFmean^2 + dfrac_IPQmean^2)^(0.5))
	sigmaOverRange	= ABS(signalOverRange)/errorOverRange

	dR = data[0].Rmax - data[0].Rmin
	R_low  = data[lower_index].Rmin
	R_high = data[upper_index].Rmax

;	result =  decimal(R_low, 2) + '  ' + decimal(R_high, 2) + '  ' + decimal(signalOverRange, 3) + '  ' + decimal(sigmaOverRange, 1)
	result =  decimal(R_low, 2) + '  ' + decimal(R_high, 2) + '  ' $
		+ decimal(signalOverRange, 3) + '  ' + decimal(errorOverRange, 3) + '  ' $
		+ decimal(sigmaOverRange, 1)

	RETURN, result
END

FUNCTION decimal, input, places
	format = '(f20.' + strtrim(places, 1) + ')'
	RETURN, strtrim(string(input, format=format),1)
END
