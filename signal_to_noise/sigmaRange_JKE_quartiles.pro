; dR = 1 Mpc unless otherwise indicated

PRO sigmaRange_JKE_quartiles
	int1Mpc   = [ [0,0], [1,2], [3,4] ]

	fields	  = ['cdfs', 'cfhtls_xmm', 'cosmos', 'es1', 'xmm']
	intervals = int1Mpc ; default

	datapath = '~/conformity/results/match_IP_sample_rigorous/jackknife_error/'
	datafiles = [$;datapath + 'normsig_allz_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE.fits', $
		datapath + 'SFR_quartiles/normsig_matchedIPsampleFBF_highLowQuarts.fits']
;		datapath + 'SFR_quartiles/normsig_matchedIPsampleFBF_allIPlowuarts.fits', $
;		datapath + 'SFR_quartiles/normsig_matchedIPsampleFBF_allIPlowuarts_inner.fits']
;	ranges 	  = []

	allResults = []
	FOR f=0,n_elements(datafiles)-1 DO BEGIN
		data = MRDFITS(datafiles[f], 1)

		datafileResults = []
		FOR i=0,n_elements(intervals[0,*])-1 DO $
		  datafileResults = [datafileResults, getSigmaRange_JKE(data, intervals[0,i], intervals[1,i])]

		allResults = [ [allResults], [datafileResults] ]
	ENDFOR

	PRINT, datapath
	FOR i=0,n_elements(datafiles)-1 DO PRINT, datafiles[i], allResults[*,i] 
END

FUNCTION getSigmaRange_JKE, data, lower_index, upper_index
	n_late_IPhigh = data.n_late_IPhigh
	n_tot_IPhigh  = data.n_tot_IPhigh
	n_late_IPlow  = data.n_late_IPlow
	n_tot_IPlow   = data.n_tot_IPlow

	normsig 	= data.normsig
	normsig_errors 	= data.normsig_errors

;	PRINT, 'normsig', normsig_errors

; wmean = weighted mean
	wmean_frac_IPhigh = TOTAL(n_late_IPhigh[lower_index:upper_index])/TOTAL(n_tot_IPhigh[lower_index:upper_index])
	wmean_frac_IPlow  = TOTAL(n_late_IPlow[lower_index:upper_index])/TOTAL(n_tot_IPlow[lower_index:upper_index])

; weighted mean conformity signal over a given range of projected radii
	signalOverRange = (wmean_frac_IPhigh - wmean_frac_IPlow)/((wmean_frac_IPhigh + wmean_frac_IPlow)/2)
;	PRINT, signalOverRange

; error of signalOverRange
	normsig_errorOfRange = MEAN(normsig_errors[lower_index:upper_index])/SQRT(FLOAT(n_elements(normsig_errors[lower_index:upper_index])))
;	PRINT, 'normsig_errorOfRange: ', normsig_errorOfRange

	sigmaOverRange	= ABS(signalOverRange)/normsig_errorOfRange
;	PRINT, sigmaOverRange
	dR = data[0].Rmax - data[0].Rmin
	R_low  = data[lower_index].Rmin
	R_high = data[upper_index].Rmax

; || MIN RADIUS | MAX RADIUS | NORMALIZED SIGNAL | ERROR | SIGMA ||
	result = decimal(R_low, 2) + '  ' + decimal(R_high, 2) + '  ' + $
		 decimal(signalOverRange, 3) + '  ' + $
		 decimal(normsig_errorOfRange, 3) + '  ' + $
		 decimal(sigmaOverRange, 1)

; for table
	RETURN, result
; for figures
;	RETURN, decimal(sigmaOverRange, 2)
; for plotting
;	RETURN, sigmaOverRange
END

FUNCTION decimal, input, places
	format = '(f20.' + strtrim(places, 1) + ')'
	RETURN, strtrim(string(input, format=format),1)
END
