PRO test, lower_index, upper_index
;	lower_index=1
;	upper_index=3
	datafile = '~/results/match_IP_sample_rigorous/jackknife_error/normsig_allz_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE.fits'
	data = mrdfits(datafile, 1)

	n_late_IPSF = data.n_late_IPSF
	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPQ  = data.n_late_IPQ
	n_tot_IPQ   = data.n_tot_IPQ

	normsig 	= data.normsig
	normsig_errors 	= data.normsig_errors

; wmean = weighted mean
	wmean_frac_IPSF = TOTAL(n_late_IPSF[lower_index:upper_index])/TOTAL(n_tot_IPSF[lower_index:upper_index])
	wmean_frac_IPQ  = TOTAL(n_late_IPQ[lower_index:upper_index])/TOTAL(n_tot_IPQ[lower_index:upper_index])

; weighted mean conformity signal over a given range of projected radii
	signalOverRange = (wmean_frac_IPSF - wmean_frac_IPQ)/((wmean_frac_IPSF + wmean_frac_IPQ)/2)

; error of signalOverRange
	normsig_errorOfRange = MEAN(normsig_errors[lower_index:upper_index])/SQRT(FLOAT(n_elements(normsig_errors[lower_index:upper_index])))

	sigmaOverRange	= ABS(signalOverRange)/normsig_errorOfRange
	dR = data[0].Rmax - data[0].Rmin
	R_low  = data[lower_index].Rmin
	R_high = data[upper_index].Rmax

; || MIN RADIUS | MAX RADIUS | NORMALIZED SIGNAL | ERROR | SIGMA ||
	result = decimal(R_low, 2) + '  ' + decimal(R_high, 2) + '  ' + $
		 decimal(signalOverRange, 3) + '  ' + $
		 decimal(normsig_errorOfRange, 3) + '  ' + $
		 decimal(sigmaOverRange, 1)

; for table
	PRINT, result
; for figures
;	RETURN, decimal(sigmaOverRange, 2)
; for plotting
;	RETURN, sigmaOverRange
END
