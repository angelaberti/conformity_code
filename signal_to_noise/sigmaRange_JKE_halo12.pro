; dR = 1 Mpc unless otherwise indicated

PRO sigmaRange_JKE_halo12, dataSet
	fields	  = ['cdfs', 'cfhtls_xmm', 'cosmos', 'es1', 'xmm']

	intervals = [ [0,0], [1,3], [4,5] ]

	datapath = '~/conformity/results/match_IP_sample_rigorous/jackknife_error/halo12/'

  datafiles = [ $
	'normsig_allz_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12.fits', $

	'zBins/normsig_H1_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12.fits', $
	'zBins/normsig_H2_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12.fits', $

	'zBins/normsig_T1_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12.fits', $
	'zBins/normsig_T2_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12.fits', $
	'zBins/normsig_T3_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12.fits', $

	'massBins/normsig_H1_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12.fits', $
	'massBins/normsig_H2_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12.fits', $

	'massBins/normsig_T1_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12.fits', $
	'massBins/normsig_T2_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12.fits', $
	'massBins/normsig_T3_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12.fits', $

	'noCosmos/normsig_allz_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12_noCosmos.fits', $

	'noCosmos/zBins/normsig_H1_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12_noCosmos.fits', $
	'noCosmos/zBins/normsig_H2_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12_noCosmos.fits', $

	'noCosmos/zBins/normsig_T1_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12_noCosmos.fits', $
	'noCosmos/zBins/normsig_T2_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12_noCosmos.fits', $
	'noCosmos/zBins/normsig_T3_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE_halo12_noCosmos.fits' $
  ]


  FOREACH datafile,datafiles DO BEGIN
	data = MRDFITS(datapath + datafile, 1)

	PRINT, datafile
	PRINT, getSigmaRange_JKE(data, intervals[0,0], intervals[1,0])
	PRINT, getSigmaRange_JKE(data, intervals[0,1], intervals[1,1])
	PRINT, getSigmaRange_JKE(data, intervals[0,2], intervals[1,2])
  ENDFOREACH

IDL> OPENW, Unit, 'test.tex', /GET_LUN
IDL> Print, Unit
         101
IDL> PRINTF, Unit, 'TEST'
IDL> FREE_LUN, Unit

END

FUNCTION getSigmaRange_JKE, data, lower_index, upper_index
	n_late_IPSF = data.n_late_IPSF
	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPQ  = data.n_late_IPQ
	n_tot_IPQ   = data.n_tot_IPQ

	normsig 	= data.normsig
	normsig_errors 	= data.normsig_errors

;	PRINT, 'normsig', normsig_errors

; wmean = weighted mean
	wmean_frac_IPSF = TOTAL(n_late_IPSF[lower_index:upper_index])/TOTAL(n_tot_IPSF[lower_index:upper_index])
	wmean_frac_IPQ  = TOTAL(n_late_IPQ[lower_index:upper_index])/TOTAL(n_tot_IPQ[lower_index:upper_index])

; weighted mean conformity signal over a given range of projected radii
	signalOverRange = (wmean_frac_IPSF - wmean_frac_IPQ)/((wmean_frac_IPSF + wmean_frac_IPQ)/2)
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
