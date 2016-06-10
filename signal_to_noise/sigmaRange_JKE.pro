; dR = 1 Mpc unless otherwise indicated

PRO sigmaRange_JKE, dataSet
	int250kpc = [ [0,3], [4,7], [4,11], [12,19] ]
;	int250kpc_to10 = [ [0,3], [4,7], [4,11], [12,19], [15,19], [20,23], [24,27], [28,31], [32,35], [36,39] ]
	int1Mpc   = [ [0,0], [1,1], [1,2], [3,4] ]
	int1Mpc_to10 = [ [0,0], [1,1], [1,2], [3,4], [4,4], [5,5], [6,6], [7,7], [8,8], [9,9] ]

	fields	  = ['cdfs', 'cfhtls_xmm', 'cosmos', 'es1', 'xmm']
	intervals = int1Mpc ; default

	datafiles = []
	ranges 	  = []


; 3 equal IP redshift bins
; results for ALL FIELDS together to projected radius 5 Mpc
; dR = 1 Mpc
  IF dataSet EQ 0 THEN BEGIN
	datapath = '~/results/match_IP_sample_rigorous/jackknife_error/zThirds/'
	FOR i=1,3 DO BEGIN
		range = 'T' + strtrim(i,2)
		datafiles = [datafiles, 'normsig_' + range + '_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE.fits']
		ranges = [ranges, range]
	ENDFOR
  ENDIF

; 3 equal IP redshift bins
; results for ALL FIELDS together EXCEPT COSMOS to projected radius 5 Mpc
; dR = 1 Mpc
  IF dataSet EQ 1 THEN BEGIN
	datapath = '~/results/match_IP_sample_rigorous/jackknife_error/noCosmos/zThirds/'
	FOR i=1,3 DO BEGIN
		range = 'T' + strtrim(i,2)
		datafiles = [datafiles, 'normsig_' + range + '_targ_weight_IPmatchFBFnoCosmos_PHI3.7_dR1Mpc_JKE.fits']
		ranges = [ranges, range]
	ENDFOR
  ENDIF

; ALL redshift range
; results for ALL FIELDS together to projected radius 5 Mpc
; dR = 1 Mpc
  IF dataSet EQ 2 THEN BEGIN
	datapath = '~/results/match_IP_sample_rigorous/jackknife_error/'
	range = 'allz'
	datafiles = [datafiles, 'normsig_' + range + '_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE.fits']
	ranges = range 
  ENDIF

; ALL redshift range
; results for ALL FIELDS together EXCEPT COSMOS to projected radius 5 Mpc
; dR = 1 Mpc
  IF dataSet EQ 3 THEN BEGIN
	datapath = '~/results/match_IP_sample_rigorous/jackknife_error/noCosmos/'
	range = 'allz'
	datafiles = [datafiles, 'normsig_' + range + '_targ_weight_IPmatchFBFnoCosmos_PHI3.7_dR1Mpc_JKE.fits']
	ranges = range
  ENDIF

; ALL redshift range
; results for ALL FIELDS together to projected radius 5 Mpc
; dR = 0.25 Mpc
  IF dataSet EQ 4 THEN BEGIN
	intervals = int250kpc
	datapath = '~/results/match_IP_sample_rigorous/jackknife_error/'
	range = 'allz'
	datafiles = [datafiles, 'normsig_' + range + '_targ_weight_IPmatchFBF_PHI3.7_dR250kpc_JKE.fits']
	ranges = range 
  ENDIF

; 3 equal IP stellar mass bins
; results for ALL FIELDS together to projected radius 5 Mpc
; dR = 1 Mpc
  IF dataSet EQ 5 THEN BEGIN
	datapath = '~/results/match_IP_sample_rigorous/jackknife_error/massThirds/'
	FOR i=1,3 DO BEGIN
		range = 'M' + strtrim(i,2) + 'equalIP'
		datafiles = [datafiles, 'normsig_' + range + '_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE.fits']
		ranges = [ranges, range]
	ENDFOR
  ENDIF

; 3 equal IP stellar mass bins
; results for ALL FIELDS together EXCEPT COSMOS to projected radius 5 Mpc
; dR = 1 Mpc
  IF dataSet EQ 6 THEN BEGIN
	datapath = '~/results/match_IP_sample_rigorous/jackknife_error/noCosmos/massThirds/'
	FOR i=1,3 DO BEGIN
		range = 'M' + strtrim(i,2) + 'equalIP'
		datafiles = [datafiles, 'normsig_' + range + '_targ_weight_IPmatchFBFnoCosmos_PHI3.7_dR1Mpc_JKE.fits']
		ranges = [ranges, range]
	ENDFOR
  ENDIF


	allResults = []
	FOR f=0,n_elements(datafiles)-1 DO BEGIN
		data = MRDFITS(datapath + datafiles[f], 1)

		datafileResults = [ranges[f], '  ']
		FOR i=0,n_elements(intervals[0,*])-1 DO $
		  datafileResults = [datafileResults, getSigmaRange(data, intervals[0,i], intervals[1,i])]

		allResults = [ [allResults], [datafileResults] ]
	ENDFOR

	PRINT, datapath
;	FOR i=0,n_elements(datafiles)-1 DO PRINT, allResults[*,i] 
END

FUNCTION getSigmaRange_JKE, data, lower_index, upper_index
	data=mrdfits(data, 1)

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
