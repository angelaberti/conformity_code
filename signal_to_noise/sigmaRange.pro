PRO sigmaRange, dataSet
	int250kpc = [ [0,3], [4,7], [4,11], [12,19] ]
	int250kpc_to10 = [ [0,3], [4,7], [4,11], [12,19], [15,19], [20,23], [24,27], [28,31], [32,35], [36,39] ]
;	int500kpc = [ [0,1], [2,3], [2,5], [6,9] ]
	int1Mpc   = [ [0,0], [1,1], [1,2], [3,4] ]
	int1Mpc_to10 = [ [0,0], [1,1], [1,2], [3,4], [4,4], [5,5], [6,6], [7,7], [8,8], [9,9] ]

	fields = ['cdfs', 'cfhtls_xmm', 'cosmos', 'es1', 'xmm']

	intervals = int1Mpc

	datafiles = []
	ranges = []

; 3 equal IP z bins; results for all fields together to projected radius 10 Mpc; dR = 1 Mpc
  IF dataSet EQ 6 THEN BEGIN
	intervals = int1Mpc_to10
	datapath = '~/conformity/results/match_IP_sample_rigorous/'
	FOR i=1,3 DO BEGIN
		range = 'T' + strtrim(i,2)
;		range = 'M' + strtrim(i,2) + 'halfDex'
;		range = 'allz'
		datafiles = [datafiles, 'latefrac_' + range + '_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits']
;		datafiles = ['latefrac_' + range + '_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits']
		ranges = [ranges, range]
	ENDFOR
  ENDIF

; 3 equal IP redshift bins; results for individual fields to projected radius 10 Mpc; dR = 1 Mpc
  IF dataSet EQ 5 THEN BEGIN
	intervals = int1Mpc_to10
;	intervals = int250kpc_to10
	datapath = '~/conformity/results/match_IP_sample_rigorous/fields/'
;	datapath = '~/conformity/results/match_IP_sample_rigorous/fields/dR_250kpc/'
    FOR f=0,4 DO BEGIN
	FOR i=1,3 DO BEGIN
		range = 'T' + strtrim(i,2) + '_' + strtrim(fields[f], 2)
		datafiles = [datafiles, 'latefrac_' + range + '_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits']
;		datafiles = [datafiles, 'latefrac_' + range + '_targ_weight_IPmatchFBF_PHI3.7_dR250kpc_BSE.fits']
		ranges = [ranges, range]
	ENDFOR
    ENDFOR
  ENDIF

; 3 0.5 dex range mass bins; results for all fields together to projected radius 5 Mpc; dR = 1 Mpc
  IF dataSet EQ 4 THEN BEGIN
	datapath = '~/conformity/results/match_IP_sample_rigorous/massThirds/'
	FOR i=1,3 DO BEGIN
		range = 'M' + strtrim(i,2) + 'halfDex'
		datafiles = [datafiles, 'latefrac_' + range + '_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits']
		ranges = [ranges, range]
	ENDFOR
  ENDIF

; 3 equal IP mass bins; results for all fields together to projected radius 5 Mpc; dR = 1 Mpc
  IF dataSet EQ 3 THEN BEGIN
	datapath = '~/conformity/results/match_IP_sample_rigorous/massThirds/'
	FOR i=1,3 DO BEGIN
		range = 'M' + strtrim(i,2) + 'equalIP'
		datafiles = [datafiles, 'latefrac_' + range + '_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits']
		ranges = [ranges, range]
	ENDFOR
  ENDIF

; 3 equal mass range bins; results for all fields together to projected radius 5 Mpc; dR = 1 Mpc
  IF dataSet EQ 2 THEN BEGIN
	datapath = '~/conformity/results/match_IP_sample_rigorous/massThirds/'
	FOR i=1,3 DO BEGIN
		range = 'M' + strtrim(i,2) + 'equalRange'
		datafiles = [datafiles, 'latefrac_' + range + '_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits']
		ranges = [ranges, range]
	ENDFOR
  ENDIF

; all redshift range; results for all fields together to projected radius 5 Mpc; dR = 1 Mpc
  IF dataSet EQ 1 THEN BEGIN
	datapath = '~/conformity/results/match_IP_sample_rigorous/'
	FOR i=1,3 DO BEGIN
		range = 'T' + strtrim(i,2)
		datafiles = [datafiles, 'latefrac_' + range + '_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits']
		ranges = [ranges, range]
	ENDFOR
  ENDIF

; all redshift range; results for all fields together to projected radius 5 Mpc; dR = 0.25 Mpc
  IF dataSet EQ 0 THEN BEGIN
	intervals = int250kpc
	datapath = '~/conformity/results/match_IP_sample_rigorous/dR_250kpc/'
	FOR i=1,3 DO BEGIN
		range = 'T' + strtrim(i,2)
		datafiles = [datafiles, 'latefrac_' + range + '_targ_weight_IPmatchFBF_PHI3.7_dR250kpc_BSE.fits']
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
	FOR i=0,n_elements(datafiles)-1 DO PRINT, allResults[*,i] 
;  ENDIF

;	FOREACH range,zranges DO BEGIN
;		datafile = 'latefrac_' + range + '_targ_weight_matchedIPsample_BSE.fits'
;		datafile = 'latefrac_' + range + '_targ_weight_matchedIPsample_BSE_FBF.fits'

;	ENDFOREACH

;	PRINT, datapath
;	FOR i=0,n_elements(zranges)-1 DO PRINT, allResults[*,i]
;	ENDIF	
END

FUNCTION getSigmaRange, data, lower_index, upper_index
	n_late_IPSF = data.n_late_IPSF
	n_tot_IPSF  = data.n_tot_IPSF
	n_late_IPQ  = data.n_late_IPQ
	n_tot_IPQ   = data.n_tot_IPQ
 
;	frac_IPSF = n_late_IPSF/n_tot_IPSF
;	frac_IPQ  = n_late_IPQ/n_tot_IPQ

	dfrac_IPSF = data.errors_IPSF
	dfrac_IPQ  = data.errors_IPQ

;	frac_IPSFmean   = MEAN(frac_IPSF[lower_index:upper_index])
	; wmean = weighted mean
	wmean_frac_IPSF = TOTAL(n_late_IPSF[lower_index:upper_index])/TOTAL(n_tot_IPSF[lower_index:upper_index])

;	PRINT, ''
;	PRINT, data[lower_index:upper_index]
	
	dfrac_IPSFmean  = MEAN(dfrac_IPSF[lower_index:upper_index])/(float(n_elements(dfrac_IPSF[lower_index:upper_index])))^0.5
;	PRINT, 'dfrac_IPSFmean: ', dfrac_IPSFmean
;	frac_IPQmean    = MEAN(frac_IPQ[lower_index:upper_index])
	wmean_frac_IPQ  = TOTAL(n_late_IPQ[lower_index:upper_index])/TOTAL(n_tot_IPQ[lower_index:upper_index])

	dfrac_IPQmean   = MEAN(dfrac_IPQ[lower_index:upper_index])/(float(n_elements(dfrac_IPQ[lower_index:upper_index])))^0.5

	signalOverRange = (wmean_frac_IPSF - wmean_frac_IPQ)/((wmean_frac_IPSF + wmean_frac_IPQ)/2)
;	PRINT, signalOverRange	
	sigmaOverRange	= ABS(signalOverRange)/((dfrac_IPSFmean^2 + dfrac_IPQmean^2)^(0.5))
;	PRINT, sigmaOverRange
	dR = data[0].Rmax - data[0].Rmin
	R_low  = data[lower_index].Rmin
	R_high = data[upper_index].Rmax

	result = decimal(R_low, 2) + '  ' + decimal(R_high, 2) + '  ' + decimal(signalOverRange, 4) + '  ' + decimal(sigmaOverRange, 4)

;	RETURN, result
	RETURN, decimal(sigmaOverRange, 2)
END

FUNCTION decimal, input, places
	format = '(f20.' + strtrim(places, 1) + ')'
	RETURN, strtrim(string(input, format=format),1)
END
