PRO sigmaRange_BSE, dataSet
	int250kpc = [ [0,3], [4,7], [4,11], [12,19] ]
	int500kpc = [ [0,1], [2,3], [2,5], [6,9] ]
	int1Mpc   = [ [0,0], [1,1], [1,2], [3,4] ]

	; equal dz = 0.2; 4 z bins; dR = 0.25 Mpc
	IF dataSet EQ 0 THEN BEGIN
		datapath = '~/conformity/results/match_IP_sample/'
		zranges = ['0.2_1.0']
		FOR i=0,3 DO zranges = [zranges, strtrim(string(0.2+0.2*i, format='(f20.1)'),1) + '_' + strtrim(string(0.2+0.2*(i+1), format='(f20.1)'),1)]
		; RADIAL RANGES OVER WHICH TO COMPUTE SIGMA
		intervals = int250kpc
	ENDIF

	; equal # IP per redshift bin; 4 z bins; dR = 0.25 Mpc
	IF dataSet EQ 1 THEN BEGIN
		datapath = '~/conformity/results/match_IP_sample/'
		zranges = []
		FOR i=1,4 DO zranges = [zranges, 'Q' + strtrim(i,1)]
		; RADIAL RANGES OVER WHICH TO COMPUTE SIGMA
		intervals = int250kpc
	ENDIF

;	; equal # IP per redshift bin; 3 z bins; dR = 0.25 Mpc
;	IF dataSet EQ 2 THEN BEGIN
;		datapath = '~/conformity/results/match_IP_sample/'
;		zranges = []
;		FOR i=1,3 DO zranges = [zranges, 'T' + strtrim(i,1)]
;		; RADIAL RANGES OVER WHICH TO COMPUTE SIGMA
;		intervals = int250kpc
;	ENDIF

	; equal # IP per redshift bin; 2 z bins; dR = 0.25 Mpc
	IF dataSet EQ 2 THEN BEGIN
		datapath = '~/conformity/results/match_IP_sample/'
		zranges = []
		FOR i=1,2 DO zranges = [zranges, 'H' + strtrim(i,1)]
		; RADIAL RANGES OVER WHICH TO COMPUTE SIGMA
		intervals = int250kpc
	ENDIF

	; equal # IP per redshift bin; 4 z bins; dR = 0.50 Mpc
	IF dataSet EQ 3 THEN BEGIN
		datapath = '~/conformity/results/match_IP_sample/dR_500kpc/'
		zranges = ['0.2_1.0']
;		zranges = []
		FOR i=1,4 DO zranges = [zranges, 'Q' + strtrim(i,1)]
		; RADIAL RANGES OVER WHICH TO COMPUTE SIGMA
		intervals = int500kpc
	ENDIF

	; equal # IP per redshift bin; 3 z bins; dR = 0.50 Mpc
	IF dataSet EQ 4 THEN BEGIN
		datapath = '~/conformity/results/match_IP_sample/dR_500kpc/'
		zranges = []
		FOR i=1,3 DO zranges = [zranges, 'T' + strtrim(i,1)]
		; RADIAL RANGES OVER WHICH TO COMPUTE SIGMA
		intervals = int500kpc
	ENDIF

	; all z range; dR = 1 Mpc
	IF dataSet EQ 5 THEN BEGIN
		datapath = '~/conformity/results/match_IP_sample/dR_1Mpc/'
		zranges = ['0.2_1.0']
;		FOR i=1,3 DO zranges = [zranges, 'T' + strtrim(i,1)]
		; RADIAL RANGES OVER WHICH TO COMPUTE SIGMA
		intervals = int1Mpc
	ENDIF

	allResults = []

	; all z range; dR = 0.25 Mpc
	IF dataSet EQ 6 THEN BEGIN
		datapath = '~/conformity/results/match_IP_sample/byField/'
		zrange = '0.20_1.00'
		; RADIAL RANGES OVER WHICH TO COMPUTE SIGMA
		intervals = int250kpc
		fields = ['cdfs', 'cfhtls_xmm', 'cosmos', 'es1', 'xmm']

		FOR j=0,4 DO BEGIN
			data = MRDFITS('~/conformity/results/match_IP_sample/byField/latefrac_' + zrange + '_targ_weight_matchedIPsample_' + fields[j] + '_BSE.fits',1)
			fieldResults = fields[j]
			FOR i=0,3 DO BEGIN
				fieldResults = [fieldResults, getSigmaRange_BSE(data, intervals[0,i], intervals[1,i])]
			ENDFOR	
			allResults = [ [allResults], [fieldResults] ]
		ENDFOR

		FOR i=0,n_elements(fields)-1 DO PRINT, allResults[*,i] 
	ENDIF

	; highest redshift third; dR = 1 Mpc
	IF dataSet EQ 7 THEN BEGIN
;		datapath = '~/conformity/results/match_IP_sample/byField/'
		datapath = '~/conformity/results/match_IP_sample/matchedFBF/'
		zrange = '0.70_1.00'
		; RADIAL RANGES OVER WHICH TO COMPUTE SIGMA
		intervals = int1Mpc
		fields = ['cdfs', 'cfhtls_xmm', 'cosmos', 'es1', 'xmm']

		FOR j=0,4 DO BEGIN
;			data = MRDFITS('~/conformity/results/match_IP_sample/byField/latefrac_' + zrange + '_targ_weight_matchedIPsample_' + fields[j] + '_dR1Mpc_BSE.fits',1)
			data = MRDFITS('~/conformity/results/match_IP_sample/matchedFBF/latefrac_' + zrange + '_targ_weight_matchedIPsample_' + fields[j] + '_dR1Mpc_BSE_FBF.fits',1)
			fieldResults = fields[j]
			FOR i=0,3 DO BEGIN
				fieldResults = [fieldResults, getSigmaRange_BSE(data, intervals[0,i], intervals[1,i])]
			ENDFOR	
			allResults = [ [allResults], [fieldResults] ]
		ENDFOR

		FOR i=0,n_elements(fields)-1 DO PRINT, allResults[*,i] 
	ENDIF

	; all z range; 3 mass bins; dR = 1 Mpc
	IF dataSet EQ 8 THEN BEGIN
		datapath = '~/conformity/results/match_IP_sample/massBins/'
		; RADIAL RANGES OVER WHICH TO COMPUTE SIGMA
		intervals = int1Mpc

		FOR j=0,2 DO BEGIN
			massbin = 'M' + strtrim(j+1, 1)
			data = MRDFITS('~/conformity/results/match_IP_sample/massBins/latefrac_' + massbin + '_targ_weight_matchedIPsample_dR1Mpc_BSE.fits',1)
			massResults = [massbin]
			FOR i=0,3 DO BEGIN
				massResults = [massResults, getSigmaRange_BSE(data, intervals[0,i], intervals[1,i])]
			ENDFOR	
			allResults = [ [allResults], [massResults] ]
		ENDFOR

		FOR i=0,2 DO PRINT, allResults[*,i] 
	ENDIF

;latefrac_M3_targ_weight_matchedIPsample_dR1Mpc_BSE.fits

      IF dataSet LE 5 THEN BEGIN	
	FOREACH range,zranges DO BEGIN
		datafile = 'latefrac_' + range + '_targ_weight_matchedIPsample_BSE.fits'
;		datafile = 'latefrac_' + range + '_targ_weight_matchedIPsample_BSE_FBF.fits'
		data = MRDFITS(datapath + datafile, 1)

		zrangeResults = [range, '  ']
		FOR i=0,n_elements(intervals[0,*])-1 DO $
		  zrangeResults = [zrangeResults, getSigmaRange_BSE(data, intervals[0,i], intervals[1,i])]

		allResults = [ [allResults], [zrangeResults] ]
	ENDFOREACH

	PRINT, datapath
	FOR i=0,n_elements(zranges)-1 DO PRINT, allResults[*,i]
      ENDIF	
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
