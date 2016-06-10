; dR = 1 Mpc unless otherwise indicated

PRO make_ubertable
	fields	  = ['cdfs', 'cfhtls_xmm', 'cosmos', 'es1', 'xmm']

	int = [ [0,0], [1,3], [4,5] ]

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

OPENW, unit, '~/latex/tables/signal_uber.tex', /GET_LUN

;==========HEADER==========;

PRINTF, unit, '\setlength{\tabcolsep}{0.03in}'
PRINTF, unit, '%\begin{deluxetable*}{cccccccccccc}[!h]'
PRINTF, unit, '\begin{deluxetable*}{ccrrccrccrcc}'
PRINTF, unit, '\tabletypesize{large}'
PRINTF, unit, '\tablecaption{Conformity Signal\label{table:}}'
PRINTF, unit, '\tablewidth{0pt}'
PRINTF, unit, '\tablehead{'
PRINTF, unit, '\colhead{} & \colhead{} & \colhead{} & \colhead{} &'
PRINTF, unit, '\multicolumn{2}{c}{$0.0 < R < 0.5$~Mpc} & {} &'
PRINTF, unit, '\multicolumn{2}{c}{$0.5 < R < 3.0$~Mpc} & {} &
PRINTF, unit, '\multicolumn{2}{c}{$3.0 < R < 5.0$~Mpc} \\'
PRINTF, unit, '\cline{5-6}\cline{8-9}\cline{11-12} \\'
PRINTF, unit, '\colhead{SF IP} & \colhead{Q IP} &'
PRINTF, unit, '\multicolumn{1}{c}{$z$} &'
PRINTF, unit, '\multicolumn{1}{c}{$\log\,(\mstar/\msun)$} &'
PRINTF, unit, '\colhead{Signal ($\%$)} & \colhead{$\sigma$} & {} &'
PRINTF, unit, '\colhead{Signal ($\%$)} & \colhead{$\sigma$} & {} &'
PRINTF, unit, '\colhead{Signal ($\%$)} & \colhead{$\sigma$} \\'
PRINTF, unit, '\cline{1-12} \\'
PRINTF, unit, '\multicolumn{12}{c}{Full Sample}'
PRINTF, unit, '}'

PRINTF, unit, '\startdata'

	IPdata = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)

	data = MRDFITS(datapath + datafiles[0], 1)
	printLine, unit, data, int, IPdata, 0, '', 0

PRINTF, unit, '\cutinhead{Redshift Bins}'
	n_bins	= 2
	binType = 'z'
	data = MRDFITS(datapath + datafiles[1], 1)
	printLine, unit, data, int, IPdata, n_bins, binType, 1
	data = MRDFITS(datapath + datafiles[2], 1)
	printLine, unit, data, int, IPdata, n_bins, binType, 2

PRINTF, unit, '\cline{1-12} \\'

	n_bins	= 3
	data = MRDFITS(datapath + datafiles[3], 1)
	printLine, unit, data, int, IPdata, n_bins, binType, 1
	data = MRDFITS(datapath + datafiles[4], 1)
	printLine, unit, data, int, IPdata, n_bins, binType, 2
	data = MRDFITS(datapath + datafiles[5], 1)
	printLine, unit, data, int, IPdata, n_bins, binType, 3

PRINTF, unit, '\cutinhead{Stellar Mass Bins}'
	n_bins	= 2
	binType = 'mass'
	data = MRDFITS(datapath + datafiles[6], 1)
	printLine, unit, data, int, IPdata, n_bins, binType, 1
	data = MRDFITS(datapath + datafiles[7], 1)
	printLine, unit, data, int, IPdata, n_bins, binType, 2

PRINTF, unit, '\cline{1-12} \\'

	n_bins	= 3
	data = MRDFITS(datapath + datafiles[8], 1)
	printLine, unit, data, int, IPdata, n_bins, binType, 1
	data = MRDFITS(datapath + datafiles[9], 1)
	printLine, unit, data, int, IPdata, n_bins, binType, 2
	data = MRDFITS(datapath + datafiles[10], 1)
	printLine, unit, data, int, IPdata, n_bins, binType, 3

;==========FOOTER==========;

PRINTF, unit, '\enddata'
PRINTF, unit, '\end{deluxetable*}'

FREE_LUN, unit

END

PRO printLine, unit, data, int, IPdata, n_bins, binType, binNumber
	IPdataAllz = IPdata
	getIPsampleStats, IPdataAllz, n_bins, binType, binNumber, N_IPSF, N_IPQ, zmin, zmax, Mmin, Mmax

	PRINTF, unit, '$' + bigint(N_IPSF) + '$ &'
	PRINTF, unit, '$' + bigint(N_IPQ) + '$ &'
	PRINTF, unit, '[' + decimal(zmin,2) + ', ' + decimal(zmax,2) + '] &'
	PRINTF, unit, '[' + decimal(Mmin,2) + ', ' + decimal(Mmax,2) + '] &'
  FOR i=0,2 DO BEGIN
	IF (i LT n_elements(int[0,*])-1) THEN (lineEnd = ' & {} &') ELSE (lineEnd = ' \\')
	getSigmaRange_JKE, data, int[0,i], int[1,i], sig, error, sigma

	PRINTF, unit, '$' + decimal(sig,3) + '\pm' + decimal(error,3) + '$ & ' + decimal(sigma,1) + lineEnd
	PRINT, '$' + decimal(sig,3) + '\pm' + decimal(error,3) + '$ & ' + decimal(sigma,1) + lineEnd
  ENDFOR

END


PRO getSigmaRange_JKE, data, lower_index, upper_index, sig, error, sigma
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

	sig	= signalOverRange
	error	= normsig_errorOfRange
	sigma	= sigmaOverRange

	RETURN
END


PRO getIPsampleStats, IPdataAllz, n_bins, binType, binNumber, N_IPSF, N_IPQ, zmin, zmax, Mmin, Mmax

	IPdataAllz = IPdataAllz[WHERE(IPdataAllz.targ_weight GE 1.)]

; DEFAULT IS SINGLE BIN
	dataIP	= IPdataAllz
	zmin	= 0.2
	zmax	= 1.0
	Mmin	= MIN(dataIP.mstar)
	Mmax	= MAX(dataIP.mstar)

; 2 Z OR MASS BINS
  IF (n_bins EQ 2) THEN BEGIN
    IF (binType EQ 'z') THEN BEGIN	
	IF (binNumber EQ 1) THEN BEGIN
		zmax	= MEDIAN(dataIP.zprimus)
		dataIP	= dataIP[WHERE(dataIP.zprimus LE zmax)]
		Mmax	= MAX(dataIP.mstar) 
	ENDIF
	IF (binNumber EQ 2) THEN BEGIN
		zmin	= MEDIAN(dataIP.zprimus)
		dataIP	= dataIP[WHERE(dataIP.zprimus GE zmin)]
		Mmin	= MIN(dataIP.mstar) 
	ENDIF
    ENDIF
    IF (binType EQ 'mass') THEN BEGIN
	IF (binNumber EQ 1) THEN BEGIN
		Mmax	= MEDIAN(dataIP.mstar)
		dataIP	= dataIP[WHERE(dataIP.mstar LE Mmax)]
		zmax	= MAX(dataIP.zprimus)
	ENDIF
	IF (binNumber EQ 2) THEN BEGIN
		Mmin	= MEDIAN(dataIP.mstar)
		dataIP	= dataIP[WHERE(dataIP.mstar GE Mmin)]
		zmin	= MIN(dataIP.zprimus)
	ENDIF
    ENDIF
  ENDIF

; 3 Z OR MASS BINS
  IF (n_bins EQ 3) THEN BEGIN
    IF (binType EQ 'z') THEN BEGIN
	orderedRedshifts = dataIP[SORT(dataIP.zprimus)].zprimus

	lower_index = CEIL(n_elements(dataIP)/3.)
	upper_index = 2*FLOOR(n_elements(dataIP)/3.)

	divArray = [0.2, orderedRedshifts[lower_index], orderedRedshifts[upper_index], 1.0]
;	PRINT, divArray
		
	IF (binNumber EQ 1) THEN BEGIN
		zmin = divArray[0]
		zmax = divArray[1]
;		PRINT, zmin, zmax
	ENDIF
	IF (binNumber EQ 2) THEN BEGIN
		zmin = divArray[1]
		zmax = divArray[2]
	ENDIF
	IF (binNumber EQ 3) THEN BEGIN
		zmin = divArray[2]
		zmax = divArray[3]
	ENDIF

	dataIP	= dataIP[WHERE( (dataIP.zprimus GE zmin) AND (dataIP.zprimus LE zmax) )]
	Mmin	= MIN(dataIP.mstar)
	Mmax	= MAX(dataIP.mstar)
    ENDIF
    IF (binType EQ 'mass') THEN BEGIN
	orderedMasses = dataIP[SORT(dataIP.mstar)].mstar

	lower_index = CEIL(n_elements(dataIP)/3.)
	upper_index = 2*FLOOR(n_elements(dataIP)/3.)

	divArray = [Mmin, orderedMasses[lower_index], orderedMasses[upper_index], Mmax]
;	PRINT, divArray

	IF (binNumber EQ 1) THEN BEGIN
		Mmin = divArray[0]
		Mmax = divArray[1]
;		PRINT, Mmin, Mmax
	ENDIF
	IF (binNumber EQ 2) THEN BEGIN
		Mmin = divArray[1]
		Mmax = divArray[2]
	ENDIF
	IF (binNumber EQ 3) THEN BEGIN
		Mmin = divArray[2]
		Mmax = divArray[3]
	ENDIF

	dataIP	= dataIP[WHERE( (dataIP.mstar GE Mmin) AND (dataIP.mstar LT Mmax) )]
	zmin	= MIN(dataIP.zprimus)
	zmax	= MAX(dataIP.zprimus)
    ENDIF
  ENDIF

	dataIPSF = dataIP[WHERE(dataIP.SFQ EQ 1)]
	dataIPQ  = dataIP[WHERE(dataIP.SFQ EQ 0)]
	N_IPSF	= n_elements(UNIQ(dataIPSF.objname))
	N_IPQ 	= n_elements(dataIPQ)

	PRINT, n_bins, '  ', binType, binNumber, zmin, zmax, Mmin, Mmax, N_IPSF, N_IPQ

END


FUNCTION decimal, input, places
	format = '(f20.' + strtrim(places, 1) + ')'
	RETURN, strtrim(string(input, format=format),1)
END
