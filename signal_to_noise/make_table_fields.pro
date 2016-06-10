; dR = 1 Mpc unless otherwise indicated

PRO make_table_fields
	fields = ['cdfs', 'cosmos', 'es1', 'cfhtls_xmm', 'xmm']
	stringFields = ['CDFS', 'COSMOS', 'ES1', 'XMM-CFHTLS', 'XMM-SXDS']

OPENW, unit, '~/latex/tables/signal_fields.tex', /GET_LUN

;==========HEADER==========;
; Field N_IPSF N_IPQ 1-halo_sigma

;PRINTF, unit, '\setlength{\tabcolsep}{0.03in}'
;PRINTF, unit, '%\begin{deluxetable}{cccc}[!h]'
PRINTF, unit, '\begin{deluxetable}{lrrr}'
PRINTF, unit, '\tabletypesize{large}'
PRINTF, unit, '\tablecaption{Signifcance of 1-halo signal ($\Rproj<1$~Mpc) by field computed with bootstrap errors\label{table:}}'
PRINTF, unit, '\tablewidth{0pt}'

PRINTF, unit, '\tablehead{'
PRINTF, unit, '\colhead{Field} & \colhead{$N_{\textrm{IP,SF}}$} & \colhead{$N_{\textrm{IP,Q}}$} & \colhead{$\sigma_{\textrm{BS}}$} \\'
PRINTF, unit, '}'

PRINTF, unit, '\startdata'

	dataJKE = MRDFITS('~/conformity/results/match_IP_sample_rigorous/jackknife_error/normsig_allz_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_JKE.fits',1)

	IPdataAll = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)

  FOR f=0,n_elements(fields)-1 DO BEGIN
	datafileField	= '~/conformity/results/match_IP_sample_rigorous/fields/latefrac_allz_' $
			+ fields[f] + '_targ_weight_IPmatchFBF_PHI3.7_dR1Mpc_BSE.fits'
	IPdata = IPdataAll[WHERE( (STRTRIM(IPdataAll.field,2) EQ fields[f]) AND (IPdataAll.targ_weight GE 1.) )]
	PRINT, f, n_elements(IPdata)

	dataBSE = MRDFITS(datafileField, 1)
	stringField = stringFields[f]
	printFieldLine, unit, dataBSE, IPdata, stringField
  ENDFOR
PRINTF, unit, '\cline{1-4} \\'
	IPSFall = IPdataAll[WHERE(IPdataAll.SFQ EQ 1)]
	N_IPSF	= n_elements(IPSFall[UNIQ(IPSFall.objname, SORT(IPSFall.objname))])
	N_IPQ	= n_elements(IPdataAll[WHERE(IPdataAll.SFQ EQ 0)])

	PRINTF, unit, 'Full Sample ($\sigma_{\textrm{JK}}$) &'	
	PRINTF, unit, '$' + bigint(N_IPSF) + '$ &'
	PRINTF, unit, '$' + bigint(N_IPQ) + '$ &'
	getSigmaRange_JKE, dataJKE, 0, 0, sig, error, sigmaJKE

	PRINTF, unit, '$' + decimal(sigmaJKE,1) + '$ \\'

;==========FOOTER==========;

PRINTF, unit, '\enddata'
PRINTF, unit, '\end{deluxetable}'

FREE_LUN, unit

END

PRO printFieldLine, unit, dataBSE, IPdata, stringField
	IPSFall = IPdata[WHERE(IPdata.SFQ EQ 1)]
	N_IPSF	= n_elements(IPSFall[UNIQ(IPSFall.objname, SORT(IPSFall.objname))])
	N_IPQ	= n_elements(IPdata[WHERE(IPdata.SFQ EQ 0, /NULL)])

	PRINTF, unit, stringField + ' &'	
	PRINTF, unit, '$' + bigint(N_IPSF) + '$ &'
	PRINTF, unit, '$' + bigint(N_IPQ) + '$ &'
;	getSigmaRange_JKE, data, int[0,i], int[1,i], sig, error, sigma
	getSigmaRange_BSE, dataBSE, 0, 0, sig, error, sigmaBSE

	PRINTF, unit, '$' + decimal(sigmaBSE,1) + '$ \\'
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


FUNCTION decimal, input, places
	format = '(f20.' + strtrim(places, 1) + ')'
	RETURN, strtrim(string(input, format=format),1)
END
