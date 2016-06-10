PRO get_neigh_data, minIPmstar
	dR = 0.1 ; Mpc
	dz_coeff = 2.0
	targPhi = -3.7
	stringPHI = strtrim(string(-1.*targPhi, format='(f20.1)'), 2)

	dataAll = mrdfits('~/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits',1)
	dataAll	= dataAll[WHERE(dataAll.targ_weight GE 1.)]
	dataIP	= dataAll[WHERE(dataAll.IP EQ 1)]

        fields = dataAll[uniq(dataAll.field, sort(dataAll.field))].field

; IP stellar mass cut
;	minIPmstar = 10.
	maxIPmstar = minIPmstar + 0.5
	dataIP	= dataIP[WHERE((dataIP.mstar GE minIPmstar) AND (dataIP.mstar LE maxIPmstar))]
	PRINT, 'IP: ', n_elements(dataIP)
; Need median SSFR as function of Mstar over [minIPmstar-0.5,maxIPmstar]
	mbin = 0.1
	mArray = minIPmstar-0.5+FINDGEN((maxIPmstar-minIPmstar+0.5)/mbin+1)*mbin
	PRINT, 'mArray: ', mArray+0.5*mbin

	medSSFRarray = []
	FOR i=0,n_elements(mArray)-2 DO BEGIN
		mlow	= mArray[i]
		mhigh	= mArray[i+1]

		subdata = dataAll[WHERE((dataAll.mstar GE mlow) AND (dataAll.mstar LT mhigh))]
		medSSFRarray = [medSSFRArray, MEDIAN(subdata.SFR-subdata.mstar)]
	ENDFOR

;	PLOT, mArray[0:-2]+0.5*mbin, medSSFRarray, PSYM=4, $
;		xrange=[minIPmstar-0.5,maxIPmstar], yrange=[MIN(medSSFRarray)-0.1,MAX(medSSFRarray)+0.1], xtitle='log(mstar)', ytitle='SSFR' 
;	OPLOT, dataIP.mstar, INTERPOL(medSSFRarray, mArray[0:-2]+0.5*mbin, dataIP.mstar), PSYM=3

	PLOT, HISTOGRAM(dataIP.SFR-dataIP.mstar, bin=0.1)

	PRINT, 'medSSFRarray: ', medSSFRarray
	PRINT, n_elements(mArray), n_elements(medSSFRarray)
	
; For each IP in Mstar range, find
; (i)	distance to near neighbor within 0.5 dex in Mstar and 2.0-sigma in z-space	
; (ii)	(SSFR of central) - median(SSFR @ central mass)
; (iii)	(SSFR of neighbor) - median(SSFR @ neighbor mass)
; (iv)	correlation coefficient between (ii) and (iii)

; Plot (iv) as function of (i)

	templateRow  = CREATE_STRUCT('dist', 0.0d, 'IP_SSFR', 0.0, 'IP_SSFRcorr', 0.0, 'IP_type', 0, $
			'neigh_SSFR', 0.0, 'neigh_SSFRcorr', 0.0, 'neigh_type', 0)
	neighData    = REPLICATE(templateRow, n_elements(dataIP))
;	neighData    = REPLICATE(templateRow, 100)

	FOR i=0,n_elements(dataIP)-1 DO BEGIN
;	FOR i=0,99 DO BEGIN
	  currentIP = dataIP[i]

	  currentIP_SSFR	= currentIP.SFR-currentIP.mstar
	  corrIP		= INTERPOL(medSSFRarray, mArray[0:-2]+0.5*mbin, currentIP.mstar)
	  IP_SSFRcorr		= currentIP_SSFR - corrIP

	  neigh = dataAll[ WHERE( (dataAll.field EQ currentIP.field) AND (dataAll.objname NE currentIP.objname) AND $
		(ABS(dataAll.zprimus - currentIP.zprimus) LE dz_coeff*0.005*(1.+currentIP.zprimus)) AND $
		((currentIP.mstar - dataAll.mstar) GT 0.0) AND ((currentIP.mstar - dataAll.mstar) LE 0.5) ) ]
	  neigh_dists	= [SQRT((neigh.xprop-currentIP.xprop)^2 + (neigh.yprop-currentIP.yprop)^2)]

	  nearNeigh	= neigh[WHERE(neigh_dists EQ MIN(neigh_dists))]

	  nearNeigh_SSFR 	= nearNeigh.SFR-nearNeigh.mstar
	  corrNeigh		= INTERPOL(medSSFRarray, mArray[0:-2]+0.5*mbin, nearNeigh.mstar)
	  nearNeigh_SSFRcorr	= nearNeigh_SSFR - corrNeigh	  

;	  PRINT, 'IP   ', ' SSFR: ', currentIP_SSFR, ' Mstar: ', currentIP.mstar, ' SSFR correction: ', corrIP, ' Corrected SSFR: ', IP_SSFRcorr
;	  PRINT, 'neigh', ' SSFR: ', nearNeigh_SSFR, ' Mstar: ', nearNeigh.mstar, ' SSFR correction: ', corrNeigh, ' Corrected SSFR: ', nearNeigh_SSFRcorr
;	  PRINT, ''

;	  PRINT, 'dist', MIN(neigh_dists), 'IP_SSFRcorr', IP_SSFRcorr, 'neigh_SSFRcorr', nearNeigh_SSFRcorr
	  newRow  = CREATE_STRUCT('dist', MIN(neigh_dists), 'IP_SSFR', currentIP_SSFR, 'IP_SSFRcorr', IP_SSFRcorr, 'IP_type', currentIP.SFQ, $
	  		'neigh_SSFR', nearNeigh_SSFR, 'neigh_SSFRcorr', nearNeigh_SSFRcorr, 'neigh_type', nearNeigh.SFQ)
;	  PRINT, newRow
	  neighData[i] = newRow

;	  IF (i mod 100 EQ 0) THEN PRINT, i
	ENDFOR
	
	IPmassRange = decimal(minIPmstar,1) + '_' + decimal(maxIPmstar,1)
	PRINT, IPmassRange
	MWRFITS, neighData, '~/results/match_IP_sample_rigorous/correlation_test/neighborData_' + IPmassRange + '.fits', /CREATE
END
