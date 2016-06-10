PRO temp, zbin
	dd = MRDFITS('matchedIPsampleFBF_PHI3.7.fits', 1)

        SFRperM_med     = MEDIAN(dd.SFR-dd.mstar)
        dd_lowHalf	= dd[WHERE((dd.SFR-dd.mstar) LE SFRperM_med)]
        dd_highHalf     = dd[WHERE((dd.SFR-dd.mstar) GE SFRperM_med)]
        SFRperM_low     = MEDIAN(dd_lowHalf.SFR-dd_lowHalf.mstar)
        SFRperM_high    = MEDIAN(dd_highHalf.SFR-dd_highHalf.mstar)

	data = dd[WHERE( ((dd.SFR-dd.mstar) LE SFRperM_low) OR ((dd.SFR-dd.mstar) GE SFRperM_high) )]
	PRINT, minmax(data.zprimus)
	PRINT, minmax(data.mstar)

	dataQall  = data[WHERE(data.SFQ EQ 0)]
	dataSFall = data[WHERE(data.SFQ EQ 1)]

;	zbin = 0.1
	zarray = 0.2+FINDGEN(CEIL(0.8/zbin)+1)*zbin

	zerodSFQ_cart = MRDFITS('~/conformity/results/zerodSFQ_all_cart.fits', 1)
	fields = zerodSFQ_cart[uniq(zerodSFQ_cart.field, sort(zerodSFQ_cart.field))].field

	ERASE
	!P.MULTI=[6,2,3]

	dataTrimmed = []
  FOR f=0,n_elements(fields)-1 DO BEGIN
	dataTrimmedField = []
	dataQ	= dataQall[WHERE(dataQall.field EQ fields[f])]
	dataSF	= dataSFall[WHERE(dataSFall.field EQ fields[f])]

	FOR i=0,n_elements(zarray)-2 DO BEGIN
		z_low 	= zarray[i]
		z_high	= zarray[i+1]

		Qbin	= dataQ[WHERE( (dataQ.zprimus GE z_low) AND (dataQ.zprimus LT z_high) )]
		SFbin	= dataSF[WHERE( (dataSF.zprimus GE z_low) AND (dataSF.zprimus LT z_high) )]
	
		seed = 1		
;		PRINT, fields[f], z_low, z_high, n_elements(SFbin), n_elements(Qbin)

		IF (n_elements(Qbin) GT n_elements(SFbin)) AND (SFbin NE !NULL) THEN BEGIN
			Qselect = Qbin[ROUND(n_elements(Qbin)*RANDOMU(seed,n_elements(SFbin)))]
			dataTrimmedField = [dataTrimmedField, SFbin, Qselect]
;			PRINT, n_elements(Qselect)
		ENDIF
		IF (n_elements(SFbin) GT n_elements(Qbin)) AND (Qbin NE !NULL) THEN BEGIN
			SFselect = SFbin[ROUND(n_elements(SFbin)*RANDOMU(seed,n_elements(Qbin)))]
			dataTrimmedField = [dataTrimmedField, Qbin, SFselect]
;			PRINT, n_elements(SFselect)
		ENDIF
	ENDFOR
	dataTrimmed = [dataTrimmed, dataTrimmedField]

	PLOT, zarray, HISTOGRAM(data.zprimus, bin=zbin), yrange=[0,150], xrange=[0.2,1.0], /NODATA
	OPLOT, zarray, HISTOGRAM(dataSF.zprimus, bin=zbin), color=cgColor('cyan'), LINESTYLE=2
	OPLOT, zarray, HISTOGRAM(dataQ.zprimus, bin=zbin), color=cgColor('orange'), LINESTYLE=2

	OPLOT, zarray, HISTOGRAM(dataTrimmedField[WHERE(dataTrimmedField.SFQ EQ 1)].zprimus, bin=zbin), color=cgColor('blue')
	OPLOT, zarray, HISTOGRAM(dataTrimmedField[WHERE(dataTrimmedField.SFQ EQ 0)].zprimus, bin=zbin), color=cgColor('red')

	PRINT, fields[f]
	PRINT, MEDIAN(dataTrimmedField[WHERE(dataTrimmedField.SFQ EQ 1)].zprimus)
	PRINT, MEDIAN(dataTrimmedField[WHERE(dataTrimmedField.SFQ EQ 0)].zprimus)
	PRINT, MEDIAN(dataTrimmedField[WHERE(dataTrimmedField.SFQ EQ 1)].mstar)
	PRINT, MEDIAN(dataTrimmedField[WHERE(dataTrimmedField.SFQ EQ 0)].mstar)
	PRINT, ''
  ENDFOR

	PRINT, n_elements(data), n_elements(dataTrimmed)
	PRINT, MEDIAN(dataTrimmed[WHERE(dataTrimmed.SFQ EQ 1)].zprimus)
	PRINT, MEDIAN(dataTrimmed[WHERE(dataTrimmed.SFQ EQ 0)].zprimus)
	PRINT, MEDIAN(dataTrimmed[WHERE(dataTrimmed.SFQ EQ 1)].mstar)
	PRINT, MEDIAN(dataTrimmed[WHERE(dataTrimmed.SFQ EQ 0)].mstar)

	MWRFITS, dataTrimmed, '~/conformity/results/match_IP_sample_rigorous/jackknife_error/SFR_quartiles/matchedIPsampleFBF.fits'
END
