PRO num_IP

	data = mrdfits('~/results/variable_mass+1.0dex/IP_data/zerodSFQ_IP_dz2.0.fits',1)
	fields = data[uniq(data.field, sort(data.field))].field
	dataFD = mrdfits('~/field-dependent_SFR_vs_mass/zerodSFQ_IP_dz2.0_1.0dex.fits',1)
 
	FOR f=0,4 DO BEGIN
		dataf = data[where(data.field EQ fields[f])]
		dataFDf = dataFD[where(dataFD.field EQ fields[f])]

		PRINT, fields[f], n_elements(dataf[where(dataf.IP EQ 1)]), n_elements(dataFDf[where(dataFDf.IP EQ 1)])
	ENDFOR	
END
