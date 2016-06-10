PRO temp
	data = MRDFITS('~/conformity/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits', 1)
	data_zlow	= data[WHERE(data.zprimus LE median(data.zprimus))]
	data_zhigh	= data[WHERE(data.zprimus GT median(data.zprimus))]
	data_mlow	= data[WHERE(data.mstar LE median(data.mstar))]
	data_mhigh	= data[WHERE(data.mstar GT median(data.mstar))]

	dataSF_zlow	= data_zlow[WHERE(data_zlow.SFQ EQ 1)]
	dataQ_zlow	= data_zlow[WHERE(data_zlow.SFQ EQ 0)]
	dataSF_zhigh	= data_zhigh[WHERE(data_zhigh.SFQ EQ 1)]
	dataQ_zhigh	= data_zhigh[WHERE(data_zhigh.SFQ EQ 0)]

	dataSF_mlow	= data_mlow[WHERE(data_mlow.SFQ EQ 1)]
	dataQ_mlow	= data_mlow[WHERE(data_mlow.SFQ EQ 0)]
	dataSF_mhigh	= data_mhigh[WHERE(data_mhigh.SFQ EQ 1)]
	dataQ_mhigh	= data_mhigh[WHERE(data_mhigh.SFQ EQ 0)]

	PRINT, MEDIAN(dataSF_zlow.zprimus), MEDIAN(dataQ_zlow.zprimus)
	PRINT, MEDIAN(dataSF_zhigh.zprimus), MEDIAN(dataQ_zhigh.zprimus)
	PRINT, MEDIAN(dataSF_mlow.mstar), MEDIAN(dataQ_mlow.mstar)
	PRINT, MEDIAN(dataSF_mhigh.mstar), MEDIAN(dataQ_mhigh.mstar)
END

