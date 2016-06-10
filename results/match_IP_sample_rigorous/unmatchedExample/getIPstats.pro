PRO getIPstats
	path = '~/conformity/results/'

	IPdatafiles = [ $
	path+'default_parameters/IP_data/zerodSFQ_IP_dz2.0_dm0.5.fits', $
	path+'conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits', $
	path+'single_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_singleMass.fits', $
	path+'conservative+0.2dex/IP_data/zerodSFQ_IP_dz2.0_dm0.2.fits', $
	path+'conservative+0.3dex/IP_data/zerodSFQ_IP_dz2.0_dm0.3.fits', $
	path+'variable_mass+0.5dex/IP_data/zerodSFQ_IP_dz2.0.fits', $
	path+'variable_mass+0.8dex/IP_data/zerodSFQ_IP_dz2.0.fits', $
	path+'variable_mass+1.0dex/IP_data/zerodSFQ_IP_dz2.0.fits']

	FOREACH file,IPdatafiles DO IPstats, file
END

PRO IPstats, datafile
	
	data	= mrdfits(datafile,1)
	dataIP	= data[where((data.IP EQ 1) AND (data.targ_weight GE 1))]
	dataIPSF= dataIP[where(dataIP.SFQ EQ 1)]
	dataIPQ	= dataIP[where(dataIP.SFQ EQ 0)]

	PRINT, datafile
	PRINT, 'N: ', n_elements(dataIPSF), n_elements(dataIPQ)
	PRINT, 'Median z: ', median(dataIPSF.zprimus), median(dataIPQ.zprimus)
	PRINT, 'Median mass: ', median(dataIPSF.mstar), median(dataIPQ.mstar)
	PRINT, ''
END
