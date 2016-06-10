PRO field_stats

	data	= MRDFITS('~/results/zerodSFQ_all_cart.fits',1)
	dataIP	= MRDFITS('~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits',1)

	fields = data[UNIQ(data.field, SORT(data.field))].field

	data	= data[WHERE(data.targ_weight GE 1.)]
	dataIP	= dataIP[WHERE(dataIP.targ_weight GE 1.)]

	FOR f=0,n_elements(fields)-1 DO BEGIN
		field = fields[f]
		dataF	= data[WHERE(data.field EQ field)]
		dataFSF = dataF[WHERE(dataF.SFQ EQ 1)]
		dataFQ  = dataF[WHERE(dataF.SFQ EQ 0)]

		dataIPF	 = dataIP[WHERE(dataIP.field EQ field)]
		dataIPFSF = dataIPF[WHERE(dataIPF.SFQ EQ 1)]
		dataIPFSFuniq = dataIPFSF[UNIQ(dataIPFSF.objname, SORT(dataIPFSF.objname))]
		dataIPFQ  = dataIPF[WHERE(dataIPF.SFQ EQ 0)]

		PRINT, field, n_elements(dataFQ), n_elements(dataFSF), n_elements(dataIPFQ), $
			n_elements(dataIPFSF), n_elements(dataIPFSFuniq)
	ENDFOR
END
