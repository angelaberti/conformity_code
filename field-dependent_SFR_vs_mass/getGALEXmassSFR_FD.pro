PRO getGALEXmassSFR_FD

	uber_fields = ['cdfs', 'cosmos', 'es1', 'cfhtls_xmm', 'xmm_swire']

	basedir = '/raid/primus/science/mf/'

	j=0
	; begin with first field
	data = mrdfits(basedir+'2165/ubersample/'+uber_fields[j]+'_zerod_v20.fits.gz', 1)
	mass = mrdfits(basedir+'isedfit/'+uber_fields[j]+'_fsps_chab_charlot_sfhgrid01.fits.gz', 1)

	mass = struct_trimtags(mass, except=['MAGGIES','IVARMAGGIES','BESTMAGGIES'])

	; add remaining fields to data and mass structures
	FOR j=1,n_elements(uber_fields)-1 DO BEGIN
		data1 = mrdfits(basedir+'2165/ubersample/'+uber_fields[j]+'_zerod_v20.fits.gz', 1)
		mass1 = mrdfits(basedir+'isedfit/'+uber_fields[j]+'_fsps_chab_charlot_sfhgrid01.fits.gz', 1)

		mass1 = struct_trimtags(mass1, except=['MAGGIES','IVARMAGGIES','BESTMAGGIES'])

		data = [data,data1]
		mass = [mass,mass1]
	ENDFOR

	PRINT, 'data:', n_elements(data), '  mass:', n_elements(mass)

	; Alison's (simpler) way
	temp_struct 	= mrd_struct(['SFR', 'mstar', 'SFQ'], ['0.0', '0.0', '0'], n_elements(data))

	; combine data and temp_struct and call the result data_massSFR
	combine_structs, data, temp_struct, data_massSFR

	; fill in actual stellar masses and SFRs
	data_massSFR.mstar = mass.mass_50
	data_massSFR.SFR   = mass.sfr100_50

	fields = ['cdfs', 'cosmos', 'es1', 'cfhtls_xmm', 'xmm']

	allData = []
	FOR f=0,4 DO BEGIN
		currentField = fields[f]

		; define threshold for SF (vs Q) status
		IF currentField EQ 'cdfs' THEN offset = 0.2 ELSE $
		IF currentFIeld EQ 'cfhtls_xmm' THEN offset = -0.2 ELSE offset = 0.0

		fieldData = data_massSFR[where(strcompress(data_massSFR.field) EQ currentField)]

		SFRcut = -1.29 + offset + 0.65*(fieldData.mstar - 10) + 1.33*(fieldData.zprimus - 0.1)
	
		; define galaxies with SFRs above the cut as SF (default is Q)
		fieldData[where(fieldData.SFR gt SFRcut)].SFQ = 1
		
		allData = [allData, fieldData]
	ENDFOR

	data_final = allData[where((allData.zprimus_zconf ge 3) AND (allData.zprimus_class eq 'GALAXY'))]
	PRINT, n_elements(data_final)

	mwrfits, data_final, '~/field-dependent_SFR_vs_mass/zerodSFQ_all.fits', /create
END
