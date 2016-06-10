PRO getGALEXmassSFR

	uber_fields = ['cdfs', 'cosmos', 'es1', 'cfhtls_xmm', 'xmm_swire']

; read in zerod for specific field
;	zerod_all = primus_read_zerod(rerun=2165)
;	zerod = zerod_all[where( (zerod_all.zprimus_zconf ge 3) and (zerod_all.zprimus_class eq 'GALAXY') )]
;	zerodSFQ_all = []

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

	; read in M* and SFR data for specific field
;	load_data, field, data, mass
	; match dataset entries by coordinates
;	spherematch, zerod.ra, zerod.dec, data.ra, data.dec, 1/3600d, zz, dd

	; zz : indices of zerod with matches in data
	; dd : indices of data with matches in zerod
	; mass[dd] is M* corresponding to data[dd]

	; uncertainty contributions from M* and z

;	n = size(dd, /n_elements)

	; add SFR, stellar mass, and SF or Q tags to data structure
	; my first way of doing this
;	templateRow = create_struct(zerod[0], 'SFR', 0.0, 'mstar', 0.0, 'SFQ', 0)
;	zerodSFQ = replicate(templateRow, n)

	; Alison's (simpler) way
	temp_struct 	= mrd_struct(['SFR', 'mstar', 'SFQ'], ['0.0', '0.0', '0'], n_elements(data))

	; combine data and temp_struct and call the result data_massSFR
	combine_structs, data, temp_struct, data_massSFR

	; fill in actual stellar masses and SFRs
	data_massSFR.mstar = mass.mass_50
	data_massSFR.SFR   = mass.sfr100_50

	; define threshold for SF (vs Q) status
	SFRcut		= -1.29 + 0.65*(data_massSFR.mstar - 10) + 1.33*(data_massSFR.zprimus - 0.1)
	
	; define galaxies with SFRs above the cut as SF (default is Q)
	data_massSFR[where(data_massSFR.SFR gt SFRcut)].SFQ = 1

	data_final 	= data_massSFR[where((data_massSFR.zprimus_zconf ge 3) AND $

					     (data_massSFR.zprimus_class eq 'GALAXY'))]

	PRINT, n_elements(data_final)

;	FOR i = 0, n - 1 DO BEGIN

	; loop through all elements of zerod with matches in data
;		z  	= zerod[zz[i]].zprimus
;		SFR	= mass[dd[i]].sfr100_50
;		SFR	= mass[dd[i]].sfr_50
;		Mstar	= mass[dd[i]].mass_50
;	 	logSFRmin = -0.49 + 0.65*( Mstar-10 )+1.07*(z - 0.1) ; Moustakas 2013
;	 	logSFRmin = -1.29 + 0.65*( Mstar-10 )+1.33*(z - 0.1)
;		IF (SFR ge logSFRmin) THEN (SFQ = 1) ELSE (SFQ = 0)

;		newRow = create_struct(zerod[zz[i]], 'SFR', SFR, 'Mstar', Mstar, 'SFQ', SFQ)
;		zerodSFQ[i] = newRow
;	ENDFOR

;	IF (j=0) THEN (zerodSFQ_all = [zerodSFQ]) ELSE (zerodSFQ_all = [zerodSFQ_all, zerodSFQ])
;	zerodSFQ_all = [zerodSFQ_all, zerodSFQ]

	mwrfits, data_final, "zerodSFQ_all.fits", /create

END
