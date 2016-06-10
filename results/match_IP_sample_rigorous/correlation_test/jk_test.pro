PRO jk_test;, minIPmstar
	dR = 0.1 ; Mpc
	dz_coeff = 2.0
	targPhi = -3.7
	stringPHI = strtrim(string(-1.*targPhi, format='(f20.1)'), 2)

	dataAll = MRDFITS('~/results/conservative_mass_cutoff/allAboveMassCompLim-0.0.fits',1)
	dataAll = dataAll[WHERE(dataAll.targ_weight GE 1.)]

	dataIPcomp	= MRDFITS('~/results/match_IP_sample_rigorous/correlation_test/allAboveMassCompLim+0.5_IP.fits',1)
	dataIPcomp 	= dataIPcomp[WHERE((dataIPcomp.IP EQ 1) AND (dataIPcomp.targ_weight GE 1.))]

	fields = dataAll[uniq(dataAll.field, sort(dataAll.field))].field

	mins = [10.1, 10.4, 10.7]
	zmax = [0.65, 0.8, 1.0]

  FOR i=0,2 DO BEGIN
	minIPmstar = mins[i]
	maxIPmstar = minIPmstar + 0.3
	dataIP	= dataIPcomp[WHERE((dataIPcomp.mstar GE minIPmstar) AND (dataIPcomp.mstar LE maxIPmstar) AND $
		(dataIPcomp.zprimus LE zmax[i]))]
	PRINT, minIPmstar, maxIPmstar, zmax[i]

	FOR f=0,n_elements(fields)-1 DO BEGIN
		dataIPfield = dataIP[WHERE(dataIP.field EQ fields[f])]
		PRINT, fields[f], n_elements(dataIPfield), '    ', decimal(n_elements(dataIPfield)/FLOAT(n_elements(dataIP)),2)
	ENDFOR
	PriNT, ''
  ENDFOR
END
