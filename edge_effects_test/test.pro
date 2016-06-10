PRO test, outputFormat
	IPcolor = cgColor('Black')
;	IPcolor = cgColor('White')

	path = '../results/default_parameters/IP_data/'
	inputFile = 'zerodSFQ_IP_dz1.0_dm0.5.fits'

	data = mrdfits(path+inputFile,1)

	fields = data[uniq(data.field, sort(data.field))].field

	!p.multi=0
	!p.thick=5
	!p.charsize=1.2

	A = FINDGEN(33)*(!PI*2/32.)
	USERSYM, 0.25*COS(A), 0.25*SIN(A), /FILL

	FOR f=0,n_elements(fields)-2 DO BEGIN
		dataf 	= data[where(data.field eq fields[f])]

		IF (string(outputFormat) eq 'ps') THEN BEGIN
;			SET_PLOT, 'ps'
;			DEVICE, file='../figures/edge_effects_test.ps', /LANDSCAPE
		        PS_OPEN, '../figures/edge_effects_test_'+strtrim(fields[f]), THICK=5, /ENCAP
		        DEVICE, /INCH, XS=8, YS=6
		ENDIF ELSE BEGIN
		        SET_PLOT, 'X'
		ENDELSE

		rand	= mrdfits("random/" + strtrim(strcompress(fields[f]),2) + "_random.fits",1)
		title 	= strtrim(string(fields(f)),1)
		
		IF (fields(f) eq 'cfhtls_xmm') OR (fields(f) eq 'xmm       ') THEN BEGIN
			dataf 	= data[where((data.field eq 'cfhtls_xmm') OR (data.field eq 'xmm       '))]
			rand 	= mrdfits("random/xmm_random.fits",1)
			title 	= 'xmm'
		ENDIF

		RAmin = min(dataf.RA)
		RAmax = max(dataf.RA)
		RArange = ABS(RAmax-RAmin)
		RAmid = mean(minmax(dataf.RA))

		DECmin = min(dataf.DEC)
		DECmax = max(dataf.DEC)
		DECrange = ABS(DECmax-DECmin)
		DECmid = mean(minmax(dataf.DEC))

		ar = DECrange/RArange

		h = 0.52

		PLOT, dataf.RA, dataf.DEC, psym=8, xtitle="RA", ytitle="DEC", title=strtrim(strcompress(fields[f]),1), $
			xrange=[RAmid-h*RArange,RAmid+h*RArange], yrange=[DECmid-h*DECrange,DECmid+h*DECrange], $
			position=aspect(ar,margin=0.15), /NODATA
;		XYOUTS, RAmid-h*RArange, 1.05*(DECmid+h*DECrange), 'Black = IP'
		OPLOT, rand.RA, rand.DEC, psym=3, color=cgColor('Green')
		OPLOT, dataf[where(dataf.IP eq 1)].RA, dataf[where(dataf.IP eq 1)].DEC, psym=8, color=IPcolor

		IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
		SET_PLOT, 'X'
	END
END
