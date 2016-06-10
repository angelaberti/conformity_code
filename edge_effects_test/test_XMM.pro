PRO test_XMM, outputFormat

	data = mrdfits('~/results/match_IP_sample_rigorous/matchedIPsampleFBF_PHI3.7.fits',1)

	dataSXDS = data[WHERE( (data.field EQ 'xmm       ') AND (data.targ_weight GE 1.) )]
	PRINT, n_elements(dataSXDS)
	dataCFH	 = data[WHERE( (data.field EQ 'cfhtls_xmm') AND (data.targ_weight GE 1.) )]
	PRINT, n_elements(dataCFH)

	dataSXDS_IPSF	= dataSXDS[WHERE(dataSXDS.SFQ EQ 1)]
	dataCFH_IPSF	= dataCFH[WHERE(dataCFH.SFQ EQ 1)]

	SXDS_IPSF	= n_elements(dataSXDS_IPSF[UNIQ(dataSXDS_IPSF.objname, SORT(dataSXDS_IPSF.objname))].objname)
	CFH_IPSF	= n_elements(dataCFH_IPSF[UNIQ(dataCFH_IPSF.objname, SORT(dataCFH_IPSF.objname))].objname)

;	SXDS_IPSF 	= n_elements(SXDS_IPSF)
;	CFH_IPSF	= n_elements(CFH_IPSF)

	SXDS_IPQ	= n_elements(dataSXDS[WHERE(dataSXDS.SFQ EQ 0)])
	CFH_IPQ		= n_elements(dataCFH[WHERE(dataCFH.SFQ EQ 0)])

	PRINT, 'SF & Q IP per pointing'
	PRINT, 'XMM-SXDS: ', SXDS_IPSF/5., SXDS_IPQ/5., FLOAT(SXDS_IPSF)/SXDS_IPQ
	PRINT, 'XMM-CFHTLS: ', CFH_IPSF/15., CFH_IPQ/15., FLOAT(CFH_IPSF)/CFH_IPQ

	!p.multi=0
	!p.thick=5
	!p.charsize=1.2

	W = 0.5
	A = FINDGEN(33)*(!PI*2/32.)
	USERSYM, W*COS(A), W*SIN(A), /FILL
	SYM = 8

	IF (string(outputFormat) EQ 'ps') THEN BEGIN
;		SET_PLOT, 'ps'
;		DEVICE, file='~/figures/edge_effects_test.ps', /LANDSCAPE
	        PS_OPEN, '~/latex/figures/edge_effects_test_allXMM, THICK=5, /ENCAP
	        DEVICE, /INCH, XS=8, YS=8
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

;	all = MRDFITS('~/zerodSFQ_all_cart.fits',1)

	

	rand	= mrdfits('random/xmm_random.fits',1)
		
	RAmin = min(rand.RA)
	RAmax = max(rand.RA)
	RArange = ABS(RAmax-RAmin)
	RAmid = mean(minmax(rand.RA))

	DECmin = min(rand.DEC)
	DECmax = max(rand.DEC)
	DECrange = ABS(DECmax-DECmin)
	DECmid = mean(minmax(rand.DEC))

	ar = DECrange/RArange

	h = 0.52

;	PLOT, rand.RA, rand.DEC, xtitle='RA', ytitle='DEC', $
;		xrange=[RAmid-h*RArange,RAmid+h*RArange], yrange=[DECmid-h*DECrange,DECmid+h*DECrange], /NODATA
;		position=aspect(ar,margin=0.15), /NODATA
;	OPLOT, rand.RA, rand.DEC, psym=3, color=cgColor('yellow')
;	OPLOT, dataXMMSXDS[where(dataXMMSXDS.IP EQ 1)].RA, dataXMMSXDS[where(dataXMMSXDS.IP EQ 1)].DEC, psym=SYM, COLOR=cgColor('red')
;	OPLOT, dataXMMCFHTLS[where(dataXMMCFHTLS.IP EQ 1)].RA, dataXMMCFHTLS[where(dataXMMCFHTLS.IP EQ 1)].DEC, psym=SYM, color=cgColor('blue')

	IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END
