PRO plot_cc, ptsPerBin, dR
	ytitle=''
	ymin = -0.2
	ymax = 0.8
	yr = ymax-ymin

	ERASE
	!P.MULTI=2
	posArray = grid_array(1,2)
;	!P.MULTI=[2,1,2]
;	dR = 0.1 ; Mpc
	dz_coeff = 2.0
	targPhi = -3.7
	stringPHI = strtrim(string(-1.*targPhi, format='(f20.1)'), 2)

	dataAll = mrdfits('~/conformity/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits',1)
	dataAll	= dataAll[WHERE(dataAll.targ_weight GE 1.)]
	dataIP	= dataAll[WHERE(dataAll.IP EQ 1)]

        fields = dataAll[uniq(dataAll.field, sort(dataAll.field))].field

; IP stellar mass cut
	dataIP	= dataIP[WHERE((dataIP.mstar GE 10.) AND (dataIP.mstar LE 10.5))]
;	PRINT, 'IP: ', n_elements(dataIP)

;	'dist', 'IP_SSFRcorr', 'neigh_SSFRcorr'

	neighData = MRDFITS('~/conformity/results/match_IP_sample_rigorous/correlation_test/neighborData.fits', 1)

;	PRINT, MINMAX(neighData.dist)
	Rmin = 0.
;	Rmax = CEIL(100*MAX(neighData.dist))/100
	Rmax = 5.
	Rarray = dR*FINDGEN(CEIL(Rmax/dR))
	PRINT, Rarray

; EQUAL # OF DATA POINTS PER BIN (VARIABLE WIDTH RADIAL BINS)
	neighData_sort = neighData[SORT(neighData.dist)]
;	ptsPerBin = 500.
	maxRadius = 10.
	IPwithinMaxRadius = n_elements(neighData[WHERE(neighData.dist LE maxRadius)])

	nBins = FLOOR(IPwithinMaxRadius/ptsPerBin)

	RarrayUnequal	= []
	cc		= []
	ccErrorArray	= []
	FOR i=0,nBins-1 DO BEGIN
;	  PRINT, i, i*ptsPerBin, (i+1)*ptsPerBin-1
	  binData = neighData_sort[i*ptsPerBin:(i+1)*ptsPerBin-1]
;	  PRINT, 'IP in bin: ', n_elements(binData)

	  binR = MEDIAN(binData.dist)
	  ccPoint = CORRELATE(binData.IP_SSFRcorr, binData.neigh_SSFRcorr)

	; ERROR: for each radial bin, resample 200 times and compute stddev of 200 results (cc)
	  binData_90array = []
	  FOR s=1,200 DO BEGIN
	    seed = s
	    random_indices = ROUND( RANDOMU( seed, ROUND( 0.9*n_elements(binData) ) )*n_elements(binData) )
	    binData_90 = binData[random_indices]
	    binData_90array = [binData_90array, CORRELATE(binData_90.IP_SSFRcorr, binData_90.neigh_SSFRcorr)]
	  ENDFOR

	  binError	= STDDEV(binData_90array)
;	  PRINT, i, ' ', binError
	  RarrayUnequal	= [RarrayUnequal, binR]
	  cc		= [cc, ccPoint]
	  ccErrorArray	= [ccErrorArray, binError]
	ENDFOR

;	PRINT, RarrayUnequal[0:-2], cc
	PLOT, [0,5], [0,0], LINESTYLE=1, xrange=[0,5], yrange=[ymin,ymax], POSITION=posArray[*,1], $
	  xtickformat='(A1)', ytitle=ytitle
;	LEGEND, [textoidl('10<M_{*IP}<10.5'),textoidl('M_{*IP}-M_{*neigh}<0.5')], PSYM=[3,3], COLOR=cGColor(['white','white']), /TOP, /LEFT, BOX=0
	LEGEND, [textoidl('10<M_{*IP}<10.5'),textoidl('M_{*IP}-M_{*neigh}<0.5')], PSYM=[3,3], COLOR=cGColor(['white','white']), /TOP, /RIGHT, BOX=0
	XYOUTS, 0.25, ymin+0.1*yr, 'Equal IP per bin'
	OPLOT, RarrayUnequal[0:-2], cc, PSYM=4
	ERRPLOT, RarrayUnequal[0:-2], cc+ccErrorArray, cc-ccErrorArray

	cc		= []
	ccErrorArray	= []
; EQUAL WIDTH RADIAL BINS
	FOR i=0,n_elements(Rarray)-2 DO BEGIN
	  Rlow	= Rarray[i]
	  Rhigh	= Rarray[i+1]
	
;	  PRINT, Rlow, Rhigh

	  Rdata = neighData[WHERE((neighData.dist GE Rlow) AND (neighData.dist LT Rhigh), /NULL)]
;	  PRINT, Rlow, n_elements(Rdata)
	  IF (Rdata NE !NULL) THEN $
	    ccPoint =  CORRELATE(Rdata.IP_SSFRcorr, Rdata.neigh_SSFRcorr) $
	  ELSE ccPoint=-99

	  Rdata_90array = []
	  IF (ccPoint NE -99) THEN BEGIN
	    FOR s=1,200 DO BEGIN
	      seed = s
	      random_indices = ROUND( RANDOMU( seed, ROUND( 0.9*n_elements(Rdata) ) )*n_elements(Rdata) )
	      Rdata_90 = Rdata[random_indices]
	      Rdata_90array = [Rdata_90array, CORRELATE(Rdata_90.IP_SSFRcorr, Rdata_90.neigh_SSFRcorr)]
	    ENDFOR
	    binError = STDDEV(Rdata_90array)
	  ENDIF ELSE BEGIN
	    binError = 0
	  ENDELSE

	  cc = [cc, ccPoint]
	  ccErrorArray = [ccErrorArray, binError]
	END
	
	PRINT, TRANSPOSE([[Rarray[0:-2]+0.5*dR], [cc], [ccErrorArray]])

	PLOT, [0,5], [0,0], LINESTYLE=1, xrange=[0,5], yrange=[ymin,ymax], POSITION=posArray[*,0], $
	  xtitle = 'Proj. dist. to nearest neighbor (Mpc)', ytitle=ytitle
	XYOUTS, 0.25, ymin+0.1*yr, 'Equal radial binwidth'
	OPLOT, Rarray[0:-2]+0.5*dR, cc, PSYM=4
	ERRPLOT, Rarray[0:-2]+0.5*dR, cc+ccErrorArray, cc-ccErrorArray
END
