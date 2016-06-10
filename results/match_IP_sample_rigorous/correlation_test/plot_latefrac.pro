PRO plot_latefrac, ptsPerBin, dR
	ytitle=textoidl('IP f_{late}')
	ymin = 0
	ymax = 1
	yr = ymax-ymin

	ERASE
	!P.MULTI=2
	posArray = grid_array(1,2)
;	!P.MULTI=[2,1,2]
;	dR = 0.1 ; Mpc
	dz_coeff = 2.0
	targPhi = -3.7
	stringPHI = strtrim(string(-1.*targPhi, format='(f20.1)'), 2)

	dataAll = mrdfits('~/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits',1)
	dataAll	= dataAll[WHERE(dataAll.targ_weight GE 1.)]
	dataIP	= dataAll[WHERE(dataAll.IP EQ 1)]

        fields = dataAll[uniq(dataAll.field, sort(dataAll.field))].field

; IP stellar mass cut
	dataIP	= dataIP[WHERE((dataIP.mstar GE 10.) AND (dataIP.mstar LE 10.5))]
;	PRINT, 'IP: ', n_elements(dataIP)

;	'dist', 'IP_SSFRcorr', 'neigh_SSFRcorr'

	neighData = MRDFITS('~/results/match_IP_sample_rigorous/correlation_test/neighborData_10.0_10.5.fits', 1)

;	PRINT, MINMAX(neighData.dist)
	Rmin = 0.
;	Rmax = CEIL(100*MAX(neighData.dist))/100
	Rmax = 15.
	Rarray = dR*FINDGEN(CEIL(Rmax/dR))
;	PRINT, Rarray

; EQUAL # OF DATA POINTS PER BIN (VARIABLE WIDTH RADIAL BINS)
	neighData_sort = neighData[SORT(neighData.dist)]
;	ptsPerBin = 500.
	maxRadius = 10.
	IPwithinMaxRadius = n_elements(neighData[WHERE(neighData.dist LE maxRadius)])

	nBins = FLOOR(IPwithinMaxRadius/ptsPerBin)

	RarrayUnequal	= []
	lf		= []
	lfErrorArray	= []
	FOR i=0,nBins-1 DO BEGIN
;	  PRINT, i, i*ptsPerBin, (i+1)*ptsPerBin-1
	  binData = neighData_sort[i*ptsPerBin:(i+1)*ptsPerBin-1]
;	  PRINT, 'IP in bin: ', n_elements(binData)

	  binR = MEDIAN(binData.dist)
	  lfPoint = n_elements(binData[WHERE(binData.IP_type EQ 1)])/FLOAT(n_elements(binData))

	; ERROR: for each radial bin, resample 200 times and compute stddev of 200 results (lf)
	  binData_90array = []
	  FOR s=1,200 DO BEGIN
	    seed = s
	    random_indices = ROUND( RANDOMU( seed, ROUND( 0.9*n_elements(binData) ) )*n_elements(binData) )
	    binData_90 = binData[random_indices]
	    binData_90array = [binData_90array, n_elements(binData_90[WHERE(binData_90.IP_type EQ 1)])/FLOAT(n_elements(binData_90))]
	  ENDFOR

	  binError	= STDDEV(binData_90array)
	  PRINT, i, ' ', binError
	  RarrayUnequal	= [RarrayUnequal, binR]
	  lf		= [lf, lfPoint]
	  lfErrorArray	= [lfErrorArray, binError]
	ENDFOR

;	PRINT, RarrayUnequal[0:-2], lf
	PLOT, [0,5], [0,0], LINESTYLE=1, xrange=[0,Rmax], yrange=[ymin,ymax], POSITION=posArray[*,1], $
	  xtickformat='(A1)', ytitle=ytitle
;	LEGEND, [textoidl('10<M_{*IP}<10.5'),textoidl('M_{*IP}-M_{*neigh}<0.5')], PSYM=[3,3], COLOR=cGColor(['white','white']), /TOP, /LEFT, BOX=0
	LEGEND, [textoidl('10<M_{*IP}<10.5'),textoidl('M_{*IP}-M_{*neigh}<0.5')], PSYM=[3,3], COLOR=cGColor(['white','white']), /TOP, /RIGHT, BOX=0
	XYOUTS, 0.25, ymin+0.1*yr, 'Equal IP per bin'
	OPLOT, RarrayUnequal[0:-2], lf, PSYM=4
	ERRPLOT, RarrayUnequal[0:-2], lf+lfErrorArray, lf-lfErrorArray

	lf		= []
	lfErrorArray	= []
; EQUAL WIDTH RADIAL BINS
	FOR i=0,n_elements(Rarray)-2 DO BEGIN
	  Rlow	= Rarray[i]
	  Rhigh	= Rarray[i+1]
	
;	  PRINT, Rlow, Rhigh

	  Rdata = neighData[WHERE((neighData.dist GE Rlow) AND (neighData.dist LT Rhigh), /NULL)]
;	  PRINT, Rlow, n_elements(Rdata)
	  IF (Rdata NE !NULL) THEN $
	    lfPoint =  n_elements(Rdata[WHERE(Rdata.IP_type EQ 1)])/FLOAT(n_elements(Rdata)) $
	  ELSE lfPoint=-99

	  Rdata_90array = []
	  IF (lfPoint NE -99) THEN BEGIN
	    FOR s=1,200 DO BEGIN
	      seed = s
	      random_indices = ROUND( RANDOMU( seed, ROUND( 0.9*n_elements(Rdata) ) )*n_elements(Rdata) )
	      Rdata_90 = Rdata[random_indices]
	      Rdata_90array = [Rdata_90array, n_elements(Rdata_90[WHERE(Rdata_90.IP_type EQ 1)])/FLOAT(n_elements(Rdata_90))]
	    ENDFOR
	    binError = STDDEV(Rdata_90array)
	  ENDIF ELSE BEGIN
	    binError = 0
	  ENDELSE

	  lf = [lf, lfPoint]
	  lfErrorArray = [lfErrorArray, binError]
	END
	
	PRINT, TRANSPOSE([[Rarray[0:-2]+0.5*dR], [lf], [lfErrorArray]])

	PLOT, [0,5], [0,0], LINESTYLE=1, xrange=[0,Rmax], yrange=[ymin,ymax], POSITION=posArray[*,0], $
	  xtitle = 'Proj. dist. to nearest neighbor (Mpc)', ytitle=ytitle
	XYOUTS, 0.25, ymin+0.1*yr, 'Equal radial binwidth'
	OPLOT, Rarray[0:-2]+0.5*dR, lf, PSYM=4
	ERRPLOT, Rarray[0:-2]+0.5*dR, lf+lfErrorArray, lf-lfErrorArray
END
