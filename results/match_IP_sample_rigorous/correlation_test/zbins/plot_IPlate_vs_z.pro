PRO plot_IPlate_vs_z, outputFormat;, ptsPerBin, dR
	dN = 2

	xmin = 0.2
	xmax = 1.0
	xr = xmax-xmin

	ymin = 0
	ymax = 1
	yr = ymax-ymin

	ERASE
	!P.MULTI=1
	!P.CHARSIZE = 1.25
	dz_coeff = 2.0

	datapath = '~/conformity/results/match_IP_sample_rigorous/correlation_test/zbins/complete/neighborDataAll_'
;	massRanges = ['11.0_11.5', '10.5_11.0', '10.0_10.5', '9.5_10.0']
	massRanges = ['10.5_11.0', '10.0_10.5']

  !P.FONT=0
  IF (STRING(outputFormat) EQ 'ps') THEN BEGIN
	PS_OPEN, '~/latex/figures/number_neigh_vs_z', /ENCAP, THICK=5
	DEVICE, /INCH, XS=6, YS=5, SET_FONT='Palatino-Roman'
  ENDIF ELSE BEGIN
	SET_PLOT, 'X'
  ENDELSE

;  colors = ['green', 'blue', 'magenta', 'red']
;  PSYMS = [1,4,5,6]
  colors = ['blue', 'magenta']
  PSYMS = [-4,-5]

  zArray = 0.2+FINDGEN(9)/10.
;  zArray = [0.2,0.3,0.4,0.6,0.8,1.0]
  zArrayPlot=[0.25,0.35,0.5,0.7,0.9]
  zbin = zArray[1]-zArray[0]
  PRINT, zArray

  FOR m=0,n_elements(massRanges)-1 DO BEGIN
	neighData = MRDFITS(datapath + massRanges[m] + '.fits', 1)
	
  zRegionArray = []
  FOR j=0,MAX(neighData.JK_region) DO BEGIN
	zDataArray = []
	neighDataJK = neighData[WHERE(neighData.JK_region NE j)]

; EQUAL WIDTH Z BINS BINS
	FOR i=0,n_elements(zArray)-2 DO BEGIN
	  zlow  = zArray[i]
	  zhigh = zArray[i+1]
	
	  zData = neighDataJK[WHERE( (neighDataJK.zIP GE zlow) AND (neighDataJK.zIP LT zhigh), /NULL)]

	  IF (zData NE !NULL) THEN BEGIN
;	    zPointArray = []
;	    FOREACH elem,zData DO zPointArray = [zPointArray,TOTAL(elem.neigh_all_weight_total)]
;	    zPoint = MEAN(zPointArray[WHERE(zPointArray NE 0.)])
	    zPoint = n_elements(zData[WHERE(zData.IP_type EQ 1)])/(FLOAT(n_elements(zData)))
	  ENDIF ELSE BEGIN
	    zPoint=0.
	  ENDELSE
	  zDataArray = [zDataArray, zPoint]
	ENDFOR
	zRegionArray = [[zRegionArray], [zDataArray]]
  ENDFOR
	
	errorArray = []
	FOR i=0,n_elements(zArray)-2 DO BEGIN
		errorVec = zRegionArray[i,*]
		errorArray = [errorArray, SQRT(MAX(neighData.JK_region)-1)*STDDEV(errorVec[WHERE(errorVec GE 0.)])]
	ENDFOR

	zArrayAllRegions= []
	meanNnumArray	= []
	medNnumArray	= []
	N_IParray	= []
	FOR i=0,n_elements(zArray)-2 DO BEGIN
	  zlow  = zArray[i]
	  zhigh = zArray[i+1]

	  zData = neighData[WHERE( (neighData.zIP GE zlow) AND (neighData.zIP LE zhigh), /NULL)]
	
	  IF (zData NE !NULL) THEN BEGIN
;	    zPointArray = []
;	    FOREACH elem,zData DO zPointArray = [zPointArray, TOTAL(elem.neigh_all_weight_total)]
;	    zPoint = MEAN(zPointArray[WHERE(zPointArray NE 0.)])
	    zPoint = n_elements(zData[WHERE(zData.IP_type EQ 1)])/(FLOAT(n_elements(zData)))
	  ENDIF ELSE BEGIN
	    zPoint = 0.
	    meanNnumArray = [meanNnumArray, 0.]
	    medNnumArray  = [medNnumArray, 0.]
	  ENDELSE
	  PRINT, zlow, zhigh, n_elements(zData), zPoint

  	  zArrayAllRegions = [zArrayAllRegions, zPoint]
	  N_IParray = [N_IParray, n_elements(zData)]
	ENDFOR

	PRINT, zArrayAllRegions

      IF (m EQ 0) THEN BEGIN
	PLOT, [0,0], [0,0], LINESTYLE=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $ ;, /XLOG, $
;, POSITION=posArray[*,0], $
	  xtitle = 'Redshift', ytitle='Mean # of Neighbors (weighted)'

;	LEGEND, massRanges, PSYM=4, NUMBER=2, PSPACING=1,  COLOR=cgColor(colors), BOX=0, /BOTTOM, /LEFT
;	LEGEND, [textoidl('log(M_{*IP}/M_{\odot})'),'11.0-11.5', '10.5-11.0', '10.0-10.5', '9.5-10.0'], $
	LEGEND, [textoidl('log(M_{*IP}/M_{\odot})'),'10.5-11.0', '10.0-10.5'], $
		PSYM=[3,PSYMS], NUMBER=2, PSPACING=1,  COLOR=[cgColor('white'),cgColor(colors)], BOX=0, CHARSIZE=1.15, $
		/TOP, /RIGHT

	XYOUTS, 80, ymin+0.925*yr, textoidl('log(M_{*IP}/M_{\odot}) - log(M_{*neigh}/M_{\odot}) < 0.5'), ALIGNMENT=1.0, CHARSIZE=1.15
      ENDIF

	OPLOT, zArray[0:-2]+zbin/2, zArrayAllRegions, PSYM=psyms[m], COLOR=cgColor(colors[m])
	ERRPLOT, zArray[0:-2]+zbin/2, zArrayAllRegions+errorArray, zArrayAllRegions-errorArray, COLOR=cgColor(colors[m])
;	OPLOT, zArrayPlot, zArrayAllRegions, PSYM=psyms[m], COLOR=cgColor(colors[m])
;	ERRPLOT, zArrayPlot, zArrayAllRegions+errorArray, zArrayAllRegions-errorArray, COLOR=cgColor(colors[m])

;	XYOUTS, 1.25, ymin+0.05*yr, textoidl('10<log(M_{*IP})<10.5'),ALIGNMENT=0.0

;	PRINT, TRANSPOSE([[N_IParray], [Narray[1:-1]], [meanNnumArray], [medNnumArray]])
;	PRINT, TRANSPOSE([[N_IParray], [Narray[1:-1]], [meanNnumArray], [zArrayAllRegions], [errorArray]])
;	PRINT, TOTAL(N_IParray)
  ENDFOR
	
  IF (STRING(outputFormat) EQ 'ps') THEN PS_CLOSE
  SET_PLOT, 'X'
	  
END
