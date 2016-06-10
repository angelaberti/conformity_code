PRO plot_neigh_vs_z, outputFormat;, ptsPerBin, dR
	dN = 2

	xmin = 0.2
	xmax = 1.0
	xr = xmax-xmin

	ymin = 0
	ymax = 40
	yr = ymax-ymin

	ERASE
	!P.MULTI=1
	!P.CHARSIZE = 1.25
	dz_coeff = 2.0

;	datapath = '~/conformity/results/match_IP_sample_rigorous/correlation_test/neighborDataAll_'
	datapath = '~/conformity/results/match_IP_sample_rigorous/correlation_test/conservative/neighborDataAll_'
;	massRanges = ['11.0_11.5', '10.5_11.0', '10.0_10.5', '9.5_10.0']
;	massRanges = ['10.5_11.0', '10.0_10.5']
	massRanges = ['10.5_11.0', '10.1_10.5']

  !P.FONT=0
  IF (STRING(outputFormat) EQ 'ps') THEN BEGIN
	PS_OPEN, '~/latex/figures/number_neigh_vs_z_M13massLim', /ENCAP, THICK=5
	DEVICE, /INCH, XS=6, YS=5, SET_FONT='Palatino-Roman'
  ENDIF ELSE BEGIN
	SET_PLOT, 'X'
  ENDELSE

;  colors = ['green', 'blue', 'magenta', 'red']
;  PSYMS = [1,4,5,6]
  colors = ['blue', 'magenta']
  PSYMS = [-3,-3]

  zArray = [0.2,0.3,0.4,0.5,0.6,0.7,0.8]
  zbin	= zArray[1]-zArray[0]
;  PRINT, zArray

  FOR m=0,n_elements(massRanges)-1 DO BEGIN
	neighData = MRDFITS(datapath + massRanges[m] + '.fits', 1)
	IF (m EQ 1) THEN neighData=neighData[WHERE(neighData.zIP LE 0.65)] ELSE $
		neighData=neighData[WHERE(neighData.zIP LE 0.8)]
		
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
	    zPointArray	= []
	    z25Array	= []
	    z75Array	= []
	    FOREACH elem,zData DO zPointArray = [zPointArray,TOTAL(elem.neigh_all_weight_total)]
;	    zPoint	= MEDIAN(zPointArray[WHERE(zPointArray NE 0.)])
	    zPoint	= MEDIAN(zPointArray)
	    z25		= MEDIAN(zPointArray[WHERE(zPointArray LE zPoint)])
	    z75		= MEDIAN(zPointArray[WHERE(zPointArray GE zPoint)])
	  ENDIF ELSE BEGIN
	    zPoint	= 0.
	  ENDELSE
	  zDataArray	= [zDataArray, zPoint]
	  z25Array	= [z25Array, z25]
	  z75Array	= [z75Array, z75]
	ENDFOR
	zRegionArray	= [[zRegionArray], [zDataArray]]
  ENDFOR
	
	errorArray = []
	FOR i=0,n_elements(zArray)-2 DO BEGIN
		errorVec = zRegionArray[i,*]
		errorArray = [errorArray, SQRT(MAX(neighData.JK_region)-1)*STDDEV(errorVec[WHERE(errorVec GE 0.)])]
	ENDFOR

	zArrayAllRegions= []
	z25ArrayAllRegions = []
	z75ArrayAllRegions = []
	z90ArrayAllRegions = []
	meanNnumArray	= []
	medNnumArray	= []
	N_IParray	= []
	FOR i=0,n_elements(zArray)-2 DO BEGIN
	  zlow  = zArray[i]
	  zhigh = zArray[i+1]

	  zData = neighData[WHERE( (neighData.zIP GE zlow) AND (neighData.zIP LE zhigh), /NULL)]
	
	  IF (zData NE !NULL) THEN BEGIN
	    zPointArray	= []
	    FOREACH elem,zData DO zPointArray = [zPointArray,TOTAL(elem.neigh_all_weight_total)]
;	    zPoint	= MEDIAN(zPointArray[WHERE(zPointArray NE 0.)])
	    zPoint	= MEDIAN(zPointArray)
	    z25		= MEDIAN(zPointArray[WHERE(zPointArray LE zPoint)])
	    z75		= MEDIAN(zPointArray[WHERE(zPointArray GE zPoint)])
	    zPointArraySorted = zPointArray[SORT(zPointArray)]
	    z90		= zPointArraySorted[CEIL(0.9*n_elements(zData))]
	    PRINT, z90

	  ENDIF ELSE BEGIN
	    zPoint	= 0.
	    meanNnumArray = [meanNnumArray, 0.]
	    medNnumArray  = [medNnumArray, 0.]
	  ENDELSE
;	  PRINT, zlow, zhigh, n_elements(zData), zPoint

	  zArrayAllRegions	= [zArrayAllRegions, zPoint]
	  z25ArrayAllRegions	= [z25ArrayAllRegions, z25]
	  z75ArrayAllRegions	= [z75ArrayAllRegions, z75]
	  z90ArrayAllRegions	= [z90ArrayAllRegions, z90]
	  N_IParray = [N_IParray, n_elements(zData)]
	ENDFOR

	PRINT, zArrayAllRegions

      IF (m EQ 0) THEN BEGIN
	PLOT, [0,0], [0,0], LINESTYLE=1, xrange=[xmin,xmax], yrange=[ymin,ymax], $ ;, /XLOG, $
;, POSITION=posArray[*,0], $
	  xtitle = 'Redshift', ytitle='Number of Neighbors (weighted)'

;	LEGEND, massRanges, PSYM=4, NUMBER=2, PSPACING=1,  COLOR=cgColor(colors), BOX=0, /BOTTOM, /LEFT
;	LEGEND, [textoidl('log(M_{*IP}/M_{\odot})'),'11.0-11.5', '10.5-11.0', '10.0-10.5', '9.5-10.0'], $
;	LEGEND, [textoidl('log(M_{*IP}/M_{\odot}) range'),'10.5-11.0', '10.0-10.5','quartiles',textoidl('90^{th} %ile')], $
	LEGEND, [textoidl('log(M_{*IP}/M_{\odot})=10.5-11.0 median'), '10.1-10.5 median','quartiles',textoidl('90^{th} %ile')], $
		PSYM=[PSYMS,'',''], LINESTYLE=['','',2,1],  NUMBER=2, PSPACING=1,  COLOR=[cgColor(colors),cgColor('black'),cgColor('black')], BOX=0, CHARSIZE=1.15, $
		/TOP, /RIGHT

	XYOUTS, xmin+0.05*xr, ymin+0.05*yr, textoidl('log(M_{*IP}/M_{\odot}) - log(M_{*neigh}/M_{\odot}) < 0.5'), ALIGNMENT=0.0, CHARSIZE=1.15
      ENDIF
	IF (m EQ 1) THEN BEGIN
	OPLOT, zArray[0:-3]+zbin/2, zArrayAllRegions, PSYM=psyms[m], COLOR=cgColor(colors[m])
	ERRPLOT, zArray[0:-3]+zbin/2, zArrayAllRegions+errorArray, zArrayAllRegions-errorArray, COLOR=cgColor(colors[m])

	OPLOT, zArray[0:-3]+zbin/2, z25ArrayAllRegions, LINESTYLE=2, COLOR=cgColor(colors[m])
	OPLOT, zArray[0:-3]+zbin/2, z75ArrayAllRegions, LINESTYLE=2, COLOR=cgColor(colors[m])
	OPLOT, zArray[0:-3]+zbin/2, z90ArrayAllRegions, LINESTYLE=1, COLOR=cgColor(colors[m])
	ENDIF

	IF (m EQ 0) THEN BEGIN
	OPLOT, zArray[0:-2]+zbin/2, zArrayAllRegions, PSYM=psyms[m], COLOR=cgColor(colors[m])
	ERRPLOT, zArray[0:-2]+zbin/2, zArrayAllRegions+errorArray, zArrayAllRegions-errorArray, COLOR=cgColor(colors[m])

	OPLOT, zArray[0:-2]+zbin/2, z25ArrayAllRegions, LINESTYLE=2, COLOR=cgColor(colors[m])
	OPLOT, zArray[0:-2]+zbin/2, z75ArrayAllRegions, LINESTYLE=2, COLOR=cgColor(colors[m])
	OPLOT, zArray[0:-2]+zbin/2, z90ArrayAllRegions, LINESTYLE=1, COLOR=cgColor(colors[m])
	ENDIF

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
