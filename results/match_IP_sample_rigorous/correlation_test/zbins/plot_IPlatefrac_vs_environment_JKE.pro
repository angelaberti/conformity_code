PRO plot_IPlatefrac_vs_environment_JKE, outputFormat, m_input;, zmin, zmax;, ptsPerBin, dR
	dN = 2
	ymin = 0.0
	ymax = 1.0
	yr = ymax-ymin

	ERASE
	!P.MULTI=1
	!P.CHARSIZE = 1.25
	dz_coeff = 2.0

;	datapath = '~/conformity/results/match_IP_sample_rigorous/correlation_test/neighborDataAll_'
	datapath = '~/conformity/results/match_IP_sample_rigorous/correlation_test/conservative/neighborDataAll_'
;	massRanges = ['10.5_11.0', '10.0_10.5']
	massRanges = ['10.5_11.0', '10.1_10.5']

	neighData_allz = MRDFITS(datapath + massRanges[m_input] + '.fits', 1)

	zz=neighData_allz
	zz=zz[WHERE(zz.zIP LE 0.65)]

	tt = ROUND(n_elements(zz)/3.)
	zzSorted = zz[SORT(zz.zIP)]

	lowHalf = zz[WHERE(zz.zIP LE MEDIAN(zz.zIP))]
	highHalf = zz[WHERE(zz.zIP GT MEDIAN(zz.zIP))]

	lowThird = zz[WHERE(zz.zIP LE zzSorted[tt].zIP)]
	midThird = zz[WHERE((zz.zIP GT zzSorted[tt].zIP) AND (zz.zIP LE zzSorted[2*tt].zIP))]
	highThird = zz[WHERE(zz.zIP GT zzSorted[2*tt].zIP)]

;	zRanges = [0.2,MEDIAN(lowHalf.zIP),MEDIAN(zz.zIP),MEDIAN(highHalf.zIP),1.0]
	zRanges = [0.2, zzSorted[tt].zIP, zzSorted[2*tt].zIP, MAX(zzSorted.zIP)]	
;	zRanges = [0.2,0.35,0.5,0.65]
;	zRanges = [0.20,0.40,0.6,0.8,1.0]

;	colors = ['red','magenta','purple','blue','cyan']
;	PSYMS = [-1,-4,-5,-6,-7]
	colors = ['cyan', 'purple', 'orange', 'red']
	PSYMS = [-1,-4,-5,-6]
;	colors = ['blue', 'cyan', 'purple', 'orange']
;	PSYMS = [-4,-5,-6,-7]
	colors = ['cyan', 'purple', 'orange']
	PSYMS = [-1,-4,-5]

	leg = []
	FOR l=0,n_elements(zRanges)-2 DO BEGIN
		leg = [leg, decimal(zRanges[l],2) + '-' + decimal(zRanges[l+1],2)]
	ENDFOR

  !P.FONT=0
  IF (STRING(outputFormat) EQ 'ps') THEN BEGIN
	PS_OPEN, '~/latex/figures/IPlatefrac_vs_environ_'+massRanges[m_input]+'_zcut_M13massLim', /ENCAP, THICK=5
	DEVICE, /INCH, XS=6, YS=5, SET_FONT='Palatino-Roman'
  ENDIF ELSE BEGIN
	SET_PLOT, 'X'
  ENDELSE

;  FOR m=0,n_elements(massRanges)-1 DO BEGIN
  FOR m=0,n_elements(zRanges)-2 DO BEGIN
	zmin = zRanges[m]
	zmax = zRanges[m+1]
	neighData = neighData_allz[WHERE((neighData_allz.zIP GE zmin) AND (neighData_allz.zIP LE zmax))]

;	PRINT, n_elements(neighData)
;	PRINT, n_elements(neighData[WHERE(neighData.n_neighbors EQ 0L)])	
	Nmin = 0.
	Nmax = CEIL(100*MAX(neighData.n_neighbors))/100
;	Narray = dN*FINDGEN(1+CEIL(Nmax/dN))
;	Narray = [10^((1+FINDGEN(13))/6.)]
	Narray = [10^((FINDGEN(9))/4.)]

	Narray=Narray[0:-1]
;	IF ( m EQ 2 ) THEN Narray=Narray[0:-2]
	IF ( m EQ (n_elements(zRanges)-2) ) THEN Narray=Narray[0:-2]

;	PRINT, MAX(neighData.neigh_all_weight_total)

  lfRegionArray = []
  FOR j=0,MAX(neighData.JK_region) DO BEGIN
	lfArray = []
	neighDataJK = neighData[WHERE(neighData.JK_region NE j)]

; EQUAL WIDTH # NEIGHBOR BINS
	FOR i=0,n_elements(Narray)-2 DO BEGIN
	  Nweight_low	= Narray[i]
	  Nweight_high	= Narray[i+1]
	
	  Ndata = neighDataJK[WHERE( (neighDataJK.neigh_all_weight_total GE Nweight_low) AND $
				     (neighDataJK.neigh_all_weight_total LE Nweight_high), /NULL)]
;	  PRINT, j, Nweight_low, n_elements(Ndata)
	  IF (Ndata NE !NULL) THEN $
;	    lfPoint = TOTAL(Ndata.neigh_late_weight_total)/FLOAT(TOTAL(Ndata.neigh_all_weight_total)) $
	    lfPoint = TOTAL(Ndata[WHERE(Ndata.IP_type EQ 1)].IP_weight)/TOTAL(Ndata.IP_weight) $
	  ELSE lfPoint=0.

	  lfArray = [lfArray, lfPoint]
	ENDFOR
	lfRegionArray = [[lfRegionArray], [lfArray]]
  ENDFOR
	
;	PRINT, Narray[0:-2]+0.5*dN
;	PRINT, lfRegionArray

	errorArray = []
	FOR i=0,n_elements(Narray)-2 DO BEGIN
;		PRINT, TRANSPOSE(lfRegionArray[i,*])
		errorVec = lfRegionArray[i,*]
		errorArray = [errorArray, SQRT(MAX(neighData.JK_region)-1)*STDDEV(errorVec[WHERE(errorVec GE 0.)])]
	ENDFOR
;	PRINT, errorArray

	lfArrayAllRegions= []
	meanNweightArray= []
	medNweightArray = []
	N_IParray	= []
	FOR i=0,n_elements(Narray)-2 DO BEGIN
	  Nweight_low	= Narray[i]
	  Nweight_high	= Narray[i+1]
	
	  Ndata = neighData[WHERE( (neighData.neigh_all_weight_total GE Nweight_low) AND $
				   (neighData.neigh_all_weight_total LT Nweight_high), /NULL)]
	  IF (Ndata NE !NULL) THEN BEGIN
;	    lfPointAllRegions = TOTAL(Ndata.neigh_late_weight_total)/FLOAT(TOTAL(Ndata.neigh_all_weight_total))
	    lfPointAllRegions = TOTAL(Ndata[WHERE(Ndata.IP_type EQ 1)].IP_weight)/TOTAL(Ndata.IP_weight)
	    meanNweightArray = [meanNweightArray, MEAN(Ndata.neigh_all_weight_total)]
	    medNweightArray  = [medNweightArray, MEDIAN(Ndata.neigh_all_weight_total)]
	  ENDIF ELSE BEGIN
	    lfPointAllRegions=0.
	    meanNweightArray = [meanNweightArray, 0.]
	    medNweightArray  = [medNweightArray, 0.]
	  ENDELSE
	  lfArrayAllRegions = [lfArrayAllRegions, lfPointAllRegions]
	  N_IParray = [N_IParray, n_elements(Ndata)]
	ENDFOR

;	PLOT, [0,MAX(neighData.n_neighbors)], [0,0], LINESTYLE=1, xrange=[1,MAX(neighData.n_neighbors)], yrange=[ymin,ymax], /XLOG, $

      IF (m EQ 0) THEN BEGIN
	PLOT, [0,MAX(neighData.n_neighbors)], [0,0], LINESTYLE=1, xrange=[0.9,100], yrange=[ymin,ymax], /XLOG, $
;, POSITION=posArray[*,0], $
	  xtitle = textoidl('Neighbors [0.3 < R_{proj} < 4 Mpc]'), ytitle=textoidl('f_{late} (weighted)')
      ENDIF
;	LEGEND, massRanges, PSYM=4, NUMBER=2, PSPACING=1,  COLOR=cgColor(colors), BOX=0, /BOTTOM, /LEFT
;	LEGEND, [textoidl('log(M_{*IP}/M_{\odot})'),'11.0-11.5', '10.5-11.0', '10.0-10.5', '9.5-10.0'], $
;	LEGEND, [textoidl('IP z range'),'0.2-0.4','0.4-0.6','0.6-0.8','0.8-1.0'], $
	LEGEND, [textoidl('IP z range'),leg], $
		PSYM=[3,PSYMS], NUMBER=2, PSPACING=1,  COLOR=[cgColor('white'),cgColor(colors)], BOX=0, $; CHARSIZE=1.15, position=[1.1,0.475+0.25]
		/BOTTOM, /LEFT, CHARSIZE=1.15

;	IF (massRanges[m_input] EQ '9.5_10.0') THEN $
;	  XYOUTS, 80, ymin+0.05*yr, textoidl('9.5 < log(M_{*IP}/M_{\odot}) < 10.0'), ALIGNMENT=1.0

;	IF (massRanges[m_input] EQ '10.0_10.5') THEN BEGIN
	IF (massRanges[m_input] EQ '10.1_10.5') THEN BEGIN
;	  XYOUTS, 80, ymin+0.05*yr, textoidl('10.0 < log(M_{*IP}/M_{\odot}) < 10.5'), ALIGNMENT=1.0
	  XYOUTS, 80, ymin+0.05*yr, textoidl('10.1 < log(M_{*IP}/M_{\odot}) < 10.5'), ALIGNMENT=1.0
	  XYOUTS, 1.2, ymin+0.35*yr, 'IP galaxies only', CHARSIZE=1.15, ALIGNMENT=0.0
	  XYOUTS, 1.2, ymin+0.425*yr, textoidl('log(M_{*IP}/M_{\odot}) - log(M_{*neigh}/M_{\odot}) < 0.5'), ALIGNMENT=0.0, CHARSIZE=1.15
	ENDIF
	IF (massRanges[m_input] EQ '10.5_11.0') THEN BEGIN
	  XYOUTS, 80, ymin+0.05*yr, textoidl('10.5 < log(M_{*IP}/M_{\odot}) < 11.0'), ALIGNMENT=1.0
	  XYOUTS, 80, ymin+0.85*yr, 'IP galaxies only', CHARSIZE=1.15, ALIGNMENT=1.0
	  XYOUTS, 80, ymin+0.925*yr, textoidl('log(M_{*IP}/M_{\odot}) - log(M_{*neigh}/M_{\odot}) < 0.5'), ALIGNMENT=1.0, CHARSIZE=1.15
	ENDIF

;	IF (massRanges[m_input] EQ '11.0_11.5') THEN $
;	  XYOUTS, 80, ymin+0.05*yr, textoidl('11.0 < log(M_{*IP}/M_{\odot}) < 11.5'), ALIGNMENT=1.0

;	OPLOT, Narray[1:-1], lfArrayAllRegions, PSYM=4
;	OPLOT, Narray[1:-1], lfArrayAllRegions, PSYM=4
;	ERRPLOT, Narray[1:-1], lfArrayAllRegions+errorArray, lfArrayAllRegions-errorArray
	OPLOT, meanNweightArray, lfArrayAllRegions, PSYM=psyms[m], COLOR=cgColor(colors[m])
	ERRPLOT, meanNweightArray, lfArrayAllRegions+errorArray, lfArrayAllRegions-errorArray, COLOR=cgColor(colors[m])
;	OPLOT, medNweightArray, lfArrayAllRegions, PSYM=4;, COLOR=cgColor('green')
;	ERRPLOT, medNweightArray, lfArrayAllRegions+errorArray, lfArrayAllRegions-errorArray

;	XYOUTS, 1.25, ymin+0.05*yr, 'IP galaxies only', ALIGNMENT=0.0

;	PRINT, TRANSPOSE([[N_IParray], [Narray[1:-1]], [meanNweightArray], [medNweightArray]])
	PRINT, TRANSPOSE([[N_IParray], [Narray[1:-1]], [meanNweightArray], [lfArrayAllRegions], [errorArray]])
	PRINT, TOTAL(N_IParray)
  ENDFOR
	
  IF (STRING(outputFormat) EQ 'ps') THEN PS_CLOSE
  SET_PLOT, 'X'
	  
END
