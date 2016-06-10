PRO plot_meanSSFR_vs_environment_JKE, outputFormat;, ptsPerBin, dR
	dN = 2
	ymin = -10.5
	ymax = -9.5
	yr = ymax-ymin

	ERASE
	!P.MULTI=1
	!P.CHARSIZE = 1.25
	dz_coeff = 2.0

	datapath = '~/conformity/results/match_IP_sample_rigorous/correlation_test/neighborDataAll_'
;	massRanges = ['11.0_11.5', '10.5_11.0', '10.0_10.5', '9.5_10.0']
	massRanges = ['10.5_11.0', '10.0_10.5']

  !P.FONT=0
  IF (STRING(outputFormat) EQ 'ps') THEN BEGIN
	PS_OPEN, '~/latex/figures/meanSSFR_vs_environ', /ENCAP, THICK=5

	DEVICE, /INCH, XS=6, YS=5, SET_FONT='Palatino-Roman'
  ENDIF ELSE BEGIN
	SET_PLOT, 'X'
  ENDELSE

;  colors = ['green', 'blue', 'magenta', 'red']
;  PSYMS = [1,4,5,6]
  colors = ['blue', 'magenta']
  PSYMS = [-4,-5]

  FOR m=0,n_elements(massRanges)-1 DO BEGIN
	neighData = MRDFITS(datapath + massRanges[m] + '.fits', 1)
	neighData = neighData[WHERE(neighData.IP_type EQ 1)] ; SF IP only
;	PRINT, n_elements(neighData)
;	PRINT, n_elements(neighData[WHERE(neighData.n_neighbors EQ 0L)])	
	Nmin = 0.
	Nmax = CEIL(100*MAX(neighData.n_neighbors))/100
;	Narray = dN*FINDGEN(1+CEIL(Nmax/dN))
;	Narray = [10^((1+FINDGEN(13))/6.)]
	Narray = [10^((FINDGEN(9))/4.)]

;	IF ( m EQ (n_elements(massRanges)-1) ) THEN Narray=Narray[0:-4]
;	Narray=Narray[0:-4]

  sSFR_RegionArray = []
  FOR j=0,MAX(neighData.JK_region) DO BEGIN
	sSFR_Array = []
	neighDataJK = neighData[WHERE(neighData.JK_region NE j)]

	FOR i=0,n_elements(Narray)-2 DO BEGIN
	  Nweight_low	= Narray[i]
	  Nweight_high	= Narray[i+1]
	
	  Ndata = neighDataJK[WHERE( (neighDataJK.neigh_all_weight_total GE Nweight_low) AND $
				     (neighDataJK.neigh_all_weight_total LE Nweight_high), /NULL)]
;	  PRINT, j, Nweight_low, n_elements(Ndata)
;	  IF (Ndata NE !NULL) THEN sSFR_Point = MEAN(Ndata.IP_SFR-Ndata.IP_Mstar) ELSE sSFR_Point=0.
	  IF (Ndata NE !NULL) THEN sSFR_Point = TOTAL((Ndata.IP_weight)*(Ndata.IP_SFR-Ndata.IP_Mstar))/TOTAL(Ndata.IP_weight) ELSE sSFR_Point=0.
	  sSFR_Array = [sSFR_Array, sSFR_Point]
	ENDFOR
	sSFR_RegionArray = [[sSFR_RegionArray], [sSFR_Array]]
  ENDFOR
	
;	PRINT, Narray[0:-2]+0.5*dN
	PRINT, sSFR_RegionArray

	errorArray = []
	FOR i=0,n_elements(Narray)-2 DO BEGIN
;		PRINT, TRANSPOSE(sSFR_RegionArray[i,*])
		errorVec = sSFR_RegionArray[i,*]
;		errorArray = [errorArray, SQRT(MAX(neighData.JK_region)-1)*STDDEV(errorVec[WHERE(errorVec GE 0.)])]
		errorArray = [errorArray, SQRT(MAX(neighData.JK_region)-1)*STDDEV(errorVec)]
	ENDFOR
;	PRINT, errorArray

	sSFR_ArrayAllRegions= []
	meanNweightArray= []
	medNweightArray = []
	N_IParray	= []
	FOR i=0,n_elements(Narray)-2 DO BEGIN
	  Nweight_low	= Narray[i]
	  Nweight_high	= Narray[i+1]
	
	  Ndata = neighData[WHERE( (neighData.neigh_all_weight_total GE Nweight_low) AND $
				   (neighData.neigh_all_weight_total LT Nweight_high), /NULL)]
	  IF (Ndata NE !NULL) THEN BEGIN
;	    sSFR_PointAllRegions = MEAN(Ndata.IP_SFR-Ndata.IP_Mstar)
	    sSFR_PointAllRegions = TOTAL((Ndata.IP_weight)*(Ndata.IP_SFR-Ndata.IP_Mstar))/TOTAL(Ndata.IP_weight)
	    meanNweightArray = [meanNweightArray, MEAN(Ndata.neigh_all_weight_total)]
	    medNweightArray  = [medNweightArray, MEDIAN(Ndata.neigh_all_weight_total)]
	  ENDIF ELSE BEGIN
	    sSFR_PointAllRegions=0.
	    meanNweightArray = [meanNweightArray, 0.]
	    medNweightArray  = [medNweightArray, 0.]
	  ENDELSE
	  sSFR_ArrayAllRegions = [sSFR_ArrayAllRegions, sSFR_PointAllRegions]
	  N_IParray = [N_IParray, n_elements(Ndata)]
	ENDFOR

;	PLOT, [0,MAX(neighData.n_neighbors)], [0,0], LINESTYLE=1, xrange=[1,MAX(neighData.n_neighbors)+1], yrange=[ymin,ymax], /XLOG, $

      IF (m EQ 0) THEN BEGIN
	PLOT, [0,MAX(neighData.n_neighbors)], [0,0], LINESTYLE=1, xrange=[0.9,100], yrange=[ymin,ymax], /XLOG, $
	  xtitle = textoidl('Neighbors [0.3 < R_{proj} < 4 Mpc]'), ytitle=textoidl('Mean sSFR for SF IP Galaxies')

;	LEGEND, massRanges, PSYM=4, NUMBER=2, PSPACING=1,  COLOR=cgColor(colors), BOX=0, /BOTTOM, /LEFT
;	LEGEND, [textoidl('log(M_{*IP}/M_{\odot})'),'11.0-11.5', '10.5-11.0', '10.0-10.5', '9.5-10.0'], $
	LEGEND, [textoidl('log(M_{*IP}/M_{\odot})'), '10.5-11.0', '10.0-10.5'], $
		PSYM=[3,PSYMS], NUMBER=2, PSPACING=1, COLOR=[cgColor('white'),cgColor(colors)], BOX=0, /BOTTOM, /LEFT, CHARSIZE=1.15

	XYOUTS, 80, ymin+0.925*yr, textoidl('log(M_{*IP}/M_{\odot}) - log(M_{*neigh}/M_{\odot}) < 0.5'), ALIGNMENT=1.0, CHARSIZE=1.15
	XYOUTS, 80, ymin+0.85*yr, 'IP galaxies only', CHARSIZE=1.15, ALIGNMENT=1.0
      ENDIF

;	OPLOT, Narray[1:-1], sSFR_ArrayAllRegions, PSYM=4
;	ERRPLOT, Narray[1:-1], sSFR_ArrayAllRegions+errorArray, sSFR_ArrayAllRegions-errorArray
	OPLOT, meanNweightArray, sSFR_ArrayAllRegions, PSYM=psyms[m], COLOR=cgColor(colors[m])
	ERRPLOT, meanNweightArray, sSFR_ArrayAllRegions+errorArray, sSFR_ArrayAllRegions-errorArray, COLOR=cgColor(colors[m])
;	OPLOT, medNweightArray, sSFR_ArrayAllRegions, PSYM=4;, COLOR=cgColor('green')
;	ERRPLOT, medNweightArray, sSFR_ArrayAllRegions+errorArray, sSFR_ArrayAllRegions-errorArray

;	XYOUTS, 1.25, ymin+0.05*yr, textoidl('10<log(M_{*IP})<10.5'),ALIGNMENT=0.0

;	PRINT, TRANSPOSE([[N_IParray], [Narray[1:-1]], [meanNweightArray], [medNweightArray]])
	PRINT, TRANSPOSE([[N_IParray], [Narray[1:-1]], [meanNweightArray], [sSFR_ArrayAllRegions], [errorArray]])
	PRINT, TOTAL(N_IParray)
  ENDFOR
	
  IF (STRING(outputFormat) EQ 'ps') THEN PS_CLOSE
  SET_PLOT, 'X'

  PRINT, Narray
	  
END

