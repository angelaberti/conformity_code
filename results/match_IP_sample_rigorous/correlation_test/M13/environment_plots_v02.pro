PRO environment_plots_v02, outputFormat;, ptsPerBin, dR
	dN = 2
	ymin = 0.2
	ymax = 1.0
	yr = ymax-ymin

	ERASE
	!P.MULTI=4
	!P.CHARSIZE = 2.25
	charsz = 1.15
	dz_coeff = 2.0

	datapath = '~/results/match_IP_sample_rigorous/correlation_test/M13/neighborDataAll_'
	massRanges = ['10.7_11.0', '10.4_10.7', '10.1_10.4']
	masses = [11., 10.7, 10.4, 10.1]
	zmax = [1.0, 0.8, 0.65]

  !P.FONT=0
  IF (STRING(outputFormat) EQ 'ps') THEN BEGIN
	PS_OPEN, '~/latex/figures/environment_plots', /ENCAP, THICK=5
	DEVICE, /INCH, XS=10, YS=8, SET_FONT='Palatino-Roman'
  ENDIF ELSE BEGIN
	SET_PLOT, 'X'
  ENDELSE

;  colors = ['green', 'blue', 'magenta', 'red']
;  PSYMS = [1,4,5,6]
  colors = ['blue', 'magenta', 'red']
  PSYMS=[2,0,4]

  posArray = grid_array_deluxe(2,2,0,0)

;=============================================================================================================================================
; IP late-type fraction				
;=============================================================================================================================================

  signalLists = []
  errorLists  = []

  FOR m=0,n_elements(massRanges)-1 DO BEGIN
;	PRINT, massRanges[m], zmax[m]

	neighData = MRDFITS(datapath + massRanges[m] + '.fits', 1)
	neighData = neighData[WHERE(neighData.zIP LE zmax[m])]

	Nmin = 0.
	Nmax = CEIL(100*MAX(neighData.n_neighbors))/100
	Narray = [10^((FINDGEN(9))/4.)]

    lfRegionArray = []
    FOR j=0,MAX(neighData.JK_region) DO BEGIN
	lfArray = []
	neighDataJK = neighData[WHERE(neighData.JK_region NE j, /NULL)]

	FOR i=0,n_elements(Narray)-2 DO BEGIN
	  Nweight_low	= Narray[i]
	  Nweight_high	= Narray[i+1]
	
	  Ndata = neighDataJK[WHERE( (neighDataJK.neigh_all_weight_total GE Nweight_low) AND $
				     (neighDataJK.neigh_all_weight_total LE Nweight_high), /NULL)]
	  IF (Ndata NE !NULL) THEN $
	    lfPoint = TOTAL(Ndata[WHERE(Ndata.IP_type EQ 1)].IP_weight)/TOTAL(Ndata.IP_weight) $
	  ELSE lfPoint=0.

	  lfArray = [lfArray, lfPoint]
	ENDFOR
	lfRegionArray = [[lfRegionArray], [lfArray]]
    ENDFOR
	
	errorArray = []
	FOR i=0,n_elements(Narray)-2 DO BEGIN
		errorVec = lfRegionArray[i,*]
		errorArray = [errorArray, SQRT(MAX(neighData.JK_region)-1)*STDDEV(errorVec[WHERE(errorVec GE 0.)])]
	ENDFOR

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

	  xtickformat = '(A1)'
	  xtitle = ''
;	  xtitle = textoidl('Neighbors [0.3 < R_{proj} < 4 Mpc]')
;	  xtickformat = ''

    IF (m EQ 0) THEN BEGIN
	PLOT, meanNweightArray, lfArrayAllRegions, xtickformat=xtickformat, xtitle=xtitle, /NODATA, $
	  xrange=[0.9,100], yrange=[ymin,ymax], /XLOG, POSITION=posArray[*,2], ytitle=textoidl('Late-type Fraction of IPs')
;	LEGEND, [textoidl('log(M_{*IP}/M_{\odot} ) = [10.7, 11.0]; z = [0.2, 1.0]'), $
;		 textoidl('log(M_{*IP}/M_{\odot} ) = [10.4, 10.7]; z = [0.2, 0.8]'), $
;		 textoidl('log(M_{*IP}/M_{\odot} ) = [10.1, 10.4]; z = [0.2, 0.65]')], $
    ENDIF
	OPLOT, meanNweightArray, lfArrayAllRegions, LINESTYLE=PSYMS[m], SYMSIZE=2, COLOR=cgColor(colors[m])
	ERRPLOT, meanNweightArray, lfArrayAllRegions+errorArray, lfArrayAllRegions-errorArray, COLOR=cgColor(colors[m])

	signalLists = [[signalLists], [lfArrayAllRegions]]
	errorLists  = [[errorLists], [errorArray]]
  ENDFOR

;	XYOUTS, 1.25, 0.485, textoidl('IP late-type fraction decrease'), ALIGNMENT=0.0, CHARSIZE=charsz
;	XYOUTS, 1.25, 0.435, textoidl('from N_{neigh} < 3 to N_{neigh} > 30'), ALIGNMENT=0.0, CHARSIZE=charsz
;	LEGEND, [textoidl('0.17 (2.6\sigma)'), $
;		 textoidl('0.23 (3.3\sigma)'), $
;		 textoidl('0.14 (1.8\sigma)')], $
;		LINESTYLE=[PSYMS], NUMBER=2, PSPACING=2, COLOR=cgColor(colors), BOX=0, /BOTTOM, /LEFT, CHARSIZE=charsz
;		LINESTYLE=[2,2,0,0,4,4], NUMBER=2, PSPACING=2, COLOR=cgColor(['blue','white','magenta','white','red','white']), BOX=0, POSITION=[0.95,0.525], CHARSIZE=charsz

  FOR i=0,1 DO BEGIN
    signal = signalLists[*,i]-signalLists[*,(i+1)]
    noise = SQRT(errorLists[*,i]^2 + errorLists[*,(i+1)]^2)
;    PRINT, signalLists[*,i]
;    PRINT, signalLists[*,(i+1)]
    PRINT, ABS(signal/noise)
  ENDFOR

  PLOT, meanNweightArray, lfArrayAllRegions, xtickformat=xtickformat, ytickformat='(A1)', /NODATA, POSITION=posArray[*,1], COLOR=cgColor('white')
  LEGEND, [textoidl('10^{10.7} < M_{*IP} < 10^{11} M_{\odot}'), '0.2 < z < 1.0', '', $
	textoidl('10^{10.4} < M_{*IP} < 10^{10.7} M_{\odot}'), '0.2 < z < 0.8', '', $
	textoidl('10^{10.1} < M_{*IP} < 10^{10.4} M_{\odot}'), '0.2 < z < 0.65'], $
;	LINESTYLE=[PSYMS], NUMBER=2, PSPACING=3, COLOR=cgColor(colors), BOX=0, POSITION=[0.0, 0.75], CHARSIZE=charsz
;	LINESTYLE=[2,2,0,0,4], NUMBER=2,PSPACING=3, COLOR=cgColor(['blue','white','magenta','white','red']), BOX=0, POSITION=[0.0,0.75], CHARSIZE=charsz
;	LINESTYLE=[2,2,0,0,4,4], NUMBER=2, PSPACING=3, COLOR=cgColor(['blue','white','magenta','white','red','white']), /BOTTOM, /LEFT, CHARSIZE=charsz
	LINESTYLE=[2,2,2,0,0,0,4,4], NUMBER=2, PSPACING=3, COLOR=cgColor(['blue','white','white','magenta','white','white','red','white']), BOX=0, position=[0,0.8], CHARSIZE=1.25

;============================================================================================================================================
; Median IP Stellar Mass			
;=============================================================================================================================================

  ymin = 10.1
  ymax = 11.
  yr = ymax-ymin

  FOR m=0,n_elements(massRanges)-1 DO BEGIN

	neighData = MRDFITS(datapath + massRanges[m] + '.fits', 1)
	neighData = neighData[WHERE(neighData.zIP LE zmax[m])]

	Nmin = 0.
	Nmax = CEIL(100*MAX(neighData.n_neighbors))/100
	Narray = [10^((FINDGEN(9))/4.)]

    mstarRegionArray = []
    FOR j=0,MAX(neighData.JK_region) DO BEGIN
	mstarArray = []
	neighDataJK = neighData[WHERE(neighData.JK_region NE j, /NULL)]

	FOR i=0,n_elements(Narray)-2 DO BEGIN
	  Nweight_low	= Narray[i]
	  Nweight_high	= Narray[i+1]
	
	  Ndata = neighDataJK[WHERE( (neighDataJK.neigh_all_weight_total GE Nweight_low) AND $
				     (neighDataJK.neigh_all_weight_total LE Nweight_high), /NULL)]
	  IF (Ndata NE !NULL) THEN $
	    mstarPoint = MEDIAN(Ndata.IP_mstar) $
	  ELSE mstarPoint=0.

	  mstarArray = [mstarArray, mstarPoint]
	ENDFOR
	mstarRegionArray = [[mstarRegionArray], [mstarArray]]
    ENDFOR
	
	errorArray = []
	FOR i=0,n_elements(Narray)-2 DO BEGIN
		errorVec = mstarRegionArray[i,*]
		errorArray = [errorArray, SQRT(MAX(neighData.JK_region)-1)*STDDEV(errorVec[WHERE(errorVec GE 0.)])]
	ENDFOR

	mstarArrayAllRegions= []
	meanNweightArray= []
	medNweightArray = []
	N_IParray	= []
	FOR i=0,n_elements(Narray)-2 DO BEGIN
	  Nweight_low	= Narray[i]
	  Nweight_high	= Narray[i+1]
	
	  Ndata = neighData[WHERE( (neighData.neigh_all_weight_total GE Nweight_low) AND $
				   (neighData.neigh_all_weight_total LT Nweight_high), /NULL)]
	  IF (Ndata NE !NULL) THEN BEGIN
	    mstarPointAllRegions = MEDIAN(Ndata.IP_mstar)
	    meanNweightArray = [meanNweightArray, MEAN(Ndata.neigh_all_weight_total)]
	    medNweightArray  = [medNweightArray, MEDIAN(Ndata.neigh_all_weight_total)]
	  ENDIF ELSE BEGIN
	    mstarPointAllRegions=0.
	    meanNweightArray = [meanNweightArray, 0.]
	    medNweightArray  = [medNweightArray, 0.]
	  ENDELSE
	  mstarArrayAllRegions = [mstarArrayAllRegions, mstarPointAllRegions]
	  N_IParray = [N_IParray, n_elements(Ndata)]
	ENDFOR

;	  xtickformat = '(A1)'
;	  xtitle = ''
	  xtitle = textoidl('Neighbors [0.3 < R_{proj} < 4 Mpc]')
	  xtickformat = ''

    IF (m EQ 0) THEN BEGIN
	PLOT, meanNweightArray, mstarArrayAllRegions, xtickformat=xtickformat, xtitle=xtitle, /NODATA, $
	  xrange=[0.9,100], yrange=[ymin,ymax], /XLOG, POSITION=posArray[*,3], ytitle=textoidl('Median IP Stellar Mass')
    ENDIF
	OPLOT, meanNweightArray, mstarArrayAllRegions, LINESTYLE=PSYMS[m], SYMSIZE=2, COLOR=cgColor(colors[m])
	ERRPLOT, meanNweightArray, mstarArrayAllRegions+errorArray, mstarArrayAllRegions-errorArray, COLOR=cgColor(colors[m])

;	PRINT, TRANSPOSE([[N_IParray], [Narray[1:-1]], [meanNweightArray], [mstarArrayAllRegions], [errorArray]])

;    IF (m EQ n_elements(massRanges)-1) THEN BEGIN
;	XYOUTS, 80, ymin+0.90*yr, textoidl('log[M_{*IP}/M_{\odot} ] - log[M_{*neigh}/M_{\odot} ] < 0.5'), ALIGNMENT=1.0, CHARSIZE=charsz
;	XYOUTS, 80, ymin+0.825*yr, 'IP galaxies only', CHARSIZE=charsz, ALIGNMENT=1.0
;    ENDIF
;	XYOUTS, 1.2, ymin+0.125*yr, '0.2 < z < ' + decimal(zmax[m],2), CHARSIZE=charsz, ALIGNMENT=0.0
;	XYOUTS, 1.2, ymin+0.05*yr, textoidl(decimal(masses[m+1],1) + ' < log[M_{*IP}/M_{\odot} ] < ') + decimal(masses[m],1) , ALIGNMENT=0.0, CHARSIZE=charsz
  ENDFOR

;=============================================================================================================================================
; Mean IP sSFR					
;=============================================================================================================================================

  ymin = -10.5
  ymax = -9.5
  yr = ymax-ymin

  FOR m=0,n_elements(massRanges)-1 DO BEGIN

	neighData = MRDFITS(datapath + massRanges[m] + '.fits', 1)
	neighData = neighData[WHERE((neighData.zIP LE zmax[m]) AND (neighData.IP_type EQ 1))] ; SF IP only

	Nmin = 0.
	Nmax = CEIL(100*MAX(neighData.n_neighbors))/100
	Narray = [10^((FINDGEN(9))/4.)]

    sSFR_RegionArray = []
    FOR j=0,MAX(neighData.JK_region) DO BEGIN
	sSFR_Array = []
	neighDataJK = neighData[WHERE(neighData.JK_region NE j)]

	FOR i=0,n_elements(Narray)-2 DO BEGIN
	  Nweight_low	= Narray[i]
	  Nweight_high	= Narray[i+1]
	
	  Ndata = neighDataJK[WHERE( (neighDataJK.neigh_all_weight_total GE Nweight_low) AND $
				     (neighDataJK.neigh_all_weight_total LE Nweight_high), /NULL)]
	  IF (Ndata NE !NULL) THEN sSFR_Point = TOTAL((Ndata.IP_weight)*(Ndata.IP_SFR-Ndata.IP_mstar))/TOTAL(Ndata.IP_weight) ELSE sSFR_Point=0.
	  sSFR_Array = [sSFR_Array, sSFR_Point]
	ENDFOR
	sSFR_RegionArray = [[sSFR_RegionArray], [sSFR_Array]]
    ENDFOR
	
	errorArray = []
	FOR i=0,n_elements(Narray)-2 DO BEGIN
		errorVec = sSFR_RegionArray[i,*]
;		PRINT, errorVec
		errorArray = [errorArray, SQRT(MAX(neighData.JK_region)-1)*STDDEV(errorVec[WHERE(errorVec)])]
	ENDFOR

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
	    sSFR_PointAllRegions = TOTAL((Ndata.IP_weight)*(Ndata.IP_SFR-Ndata.IP_mstar))/TOTAL(Ndata.IP_weight)
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

;	  xtickformat = '(A1)'
;	  xtitle = ''
	  xtitle = textoidl('Neighbors [0.3 < R_{proj} < 4 Mpc]')
	  xtickformat = ''

    IF (m EQ 0) THEN BEGIN
	PLOT, meanNweightArray, sSFR_ArrayAllRegions, xtickformat=xtickformat, xtitle=xtitle, /NODATA, $
	  xrange=[0.9,100], yrange=[ymin,ymax], /XLOG, POSITION=posArray[*,0], ytitle=textoidl('Mean sSFR [yr^{-1}] for SF IPs')
    ENDIF
	OPLOT, meanNweightArray, sSFR_ArrayAllRegions, LINESTYLE=PSYMS[m], SYMSIZE=2, COLOR=cgColor(colors[m])
	ERRPLOT, meanNweightArray, sSFR_ArrayAllRegions+errorArray, sSFR_ArrayAllRegions-errorArray, COLOR=cgColor(colors[m])

;	PRINT, TRANSPOSE([[N_IParray], [Narray[1:-1]], [meanNweightArray], [sSFR_ArrayAllRegions], [errorArray]])
;	PRINT, errorArray

;    IF (m EQ n_elements(massRanges)-1) THEN BEGIN
;	XYOUTS, 80, ymin+0.90*yr, textoidl('log[M_{*IP}/M_{\odot} ] - log[M_{*neigh}/M_{\odot} ] < 0.5'), ALIGNMENT=1.0, CHARSIZE=charsz
;	XYOUTS, 80, ymin+0.825*yr, 'IP galaxies only', CHARSIZE=charsz, ALIGNMENT=1.0
;    ENDIF
;	XYOUTS, 1.2, ymin+0.125*yr, '0.2 < z < ' + decimal(zmax[m],2), CHARSIZE=charsz, ALIGNMENT=0.0
;	XYOUTS, 1.2, ymin+0.05*yr, textoidl(decimal(masses[m+1],1) + ' < log[M_{*IP}/M_{\odot} ] < ') + decimal(masses[m],1) , ALIGNMENT=0.0, CHARSIZE=charsz
  ENDFOR

;=============================================================================================================================================
	
  IF (STRING(outputFormat) EQ 'ps') THEN PS_CLOSE
  SET_PLOT, 'X'	  
END
