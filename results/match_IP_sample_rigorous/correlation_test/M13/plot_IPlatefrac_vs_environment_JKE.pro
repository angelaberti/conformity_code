PRO plot_IPlatefrac_vs_environment_JKE, outputFormat;, ptsPerBin, dR
	dN = 2
	ymin = 0.3
	ymax = 1.0
	yr = ymax-ymin

	ERASE
	!P.MULTI=3
	!P.CHARSIZE = 2.25
	dz_coeff = 2.0

	datapath = '~/conformity/results/match_IP_sample_rigorous/correlation_test/M13/neighborDataAll_'
	massRanges = ['10.7_11.0', '10.4_10.7', '10.1_10.4']
	masses = [11., 10.7, 10.4, 10.1]
	zmax = [1.0, 0.8, 0.65]

  !P.FONT=0
  IF (STRING(outputFormat) EQ 'ps') THEN BEGIN
	PS_OPEN, '~/latex/figures/IPlatefrac_vs_environ_M13massLim', /ENCAP, THICK=5

	DEVICE, /INCH, XS=6, YS=9, SET_FONT='Palatino-Roman'
  ENDIF ELSE BEGIN
	SET_PLOT, 'X'
  ENDELSE

;  colors = ['green', 'blue', 'magenta', 'red']
;  PSYMS = [1,4,5,6]
  colors = ['blue', 'magenta', 'red']
  PSYMS=[-3,-3,-3]
;  PSYMS = [4,5,6]
;  colors = ['blue', 'magenta']
;  PSYMS = [-4,-5]

  posArray = grid_array(1,3)
;  PRINT, posArray

  FOR m=0,n_elements(massRanges)-1 DO BEGIN
	PRINT, massRanges[m], zmax[m]

	neighData = MRDFITS(datapath + massRanges[m] + '.fits', 1)
	neighData = neighData[WHERE(neighData.zIP LE zmax[m])]

;	PRINT, n_elements(neighData)
;	PRINT, n_elements(neighData[WHERE(neighData.n_neighbors EQ 0L)])	
	Nmin = 0.
	Nmax = CEIL(100*MAX(neighData.n_neighbors))/100
;	Narray = dN*FINDGEN(1+CEIL(Nmax/dN))
;	Narray = [0,0.5,1.0,10^((1+FINDGEN(13))/6.)]
;	Narray = [10^((FINDGEN(9))/4.)]
	Narray = [0,10,10^1.5,100]
	PRINT, Narray

  lfRegionArray = []
  FOR j=0,MAX(neighData.JK_region) DO BEGIN
	lfArray = []
	neighDataJK = neighData[WHERE(neighData.JK_region NE j, /NULL)]

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

;	PLOT, [0,MAX(neighData.n_neighbors)], [0,0], LINESTYLE=1, xrange=[1,MAX(neighData.n_neighbors)+1], yrange=[ymin,ymax], /XLOG, $

;	LEGEND, massRanges, PSYM=4, NUMBER=2, PSPACING=1,  COLOR=cgColor(colors), BOX=0, /BOTTOM, /LEFT
;	LEGEND, [textoidl('log[M_{*IP}/M_{\odot} ]'),'11.0-11.5', '10.5-11.0', '10.0-10.5', '9.5-10.0'], $
;	LEGEND, [textoidl('log[M_{*IP}/M_{\odot} ] range'), '10.7-11.0', '10.4-10.7', '10.1-10.4'], $
;		PSYM=[3,PSYMS], SYMSIZE=[1,1,1], NUMBER=2, PSPACING=2, COLOR=[cgColor('white'),cgColor(colors)], BOX=0, /BOTTOM, /LEFT, CHARSIZE=1.1
;		PSYM=[3,PSYMS], NUMBER=2, PSPACING=2, COLOR=[cgColor('white'),cgColor(colors)], BOX=0, /BOTTOM, /LEFT, CHARSIZE=1.1

	IF (m GT 0) THEN BEGIN
	  xtickformat = '(A1)'
	  xtitle = ''
	ENDIF ELSE BEGIN
	  xtitle = textoidl('Neighbors {0.3 < R_{proj} < 4 Mpc}')
	ENDELSE

	PLOT, meanNweightArray, lfArrayAllRegions, PSYM=psyms[m], SYMSIZE=2, xtickformat=xtickformat, xtitle=xtitle, $
	  xrange=[0.9,100], yrange=[ymin,ymax], /XLOG, POSITION=posArray[*,m], ytitle=textoidl('Late-type Fraction of IPs')
;	  xtitle = textoidl('Neighbors {0.3 < R_{proj} < 4 Mpc}'), ytitle=textoidl('Late-type Fraction of IPs')
	ERRPLOT, meanNweightArray, lfArrayAllRegions+errorArray, lfArrayAllRegions-errorArray

	PRINT, lfArrayAllRegions, errorArray
	PRINT, ABS(lfArrayAllRegions[0] - lfArrayAllRegions[-1])
	PRINT, ABS(lfArrayAllRegions[0] - lfArrayAllRegions[-1])/SQRT(errorArray[0]^2+errorArray[-1]^2)	

    IF (m EQ n_elements(massRanges)-1) THEN BEGIN
	XYOUTS, 80, ymin+0.90*yr, textoidl('log[M_{*IP}/M_{\odot} ] - log[M_{*neigh}/M_{\odot} ] < 0.5'), ALIGNMENT=1.0, CHARSIZE=1
	XYOUTS, 80, ymin+0.825*yr, 'IP galaxies only', CHARSIZE=1, ALIGNMENT=1.0
    ENDIF
	XYOUTS, 1.2, ymin+0.125*yr, '0.2 < z < ' + decimal(zmax[m],2), CHARSIZE=1.1, ALIGNMENT=0.0
	XYOUTS, 1.2, ymin+0.05*yr, textoidl(decimal(masses[m+1],1) + ' < log[M_{*IP}/M_{\odot} ] < ') + decimal(masses[m],1) , ALIGNMENT=0.0, CHARSIZE=1.1
  ENDFOR
	
  IF (STRING(outputFormat) EQ 'ps') THEN PS_CLOSE
  SET_PLOT, 'X'

;  PRINT, Narray
;  PRINT, meanNweightArray
	  
END
