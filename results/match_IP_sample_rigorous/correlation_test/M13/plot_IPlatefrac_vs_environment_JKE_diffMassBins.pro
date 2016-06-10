PRO plot_IPlatefrac_vs_environment_JKE_diffMassBins, outputFormat;, ptsPerBin, dR
	dN = 2
	ymin = 0.3
	ymax = 1.0
	yr = ymax-ymin

	ERASE
	!P.MULTI=2
	!P.CHARSIZE = 2.25
	dz_coeff = 2.0

	datapath = '~/conformity/results/match_IP_sample_rigorous/correlation_test/M13/neighborDataAll_'
	massRanges = ['10.7_11.0', '10.4_10.7', '10.1_10.4']
	masses = [11., 10.7, 10.4, 10.1]
	zmax = [0.8, 0.65]

	low  = MRDFITS(datapath + massRanges[2] + '.fits', 1)
	med  = MRDFITS(datapath + massRanges[1] + '.fits' ,1)
	high = MRDFITS(datapath + massRanges[0] + '.fits', 1)

	low_med  = [low, med]
	med_high = [med, high]

	PRINT, n_elements(low_med), n_elements(med_high)

	low_med  = low_med[WHERE(low_med.zIP LE 0.65)]
	med_high = med_high[WHERE(med_high.zIP LE 0.8)]

	PRINT, n_elements(low_med), n_elements(med_high)

	massRanges = ['10.4_11.0', '10.1_10.7']

  !P.FONT=0
  IF (STRING(outputFormat) EQ 'ps') THEN BEGIN
	PS_OPEN, '~/latex/figures/IPlatefrac_vs_environ_M13massLim_diffMassBins', /ENCAP, THICK=5

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

  posArray = grid_array(1,2)

  FOR m=0,n_elements(massRanges)-1 DO BEGIN
	PRINT, massRanges[m]

	IF (m EQ 0) THEN neighData = med_high
	IF (m EQ 1) THEN neighData = low_med

	Nmin = 0.
	Nmax = CEIL(100*MAX(neighData.n_neighbors))/100
;	Narray = dN*FINDGEN(1+CEIL(Nmax/dN))
;	Narray = [0,0.5,1.0,10^((1+FINDGEN(13))/6.)]
	Narray = [10^((FINDGEN(9))/4.)]
;	Narray = [0,10,10^1.5,100]
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
	XYOUTS, 1.2, ymin+0.05*yr, textoidl(decimal(masses[m+2],1) + ' < log[M_{*IP}/M_{\odot} ] < ') + decimal(masses[m],1) , ALIGNMENT=0.0, CHARSIZE=1.1
  ENDFOR
	
  IF (STRING(outputFormat) EQ 'ps') THEN PS_CLOSE
  SET_PLOT, 'X'
	  
END
