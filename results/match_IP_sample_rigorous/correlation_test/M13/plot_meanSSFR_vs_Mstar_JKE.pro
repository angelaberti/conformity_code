PRO plot_meanSSFR_vs_Mstar_JKE, outputFormat;, ptsPerBin, dR
	dm = 0.2
	zmax = 0.65
	dN = 2

	xmin = 10.2
	xmax = 11
	xr = xmax-xmin

	ymin = -11.1
	ymax = -9.8
	yr = ymax-ymin

	ERASE
	!P.MULTI=1
	!P.CHARSIZE = 1.5
	charsz = 1.25

	datapath = '~/results/match_IP_sample_rigorous/correlation_test/M13/neighborDataAll_'
	massRanges = ['10.7_11.0', '10.4_10.7', '10.1_10.4']

	neighData = []
	FOR m=0,2 DO BEGIN
	 neighDataPart = MRDFITS(datapath + massRanges[m] + '.fits', 1)
	 neighData = [neighData, neighDataPart]
	ENDFOR

;	neighData = neighData[WHERE((neighData.IP_type EQ 1) AND (neighData.zIP LE zmax))] ; SF IP only
	neighData = neighData[WHERE(neighData.zIP LE zmax)]

  !P.FONT=0
  IF (STRING(outputFormat) EQ 'ps') THEN BEGIN
	PS_OPEN, '~/latex/figures/meanSSFR_vs_Mstar_M13massLim_allIP', /ENCAP, THICK=5
	DEVICE, /INCH, XS=6, YS=5, SET_FONT='Palatino-Roman'
  ENDIF ELSE BEGIN
	SET_PLOT, 'X'
  ENDELSE

  colors = ['black', 'purple', 'orange']
;  colors = ['white', 'purple', 'orange']
;  colors = ['green', 'cyan', 'blue', 'purple', 'magenta', 'red', 'orange']
;  colors = ['blue', 'magenta', 'red']
  PSYMS = [0,2,3,4,5,6]

;	Narray = [1,3,6,12,24,48,100]
;	Narray = [0,4,10,30,100]
;	Narray = [0,5,25,100]
	Narray = [0,10, 10^1.5,100]

  legend = [textoidl('N_{neigh} < ') + decimal(Narray[1],1), decimal(Narray[1],1) + textoidl(' < N_{neigh} < ') + decimal(Narray[2],1), textoidl('N_{neigh} > ') + decimal(Narray[2],1)]
;  legend = [textoidl('N_{neigh} < ') + decimal(Narray[1],1), textoidl('N_{neigh} > ') + decimal(Narray[1],1)]

  yLists = []
  errorLists = []

  FOR i=0,n_elements(Narray)-2 DO BEGIN
    Nweight_low	= Narray[i]
    Nweight_high = Narray[i+1]
    Ndata = neighData[WHERE( (neighData.neigh_all_weight_total GE Nweight_low) AND $
			   (neighData.neigh_all_weight_total LE Nweight_high), /NULL)]
    mmin = MIN(Ndata.IP_weight*Ndata.IP_mstar)
    mmax = MAX(Ndata.IP_weight*Ndata.IP_mstar)

;    massArray = 0.1*FLOOR(10*mmin) + dm*INDGEN((mmax-mmin)/dm)
;    massArray = massArray[WHERE(massArray LE 11.)]
    massArray = [10.1, 10.4, 10.7, 11.]
;    PRINT, Narray[i], 'massArray: ', massArray	

; ACTUAL DATA POINTS
    xvec = []
    yvec = []
    FOR k=0,n_elements(massArray)-2 DO BEGIN
      Mdata = Ndata[WHERE((Ndata.IP_weight*Ndata.IP_mstar GT massArray[k]) AND (Ndata.IP_weight*Ndata.IP_mstar LE massArray[k+1]), /NULL)]
      PRINT, k, massArray[k], MINMAX(Mdata.IP_weight*Mdata.IP_mstar), n_elements(Mdata)
      xvec = [xvec, MEDIAN([Mdata.IP_weight*Mdata.IP_mstar])]
      yvec = [yvec, MEDIAN([Mdata.IP_weight*(Mdata.IP_SFR-Mdata.IP_mstar)])]
;      xvec = [xvec, TOTAL(Mdata.IP_weight*Mdata.IP_mstar)/TOTAL(Mdata.IP_weight)]
;      yvec = [yvec, TOTAL(Mdata.IP_weight*(Mdata.IP_SFR-Mdata.IP_mstar))/TOTAL(Mdata.IP_weight)]
    ENDFOR

; JACKKNIFE ERRORS ON DATA POINTS
    yvec_RegionArray = []
    FOR j=0,MAX(neighData.JK_region) DO BEGIN
      NdataJK = neighData[WHERE(neighData.JK_region NE j)]
;      PRINT, j, n_elements(NdataJK)
 
      yvec_Region = []
      FOR k=0,n_elements(massArray)-2 DO BEGIN
	MdataJK = NdataJK[WHERE((NdataJK.IP_weight*NdataJK.IP_mstar GT massArray[k]) AND (NdataJK.IP_weight*NdataJK.IP_mstar LE massArray[k+1]), /NULL)]
;	PRINT, k, massArray[k], MINMAX(MdataJK.IP_weight*MdataJK.IP_mstar), n_elements(MdataJK)
 	yvec_Region = [yvec_Region, MEDIAN([MdataJK.IP_weight*(MdataJK.IP_SFR-MdataJK.IP_mstar)])]
      ENDFOR
      yvec_RegionArray = [[yvec_RegionArray], [yvec_Region]]
    ENDFOR
	
    errorArray = []
    FOR r=0,n_elements(massArray)-2 DO BEGIN
;      PRINT, TRANSPOSE(yvec_RegionArray[r,*])
      errorVec = yvec_RegionArray[r,*]
      errorArray = [errorArray, SQRT(MAX(neighData.JK_region)-1)*STDDEV(errorVec)]
    ENDFOR
;    PRINT, errorArray

    IF (i EQ 0) THEN PLOT, xvec, yvec, xrange=[xmin,xmax], yrange=[ymin,ymax], /NODATA, XMINOR=2, $
	xtitle = textoidl('log (M_{*IP}/M_{\odot} )'), ytitle = textoidl('Median sSFR [yr^{-1}] for all IPs')
    OPLOT, xvec, yvec, LINESTYLE=PSYMS[i], COLOR=cgColor(colors[i])
    ERRPLOT, xvec, yvec+errorArray, yvec-errorArray, COLOR=cgColor(colors[i])

    PRINT, yvec[0]-yvec[-1]
    PRINT, ABS(yvec[0]-yvec[-1])/SQRT(errorArray[0]^2+errorArray[-1]^2)

    yLists = [[yLists], [yvec]]
    errorLists = [[errorLists], [errorArray]]
  ENDFOR

  LEGEND, legend, LINESTYLE=PSYMS[0:n_elements(Narray)-2], COLOR=cgColor(colors[0:n_elements(Narray)-2]), /TOP, /RIGHT, BOX=0, NUMBER=2, PSPACING=2.5, CHARSIZE=charsz

;  XYOUTS, 80, ymin+0.925*yr, textoidl('log(M_{*IP}/M_{\odot} ) - log(M_{*neigh}/M_{\odot} ) < 0.5'), ALIGNMENT=1.0;, CHARSIZE=1.15
;  XYOUTS, 80, ymin+0.875*yr, 'IP galaxies only', CHARSIZE=1.15, ALIGNMENT=1.0
  XYOUTS, xmin+0.05*xr, ymin+0.05*yr, '0.2 < z < 0.65', ALIGNMENT=0.0

  FOR i=0,0 DO BEGIN
    signal = yLists[*,i]-yLists[*,(i+1)]
    noise = SQRT(errorLists[*,i]^2 + errorLists[*,(i+1)]^2)
    PRINT, signal
    PRINT, ABS(signal/noise)
  ENDFOR

  IF (STRING(outputFormat) EQ 'ps') THEN PS_CLOSE
  SET_PLOT, 'X'	  
END
