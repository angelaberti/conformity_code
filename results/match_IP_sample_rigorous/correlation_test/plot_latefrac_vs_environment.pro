PRO plot_latefrac_vs_environment, dN, outputFormat;, ptsPerBin, dR
	ymin = 0.5
	ymax = 1
	yr = ymax-ymin

	ERASE
	!P.MULTI=1
;	posArray = grid_array(1,2)
;	!P.MULTI=[2,1,2]
;	dR = 0.1 ; Mpc
	dz_coeff = 2.0

	neighData = MRDFITS('~/conformity/results/match_IP_sample_rigorous/correlation_test/neighborDataAll_10.0_10.5.fits', 1)

;	PRINT, MINMAX(neighData.n_neighbors)
	Nmin = 0.
	Nmax = CEIL(100*MAX(neighData.n_neighbors))/100
	PRINT, Nmax
	Narray = dN*FINDGEN(1+CEIL(Nmax/dN))
;	Narray = [0,4,8,12,16,20,30]
	PRINT, Narray

	lf = []
	lfErrorArray = []
; EQUAL WIDTH # NEIGHBOR BINS
	FOR i=0,n_elements(Narray)-2 DO BEGIN
	  Nlow	= Narray[i]
	  Nhigh	= Narray[i+1]
	
	  Ndata = neighData[WHERE((neighData.n_neighbors GE Nlow) AND (neighData.n_neighbors LT Nhigh), /NULL)]
;	  PRINT, Nlow, n_elements(Ndata)
	  IF (Ndata NE !NULL) THEN $
	    lfPoint =  TOTAL(Ndata.n_neigh_late)/FLOAT(TOTAL(Ndata.n_neighbors)) $
	  ELSE lfPoint=-99

	  Ndata_90array = []
	  IF (lfPoint NE -99) THEN BEGIN
	    FOR s=1,200 DO BEGIN
	      seed = s
	      random_indices = ROUND( RANDOMU( seed, ROUND( 0.9*n_elements(Ndata) ) )*n_elements(Ndata) )
	      Ndata_90 = Ndata[random_indices]
	      Ndata_90array = [Ndata_90array, TOTAL(Ndata_90.n_neigh_late)/FLOAT(TOTAL(Ndata_90.n_neighbors))]
	    ENDFOR
	    binError = STDDEV(Ndata_90array)
	  ENDIF ELSE BEGIN
	    binError = 0
	  ENDELSE

	  lf = [lf, lfPoint]
	  lfErrorArray = [lfErrorArray, binError]
	END
	
	PRINT, TRANSPOSE([[Narray[0:-2]+0.5*dN], [lf], [lfErrorArray]])

  !P.FONT=0
  IF (STRING(outputFormat) EQ 'ps') THEN BEGIN
	PS_OPEN, '~/latex/figures/latefrac_vs_environ', /ENCAP, THICK=5
	DEVICE, /INCH, XS=6, YS=5, SET_FONT='Palatino-Roman'
  ENDIF ELSE BEGIN
	SET_PLOT, 'X'
  ENDELSE
	PLOT, [0,MAX(neighData.n_neighbors)], [0,0], LINESTYLE=1, xrange=[0,MAX(neighData.n_neighbors)], yrange=[ymin,ymax], $
;, POSITION=posArray[*,0], $
	  xtitle = textoidl('Number of Neighbors (0.3<R_{proj}<4 Mpc)'), ytitle=textoidl('f_{late} for All Neigh. in # Neigh. Range')

;	XYOUTS, 0.25, ymin+0.1*yr, 'Equal radial binwidth'
	OPLOT, Narray[0:-2]+0.5*dN, lf, PSYM=4
	ERRPLOT, Narray[0:-2]+0.5*dN, lf+lfErrorArray, lf-lfErrorArray
	XYOUTS, 2, ymin+0.05*yr, textoidl('10<log(M_{*IP})<10.5'),ALIGNMENT=0.0
	XYOUTS, 2, ymin+0.15*yr, textoidl('log(M_{*IP})-log(M_{*neigh})<0.5')

  IF (STRING(outputFormat) EQ 'ps') THEN PS_CLOSE
  SET_PLOT, 'X'
	  
END
