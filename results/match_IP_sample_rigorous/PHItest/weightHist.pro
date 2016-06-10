PRO weightHist, outputFormat
	dataAll = mrdfits('~/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits',1)
        fields = dataAll[uniq(dataAll.field, sort(dataAll.field))].field

;  FOR f=0,n_elements(fields)-1 DO BEGIN
;	data = dataIPallFields[WHERE(dataIPallFields.field EQ fields[f])]
;	PRINT, 'FIELD: ' + strcompress(fields[f])
;	PRINT, 'IPs: ' + strcompress(n_elements(data))
  
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/matchedIPsampleFBF/PHItest_weightHist, /ENCAP, THICK=3
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	!P.MULTI = 0
	!P.CHARSIZE = 1.1

	xmin = 1
	xmax = 12
	xr = xmax-xmin

	ymin = 0
	ymax = 50
	yr = ymax-ymin

	PHIarray = -1*(3.0 + 0.1*findgen(10))
        colors = ['red', 'orange', 'yellow', 'green', 'cyan', 'blue', 'purple', 'magenta', 'pink', 'gray']

  FOR p=0,n_elements(PHIarray)-1 DO BEGIN 
	stringPHI = strtrim(string(-1.*PHIarray[p], format='(f20.1)'), 2)

	data = MRDFITS('~/results/match_IP_sample_rigorous/PHItest/matchedIPsampleFBF_PHI' + stringPHI + '.fits', 1)

	SFIP = data[where(data.SFQ EQ 1)]
	objnamesUniq = SFIP[uniq(SFIP.objname, sort(SFIP.objname))].objname

	weights = []
	FOR n=0,n_elements(objnamesUniq)-1 DO weights = [weights, n_elements( SFIP[ where(SFIP.objname EQ objnamesUniq[n]) ] )]

;	freq = []
;	FOR w=2,MAX(weights) DO freq = [freq, w, n_elements(weights[where(weights EQ w)])]
;	PRINT, stringPHI, ' ', float(n_elements(weights[where(weights EQ 1)]))/n_elements(weights), freq
	PRINT, stringPHI, histogram(weights), n_elements(weights[where(weights GE 5)])/float(n_elements(weights)), $
		n_elements(weights[where(weights GE 7)])/float(n_elements(weights))
	IF (p EQ 0) THEN PLOTHIST, weights[where(weights GT 2)], color=cgcolor(colors[p]), xrange=[xmin,xmax], yrange=[ymin,ymax] $
	ELSE PLOTHIST, weights[where(weights GT 2)], color=cgcolor(colors[p]), /OVERPLOT
  ENDFOR

	LEGEND, [decimal(PHIarray, 2)], COLOR=cgcolor(colors[indgen(n_elements(PHIarray))]), LINESTYLE=0, THICK=3, BOX=0, /TOP, /RIGHT, CHARSIZE=1.0

END
