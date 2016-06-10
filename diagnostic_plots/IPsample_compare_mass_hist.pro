PRO IPsample_compare, var, outputFormat
; default var is mstar
	IF (string(outputFormat) eq 'ps') THEN BEGIN
		IF var EQ 'z' THEN BEGIN
			PS_OPEN, '../figures/IPsample_compare_z_hist', THICK=5, /ENCAP
	        ENDIF ELSE BEGIN
			PS_OPEN, '../figures/IPsample_compare_mass_hist', THICK=5, /ENCAP		
		ENDELSE
	        DEVICE, /INCH, XS=8, YS=8
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	!x.margin = [6,2.5] ; left, right
	!y.margin = [3.5,1] ; bottom, top
	!p.multi=[3,1,3]

	!p.charsize = 2

;	datapaths = [	'~/conformity/results/SF_status_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_variableMassCutoff.fits', $
;		    	'~/conformity/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits', $
;			'~/conformity/results/conservative+0.3dex/IP_data/zerodSFQ_IP_dz2.0_dm0.3.fits']
;
;	labels = ['Variable mass limit: M13 + 0.5 (0.0) dex for SF (Q) IP','Moustakas13 mass limit','Moustakas13 + 0.3 dex']

	string_dm = ['0.5', '0.8', '1.0']

	xtickformat = ['','','']
;	xtickformat = ['(A1)', '(A1)', '']
	IF var EQ 'z' THEN xtitle = ['redshift','redshift','redshift'] $
	ELSE xtitle = [textoidl('Log(M_{IP}/M_{*})'), textoidl('Log(M_{IP}/M_{*})'), textoidl('Log(M_{IP}/M_{*})')]
;	xtitle = ['', '', textoidl('Log(M_{IP}/M_{*})')]

	ERASE
;	FOR i=0,2 DO makeplot, datapaths[i], labels[i], xtickformat[i], xtitle[i]
	FOR i=0,2 DO makeplot, '~/conformity/results/variable_mass+' + string_dm[i] + 'dex/IP_data/zerodSFQ_IP_dz2.0.fits', $
		'Variable mass limit: M13 + ' + string_dm[i] + ' (0.0) dex for SF (Q) IP', $
		xtickformat[i], xtitle[i], var

        IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
       	SET_PLOT, 'X'
END

PRO makeplot, datapath, label, xtickformat, xtitle, var
	data = mrdfits(datapath, 1)

	dataIP = data[where(data.IP eq 1)]
	dataIPSF = dataIP[where(dataIP.SFQ eq 1)]
	dataIPQ  = dataIP[where(dataIP.SFQ eq 0)]

	bin = 0.05

	xmin = 8.5
	xmax = 12
	IF var EQ 'z' THEN BEGIN
		bin = 0.02

		xmin = 0.2
		xmax = 1.0
	ENDIF
	xr = xmax-xmin

	ymin = 0
	ymax = 1000
	yr = ymax-ymin

	IF var EQ 'z' THEN BEGIN
		plothist, dataIP.zprimus, bin=bin, xrange=[xmin,xmax], yrange=[ymin,ymax], xtickformat=xtickformat, xtitle=xtitle, /nodata
		plothist, dataIPSF.zprimus, bin=bin, /overplot, color=cgColor('Blue')
		plothist, dataIPQ.zprimus, bin=bin, /overplot, color=cgColor('Red')
	ENDIF ELSE BEGIN
		plothist, dataIP.mstar, bin=bin, xrange=[xmin,xmax], yrange=[ymin,ymax], xtickformat=xtickformat, xtitle=xtitle, /nodata
		plothist, dataIPSF.mstar, bin=bin, /overplot, color=cgColor('Blue')
		plothist, dataIPQ.mstar, bin=bin, /overplot, color=cgColor('Red')
	ENDELSE

	medianMasses = textoidl('Median M_{IP}: ') + strtrim(string(median(dataIPSF.mstar), format='(f20.2)'),1) + ' (SF) ' + $
	  strtrim(string(median(dataIPQ.mstar), format='(f20.2)'),1) + ' (Q)'

	XYOUTS, xmin+0.05*xr, ymin+0.90*yr, label, ALIGNMENT=0.0, charsize=1
	XYOUTS, xmin+0.05*xr, ymin+0.80*yr, medianMasses, ALIGNMENT=0.0, charsize=1
END
