PRO IPsample_compare_zhist, outputFormat
	IF (string(outputFormat) eq 'ps') THEN BEGIN
		PS_OPEN, '../figures/IPsample_compare_zhist', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=11, YS=8.5
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

;	!x.margin = [6,2.5] ; left, right
;	!y.margin = [3.5,1] ; bottom, top
;	!p.multi=[4,2,2]
	!p.multi=0

	!p.charsize = 1.3

	datapath = '~/results/conservative_mass_cutoff/IP_data/zerodSFQ_IP_dz2.0_dm0.0.fits'

        data = mrdfits(datapath, 1)

        dataIP = data[where(data.IP eq 1)]
        dataIPSF = dataIP[where(dataIP.SFQ eq 1)]
        dataIPQ  = dataIP[where(dataIP.SFQ eq 0)]

        bin = 0.02              

        xmin = 0.2
        xmax = 1.0

        ymin = 0
        ymax = 800

        xr = xmax-xmin
        yr = ymax-ymin

	PLOTHIST, dataIPQ.zprimus, bin=bin, color=cgColor('Orange'), xrange=[xmin,xmax], yrange=[ymin,ymax], xtitle='Redshift', ytitle='Number', title='Redshift Distribution of Different IP Samples'
        PLOTHIST, dataIPSF.zprimus, bin=bin, color=cgColor('Magenta'), /OVERPLOT

	medianMasses = [strtrim(string(median(dataIPQ.mstar), format='(f20.2)'),1), strtrim(string(median(dataIPSF.mstar), format='(f20.2)'),1)]
	medianRedshifts = [strtrim(string(median(dataIPQ.zprimus), format='(f20.2)'),1), strtrim(string(median(dataIPSF.zprimus), format='(f20.2)'),1)]

	string_dm = ['0.5', '0.8', '1.0']
	colors = ['Purple', 'Cyan', 'Green']
	
  FOR i=0,2 DO BEGIN	
	datapath = '~/results/variable_mass+' + string_dm[i] + 'dex/IP_data/zerodSFQ_IP_dz2.0.fits'
	data 	 = mrdfits(datapath, 1)

	dataIP 	 = data[where(data.IP eq 1)]
	dataIPSF = dataIP[where(dataIP.SFQ eq 1)]
	dataIPQ  = dataIP[where(dataIP.SFQ eq 0)]

	PLOTHIST, dataIPSF.zprimus, bin=bin, color=cgColor(colors[i]), /OVERPLOT
;	plothist, dataIPQ.zprimus, bin=bin, color=cgColor('Red')

	medianMasses = [medianMasses, strtrim(string(median(dataIPSF.mstar), format='(f20.2)'),1)]
	medianRedshifts = [medianRedshifts, strtrim(string(median(dataIPSF.zprimus), format='(f20.2)'),1)]

;	XYOUTS, xmin+0.95*xr, ymin+0.90*yr, label, ALIGNMENT=1.0, charsize=1
;	XYOUTS, xmin+0.95*xr, ymin+0.80*yr, medianMasses, ALIGNMENT=1.0, charsize=1
;	XYOUTS, xmin+0.95*xr, ymin+0.70*yr, medianRedshifts, ALIGNMENT=1.0, charsize=1
  END

	LEGEND, ['Early-type IP; M13 mass limit', $
		'Late-type IP; M13 mass limit', $
		'Late-type IP; M13 mass limit + 0.5 dex', $
		'Late-type IP; M13 mass limit + 0.8 dex', $
		'Late-type IP; M13 mass limit + 1.0 dex'] + textoidl('; Med M_{IP}: ') + medianMasses + '; Med z: ' + medianRedshifts, $
		 LINESTYLE=[0,0,0,0,0], color=[cgColor('Orange'), cgColor('Magenta'), cgColor(colors[0]), cgColor(colors[1]), cgColor(colors[2])], BOX=0, /RIGHT, /TOP

        IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
       	SET_PLOT, 'X'
END
