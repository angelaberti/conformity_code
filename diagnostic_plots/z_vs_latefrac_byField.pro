PRO z_vs_latefrac_byField, outputFormat
        IF (string(outputFormat) eq 'ps') THEN BEGIN
		PS_OPEN, '~/figures/z_vs_latefrac_byField', THICK=3, /ENCAP
		DEVICE, /INCH, XS=8, YS=6
        ENDIF ELSE BEGIN
                SET_PLOT, 'X'
        ENDELSE

	binsize = 0.1

	A = FINDGEN(17) * (!PI*2/16.)
	USERSYM, 0.5*COS(A), 0.5*SIN(A), /FILL
	ptsymbol = 8

	!p.multi=0
	!p.charsize=1.25

	zmin = 0.2
	zmax = 1.0

	ymin = 0.2
	ymax = 1.0
	
	dataAll = mrdfits('~/results/zerodSFQ_all_cart.fits', 1)

	fracSFall= []
	zbins	= []

	FOR z_index=0,(zmax-zmin)/binsize-1 DO BEGIN
		zbin = mean([zmin+z_index*binsize, zmin+(z_index+1)*binsize])

		dataAllzbin	= dataAll[where( (dataAll.zprimus gt zmin+z_index*binsize) AND (dataAll.zprimus le zmin+(z_index+1)*binsize) )]
		nAllzbin 	= float(n_elements(dataAllzbin))

		fracSFzbin	= n_elements(dataAllzbin[where(dataAllzbin.SFQ eq 1)])/nAllzbin
		fracSFall	= [fracSFall, fracSFzbin]

		zbins = [zbins, zbin]
	ENDFOR

	fields = dataAll[uniq(dataAll.field, sort(dataAll.field))].field

	colors = ['Green', 'Cyan', 'Blue', 'Magenta', 'Red']

	PLOT, zbins, fracSFall, xrange=[0.2,1.0], yrange=[ymin,ymax], xtitle='Redshift', ytitle='Late-type fraction', title='Late-type fraction vs. redshift for entire sample', LINESTYLE=2, THICK=5
;	LEGEND, [strcompress(fields), 'entire sample', 'deviation from entire sample'], LINESTYLE=[0,0,0,0,0,2,1], THICK=[3,3,3,3,3,5,3], COLOR=[cgColor(colors),cgColor('Black'),cgColor('Black')], BOX=0, POSITION=[0.2, 0.6]
	LEGEND, [strcompress(fields), 'entire sample'], LINESTYLE=[0,0,0,0,0,2], THICK=[5,5,5,5,5,5], COLOR=[cgColor(colors),cgColor('Black')], BOX=0, /BOTTOM, /LEFT
	
	FOR f=0,4 DO BEGIN
		data = dataAll[where(dataAll.field EQ fields[f])]
;		dataSF = data[where(data.SFQ eq 1)]

		fracSF	= []

		FOR z_index=0,(zmax-zmin)/binsize-1 DO BEGIN
			zbin = mean([zmin+z_index*binsize, zmin+(z_index+1)*binsize])

			dataAllzbin	= data[where( (data.zprimus gt zmin+z_index*binsize) AND (data.zprimus le zmin+(z_index+1)*binsize) )]
			nAllzbin 	= float(n_elements(dataAllzbin))

			fracSFzbin	= n_elements(dataAllzbin[where(dataAllzbin.SFQ eq 1)])/nAllzbin
			fracSF		= [fracSF, fracSFzbin]
		ENDFOR

		OPLOT, zbins, fracSF, COLOR=cgColor(colors[f]), LINESTYLE=0, THICK=5
;		OPLOT, zbins, fracSFall-fracSF, COLOR=cgColor(colors[f]), LINESTYLE=1
	END

        IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
       	SET_PLOT, 'X'
END
