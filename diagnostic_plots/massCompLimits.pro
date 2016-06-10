PRO massCompLimits, outputFormat
	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '../figures/massCompLimits', /ENCAP, THICK=5
	        DEVICE, /INCH, XS=8, YS=6, LANGUAGE_LEVEL=2
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

		data=mrdfits('~/conformity/results/zerodSFQ_all_cart.fits',1)
		
		!p.multi=0

;	        !x.margin=[6,2.5] ; left, right
;	        !y.margin=[3.5,2.5] ; bottom, top

		!p.charsize=1.1

		xmin = 0.2
		xmax = 1.0
		xr = xmax-xmin

		ymin = 8
		ymax = 12
		yr = ymax-ymin

		IF outputFormat eq 'ps' THEN axisColor = cgColor('Black') ELSE axisColor = cgColor('White')
		Qcolor  = cgColor('Red')
		SFcolor = cgColor('Blue')

	        mass_binwidth 	= 0.15
	        z_binwidth	= 0.05
		
		xtitle = 'Redshift'	
		ytitle = textoidl('Log (M / M_{\odot})')

		zmin = 0.2
		zmax = 1.0
		zArray = zmin + z_binwidth*findgen((zmax-zmin)/z_binwidth + 1)

	PLOT, data[where(data.SFQ eq 0)].zprimus, data[where(data.SFQ eq 0)].mstar, psym=3, color=axisColor, $
		xrange=[xmin, xmax], yrange=[ymin, ymax], xtitle=xtitle, ytitle=ytitle, title='Mass Completeness Limits for IP Samples', /NODATA
;	OPLOT, data[where(data.SFQ eq 0)].zprimus, data[where(data.SFQ eq 0)].mstar, psym=3, color=Qcolor
;	OPLOT, data[where(data.SFQ eq 1)].zprimus, data[where(data.SFQ eq 1)].mstar, psym=3, color=SFcolor

	FOR SFstatus=0,1 DO BEGIN
		IF SFstatus EQ 1 THEN SFQcolor=SFcolor ELSE SFQcolor=Qcolor
;		SFQcolor = cgColor('Gray')
		dataSFQ = data[where(data.SFQ eq SFstatus)]
	        mass	= dataSFQ.mstar
	        z	= dataSFQ.zprimus
		
	        n_massbins = ceil((max(mass)-min(mass))/mass_binwidth)
	        n_zbins    = floor((max(z)-min(z))/z_binwidth)

	        hist2d = []
		FOR s_index=0,n_zbins DO BEGIN
       		        z_low  = min(z)+z_binwidth*float(s_index)
	        	z_high = min(z)+z_binwidth*float(s_index+1)

			in_z_range = dataSFQ[where((dataSFQ.zprimus ge z_low) AND (dataSFQ.zprimus le z_high), /NULL)]

	                z_row = []
	                FOR m_index=0,n_massbins DO BEGIN
	                        m_low   = min(mass)+mass_binwidth*float(m_index)
	                        m_high  = min(mass)+mass_binwidth*float(m_index+1)

	                        IF in_z_range ne !NULL THEN $
	                        in_unit = n_elements(in_z_range[where((in_z_range.mstar ge m_low) AND (in_z_range.mstar le m_high), /NULL)]) $
	                        ELSE in_unit = 0

	                        z_row = [z_row, in_unit]
	                ENDFOR

	                hist2d = [[hist2d], [z_row]]
	        ENDFOR

		masspts	= min(mass)+mass_binwidth*float(findgen(n_massbins+1))
		zpts	= min(z)+z_binwidth*float(findgen(n_zbins+1))

		PRINT, minmax(hist2d)

		CONTOUR, TRANSPOSE(hist2d), zpts+0.5*z_binwidth, masspts+0.5*mass_binwidth, levels=float(max(hist2d))*findgen(8)/7., $
			color=SFQcolor, c_linestyle=1, c_thick=2, /OVERPLOT

;		mlist = 8.+findgen(5)
;		OPLOT, mlist, -1.29+0.65*(mlist-10.)+1.33*(zmean-0.1), linestyle=2, color=axisColor, thick=linethick
;		OPLOT, mlist, -0.49+0.65*(mlist-10.)+1.07*(zmean-0.1), linestyle=2, color=axisColor, thick=linethick
	ENDFOR

        zbins = [0.20, 0.30, 0.40, 0.50, 0.65, 0.80, 1.00]

        zbinsMidpts = [0.2]
        FOR i=1,n_elements(zbins)-3 DO zbinsMidpts = [zbinsMidpts, [mean([zbins[i],zbins[i+1]])]]
	zbinsMidpts = [zbinsMidpts,[1.0]]

;	default_mean = [9.172,9.492,9.792,10.138,10.508,10.886]
	ETmean = [9.440, 9.738, 10.012, 10.312, 10.610, 10.878]
	LTmean = [9.090, 9.420, 9.712, 10.026, 10.328, 10.576]
	
	colors = ['Orange', 'Magenta', 'Purple', 'Cyan', 'Green']
	dm = [0.0,0.5,0.8,1.0]

	OPLOT, zbinsMidpts, ETmean, COLOR=cgColor(colors[0]), THICK=3
	FOR i=0,3 DO OPLOT, zbinsMidpts, LTmean+dm[i], COLOR=cgColor(colors[i+1]), THICK=3

        LEGEND, ['Early-type IP; M13 mass limit', $
                'Late-type IP; M13 mass limit', $
                'Late-type IP; M13 mass limit + 0.5 dex', $
                'Late-type IP; M13 mass limit + 0.8 dex', $
                'Late-type IP; M13 mass limit + 1.0 dex'], $
                 LINESTYLE=[0,0,0,0,0], color=cgColor(colors), BOX=0, /RIGHT, /BOTTOM, THICK=3

	LEGEND, ['Complete ET sample', 'Complete LT sample'], LINESTYLE=[1,1], COLOR=[Qcolor, SFcolor], THICK=2, BOX=0, /TOP, /LEFT

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

