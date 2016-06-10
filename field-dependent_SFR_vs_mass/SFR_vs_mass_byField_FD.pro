PRO SFR_vs_mass_byField_FD, field, outputFormat
	dataAll = mrdfits('~/field-dependent_SFR_vs_mass/zerodSFQ_all.fits',1)

	string_field = strtrim(strcompress(field),2)
	data = dataAll[where(strtrim(strcompress(dataAll.field),2) EQ string_field)]		
	PRINT, string_field + ': ' + string(n_elements(data))

	IF (string(outputFormat) eq 'ps') THEN BEGIN
	        PS_OPEN, '~/figures/SFR_vs_mass_' + string_field + '_FD', THICK=3, /ENCAP
	        DEVICE, /INCH, XS=8, YS=6
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

;		!p.multi=0
		!p.multi=[0,3,2]

;	        !x.margin=[6,2.5] ; left, right
;	        !y.margin=[3.5,2.5] ; bottom, top

		!p.charsize=1.2
		sym = 8

		A = FINDGEN(17)*(!PI*2/16.)
		R = 0.25
		USERSYM, R*COS(A), R*SIN(A), /FILL
		
		xmin = 8.05
		xmax = 12
		xr = xmax-xmin

		ymin = -3
		ymax = 2.5
		yr = ymax-ymin

		IF outputFormat eq 'ps' THEN axisColor = cgColor('Black') ELSE axisColor = cgColor('White')
		Qcolor  = cgColor('Red')
		SFcolor = cgColor('Blue')

		linethick = 3

		zArray=[0.2, 0.3, 0.4, 0.5, 0.65, 0.8, 1.0]

	        mass_binwidth = 0.15
	        sfr_binwidth  = 0.15
		
		nx=3.
		ny=2.

		xm=0.08
		ym=0.08
		
		dx = (1 - 1.5*xm)/nx
        	dy = (1 - 2.*ym)/ny

		FOR i=0,n_elements(zArray)-2 DO BEGIN
			IF i le 2 THEN ii=i+3 ELSE ii=i-3
			IF (ii le (nx-1)) THEN BEGIN
				xtitle = textoidl('Log (M_{*} / M_{\odot})')
				xtickformat = ''
			ENDIF ELSE BEGIN
				xtitle = ''
				xtickformat = '(A1)'
			ENDELSE
			IF (ii mod nx eq 0) THEN BEGIN
				ytitle = textoidl('Log (SFR / M_{\odot} yr^{-1})')	
				ytickformat = ''
			ENDIF ELSE BEGIN
				ytitle = ''
				ytickformat = '(A1)'
			ENDELSE

;		        pos = [ xm + dx*float(i mod nx), $
;        		        ym + dy*float(floor(float(i)/nx)), $
;      			        xm + dx*(1. + float(i mod nx)), $
;        		        ym + dy*(1. + float(floor(float(i)/nx))) ]
		        pos = [ xm + dx*float(i mod nx), $
        		        ym + dy*float(floor(float(ii)/nx)), $
      			        xm + dx*(1. + float(i mod nx)), $
        		        ym + dy*(1. + float(floor(float(ii)/nx))) ]
			PRINT, pos

			zmean = mean([zArray[i],zArray[i+1]])

			data_zcut = data[where( (data.zprimus ge zArray[i]) AND $
						(data.zprimus le zArray[i+1]) )]

			IF i eq nx*ny-5 THEN title='Stellar Mass vs. Star Formation Rate: '+string_field ELSE title=''

			PLOT, data_zcut[where(data_zcut.SFQ eq 0)].mstar, data_zcut[where(data_zcut.SFQ eq 0)].SFR, psym=sym, $
			  color=axisColor, xrange=[xmin, xmax], yrange=[ymin, ymax], position=pos, charsize=2, $
			  xtitle=xtitle, ytitle=ytitle, xtickformat=xtickformat, ytickformat=ytickformat, title=title, /NODATA
			OPLOT, data_zcut[where(data_zcut.SFQ eq 0)].mstar, data_zcut[where(data_zcut.SFQ eq 0)].SFR, psym=sym, color=Qcolor
			OPLOT, data_zcut[where(data_zcut.SFQ eq 1)].mstar, data_zcut[where(data_zcut.SFQ eq 1)].SFR, psym=sym, color=SFcolor
			XYOUTS, xmin+0.10*xr, ymin+0.90*yr, textoidl('z=')+strtrim(string(zArray[i],format='(f20.2)'),1)+'-'+strtrim(string(zArray[i+1],format='(f20.2)'),1), color=axisColor

		        mass    = data_zcut.mstar
		        sfr     = data_zcut.sfr
		
		        n_massbins = ceil((max(mass)-min(mass))/mass_binwidth)
		        n_sfrbins  = floor((max(sfr)-min(sfr))/sfr_binwidth)

		        hist2d = []
		        FOR s_index=0,n_sfrbins DO BEGIN
        		        sfr_low  = min(sfr)+sfr_binwidth*float(s_index)
       		        	sfr_high = min(sfr)+sfr_binwidth*float(s_index+1)

	                in_sfr_range = data_zcut[where((data_zcut.sfr ge sfr_low) AND (data_zcut.sfr le sfr_high), /null)]

	                sfr_row = []
	                FOR m_index=0,n_massbins DO BEGIN
	                        m_low   = min(mass)+mass_binwidth*float(m_index)
	                        m_high  = min(mass)+mass_binwidth*float(m_index+1)

	                        IF in_sfr_range ne !null THEN $
	                        in_unit = n_elements(in_sfr_range[where((in_sfr_range.mstar ge m_low) AND (in_sfr_range.mstar le m_high), /null)]) $
	                        ELSE in_unit = 0

	                        sfr_row = [sfr_row, in_unit]
	                ENDFOR

	                hist2d = [[hist2d], [sfr_row]]
	        ENDFOR

		masspts = min(mass)+mass_binwidth*float(findgen(n_massbins+1))
		sfrpts  = min(sfr)+sfr_binwidth*float(findgen(n_sfrbins+1))

		PRINT, minmax(hist2d)

;	        CONTOUR, hist2d, masspts, sfrpts, nlevels=10, color=axisColor, c_thick=linethick, /OVERPLOT
		CONTOUR, hist2d, masspts, sfrpts, levels=float(max(hist2d))*findgen(8)/7.,color=axisColor, c_thick=linethick, /OVERPLOT

		mlist = 8.+findgen(5)
		IF field EQ 'cdfs' THEN offset=0.2 ELSE $
		IF field EQ 'cfhtls_xmm' THEN offset = -0.2 ELSE offset=0.0
		OPLOT, mlist, -1.29+offset+0.65*(mlist-10.)+1.33*(zmean-0.1), linestyle=2, color=axisColor, thick=linethick
;		OPLOT, mlist, -0.49+0.65*(mlist-10.)+1.07*(zmean-0.1), linestyle=2, color=axisColor, thick=linethick
	END

	IF (string(outputFormat) eq 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

