PRO SFR_vs_mass_pres, outputFormat
	!P.FONT=0

	IF (string(outputFormat) EQ 'ps') THEN BEGIN
	        PS_OPEN, '~/latex/figures/SFR_vs_mass_pres', THICK=5, /ENCAP
	        DEVICE, /INCH, XS=6, YS=6, SET_FONT='Palatino-Roman'
	ENDIF ELSE BEGIN
	        SET_PLOT, 'X'
	ENDELSE

	A = 0.15
	th = FINDGEN(17)*(!PI*2/16.)
	USERSYM, A*COS(th), A*SIN(th), /FILL

	SYM=8	
		data=mrdfits('~/results/zerodSFQ_all_cart.fits',1)
		
;		!P.FONT = 0
		!p.multi=[0,3,2]

;	        !x.margin=[6,2.5] ; left, right
;	        !y.margin=[3.5,2.5] ; bottom, top

		!p.charsize=2.25
		charsz=1.25

		xmin = 8.05
		xmax = 12
		xr = xmax-xmin

		ymin = -2.8
		ymax = 2.5
		yr = ymax-ymin

		IF outputFormat EQ 'ps' THEN axisColor = cgColor('Black') ELSE axisColor = cgColor('White')
		Qcolor  = cgColor('RED4')
		SFcolor = cgColor('BLU4')
;		Qcolor = cgColor('gray')
;		SFcolor = cgColor('gray')
		zArray=[0.2, 0.3, 0.4, 0.5, 0.65, 0.8, 1.0]

	        mass_binwidth = 0.15
	        sfr_binwidth  = 0.15
		
		nx=3.
		ny=2.

		POSITIONarray=grid_array_deluxe(3,2,0,0)

		xm=0.12
		ym=0.12
		
		dx = (1 - 1.5*xm)/nx
        	dy = (1 - 2.*ym)/ny

		FOR i=0,n_elements(zArray)-2 DO BEGIN
			IF i LE 2 THEN ii=i+3 ELSE ii=i-3
			IF (ii LE (nx-1)) THEN BEGIN
;				xtitle = textoidl('log (M_{stellar}/M_{\odot})')
				xtitle = textoidl('log ( M_{*}/M'+sunsymbol()+' )')
				xtickformat = ''
			ENDIF ELSE BEGIN
				xtitle = ''
				xtickformat = '(A1)'
			ENDELSE
			IF (ii mod nx EQ 0) THEN BEGIN
;				ytitle = textoidl('log (SFR / M_{\odot} yr^{-1})')	
				ytitle = textoidl('log ( SFR / M') + sunsymbol() + textoidl(' yr^{-1} )')	
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

			data_zcut = data[where( (data.zprimus GE zArray[i]) AND $
						(data.zprimus LE zArray[i+1]) )]

;			IF i EQ nx*ny-5 THEN title='Stellar Mass vs. Star Formation Rate' ELSE title=''

			PLOT, data_zcut[where(data_zcut.SFQ EQ 0)].mstar, data_zcut[where(data_zcut.SFQ EQ 0)].SFR, $
;			  color=axisColor, xrange=[xmin, xmax], yrange=[ymin, ymax], position=pos, $
			  color=axisColor, xrange=[xmin, xmax], yrange=[ymin, ymax], position=POSITIONarray[*,(i+3) MOD 6], $
			  xtitle=xtitle, ytitle=ytitle, xtickformat=xtickformat, ytickformat=ytickformat, title=title, /NODATA
			OPLOT, data_zcut[where(data_zcut.SFQ EQ 0)].mstar, data_zcut[where(data_zcut.SFQ EQ 0)].SFR, PSYM=SYM, color=Qcolor
			OPLOT, data_zcut[where(data_zcut.SFQ EQ 1)].mstar, data_zcut[where(data_zcut.SFQ EQ 1)].SFR, PSYM=SYM, color=SFcolor
			IF (i LT 3) THEN $
			  XYOUTS, xmin+0.10*xr, ymin+0.875*yr, textoidl('z = ')+strtrim(string(zArray[i],format='(f20.2)'),1)+'-'+strtrim(string(zArray[i+1],format='(f20.2)'),1), color=axisColor, CHARSIZE=charsz $
			ELSE $
			  XYOUTS, xmin+0.10*xr, ymin+0.05*yr, textoidl('z = ')+strtrim(string(zArray[i],format='(f20.2)'),1)+'-'+strtrim(string(zArray[i+1],format='(f20.2)'),1), color=axisColor, CHARSIZE=charsz
;			IF (i EQ 5) THEN $
;			  XYOUTS, 10, -2.5, textoidl('-1.29+0.65(M_{*}-10)+1.33(z-0.1)'), CHARSIZE=0.8, ALIGNMENT=0.0
;			  LEGEND, [textoidl('log(SFR/M_{\odot} yr^{-1}) ='),textoidl('-1.29+0.65 log(M_{*}-10)+1.33(z-0.1)')], PSYM=[], COLOR=cgcolor('White'), CHARSIZE=0.8, BOX=0, /BOTTOM, /LEFT

		        mass    = data_zcut.mstar
		        sfr     = data_zcut.sfr
		
		        n_massbins = ceil((max(mass)-min(mass))/mass_binwidth)
		        n_sfrbins  = floor((max(sfr)-min(sfr))/sfr_binwidth)

		        hist2d = []
		        FOR s_index=0,n_sfrbins DO BEGIN
        		        sfr_low  = min(sfr)+sfr_binwidth*float(s_index)
       		        	sfr_high = min(sfr)+sfr_binwidth*float(s_index+1)

	                in_sfr_range = data_zcut[where((data_zcut.sfr GE sfr_low) AND (data_zcut.sfr LE sfr_high), /null)]

	                sfr_row = []
	                FOR m_index=0,n_massbins DO BEGIN
	                        m_low   = min(mass)+mass_binwidth*float(m_index)
	                        m_high  = min(mass)+mass_binwidth*float(m_index+1)

	                        IF in_sfr_range ne !null THEN $
	                        in_unit = n_elements(in_sfr_range[where((in_sfr_range.mstar GE m_low) AND (in_sfr_range.mstar LE m_high), /null)]) $
	                        ELSE in_unit = 0

	                        sfr_row = [sfr_row, in_unit]
	                ENDFOR

	                hist2d = [[hist2d], [sfr_row]]
	        ENDFOR

		masspts = min(mass)+mass_binwidth*float(findgen(n_massbins+1))
		sfrpts  = min(sfr)+sfr_binwidth*float(findgen(n_sfrbins+1))

		PRINT, minmax(hist2d)

;	        CONTOUR, hist2d, masspts, sfrpts, nlevels=10, color=axisColor, c_thick=linethick, /OVERPLOT
		CONTOUR, hist2d, masspts, sfrpts, levels=float(max(hist2d))*findgen(8)/7.,color=axisColor, /OVERPLOT

		mlist = 8.+findgen(5)
		OPLOT, mlist, -1.29+0.65*(mlist-10.)+1.33*(zmean-0.1), linestyle=2, color=axisColor
;		OPLOT, mlist, -0.49+0.65*(mlist-10.)+1.07*(zmean-0.1), linestyle=2, color=axisColor
	END

	IF (string(outputFormat) EQ 'ps') THEN PS_CLOSE
	SET_PLOT, 'X'
END

