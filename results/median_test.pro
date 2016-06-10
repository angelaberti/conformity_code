PRO median_test
	string_dm = '0.5'
	topLevel = '~/conformity/results/variable_mass+' + string_dm + 'dex/'

	RESTORE, topLevel + 'latefrac_annulus_array_IPSF.sav'
	latefrac_array_IPSF = IP_latefrac_array
	RESTORE, topLevel + 'latefrac_annulus_array_IPQ.sav'
	latefrac_array_IPQ  = IP_latefrac_array

;	PLOTHIST, findgen(10), bin=0.05, yrange=[0,1], /NODATA
	PLOT, findgen(10), xrange=[0,n_elements(latefrac_array_IPSF[0,*])], yrange=[0.4,1], /NODATA

	colors = ['Red', 'Orange', 'Yellow', 'Green', 'Cyan', 'Blue', 'Purple', 'Magenta']

;	FOR c=4,11 DO PLOTHIST, getcol(c,latefrac_array_IPSF), bin=0.001, color=cgColor(colors[c-4]), /OVERPLOT
	FOR c=4,11 DO OPLOT, getcol(c,latefrac_array_IPSF), color=cgColor(colors[c-4])
	FOR c=4,11 DO OPLOT, getcol(c,latefrac_array_IPQ), LINESTYLE=1, color=cgColor(colors[c-4])
END

FUNCTION getcol, c, latefrac_array

	col = latefrac_array[c,sort(latefrac_array[c,*])] 
	colpos = col[where(col GE 0)]

	PRINT, median(colpos)

	RETURN, colpos
END
