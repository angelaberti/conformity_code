PRO SFR_mass_contour

	data=mrdfits("../results/zerodSFQ_all_cart.fits",1)

	!p.multi=0

	mass	= data.mstar
	sfr	= data.sfr

	mass_binwidth = 0.1
	sfr_binwidth  = 0.1

	n_massbins = ceil((max(mass)-min(mass))/mass_binwidth)
	n_sfrbins  = floor((max(sfr)-min(sfr))/sfr_binwidth)

;	PRINT, minmax(mass), minmax(sfr)
;	PRINT, "mass bins:", n_massbins, " sfr bins:", n_sfrbins

	hist2d = []
	FOR s_index=0,n_sfrbins DO BEGIN
		sfr_low  = min(sfr)+sfr_binwidth*float(s_index)
		sfr_high = min(sfr)+sfr_binwidth*float(s_index+1)

		in_sfr_range = data[where((data.sfr ge sfr_low) AND (data.sfr le sfr_high), /null)]
	
;		PRINT, "low SFR:", sfr_low, " high SFR:", sfr_high, " in range:", n_elements(in_sfr_range)
		
		sfr_row = []
		FOR m_index=0,n_massbins DO BEGIN
			m_low 	= min(mass)+mass_binwidth*float(m_index)
			m_high  = min(mass)+mass_binwidth*float(m_index+1)
	
			IF in_sfr_range ne !null THEN $
			in_unit = n_elements(in_sfr_range[where((in_sfr_range.mstar ge m_low) AND (in_sfr_range.mstar le m_high), /null)]) $
			ELSE in_unit = 0

;			PRINT, "low mass:", m_low, " high mass:", m_high, " in range:", in_unit		

			sfr_row = [sfr_row, in_unit]
		ENDFOR

		hist2d = [[hist2d], [sfr_row]]
	ENDFOR

;	PRINT, min(mass)+mass_binwidth*float(findgen(n_massbins+1))
;	PRINT, min(sfr)+sfr_binwidth*float(findgen(n_sfrbins+1))

;	PLOT, mass, sfr, psym=3

	CONTOUR, hist2d, min(mass)+mass_binwidth*(findgen(n_massbins+1)), $
			 min(sfr)+sfr_binwidth*(findgen(n_sfrbins+1)), $
			 nlevels=10, color=cgColor('Red'), c_thick=3;, /overplot
END
