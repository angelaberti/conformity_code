PRO num_density
	data 	= mrdfits('~/conformity/results/zerodSFQ_all_cart.fits', 1)
	fields 	= data[uniq(data.field, sort(data.field))].field
	
	total_pyr_height = dcomovinglos(1., /Mpc)/(1 + 1.)

	basetotal = 0
	FOREACH field,fields DO BEGIN
		dataf 	= data[where(data.field EQ field)]
		base 	= ABS(MAX(dataf.xprop) - MIN(dataf.xprop))*ABS(MAX(dataf.yprop) - MIN(dataf.yprop))
		basetotal += base
	ENDFOREACH
	
	total_vol	= basetotal*total_pyr_height/3.
	density 	= n_elements(data)/total_vol

	PRINT, 'Total volume (cubic Mpc): ', total_vol
	PRINT, 'Galaxies per cubic Mpc: ', density
	PRINT, 'Number density at...'
	PRINT, 'z = 0.2: ', 30./cylvol(0.2)
	PRINT, 'z = 0.6: ', 30./cylvol(0.6)
	PRINT, 'z = 1.0: ', 30./cylvol(1.0)
END

FUNCTION cylvol, z
	dz		= 2.0*0.005*(1 + z)
	cylheight	= ABS(dcomovinglos(z+dz, /Mpc)-dcomovinglos(z-dz, /Mpc))/(1 + z)
	cylbase		= !PI*4.^2
	cylvol		= cylbase*cylheight
	RETURN, cylvol
END
