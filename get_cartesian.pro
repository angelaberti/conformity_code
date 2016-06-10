PRO get_cartesian, path

	inputFile = path + 'zerodSFQ_all.fits'
	zerodInput = mrdfits(inputFile, 1)

	fields = zerodInput[uniq(zerodInput.field, sort(zerodInput.field))].field
	print, fields
	n_fields = n_elements(fields)

;	templateRow = create_struct(zerodInput[0], 'dx', 0.0d, 'dy', 0.0d, 'dz', 0.0d)
	templateRow = create_struct(zerodInput[0], 'dx', 0.0d, 'dy', 0.0d, 'dz', 0.0d, 'xprop', 0.0d, 'yprop', 0.0d)
	;help, templateRow
	;print, templateRow
	zerodOutput = replicate(templateRow, n_elements(zerodInput))
	outputIndex = 0L

	FOR f=0,(n_fields-1) DO BEGIN ; loop through all fields
		field = fields[f]
		currentField = zerodInput[ where( zerodInput.field eq field, /null ) ]
		n_gal = size(currentField, /n_elements)
; loop through all galaxies in field
		FOR i=0,(n_gal-1) DO BEGIN
			gal	= currentField[i]

			; comoving line-of-sight distance (Mpc/h)
			dz 	= redh100()*dcomovinglos(gal.zprimus, /Mpc)
			
			dy	= dz*(!pi/180.0)*( gal.DEC - mean(currentField.DEC) )
			dx	= dz*(!pi/180.0)*( gal.RA - mean(currentField.RA) )*cos((!pi/180.0)*mean(currentField.DEC))

			; proper projected x and y distances (Mpc)
			xprop	= dx/redh100()/(1 + gal.zprimus)
			yprop	= dy/redh100()/(1 + gal.zprimus)

			newRow = create_struct(gal, 'dx', dx, 'dy', dy, 'dz', dz, 'xprop', xprop, 'yprop', yprop)
;			IF (field eq 'cfhtls_xmm') THEN print, dx, dy, dz, xprop, yprop
			zerodOutput[outputIndex] = newRow
			outputIndex += 1
		
;			IF outputIndex mod 10 eq 0 THEN print, outputIndex, " ", dx, dy, dz, xprop, yprop
		ENDFOR
	ENDFOR

	print, fields
	print, "input size: ", n_elements(zerodInput)
	print, "outputIndex: ", outputIndex
	IF outputIndex ne n_elements(zerodInput) THEN message, "outputIndex mismatch"

	mwrfits, zerodOutput, path + 'zerodSFQ_all_cart.fits', /create
END
