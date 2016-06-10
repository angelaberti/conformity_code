; given a field, a z value, and SFQstatus (0 or 1), return mass limit

FUNCTION mlimit, massStruct, fields, inputGalaxy

	zbins = [0.20, 0.30, 0.40, 0.50, 0.65, 0.80, 1.00]

	mtableSF = {cdfs:	[9.60,	9.92,	10.19,	10.44,	10.63,	10.69], $
		    cfhtls_xmm:	[8.80,	9.06,	 9.30,	 9.58,	 9.89,	10.21], $
		    cosmos:	[8.68,	9.05,	 9.38,	 9.75,	10.12,	10.46], $
		    es1:	[9.58,	9.94,	10.25,	10.59,	10.90,	11.14], $
		    xmm:	[8.79,	9.13,	 9.44,	 9.77,	10.10,	10.38]}

	mtableQ =  {cdfs:	[9.65,	 9.92,	10.17,	10.44,	10.71,	10.96], $
		    cfhtls_xmm:	[9.17,	 9.52,	 9.85,	10.22,	10.60,	10.96], $
		    cosmos:	[9.23,	 9.58,	 9.89,	10.22,	10.52,	10.75], $
		    es1:	[9.80,	10.06,	10.30,	10.55,	10.79,	10.99], $
		    xmm:	[9.35,	 9.61,	 9.85,	10.13,	10.43,	10.73]}

	field	= inputGalaxy.field
	z 	= inputGalaxy.zprimus
	SFQ	= inputGalaxy.SFQ

	field_index = array_indices(fields, where(fields eq field))

;	print, "field_index: ", field_index

;	print, field, z, SFQ

	IF (SFQ eq 1) THEN (mtable = mtableSF) ELSE (mtable = mtableQ)

	templateRow =	create_struct(	"cdfs",		0.0, $
					"cfhtls_xmm",	0.0, $
					"cosmos",	0.0, $
					"es1", 		0.0, $
					"xmm", 		0.0)
	massStruct =	replicate(templateRow, n_elements(zbins)-1)

	FOR i=0,n_elements(zbins)-2 DO BEGIN
		newRow = create_struct( "cdfs", 	mtable.cdfs[i], $
					"cfhtls_xmm", 	mtable.cfhtls_xmm[i], $
					"cosmos", 	mtable.cosmos[i], $
					"es1", 		mtable.es1[i], $
					"xmm", 		mtable.xmm[i])
	        massStruct[i] = newRow
	ENDFOR

	RETURN, massStruct[max(array_indices(zbins, where(z gt zbins)))].(field_index)
END
