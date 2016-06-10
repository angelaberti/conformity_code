PRO ipfrac, inputFile
	zerodInput = mrdfits(inputFile, 1)

	fields = zerodInput[uniq(zerodInput.field, sort(zerodInput.field))].field
	n_fields = n_elements(fields)

	ipfracArray = []

	FOR j=0,n_fields-1 DO BEGIN
		data = zerodInput[where(zerodInput.field eq fields[j])]
		dataIP  = data[where(data.IP eq 1)]
		fieldStats = [strtrim(fields[j]), string(n_elements(data)), $
			      string(float(n_elements(dataIP))/float(n_elements(data)))]
		ipfracArray = [ipfracArray, fieldStats]
	ENDFOR

	print, inputFile
	print, ipfracArray
END
