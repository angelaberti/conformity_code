FUNCTION ipfrac_by_field, inputFile, fields

	zerodInput = mrdfits(inputFile, 1)

;	fields = zerodInput[uniq(zerodInput.field, sort(zerodInput.field))].field
	n_fields = n_elements(fields)

	ipfracArray = []

	FOR j=0,n_fields-1 DO BEGIN
		data = zerodInput[where(zerodInput.field eq fields[j])]
		n_ip = data[where(data.IP eq 1)]
		fieldStats = [strtrim(fields[j]), strtrim(n_elements(data)), strtrim(float(n_elements(n_ip))/float(n_elements(data)))]
	;	print, fieldStats
		ipfracArray = [ipfracArray, fieldStats]
	ENDFOR	

	RETURN, ipfracArray
END
