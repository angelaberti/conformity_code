FUNCTION bigint, num

	strnum	= strtrim(string(num, format='(i20)'),1)
	numSegs = floor(strlen(strnum)/3.)
	leadSegLength = strlen(strnum)-3.*numSegs
	leadSeg	= strmid(strnum, 0, leadSegLength)
	
	result = leadSeg
	iStart = 0

	IF strlen(leadSeg) eq 0 THEN BEGIN
		result = strmid(strnum, 0, 3)
		iStart = 1
	ENDIF

	FOR i=iStart,numSegs-1 DO result += ',' + strmid(strnum, leadSegLength+3*i, 3)
	result = result + strmid(strnum, leadSegLength+3*numSegs, 3)

	IF strlen(strnum) le 3 THEN RETURN, strnum ELSE RETURN, result
END
