; returns float as string with specified decimal places

FUNCTION decimal, input, places
	
	format = '(f20.' + strtrim(string(places, format='(i20)'),1) + ')'

	RETURN, strtrim(string(input, format=format),1)	
END
