FUNCTION poissonError, late, total

	RETURN, (late/total)*( (1./late) + (1./total) )^(0.5)
;	RETURN, ( ((late^2)/(total^3)) + (late/(total^2)) )^(0.5)

END
