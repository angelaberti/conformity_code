FUNCTION grid_array, NX, NY

; OUTER MARGINS
	XMouter=0.10
	YMouter=0.10
; INNER MARGINS
	XMinner = 0.02
	YMinner = 0.02

	IF (NX EQ 1) THEN XMinner = 0.
	IF (NY EQ 1) THEN YMinner = 0.

	DX = (1 - 1.5*XMouter)/NX
	DY = (1 - 2.0*YMouter)/NY

  posArray = []
  FOR i=0,NX*NY-1 DO BEGIN

; EQUAL SIZE GRID
	POS = [ XMouter + DX*FLOAT(i MOD NX) + XMinner/2., $
		YMouter + DY*FLOAT(FLOOR(FLOAT(i)/NX)) + YMinner/2., $
		XMouter + DX*(1. + FLOAT(i MOD NX)) - XMinner/2., $
		YMouter + DY*(1. + FLOAT(FLOOR(FLOAT(i)/NX))) - YMinner/2. ]
	
	posArray = [[posArray], [POS]]
  ENDFOR

  RETURN, posArray
END

; UNEQUAL HEIGHT UPPER AND LOWER PANELS (2 ROWS)
;	S = 0.67
;  IF (i LE NX-1) THEN BEGIN
;	POS = [ XMouter + DX*FLOAT(i MOD NX) + XMinner/2., $
;		YMouter + DY*FLOAT(FLOOR(FLOAT(i)/NX)) + YMinner/2., $
;		XMouter + DX*(1. + FLOAT(i MOD NX)) - XMinner/2., $
;		S*(YMouter + DY*(1. + FLOAT(FLOOR(FLOAT(i)/NX)))) - YMinner/2. ]
;  ENDIF ELSE BEGIN
;	POS = [ XMouter + DX*FLOAT(i MOD NX) + XMinner/2., $
;		S*(YMouter + DY*FLOAT(FLOOR(FLOAT(i)/NX))) + YMinner/2., $
;		XMouter + DX*(1. + FLOAT(i MOD NX)) - XMinner/2., $
;		YMouter + DY*(1. + FLOAT(FLOOR(FLOAT(i)/NX))) - YMinner/2. ]
;  ENDELSE
