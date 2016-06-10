PRO medianStats

	RESTORE, 'IPQ_latefrac_array_dR1Mpc.sav'
	IPQ = IP_latefrac_array
	RESTORE, 'IPSF_latefrac_array_dR1Mpc.sav'
	IPSF = IP_latefrac_array

;	intervals = [[0,1],[1,2],[1,3],[3,5]]
	intervals = [[0,5]]

;	PRINT, intervals

	FOR i=0,n_elements(intervals)-1 DO BEGIN
		pair = intervals[*,i]
		allRangeQ = [TRANSPOSE(IPQ[pair[0]:pair[1]-1,*]), TRANSPOSE(IPQ[pair[0]:pair[1]-1,*])]
		allRangeQ = allRangeQ[where(allRangeQ GT 0.)]

		allRangeSF = [TRANSPOSE(IPSF[pair[0]:pair[1]-1,*]), TRANSPOSE(IPSF[pair[0]:pair[1]-1,*])]
		allRangeSF = allRangeSF[where(allRangeSF GT 0.)]

;		PRINT, pair[0], pair[1]-1, n_elements(allRangeQ), MEDIAN(allRangeQ), n_elements(allRangeSF), MEDIAN(allRangeSF)
		PRINT, pair[0], pair[1], MEDIAN(allRangeSF)-MEDIAN(allRangeQ)
;		PRINT, pair[0], pair[1], MEAN(allRangeSF)-MEAN(allRangeQ)
	ENDFOR
END
