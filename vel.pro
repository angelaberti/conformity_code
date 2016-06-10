FUNCTION vel, z, sigfactor

	H = 100.*REDH100()

	sigz = 0.005*(1.+z)
	
	delta_v = H*(DCOMOVINGLOS(z+sigfactor*sigz,/Mpc)-DCOMOVINGLOS(z,/Mpc))

	RETURN, delta_v
END
