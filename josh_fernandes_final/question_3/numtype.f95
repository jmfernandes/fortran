module numtype
	save
	integer, parameter :: dp = kind(1.d0)			!double precision
	real(dp), parameter :: pi = 4*atan(1._dp), e = exp(1._dp)		
	complex(dp), parameter :: iic = (0._dp,1._dp)	!complex part of function

end module numtype 

