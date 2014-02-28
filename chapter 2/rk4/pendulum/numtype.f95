
module NumType

	save
	integer, parameter :: dp = kind(1.d0)
	real(dp), parameter :: pi = 4*atan(1._dp), g = 9.80665_dp 
	complex(dp), parameter :: iic = (0._dp,1._dp)
	
end module NumType
