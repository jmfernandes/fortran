!first program for the computational physics class
!assumed integers are i,j,k,l,m,n in the fortran language

module NumType

	save
	integer, parameter :: dp = kind(1.d0)				!double precision
	real(dp), parameter :: pi = 4*atan(1._dp)			!pi=3.14
	complex(dp), parameter :: iic = (0._dp,1._dp)		!complex unit

end module NumType