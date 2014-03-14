


module setup

	use NumType
	implicit none
	real(dp), parameter :: n=2,m=3,l=4, dim=2

end module setup

program matrix

	use setup
	implicit none

	real(dp), dimension(dim) :: A(1:n,1:m), B(1:m,1:l), C(1:n,1:l) 

end program matrix