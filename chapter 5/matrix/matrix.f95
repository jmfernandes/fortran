
module setup

	use NumType
	implicit none
	real(dp), parameter :: n=2,m=3,l=4
	integer(dp), parameter :: ndim=5

end module setup

program matrix

	use setup
	implicit none

	complex(dp), dimension(ndim,ndim) :: A,B,C,D,E,F

	A(1:3,1:3) = reshape((/ 7*one,zero,zero,zero,one,-iic,zero,iic,-one /), &
		(/3,3/))

	B(1:3,1:3) = reshape((/ one,zero,3*one,zero,2*iic,zero,iic,zero,-5*iic /), &
		(/3,3/))

	print *, A(1:3,1:3)
	print *, '------------'
	print *, B(1:3,1:3)


end program matrix