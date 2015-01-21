
program walking

	use NumType
	implicit none

	!===================================================

	real(dp)	:: r, r_rad, x,y, length,x_0,y_0
	integer		:: i

	!===================================================

	length = 1._dp

	x = 0.0_dp
	y = 0.0_dp

	call init_random_seed()

	do i = 0,1000000
! 		print *, i
! 		call init_random_seed()
		call random_number(r)
		r_rad = r*2*pi
		x = x + length*cos(r_rad)
		y = y + length*sin(r_rad)
! 		print *, r_rad, x, y, x**2+y**2
		write(1,*) x,y
	enddo 

	print *, x,y

end program walking