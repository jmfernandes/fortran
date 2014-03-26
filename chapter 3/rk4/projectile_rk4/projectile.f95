!header

module setup

	use NumType
	implicit none
	real(dp), parameter :: g = 9.81_dp, mass = 1._dp
	real(dp), parameter :: X0 = 0._dp, Y0 = 0._dp
	real(dp) :: drag, alpha, V0, wind(2)
	integer, parameter :: number_of_equations = 4

end module setup

program kicking

	use setup
	implicit none
	real(dp) :: y(number_of_equations), t, dt, x


	t = 0._dp
	dt = 0.1_dp
	drag = 0.05_dp
	wind = (10._dp, 0._dp)
	alpha = 75._dp*(pi/180)
	V0 = 20._dp

	y(1) = X0
	y(2) = mass*V0*cos(alpha)
	y(3) = Y0
	y(4) = mass*V0*sin(alpha)
	

	open(unit = 3, file='x.data', action = 'write', status = 'replace')
	open(unit = 4, file='y.data', action = 'write', status = 'replace')
	open(unit = 5, file='x-y.data', action = 'write', status = 'replace')

	do while( y(3) >= 0._dp )
		
		write(3, fmt = '(2f20.10)') t, y(1)
		write(4, fmt = '(2f20.10)') t, y(3)
		write(5, fmt = '(2f20.10)') y(1), y(3)
		call rk4step(t,y,dt)

	end do


end program kicking