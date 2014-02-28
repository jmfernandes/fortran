!header

module setup

	use NumType
	implicit none
	real(dp), parameter :: g = 9.81_dp, length = 10._dp, mass = 1._dp
	real(dp), parameter :: mu = 0.1_dp
	integer, parameter :: number_of_equations = 2

end module setup

program swinging

	use setup
	implicit none
	real(dp) :: y(number_of_equations), t, dt, tmax, theta, omega


	t = 0._dp
	tmax = 100._dp
	dt = 0.1_dp
	theta = 42_dp*(pi/180)
	omega = 0._dp

	y(1) = theta
	y(2) = omega
	

	open(unit = 3, file='pendulum_angle.data', action = 'write', status = 'replace')
	open(unit = 4, file='pendulum_omega.data', action = 'write', status = 'replace')


	do while( t <= tmax )
		
		write(3, fmt = '(2f20.10)') t, y(1)
		write(4, fmt = '(2f20.10)') t, y(2)
		print *, t, y
		call rk4step(t,y,dt)

	end do


end program swinging