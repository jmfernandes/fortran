
module setup

	use NumType
	implicit none
	integer, parameter :: n_eq = 2
	real(dp), parameter:: g = 10.0_dp, lenght = 10.0_dp
	
end module setup

program pendulum

	use setup
	implicit none
	real(dp), dimension(n_eq) :: y
	real(dp) :: y0, v0, t, tmax, dt

		
	t = 0._dp		      ! time to start
	tmax = 60._dp         ! time to exit
	dt = 0.1_dp		      ! time step	
	y(1) = 30*(pi/180)    ! initial angle
	y(2) = 0._dp          ! initial angular speed	
	
	open(unit = 3, file = '1_pend.data', &
	    action = 'write', status = 'replace')

	do while ( t < tmax )	
		if ( t + dt > tmax) dt = tmax -t
		call rkf45step(t,y,dt)	
	end do
	
end program pendulum
