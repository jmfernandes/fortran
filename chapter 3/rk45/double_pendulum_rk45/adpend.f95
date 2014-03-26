
module setup

	use NumType
	implicit none
	integer, parameter :: n_eq = 2
	real(dp), parameter:: g = 10.0_dp, length = 10.0_dp, mass = 1.0_dp
	
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
! 	y(2) = 0._dp          ! initial angular speed
	y(2) = 0._dp		  ! initial momentum
! 	y(4) = 0._dp		  ! initial change in momentum	
	
	open(unit = 3, file = 'aaaapend_angle.data', &
	    action = 'write', status = 'replace')

	open(unit = 4, file = 'aaaaapend_momentum.data', &
	    action = 'write', status = 'replace')

	open(unit = 5, file = 'aapend_ang_mo.data', &
	    action = 'write', status = 'replace')

! 	open(unit = 5, file = 'pend_momentum.data', &
! 	    action = 'write', status = 'replace')

	do while ( t < tmax )	
		if ( t + dt > tmax) dt = tmax -t
		call rkf45step(t,y,dt)	
	end do
	
end program pendulum
