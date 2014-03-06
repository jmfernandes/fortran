
module setup

	use NumType
	implicit none
	integer, parameter :: n_eq = 4
	real(dp), parameter:: g = 10.0_dp, length = 1.0_dp, mass1 = 1.0_dp, mass2 = 1.0_dp
	
end module setup

program pendulum

	use setup
	implicit none
	real(dp), dimension(n_eq) :: y
	real(dp) :: y0, v0, t, tmax, dt, c1, c2

		
	t = 0._dp		      ! time to start
	tmax = 60._dp         ! time to exit
	dt = 0.1_dp		      ! time step	
	y(1) = 30*(pi/180)    ! initial angle 1
	y(2) = 0._dp          ! initial momentum 1
	y(3) = 30*(pi/180) 	  ! initial angle 2
	y(4) = 0._dp		  ! initial momentum 2


! 	c1 = (y(3)*y(7)*sin(y(1)-y(5)))/ &
! 		(length*length(mass1+mass2*(sin(y(1)-y(5)))**2))
! 	c2 = (length**2*mass2*y(3)**2 + &
! 		length**2*(mass1+mass2)*y(7)**2- &
! 		length*length*mass2*y(3)*y(7)*cos(y(1)-y(5)))/ &
! 		(2*length**4(mass1+mass2*(sin(y(1)-y(5)))**2)**2)* &
! 		sin(2*(y(1)-y(5)))
	
	open(unit = 3, file = 'pend_angle1.data', &
	    action = 'write', status = 'replace')

	open(unit = 4, file = 'pend_mom1.data', &
	    action = 'write', status = 'replace')

	open(unit = 5, file = 'pend_angle2.data', &
	    action = 'write', status = 'replace')

	open(unit = 6, file = 'pend_mom2.data', &
	    action = 'write', status = 'replace')

	do while ( t < tmax )	
		if ( t + dt > tmax) dt = tmax -t
		call rkf45step(t,y,dt)	
	end do
	
end program pendulum
