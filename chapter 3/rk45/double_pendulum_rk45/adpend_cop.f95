
module setup

	use NumType
	implicit none
	integer, parameter :: n_eq = 4
	real(dp), parameter:: g = 10.0_dp, length = 1.0_dp, mass1 = 1.0_dp, mass2 = 2.0_dp
	real(dp) :: t, tmax, dt, lambda
	real(dp), dimension(n_eq) :: y
	
end module setup

program pendulum

	use setup
	implicit none

	!!initial conditions!!
		
	t = 0._dp					! time to start
	tmax = 100._dp				! time to exit
	dt = 0.001_dp				! time step
	lambda = 0._dp				!intiate a value for lyapanov exponent


	!!Below are four cases - Make sure only one set of y is uncommented!!

	!!Low Energy - No Chaos
! 	y(1) = 1._dp*(pi/180)		! initial angle 1
! 	y(2) = 0._dp				! initial momentum 1
! 	y(3) = 1._dp*(pi/180)		! initial angle 2
! 	y(4) = 0._dp				! initial momentum 2

	!!Medium Energy - Emerging Chaos But Stable Orbit
!  	y(1) = 1._dp*(pi/180)		! initial angle 1
!  	y(2) = 1._dp				! initial momentum 1
!  	y(3) = 1._dp*(pi/180)		! initial angle 2
!  	y(4) = 1.01_dp				! initial momentum 2

	!!High Energy - Full Chaos
! 	y(1) = 1._dp*(pi/180)		! initial angle 1
! 	y(2) = 1._dp				! initial momentum 1
! 	y(3) = 1._dp*(pi/180)		! initial angle 2
! 	y(4) = 1.3_dp				! initial momentum 2

	!!High Energy - Full Chaos
	y(1) = 60._dp*(pi/180)		! initial angle 1
	y(2) = 1._dp				! initial momentum 1
	y(3) = 60._dp*(pi/180)		! initial angle 2
	y(4) = 2._dp				! initial momentum 2

	
	!!open all the files that data will be written to!!

	open(unit = 3, file = 'pend_angle1.data', &
	    action = 'write', status = 'replace')
	open(unit = 4, file = 'pend_mom1.data', &
	    action = 'write', status = 'replace')
	open(unit = 5, file = 'pend_angle2.data', &
	    action = 'write', status = 'replace')
	open(unit = 6, file = 'pend_mom2.data', &
	    action = 'write', status = 'replace')
	open(unit = 7, file = 'angle_mom1.data', &
	    action = 'write', status = 'replace')
	open(unit = 8, file = 'angle_mom2.data', &
	    action = 'write', status = 'replace')
	open(unit = 9, file = 'pend_xy1.data', &
	    action = 'write', status = 'replace')
	open(unit = 10, file = 'pend_xy2.data', &
	    action = 'write', status = 'replace')
	open(unit = 11, file = 'lambda.data', &
	    action = 'write', status = 'replace')

	!!calculate the momenta and angle of the pendulums!!

	do while ( t < tmax )
		if ( t + dt > tmax) dt = tmax -t
		call rkf45step(t,y,dt)						!use the runga kutta method to calculate for theta and momentum
		lambda = lambda + log(abs(y(4)-y(2)))/t 	!update the summation of the lyapanov exponent
		write (unit = 11,fmt='(3f20.10)') lambda 	!write lamba to file
	end do

end program pendulum
