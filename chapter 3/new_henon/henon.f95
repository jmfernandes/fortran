
module setup

	use NumType
	implicit none
	real(dp), parameter :: mass = 1._dp
	integer, parameter :: n_eq = 4

end module setup

program orbit

	use setup
	implicit none

	real(dp) :: y(n_eq), t, dt, tmax
	real(dp) :: E

	!!initial conditions!!
		
	t = 0._dp
	tmax = 100._dp !years
	dt = 0.1_dp !days
	

	E = 1/200._dp

	y(1) = 0 !x
	y(2) = 0 !px
	y(3) = 0 !y
	y(4) = sqrt(2*mass*E)!py


	


	open(unit = 11, file='point1.data', action = 'write', status = 'replace')
	open(unit = 12, file='point2.data', action = 'write', status = 'replace')
	open(unit = 13, file='pointxy.data', action = 'write', status = 'replace')


	do while( t <= tmax )
		
		if (0.5_dp*(y(2)**2+y(4)**2) >= 0._dp) then
			write(11, *) y(1), y(2)
			write(12, *) y(3), y(4)
			write(13, *) y(1), y(3)
		end if 
		
		call rk4step(t,y,dt)

	end do







end program orbit
