
module setup

	use NumType
	implicit none
	real(dp), parameter :: mass = 1._dp
	integer, parameter :: n_eq = 3
	real(dp), parameter :: a(1:n_eq) = (/ 0.2_dp, 0.2_dp, 5.7_dp /) 

end module setup

program orbit

	use setup
	implicit none

	real(dp) :: y(n_eq), t, dt, tmax
	real(dp) :: E

	!!initial conditions!!
		
	t = 0._dp
	tmax = 100._dp 
	dt = 0.01_dp 
	


	y(1:n_eq) = (/ -1._dp, 0._dp, 0._dp /) 



! 	open(unit = 13, file='xt.data', action = 'write', status = 'replace')


	do while( t <= tmax )
	
		write(13,*) t,y(1),y(2),y(3)
		call rk4step(t,y,dt)

	end do







end program orbit
