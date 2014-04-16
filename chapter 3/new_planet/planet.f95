
module setup

	use NumType
	implicit none
	real(dp), parameter :: gravity = 6.673e-11_dp, mass_sun = 1.9891e+30_dp
	real(dp), parameter :: mass_earth = 5.9736e+24_dp,mass_moon = 7.3476e+22_dp
	integer, parameter :: n_eq = 6

end module setup

program orbit

	use setup
	implicit none

	real(dp) :: y(n_eq), t, dt, tmax
	real(dp) :: X0, Y0, Z0, VX0, VY0, VZ0

	!!initial conditions!!
		
	t = 0._dp
	tmax = 7*(60*60*24*365) !years
	dt = 7*(24*60*60) !days
	

	X0 = 1.596e+11_dp
	Y0 = 0._dp
	Z0 = 0._dp

	VX0 = 0._dp
	VY0 = 29.783e+3_dp
	VZ0 = 30.783e+2_dp

	y(1) = X0
	y(2) = Y0
	y(3) = Z0
	y(4) = VX0*mass_earth
	y(5) = VY0*mass_earth
	y(6) = VZ0*mass_earth


	


	open(unit = 3, file='earth.data', action = 'write', status = 'replace')


	do while( t <= tmax )
		
		write(3, fmt = '(2f20.0)') y(1), y(2)
		
		call rk4step(t,y,dt)

	end do







end program orbit
