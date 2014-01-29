!this program shoots a ball

module setup

	use NumType
	implicit none
	real(dp), parameter :: g = 9.81_dp
	real(dp) :: h0, v0, alpha

end module setup

program ball_path 

	use setup
	implicit none
	real(dp) :: t, dt, x, y, vx, vy, vx0, vy0, x0, y0

	!!define variables
	h0 = 100._dp
	v0 = 15._dp
	alpha = 36*(pi/180) !convert from degrees to radians
	
	t = 0._dp
	dt = 0.01_dp

	x0 = 0._dp
	y0 = h0

	vx0 = v0*cos(alpha)
	vy0 = v0*sin(alpha)
	

	!!calculate

	do while (y >= 0._dp)

		vx = vx0
		vy = vy0 - g*t
		x = vx*t
		y = y0 + vy0*t - g/2 * t**2

		print *, t, vx, vy
		t = t + dt

	end do






end program ball_path