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
	real(dp) :: t, dt, x, y, vx, vy

end program ball_path