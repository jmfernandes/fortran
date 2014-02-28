!this programs solves the 1 dimensional quantum harmonic oscillator

module setup

	use NumType
	implicit none
	real(dp), parameter ::  m = 10._dp, l = 5._dp, tmin = 0._dp, tmax = 10._dp
	real(dp) :: t, dt, results, steps, theta_max
	

end module setup

program potatoe

	use setup
	
	steps = 500._dp
	theta_max = .5_dp
	t = tmin
	dt = (tmax-tmin)/steps


	do while (t < tmax)

		t = t + dt
		results = pendulum_osc(theta_max,m,t)
		print *, t, results

	end do
	

	contains

		function pendulum_osc(a,m,t) result(theta)

			implicit none
			real(dp) :: a,m,t,theta

			theta = a * sin(sqrt(g/l)*t)

		end function pendulum_osc


end program potatoe
