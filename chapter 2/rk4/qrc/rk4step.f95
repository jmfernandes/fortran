!header

subroutine rk4step(t,y,dt,k)

	use setup
	implicit none
	real(dp), intent(in) :: dt, k
	real(dp), intent(inout) :: t
	real(dp), intent(out) :: y
	real(dp) :: k1, k2, k3, k4, dy

	call deriv(t,y,dt,k1)
	call deriv(t+dt/2,y+k1/2,dt,k2)
	call deriv(t+dt/2,y+k2/2,dt,k3)
	call deriv(t+dt,y+k3,dt,k4)

	dy = (k1 + 2*k2 + 2*k3 + k4)/6

	t = t + dt
	y = y + dy

	contains

		subroutine deriv(t,q,dt,k)

			real(dp), intent(in) :: t, dt, q
			real(dp) :: k

			k = dt * (emf - q/capacity)/resistance

		end subroutine deriv

end subroutine rk4step