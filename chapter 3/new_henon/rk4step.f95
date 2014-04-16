!header

subroutine rk4step(t,y,dt,k)

	use setup
	implicit none
	real(dp), intent(in) :: dt, k
	real(dp), intent(inout) :: t
	real(dp), dimension(n_eq), intent(out) :: y
	real(dp), dimension(n_eq) :: k1, k2, k3, k4, dy

	call deriv(t,y,dt,k1)
	call deriv(t+dt/2,y+k1/2,dt,k2)
	call deriv(t+dt/2,y+k2/2,dt,k3)
	call deriv(t+dt,y+k3,dt,k4)

	dy = (k1 + 2*k2 + 2*k3 + k4)/6

	t = t + dt
	y = y + dy

	!see what's going on print *, k1,k2,k3,k4,y(2)

	!never change above code

	contains

		subroutine deriv(t,y,dt,k)

			use setup

			real(dp), intent(in) :: t, dt
			real(dp), dimension(n_eq), intent(in) :: y
			real(dp), dimension(n_eq) :: k, f
			

			f(1) = y(2)/mass
			f(2) = -y(1)-2*y(1)*y(3)
			f(3) = y(4)/mass
			f(4) = -y(3)-y(1)**2+y(3)**2

			
			k = dt * f !step size times derivative

		end subroutine deriv

end subroutine rk4step