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
			real(dp) :: c1,c2


			c1 = (y(2)*y(4)*sin(y(1)-y(3)))/ &
	        	(length*length*(mass1+mass2*(sin(y(1)-y(3)))**2))
	        c2 = (length**2*mass2*y(2)**2 + &
				length**2*(mass1+mass2)*y(4)**2- &
				length*length*mass2*y(2)*y(4)*cos(y(1)-y(3)))/ &
				(2*length**4*(mass1+mass2*(sin(y(1)-y(3)))**2)**2)* &
				sin(2*(y(1)-y(3)))

			f(1) = (length*y(2) - length*y(4)*cos(y(1)-y(2)))/ &
					(length**3*(mass1+mass2*(sin(y(1)-y(2)))**2))
			f(2) = -(mass1+mass2)*g*length*sin(y(1)) - c1 + c2
			f(3) = (length*(mass1+mass2)*y(4) - length*mass2*y(2)*cos(y(1)-y(2)))/ &
					(length**3*mass2*(mass1+mass2*(sin(y(1)-y(2)))**2))
			f(4) = -mass2*g*length*sin(y(3)) + c1 - c2

			k = dt * f !step size times derivative

		end subroutine deriv

end subroutine rk4step