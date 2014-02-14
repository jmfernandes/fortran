!this program is Dr. Papp's way of solving the 1D harmonic oscillator
program stuff
	use NumType
	implicit none
	real(dp) :: x, Xmin, Xmax, dx
	integer :: steps, i

	x = 2.5_dp
	Xmin = -5._dp
	Xmax = 5._dp

	steps = 500

	dx = (Xmax-Xmin)/steps

	do i=0,steps
		x = Xmin + i*dx
		print *,x,psi_ho1d(0,x),psi_ho1d(1,x),psi_ho1d(2,x)

	end do


contains

	recursive function psi_ho1d(n,x) result(psi)

		real(dp) :: x, y, psi, alpha
		real(dp), parameter :: hbar = 1._dp, mass = 1._dp, omega = 1._dp
		integer :: n

		alpha = sqrt(mass*omega/hbar)
		y = alpha*x

		if ( n < 0 ) then
			psi = 0._dp
		else if ( n == 0) then
			psi = sqrt(alpha/sqrt(pi))*exp(-y**2/2)
		else
			psi = sqrt(2._dp/n)*y*psi_ho1d(n-1,y) &
				- sqrt((n-1._dp)/n)*psi_ho1d(n-2,y)
		end if 


	end function psi_ho1d

end program stuff