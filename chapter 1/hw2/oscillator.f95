!this programs solves the 1 dimensional quantum harmonic oscillator

module initiate_phase_one

	use NumType
	implicit none
	real(dp), parameter :: w = 1._dp, &
	m = 1._dp, Xmin = -10._dp, Xmax = 10._dp
	real(dp) :: x, k, plot1(5000),plot2(5000),plot3(5000), &
	plot4(5000),plot5(5000),plot6(5000),plot7(5000), plot8(5000)
	REAL, DIMENSION(:, :), ALLOCATABLE :: psi
	integer, Dimension(:), ALLOCATABLE :: n
	real(dp), Dimension(:), ALLOCATABLE :: curve
	integer :: i, j, DeAllocateStatus, AllocateStatus, steps
	

end module initiate_phase_one

program potatoe

	use initiate_phase_one
	
	n = (/0,1, 2, 3, 4, 5, 6, 7/)

	j = 0._dp
	k = 0.5_dp

	steps = 5000
	x = Xmin
	dx = (Xmax-Xmin)/steps



	ALLOCATE (curve(steps), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

	ALLOCATE ( psi(size(n), steps), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

	do while(x < Xmax)
		j = j + 1
		curve(j) = 0.5_dp*k*x**2
		plot1(j) = 1
		plot2(j) = 3
		plot3(j) = 5
		plot4(j) = 7
		plot5(j) = 9
		plot6(j) = 11
		plot7(j) = 13
		plot8(j) = 15
		do i=1,size(n)
			psi(i,j) = quantum_oscillator(n(i),m,w,x)
		end do
		print *, x, 1+psi(1,j), 3+psi(2,j), &
		5+psi(3,j), 7+psi(4,j), 9+psi(5,j), &
		11+psi(6,j), 13+psi(7,j), 15+psi(8,j), &
		curve(j), plot1(j), plot2(j), plot3(j), &
		plot4(j), plot5(j), plot6(j), plot7(j), plot8(j)
		!iterate x
		x = x + dx
	end do

	DEALLOCATE (curve, STAT = DeAllocateStatus)
	DEALLOCATE (psi, STAT = DeAllocateStatus)


	

	contains

		function quantum_oscillator(n,m,w,x) result(psi)

			implicit none
			real(dp), parameter :: hbar = 3._dp
			real(dp) :: m,w,x,psi
			integer :: n

			psi = 1/sqrt(real((2**n)*factorial(n)))*((m*w)/(pi*hbar))**(1/4)&
			*exxp(-(m*w*x**2)/(2*hbar))*hermite(n,(sqrt(m*w/hbar)*x))

		end function quantum_oscillator

		recursive function exxp(x) result(ex)

		real(dp) :: x, ex, E0
		integer :: i, imax

		imax = 20
		if ( abs(x) < 1._dp) then
			E0 = 1._dp
			ex = 1._dp
			do i =1, imax
				E0 = E0*x/i
				ex = ex + E0
			end do
		else if (1._dp <= x) then
			ex = e * exxp(x-1)
		else if (x <= 1._dp) then
			ex = exxp(x+1) / e
		end if

		end function exxp 

		recursive function hermite(n,x) result(hpol)

			real(dp) :: x, hpol
			integer :: n

			if ( n < 0 ) then
				stop "n can't be less than zero"
			else if ( n == 0) then
				hpol = 1._dp
			else if ( n == 1) then
				hpol = 2*x
			else
				hpol = 2*x*hermite(n-1,x) - 2*(n-1)*hermite(n-2,x)
			end if 


		end function hermite


		recursive function factorial(n) result(factorial_number)

			implicit none
			integer :: n, factorial_number


			if ( n < 0 ) then
				stop "you know that's wrong"
			else if ( n == 0) then
				factorial_number = 1
			else
				factorial_number = n*factorial(n-1)
			end if 



		end function factorial


end program potatoe
