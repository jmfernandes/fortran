!this programs solves the hermite equation

module initiate_phase_one

	use NumType
	implicit none
	!real(dp), parameter :: X0 = 0._dp
	real(dp) :: x
	!REAL, DIMENSION(:, :), ALLOCATABLE :: result
 	!integer :: i, j, DeAllocateStatus, AllocateStatus, steps
	

end module initiate_phase_one

program power_up_the_bass_cannon

	use initiate_phase_one

	x = 2.5_dp

	print *, 'result', x, hermite(4,x), hermite(5,-x)

contains

	recursive function hermite(n,x) result(hpol)

		real(dp) :: x, hpol
		integer :: n

		if ( n < 0 ) then
			stop 'error'
		else if ( n == 0) then
			hpol = 1._dp
		else if ( n == 1) then
			hpol = 2*x
		else
			hpol = 2*x*hermite(n-1,x) - 2*(n-1)*hermite(n-2,x)
		end if 


	end function hermite

	
end program power_up_the_bass_cannon

! program potatoe

! 	use initiate_phase_one

! 	x = 20.5_dp

! 	print *, 'result', x, exxp(x), exp(x)

! contains
! 	recursive function exxp(x) result(ex)

! 		real(dp) :: x, ex, E0
! 		integer :: i, imax

! 		imax = 20
! 		if ( abs(x) < 1._dp) then
! 			E0 = 1._dp
! 			ex = 1._dp
! 			do i =1, imax
! 				E0 = E0*x/i
! 				ex = ex + E0
! 			end do
! 		else if (1._dp <= x) then
! 			ex = e * exxp(x-1)
! 		else if (x <= 1._dp) then
! 			ex = exxp(x+1) / e
! 		end if

! 		end function exxp 

! end program potatoe

! program potatoe

! 	use initiate_phase_one
! 	steps = 10
! 	x = 2._dp
! 	E0 = 1._dp
! 	E1 = 1._dp

! 	print *, 0, E0, E1
! 	do i = 1,steps
! 		E0 = E0*(x/i)
! 		E1 = E1 + E0
! 		print *, i, E0, E1
! 	end do

! 	print *, 'result',x,E1,exp(x)

! end program potatoe

! program heat

! 	use setup
! 	implicit none
	
! 	t  = (/ 5.0, 1.0, 0.5, 0.1, 0.01/) !set the values of Q

! 	steps = 500
! 	x = Xmin
! 	dx = (Xmax-Xmin)/steps

! 	j = 0._dp		

! 	!google told me I had to do this to allocate space
! 	ALLOCATE ( result(size(t), steps), STAT = AllocateStatus)
! 	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
	
! 	do while(x < Xmax)
! 		j = j + 1
! 		do i=1,5
! 			result(i,j) = tcoeff(x,i)
! 		end do
! 		print *, x, result(1,j), result(2,j), result(3,j), result(4,j), result(5,j)
! 		!iterate x
! 		x = x + dx
! 	end do

! 	DEALLOCATE (result, STAT = DeAllocateStatus)

! 	contains
! 		function tcoeff(x,i) result(theta)

! 			implicit none

! 			real(dp) :: x, theta
! 			integer :: i

! 			theta = Q/(2*sqrt(pi*conductivity*t(i)))*exp(-(x-X0)**2/(4*conductivity*t(i)))

! 		end function tcoeff
		
! end program heat