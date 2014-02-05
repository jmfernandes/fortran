!this programs solves the 1 dimensional quantum harmonic oscillator

module initiate_phase_one

	use NumType
	implicit none
	real(dp), parameter :: X0 = 0._dp
	real(dp) :: x, E0, E1
	REAL, DIMENSION(:, :), ALLOCATABLE :: result
	integer :: i, j, DeAllocateStatus, AllocateStatus, steps
	

end module initiate_phase_one

program potatoe

	steps = 10
	x = -2._dp
	E0 = 1._dp
	E1 = 1._dp

	do i = 1,steps
		E0 = E0*(x/i)
		E1 = E1 + E0
	end do

end program potatoe

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