!this program calculates Scattering

module setup

	use NumType
	implicit none
	real(dp), parameter :: e0 = 5._dp, emax = 10._dp
	integer, parameter :: imax = 500

end module setup

program transmit

	use setup
	implicit none
	real(dp) :: e, de, gam, t1, t2, t3, t4
	integer :: i 
	real(dp), external :: tcoeff

	de = emax/imax

	do i = 0, imax
		e = de*i

		gam = 1._dp
		t1 = tcoeff(e,gam)

		gam = 0.5_dp
		t2 = tcoeff(e,gam)

		gam = 0.1_dp
		t3 = tcoeff(e,gam)

		gam = 0.01_dp
		t4 = tcoeff(e,gam)

		print *, e, t1, t2, t3, t4

	end do


end program transmit


function tcoeff(e,gam) result(t)

	use NumType
	use setup
	implicit none

	real(dp) :: e, gam, t

	t = (gam/2)**2 / ( (e-e0)**2 + (gam/2)**2 )

end function tcoeff




















! 	real(dp) :: E, dE, E0, gamma(4)
! 	REAL, DIMENSION(:, :), ALLOCATABLE :: T
! 	integer :: i, j, N, M, DeAllocateStatus, AllocateStatus, Emax
	

! end module setup

! program ball_path 

! 	use setup
! 	implicit none
	
! 	gamma  = (/ 1.0, 0.5, 0.1, 0.01/) !set the values of gamma

! 	E = 0._dp		!this is min E value	
! 	dE = 0.01_dp	!this is the step size
! 	Emax = 10_dp	!this is max E value

! 	E0 = 5_dp		!this is a constant
! 	j = 0._dp		

! 	!calculate the size of the array based on how many gammas, and how many steps I am going to take
! 	N = size(gamma)
! 	M = (Emax/dE)+1

! 	!google told me I had to do this to allocate space
! 	ALLOCATE ( T(N, M), STAT = AllocateStatus)
! 	IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		
! 	do while (E < Emax)
! 		j = j + 1
! 		!this creates a 4 column array of the T data. Each column is for a different value of gamma
! 		do i=1,4
! 			T(i,j) = (.25*(gamma(i)**2))/((E-E0)**2 + .25*(gamma(i)**2))
! 		end do
! 		!print the data
! 		print *, E, T(1,j), T(2,j), T(3,j), T(4,j)
! 		!iterate E
! 		E = E + DE
! 	end do 
		
! 	DEALLOCATE (T, STAT = DeAllocateStatus)
! end program ball_path