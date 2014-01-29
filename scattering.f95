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
	real(dp) :: E, dE, E0, gamma(4), T(4,1001)
	integer :: i, j


	gamma  = (/ 1.0, 0.5, 0.1, 0.01/)
	!print *, a

	E = 0._dp
	dE = 0.1_dp

	E0 = 5_dp
	!!a = (/ 2, 3, 5, 7, 11 /)
	j = 0._dp

	
		
		do while (E < 10)
		j = j + 1
		!print *, a(i)
		do i=1,4
			T(i,j) = (.25*(gamma(i)**2))/((E-E0)**2 + .25*(gamma(i)**2))
		end do
		print *, E, T(1,j), T(2,j), T(3,j), T(4,j)
		E = E + DE
		end do 
		E = 0
		j = 0
		!print *, a(i)
	


	






end program ball_path