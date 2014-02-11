!this programs solves the 1 dimensional quantum harmonic oscillator

module initiate_phase_one

	use NumType
	implicit none
! 	real(dp), parameter :: X0 = 0._dp
! 	real(dp) :: x, E0, E1
! 	REAL, DIMENSION(:, :), ALLOCATABLE :: result
	complex(dp) :: z, expcc
	

end module initiate_phase_one

program spud

	use initiate_phase_one
	
	z = (2.5_dp, 1.5_dp)
	expcc = 1/(1-ss(1,z))

	print *, 'result', z, exp(z)
	print *, 'result', z, expcc

	contains

		recursive function ss(i,z) result(scf)

			implicit none
			integer :: i, imax
			complex(dp) :: z, scf


			imax = 5000

			if (i >= imax) then
				scf = 0._dp
			else
				scf = (z/i) / (1 + (z/i) - ss(i+1,z))
			end if



		end function ss


end program spud
