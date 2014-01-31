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

	contains

		function tcoeff(e,gam) result(t)

			!use NumType
			!use setup
			implicit none

			real(dp) :: e, gam, t

			t = (gam/2)**2 / ( (e-e0)**2 + (gam/2)**2 )

		end function tcoeff


end program transmit


