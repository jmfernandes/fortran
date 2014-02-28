!header

module setup

	use NumType
	implicit none
	real(dp), parameter :: emf = 40._dp, capacity = 1._dp, resistance = 1._dp

end module setup

program charging

	use setup
	implicit none
	real(dp) :: q, t, dt, tmax

	t = 0._dp
	tmax = 5._dp
	dt = 0.1_dp
	q = 0._dp

	open(unit = 3, file='qc.data', action = 'write', status = 'replace')

	do while( t <= tmax )
		
		write(3, fmt = '(2f20.10)') t, q
! 		print *, t, q
		call rk4step(t,q,dt)

	end do


end program charging