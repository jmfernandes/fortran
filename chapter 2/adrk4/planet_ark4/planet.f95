!header


module setup

	use NumType
	implicit none
	real(dp), parameter :: gravity = 6.673e-11_dp, mass_sun = 1.9891e+30_dp
	real(dp), parameter :: mass_earth = 5.9736e+24_dp
	integer, parameter :: number_of_equations = 6

end module setup

module gnuplot_fortran

	use NumType
	implicit none
	contains
	subroutine plot2d(x,y1,y2,y3,y4,y5,y6)
		real(dp), intent(in) :: x,y1,y2,y3,y4,y5,y6
		!integer :: size_x, size_y, i
		!size_x = size(x)
		!size_y = size(y1)
		open(unit=1, file = 'data.data')
		write(1,*) x, ' ', y1, ' ', y2, ' ', y3, ' ', y4, ' ', y5, ' ', y6

	end subroutine plot2d


end module gnuplot_fortran

program orbit

	use setup
	use gnuplot_fortran
	implicit none
	real(dp) :: y(number_of_equations), t, dt, tmax
	real(dp) :: X0, Y0, Z0, VX0, VY0, VZ0


	X0 = 1.496e+11_dp
	Y0 = 0._dp
	Z0 = 0._dp

	VX0 = 0._dp
	VY0 = 29.783e+3_dp
	VZ0 = 30.783e+2_dp

	t = 0._dp
	tmax = 7*(60*60*24*365) !years
	dt = 7*(24*60*60) !days
	
	y(1) = X0
	y(2) = Y0
	y(3) = Z0
	y(4) = VX0*mass_earth
	y(5) = VY0*mass_earth
	y(6) = VZ0*mass_earth

	!open(unit = 3, file='earth.data', action = 'write', status = 'replace')


	do while( t <= tmax )
		
		!write(3, *) y(1), y(2)
		call ark4(t,y,dt)
		call plot2d(t,y(1),y(2),y(3),y(4),y(5),y(6))

	end do


end program orbit