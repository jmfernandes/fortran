!this program shoots a ball

module gnuplot_fortran
	implicit none
	contains
	subroutine plot2d(x,y1,y2,y3)
		real, intent(in), dimension(:) :: x, y1, y2, y3
		integer :: size_x, size_y, i
		size_x = size(x)
		size_y = size(y1)
		if (size_x /= size_y) then
			print *, "Array size mismatch"
		else
			open(unit=1, file = 'data.dat')
			do i = 1, size(x)
				write(1,*) x(i), ' ', y1(i), ' ', y2(i), ' ', y3(i)
			end do
		end if

	end subroutine plot2d


end module gnuplot_fortran

program plotter 
	use gnuplot_fortran
	implicit none

	integer, parameter :: n = 100
	real, dimension(0:n) :: x, y1,y2,y3
	real :: x_start = 0.0, x_end =20, dx
	integer :: i

	!make x array
	dx = (x_end - x_start)/n
	x(0:n) = [(i*dx, i=0, n)] 

	!make y array
	y1= sin(x) / (x + 1)
	y2= 5*sin(x) / (x + 1)
	y3= 10*sin(x) / (2*x + 1)

	!generate data
	call plot2d(x,y1,y2,y3)

end program plotter


! module setup

! 	use NumType
! 	implicit none
! 	real(dp), parameter :: g = 9.81_dp
! 	real(dp) :: h0, v0, alpha

! end module setup

! program ball_path 

! 	use setup
! 	implicit none
! 	real(dp) :: t, dt, x, y, vx, vy, vx0, vy0, x0, y0

! 	!!define variables
! 	h0 = 100._dp
! 	v0 = 15._dp
! 	alpha = 36*(pi/180) !convert from degrees to radians
	
! 	t = 0._dp
! 	dt = 0.01_dp

! 	x0 = 0._dp
! 	y0 = h0

! 	vx0 = v0*cos(alpha)
! 	vy0 = v0*sin(alpha)
	

! 	!!calculate

! 	do while (y >= 0._dp)

! 		vx = vx0
! 		vy = vy0 - g*t
! 		x = vx*t
! 		y = y0 + vy0*t - g/2 * t**2

! 		print *, g
! 		t = t + dt

! 	end do


! end program ball_path