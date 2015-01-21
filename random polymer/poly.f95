
program walking

	use NumType
	implicit none

	!===================================================

	real(dp)	:: r, r_rad, r_rad_old,r_last,x_1,x_2,y_1,y_2
	real(dp), parameter :: length=1._dp
	integer, parameter :: iter = 100
	integer		:: i, j
	real(dp) :: x(0:iter),y(0:iter)

	!===================================================

	x = 0.0_dp
	y = 0.0_dp
	r = 0.0_dp

	call init_random_seed()

	r_last = 0.0_dp

	do i = 1,iter
! 		print *, i
! 		call init_random_seed()
		r_rad_old = r*2*pi
! 		print *, r_rad_old, 'old'
		call random_number(r)
		r_rad = r*2*pi
! 		print *, r_rad, 'new'
		if (r_rad >  r_last-(pi/2) .and. r_rad < r_last+(pi/2) ) then
			x(i) = x(i-1) + length*cos(r_rad)
			y(i) = y(i-1) + length*sin(r_rad)
! 			print *, r_rad, r_rad_old , 'trigger'
			r_last = r_rad
			r_rad_old = r_last
			do j=1,i 
				if (x(j-1) < x(j)) then
					x_1=x(j-1)
					x_2=x(j)
					y_1=y(j-1)
					y_2=y(j)
				else
					x_1=x(j)
					x_2=x(j-1)
					y_1=y(j)
					y_2=y(j-1)
				end if
				if (y(j-1) < y(j)) then
					y_1=y(j-1)
					y_2=y(j)
				else
					y_1=y(j)
					y_2=y(j-1)
				end if

				if (x(i-1) > x_1 .and. x(i-1) < x_2 .and. y(i-1) > y_2 .and. y(i) < y_2) then
					print *, 'collision'
					x(i) = x(i-1)
					y(i) = y(i-1)
					EXIT
				else if (y(i-1) > y_1 .and. y(i-1) < y_2 .and. x(i-1) < x_1 .and. x(i) > x_1) then
					print *, 'collision'
					x(i) = x(i-1)
					y(i) = y(i-1)
					EXIT
				else if (y(i) > y_1 .and. y(i) < y_2 .and. x(i) > x_1 .and. x(i-1) < x_1) then
					print *, 'collision'
					x(i) = x(i-1)
					y(i) = y(i-1)
					EXIT
				else if (y(i-1) > y_1 .and. y(i-1) < y_2 .and. x(i-1) > x_1 .and. x(i) < x_1) then
					print *, 'collision'
					x(i) = x(i-1)
					y(i) = y(i-1)
					EXIT
				else if (x(i-1) > x_1 .and. x(i-1) < x_2 .and. y(i-1) < y_1 .and. y(i) > y_1) then
					print *, 'collision'
					x(i) = x(i-1)
					y(i) = y(i-1)
					EXIT
				else
					continue
				end if
			enddo
		else
			x(i) = x(i-1)
			y(i) = y(i-1)
! 			print *, r_rad, r_rad_old
		end if
! 		print *, r_rad, x, y, x**2+y**2
		write(1,*) x(i),y(i)
	enddo 

	print *, x(i-1),y(i-1)

contains 

	

end program walking