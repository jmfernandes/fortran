
program walking

	use NumType
	implicit none

	!===================================================

	integer, parameter :: flips = 3, sides = 3, repeats=100
	real(dp)	:: division,r
	integer		:: i,j,x1,x2,x3,x(1:flips)

	!===================================================

	division = 1.0/sides

	x1=0
	x2=0
	x3=0

	call init_random_seed()
	'for bosons'
	do i=0,repeats
		do j=1,flips
			call random_number(r)
			r = r*sides
			r = ceiling(r)
			if (j==1) then
				x(1) = r
			else 
				do while (r < x(j-1))
					call random_number(r)
					r = r*sides
					r = ceiling(r)
! 					print *, r, x(j-1)
				end do
				x(j) = r
			end if 
! 			print *, j, r
		end do
		print *, x
	! 	print *, r
	!=========== check for fairness
! 		if (r == 1) then 
! 			x1 = x1 +1
! 		else if (r == 2 ) then
! 			x2 = x2 + 1
! 		else
! 			x3 = x3 + 1
! 		end if
	!=========================

	end do
! 	print *, x1, x2,x3, 'is it fair?'
! 	do i=1,flips
! 		do j=1,sides
! 			print *, j
! 		end do
! 	end do


end program walking