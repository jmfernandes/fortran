
program walking

	use NumType
	implicit none

!===================================================

	integer, parameter	::	steps = 1, repeat =100000
	real(dp), parameter ::	bins = 2*steps/50._dp
	real(dp)			::	histogram(0:49)
	real(dp)			::	r, r_rad, x,y, length,x_0,y_0, &
							distance, sum, counting(0:2,0:steps)
	integer				::	i, j,k

!===================================================
	!seed

	call init_random_seed()

!===================================================
 	!regular walking code

! 	length = 1._dp

! 	x = 0.0_dp
! 	y = 0.0_dp

! 	do i = 0,1000000
! 		call random_number(r)
! 		r_rad = r*2*pi
! 		x = x + length*cos(r_rad)
! 		y = y + length*sin(r_rad)
! 		write(1,*) x,y
! 	enddo 


!===================================================
!photon in sun

! 	length =5._dp
! 	distance = 0.0_dp
! 	i = 0
! 	sum = 0
! 	do j=0,100
! 		i=0
! 		distance = 0.0_dp
! 		x = 0.0_dp
! 		y = 0.0_dp

! 		do while (distance <= 500000000)
! 	   		call random_number(r)
! 			r_rad = r*2*pi
! 			x = x + length*cos(r_rad)
! 			y = y + length*sin(r_rad)
! 			distance = x**2 + y**2
! 			i = i+1
! 	! 		print *, 'the distance is ', distance
! 		end do
! 	print *, j,i
! 	! write(1,*) j,i
! 	sum = sum + i
! 	end do

! 	print *, sum/j


!===================================================

! 	length =1._dp
! 	do j=0,10000
! 		i=0
! 		distance = 0.0_dp
! 		x = 0.0_dp
! 		y = 0.0_dp

! 		do while (i <= 100000)
! 	   		call random_number(r)
! 			r_rad = r*2*pi
! 			x = x + length*cos(r_rad)
! 			y = y + length*sin(r_rad)
! 			i = i+1
! 	! 		print *, 'the distance is ', distance
! 		end do
! 	! print *, x,y, 'endpoint'
! 	write(1,*) x,y

! 	end do

!===================================================
!x vs t

! 	counting(0:2,0:steps) = 0._dp
! 	length =1._dp
! 	do j=0,repeat
! 		i=0
! 		distance = 0.0_dp
! 		x = 0.0_dp
! 		y = 0.0_dp

! 		do while (i <= steps)
! 	   		call random_number(r)
! 			r_rad = r*2*pi
! 			x = x + length*cos(r_rad)
! 			y = y + length*sin(r_rad)
! 			i = i+1
! 	! 		write(1,*) i, x
! 			counting(1,i) = i
! 			counting(2,i) = counting(2,i) + x
! 	! 		print *, 'the distance is ', distance
! 		end do
! 	! print *, x,y, 'endpoint'
! 	print *, j
! 	end do

! 	do i=0,steps
! 		write(1,*) counting(1,i), counting(2,i)
! 	end do

!===================================================
!square plot

! 	counting(0:2,0:steps) = 0._dp
! 	length =1._dp
! 	do j=0,repeat
! 		i=0
! 		distance = 0.0_dp
! 		x = 0.0_dp
! 		y = 0.0_dp

! 		do while (i <= 0)
! 	   		call random_number(r)
! ! 			r_rad = (2._dp*r-1._dp)/2._dp
! 			x = x + (2._dp*r-1._dp)/2._dp
! 			call random_number(r)
! 			y = y + (2._dp*r-1._dp)/2._dp
! 			i = i+1
! 		end do

! 		write(1,*) x,y

! 	end do

! 	length =1._dp
! 	do j=0,repeat
! 		i=0
! 		distance = 0.0_dp
! 		x = 0.0_dp
! 		y = 0.0_dp

! 		do while (i <= 9)
! 	   		call random_number(r)
! ! 			r_rad = (2._dp*r-1._dp)/2._dp
! 			x = x + (2._dp*r-1._dp)/2._dp
! 			call random_number(r)
! 			y = y + (2._dp*r-1._dp)/2._dp
! 			i = i+1
! 		end do

! 		write(2,*) x,y

! 	end do

!===================================================
!histogram plot

! 	counting(0:2,0:steps) = 0._dp
	length =1._dp
	histogram=0._dp
	do j=0,repeat
		i=0
		distance = 0.0_dp
		x = 0.0_dp
		y = 0.0_dp
! 		bins = 2*steps/4._dp

		do i=1,steps
	   		call random_number(r)
! 			r_rad = (2._dp*r-1._dp)/2._dp
			x = x + (2._dp*r-1._dp)/2._dp
			call random_number(r)
			y = y + (2._dp*r-1._dp)/2._dp

			do k=0,49
				if (x > -steps+(k*bins) .AND. x<= -steps+(k*bins)+bins) then
					histogram(k)=histogram(k)+1
				end if
			end do

		end do

		

	end do

	do i=0,49
		write(3,*) -steps+(i*bins), histogram(i)
	end do

end program walking