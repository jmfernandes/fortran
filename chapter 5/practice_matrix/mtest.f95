
module setup

	use NumType
	use matrix
	implicit none
! 	real(dp), parameter :: n=2,m=3,l=4
	integer, parameter :: ndim=2

end module setup

program mattest

	use setup
	implicit none

	real(dp), dimension(ndim,ndim) :: a,b,c,fac,det
	real(dp), dimension(3) :: e,f,g
! 	real(dp) :: det(2)
	integer :: i, j

	a(1:2,1:2) = reshape((/ 1,2,2,1 /), &
		(/2,2/))



	b(1,1) = 4
	b(2,1) = 1
	b(1,2) = 1
	b(2,2) = 4

	e = (/ 1,2,2 /)
	f = (/ -2,4,3 /)

	g = dot_product(e,f)
	print *, e
	print *, f
	print *, g(1)

	print *, '--------------------------------'

! 	do n=1,2
! 		do m=1,2
! 			b(n,m) = 1
! 		end do
! 	end do
! 	b(1:2,1:2) = reshape((/ 4,2,1,1 /), &
! 		(/2,2/))

! 	do i=1,2
! 		print '(10f12.4)', a(i,1:2)
! 	end do

! 	print *, '---------------'

! 	do i=1,2
! 		print '(10f12.4)', b(i,1:2)
! 	end do
	

! 	c(1:2,1:2) = matmul(a(1:2,1:2),b(1:2,1:2))
! 	call dsydet(2,b,2,det)

! 	print *, '---------------'

	do i=1,2
		print '(10f12.4)', a(i,1:2)
	end do


	print *, '-------------b-------------------'

    do i=1,2
		print '(10f12.4)', b(i,1:2)
	end do

! 	print *, '------------inner product--------------------'

! 	do j = 1,2
!             do i = 1,2
!                 print *, i, j, dot_product(a(1:2,i) ,b(1:2,j) ) !determines if matrix is orthogonal
            
!             end do
!         end do

! 	do i=1,2
! 		print '(10f12.4)', c(i,1:2)
! 	end do

	call dsydet(ndim,b,ndim,det)

	print *, '------------det--------------------'

	do i=1,2
		print '(10f12.4)', det(i,1:2)
	end do

	

	

! 	print *, -b(1,1)*det(1,1)


! 	call dgei(2,b,2)

! 	do i=1,2
! 		print '(10f12.4)', b(i,1:2)
! 	end do

! 	fac = b(2,2)*b(1,1)-b(2,1)*b(1,2)
! 	print *, fac

	


end program mattest