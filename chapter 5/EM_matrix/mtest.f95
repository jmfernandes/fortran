
module setup

	use NumType
	implicit none
! 	real(dp), parameter :: n=2,m=3,l=4
	integer, parameter :: ndim=3, lwork=2*ndim+1

end module setup

program matrix

	use setup
	implicit none

	real(dp), dimension(ndim,ndim) :: f(0:ndim,0:ndim), l(0:ndim,0:ndim)
	real(dp), dimension(3) :: b,e
	real(dp) :: gamma, beta
	integer :: n,m

	e = (/1,2,3/)
	b = (/4,5,6/)
	f = 0._dp
	l = 0._dp

	gamma = 1._dp
	beta = 1._dp

	!!create F matrix
	do n=0,ndim
		do m=0,ndim
			if (n==0) then
				f(n,m) = e(m)
			else
				f(n,n+1) = b(2*ndim-(2*n+1))
				f(n,n+2) = -b(2*ndim-(2*n+2))

			end if
			
			if (n==m) then
				continue
			else
				f(m,n) = -f(n,m)
			end if

 		end do
	end do

	!!create L matrix
	do n=0,ndim
        print '(5f12.3)', f(n,0:ndim)
    end do
		
	do n=0,ndim
		do m=0,ndim
			if (n==m) then
				l(n,m) = 1
			end if
			l(0,0) = gamma
			l(0,1) = iic*gamma*beta
			l(1,0) = -iic*gamma*beta
			l(1,1) = gamma
		end do
			
	end do

	print *, '--------------------------------------'

	do n=0,ndim
        print '(5f12.3)', l(n,0:ndim)
    end do

! 	real(dp), dimension(ndim) :: w
! 	integer :: i, nn, info, j, ipiv(ndim)
! 	integer, parameter :: lwork = 5*ndim
! 	complex(dp) :: work(lwork)
! 	real(dp) :: rwork(lwork)


! 	print *, A(1:3,1:3)
! 	print *, '------------'
! 	print *, B(1:3,1:3)
! 	print *, '______________'

! 	do i=1,3
! 		print '(10f12.4)',A(i,1:3)
! 	end do
! 	print *, '_____________'
! 	do i=1,3
! 		print '(10f12.4)', B(i,1:3)
! 	end do

! 	c(1:3,1:3) = conjg(transpose(a(1:3,1:3)))
! 	d(1:3,1:3) = conjg(transpose(b(1:3,1:3)))

! 	print *, '_____hermiticity test________'
! 	do i = 1,3
! 		print '(a,10f12.4)','a matrix',C(i,1:3)-A(i,1:3)
! 	end do

! 	print *, '_____________'
! 	do i = 1,3
! 		print '(a,10f12.4)','b matrix',D(i,1:3)-B(i,1:3)
! 	end do

! 	print *, '______commutator test_______'

! 	c(1:3,1:3) = matmul(a(1:3,1:3),b(1:3,1:3))
! 	d(1:3,1:3) = matmul(b(1:3,1:3),a(1:3,1:3))
! 	e(1:3,1:3) = c(1:3,1:3)-d(1:3,1:3)

! 	do i = 1,3
! 		print '(a,10f12.4)','[a,b]',c(i,1:3)-d(i,1:3)
! 	end do

! 	s = 0._dp
! 	do i=1,3
! 		s = s + e(i,i)
! 	end do

! 	print *, '_________________'

! 	print '(a,10f12.4)','Tr([a,b])',s

! 	print *, '__________eigenvalues____________'

! 	nn = 3
! 	info = 0
! 	e(1:nn,1:nn) = a(1:nn,1:nn)

! 	call zheev('v','u',nn,e,ndim,w,work,lwork,rwork,info)

! 	print *, 'info=', info

! 	do i = 1,3
! 		print '(a,i5,f15.8,a,6f12.4)','eigenvalues',i,w(i),'  vector  ', e(1:nn,i) 
! 	end do

! 	print *, '------- Orthogonality ------- '
        
!         do j = 1,3
!             do i = 1,3
!                 print *, i, j, dot_product(e(1:3,i) ,e(1:3,j) ) !determines if matrix is orthogonal
            
!             end do
!         end do
        
        
!           print *, '------- completeness & spectral decomp------- '
        
!         do j = 1,3
!             do i = 1,3
!                 f(i,j)= dot_product(e(i,1:3) ,e(j,1:3)*w(1:3))  !determines if matrix is complete
!             end do
!         end do 
        
!         Do i= 1,3
!             print '(10f7.2)', f(i,1:3)
!         end do 
        
!         print *, '+++++++++++++++++++++++++++++++++++'
        
!         c(1:3,1:3) = conjg(transpose(e(1:3,1:3)))
        
!         f(1:3,1:3) = matmul(c(1:3,1:3),e(1:3,1:3))
!         f(1:3,1:3) = matmul(e(1:3,1:3),c(1:3,1:3))
        
!         Do i= 1,3
!             print '(10f7.2)', f(i,1:3)
!         end do 
        
!         print *, '++++++++++++++++++++++++++++++++++++++++++++'
        
!         f(1:3,1:3) = matmul(a(1:3,1:3),e(1:3,1:3))
!         d(1:3,1:3) = matmul(c(1:3,1:3),f(1:3,1:3))
        
!         Do i= 1,3
!             print '(10f7.2)', d(i,1:3)
!         end do 
        
        
!         print *, '++++++++ inversions of e++++++++++++++ '
!         info = 0
!         call zgetrf(nn,nn,e, ndim, ipiv,info)
        
!         if(info/= 0) stop 'stop'
        
!         call zgetri(nn,e,ndim,ipiv,work, lwork,info)
        
!         if(info/= 0) stop 'stop'
        
!         Do i= 1,3
!             print '(6f7.2,5x,6f7.2)', c(i,1:3)-e(i,1:3)
!         end do 


end program matrix