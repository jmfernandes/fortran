
program cf_test

	use numtype
	use cf_approx
	implicit none
	real(dp), dimension(0:ncf) :: taylor, cf
	integer :: n, i
	real(dp) :: sign, x
	
	n = 50
	x = 2._dp

!	Taylor coefficients of logarithmic function

! 	taylor(0) = 0._dp
! 	sign = 1._dp
! 	do i = 1, n
! 	   taylor(i) = sign/i
! 	   print*, i,taylor(i)
! 	   sign = -sign
! 	end do
	
! 	print *,'  taylor sum ',x,'=', horner(taylor,n,x)
	
! 	call taylor_cfrac(taylor,n,cf)	

! 	print *,'          cf ',x,' =',evalcf(cf,n,x)
!     print *,'         log ',x,'=',log(1+x)	

    print *,'--------------------------'
!	Taylor coefficients of exponential function
	
! 	taylor(0) = 1._dp	
! 	taylor(1) = 1_dp
! 	do i = 2 , n
! 	   taylor(i) = taylor(i-1)/i
! 	   print*, i,taylor(i)
! 	end do   
	
! 	print *,'  taylor sum ',x,'=', horner(taylor,n,x)
	
! 	call taylor_cfrac(taylor,n,cf)	

! 	print *,'          cf ',x,'=',evalcf(cf,n,x)
!     print *,'         epx ',x,'=',exp(x)	

	taylor(0) = 0._dp	
	taylor(1) = 1_dp
	sign = -1._dp
	do i = 2 , n
	   taylor(i) = sign*taylor(i-1)/((2*i-1)*(2*i-2))
	   print *, i,taylor(i)
	   sign = -sign
	end do   
	
	print *,'  taylor sum ',x,'=', horner(taylor,n,x)
	
	call taylor_cfrac(taylor,n,cf)	

	print *,'          cf ',x,'=',evalcf(cf,n,x)
    print *,'         epx ',x,'=',sin(x)

    contains
    
        function horner(f,n,x) result(y)

            implicit none
            real(dp), dimension(0:ncf) :: f
            integer :: n, i
            real(dp) :: x, y
            
            y = f(n)
            do i = n-1, 0, -1
                y = f(i) + x*y !if the x increses by factor of 2, need x*x
            end do

        end function horner

end program cf_test

