
program extrap_test

	use numtype
	use thiele_approx
	implicit none
    real(dp), dimension(maxpt) :: zn, fn, an
    integer :: n, i
    real(dp), external :: func
    real(dp) :: x, dx, yy
    
    n = 20
    dx = 0.5_dp
    do i = n, 1, -1
        zn(i) = i*dx
        fn(i) = func(zn(i))
        print *,zn(i),fn(i)
        write(1,*) zn(i),fn(i)
    end do
    print *,'------------'
    
    call thiele_coef( n, zn, fn, an) 
    
    do i = n/2, -n/2, -1
        x = dx/2 + i*dx + 1.e-30_dp
        yy =  thiele_cf( x, n, zn, an)
        fn(i) = func(zn(i))
        print *,x, func(x), yy
        write(2,*) x, func(x)
        write(3,*) x, yy
    end do
    
    print *,
    
    print *,'------------------------------------------'

end program extrap_test

function func(x) result(f)
	    
	use numtype, only: dp,pi,iic
	implicit none
	real(dp) :: x, f, z, alpha
	    
!	f = x*cos(x) ! f = 1/(1+x**2) - x**3

    !f = sin(x)/x
!     f = cos(x)/x
    alpha = pi/2._dp
    z = x*(cos(alpha)+iic*sin(alpha))
    f = pi*(1/(tan(pi*x)))

	    
end function func

