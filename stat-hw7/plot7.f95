
program extrap_test

	use numtype
    use chebyshev
	implicit none
    integer :: n, i, max, min,nch
    real(dp) :: x, dx, yy,j, min,ymax offset
    
    max=3
    min=-3
    n = 100
    dx = abs(max-min)/(n*1._dp)
    j = min

    open(unit = 1, file='tanh.data', action = 'write', status = 'replace')
    open(unit = 2, file='mx.data', action = 'write', status = 'replace')
    open(unit = 3, file='tanh2.data', action = 'write', status = 'replace')
    open(unit = 4, file='mx2.data', action = 'write', status = 'replace')

    do while (j <= max)
        write(1,*) j, func(j)
        write(2,*) j, func2(j)
        write(3,*) j, func4(j)
        write(4,*) j, func5(j)
!         print *, j, func(j), func3(j)
        j = j + dx
    end do


    CLOSE(1)
    CLOSE(2)
    CLOSE(3)
    CLOSE(4)
    
    nch = 50

    offset = 0
    
    call chebyex(func3, nch, cheb, min, max)

    call chebyzero(nch, cheb, min, max, z0, iz0)
    
!     print *, iz0, z0(1:iz0)*180._dp/pi



    do i=1,iz0
        offset = func3(z0(i))
!         print *, offset, z0(i), 'heyeye'
        write(3,*) z0(i), func(z0(i))
    end do

contains 


    function func(x) result(f)

    use numtype , only : dp 
    implicit none
    real, parameter :: Tc = 1._dp, T = 0.5_dp

    real(dp) :: x, f

    f = tanh(x*Tc/T)

    end function func

    function func2(x) result(f)

    use numtype , only : dp 
    implicit none

    real(dp) :: x, f

    f = x

    end function func2

    function func3(x) result(f)

    use numtype , only : dp 
    implicit none

    real(dp) :: x, f

    f = func2(x) - func(x)

    end function func3

    function func4(x) result(f)

    use numtype, only: dp
    implicit none

    real(dp) :: x,f
    real, parameter :: Tc = 1._dp, T = 0.5_dp

    f = x*(Tc/T**2)*(1-tanh(x*Tc/T)**2)

    end function func4

    function func5(x) result(f)

    use numtype, only: dp
    implicit none

    real(dp) :: x,f
    real, parameter :: Tc = 1._dp, T = 0.5_dp

    f = -T*x**2

    end function func5
    


end program extrap_test


