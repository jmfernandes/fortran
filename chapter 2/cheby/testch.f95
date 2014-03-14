
program chebytest

    use numtype
    use chebyshev
    implicit none
    real(dp) :: ya, yb, dx, eps
    complex(dp) :: x, fc, fx
    integer :: n, np, i


    ya = -5
    yb = 5
    n = 10

    call chebyex(ourf,n,cheb,ya,yb)

    np = 50
    dx = (yb-ya)/np
    do i=0,np
        x = ya + i*dx
        fc = cheby(x,cheb,n,ya,yb) 
        fx = ourf(x)
        write(1,*) x, 0._dp, fx, fc
    end do

    call chebyzero(n,cheb,ya,yb,z0,iz0)
    print *, iz0
    print *, z0(1:iz0)

    contains

        function ourf(x) result(f)

            implicit none
            complex(dp) :: x,f,z
            real(dp) :: alpha


            alpha = pi/2._dp
            z = x*(cos(alpha) + iic*sin(alpha))
            f = x**2 + 1

        end function ourf



end program chebytest