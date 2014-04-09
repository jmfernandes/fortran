
program chebytest

    use numtype
    use chebyshev_c
    implicit none
    real(dp) :: ya, yb, dx, eps
    complex(dp) :: chder2(0:maxch)
    complex(dp) x, fc, fx, fz, fs, zz, dz
    integer :: n, np, i

    
    ya = -5
    yb = 5
    
    n = 10
    
    call chebyex(our_f,n,cheb,ya,yb) 

    call chebyderiv(cheb,n,chder,ya,yb)
    call chebyderiv(chder,n-1,chder2,ya,yb)


    np = 50
    dx = (yb-ya)/np
    do i = 0, np
        x = ya + i*dx
        fc = cheby(x,cheb,n,ya,yb)
        fx = our_f(x)
        fz = cheby(x,chder,n,ya,yb)
        write(1,*) dble(x), 0._dp, dble(fx), imag(fx)
        write(2,*) dble(x), 0._dp, dble(fc), imag(fc)
        write(3,*) dble(x), 0._dp, dble(fz), imag(fz)
    end do


    call chebyzero(n,cheb,ya,yb,z0,iz0)
    print *,iz0
    print *, '_____________'
    print *,z0(1:iz0)
    print *, '============='



    eps = 1.e-15_dp
    dz = 0.01_dp


    do i = 1, iz0
        zz = z0(i)
        call root_polish(dfunc,zz,dz,eps,10) 
        fs = cheby(z0(i),chder2,n-2,ya,yb)
!         print *, 'df(x)/dx = 0',i,z0(i),zz,fs
    end do
!     print *, '_____________'
!     print *, zz




    contains 
    
        function our_f (x) result(f)
        
            implicit none 
            complex(dp) :: x, f, z
            real (dp) :: alpha
            
!           alpha = pi/2._dp
            alpha = 0._dp
            z = x*(cos(alpha)+iic*sin(alpha))
            f = (z-1)**2 - 1
        
        end function our_f

        function dfunc(x) result(df)

            implicit none
            complex(dp) :: x, df, z
            real(dp) :: alpha

!           alpha = pi/2._dp !imaginary
            alpha = 0._dp
            z = x*(cos(alpha)+iic*sin(alpha))
            df = 2*(z-1)

        end function


end program chebytest

