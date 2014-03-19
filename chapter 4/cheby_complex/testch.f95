
program chebytest

    use numtype
    use chebyshev_c
    implicit none
    real(dp) :: ya, yb, dx, eps
    complex(dp) x, fc, fx
    integer :: n, np, i

    
    ya = -5
    yb = 5
    
    n = 10
    
    call chebyex(our_f,n,cheb,ya,yb) 

    


    np = 50
    dx = (yb-ya)/np
    do i = 0, np
        x = ya + i*dx
        fc = cheby(x,cheb,n,ya,yb)
        fx = our_f(x)
        write(1,*) dble(x), 0._dp, dble(fx), imag(fx)
        write(2,*) dble(x), 0._dp, dble(fc), imag(fc)
    end do


    call chebyzero(n,cheb,ya,yb,z0,iz0)
    print *,iz0
    print *,z0(1:iz0)
    
 



    contains 
    
        function our_f (x) result(f)
        
            implicit none 
            complex(dp) :: x, f, z
            real (dp) :: alpha
            
            alpha = pi/2._dp
            z = x*(cos(alpha)+iic*sin(alpha))
            f = z**2 + 1
        
        end function our_f


end program chebytest

