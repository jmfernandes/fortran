
module chebyshev

	use numtype
	implicit none
	integer, parameter :: maxch = 50
	complex(dp), dimension(0:maxch) :: cheb, chder
	complex(dp), dimension(maxch) :: z0
	integer :: iz0
	
	contains
        
        subroutine chebyex(func,n,a,ya,yb) 
        !   func([ya,yb]) = sum_{i=0}^n  a_i T_i
    
            real(dp), external :: func
            integer :: n
            complex(dp), dimension(0:maxch) :: f, a
            real(dp) :: ya, yb, aa, bb
            complex(dp) :: x, ss
            integer :: i, j
    
            if ( n > maxch ) stop '  n > maxch '
            aa = (yb-ya)/2; bb = (yb+ya)/2
            do i = 0, n
                x = cos(pi/(n+1)*(i+0.5_dp))
                f(i) = func(aa*x+bb)
            end do
            do j = 0, n
                ss = 0._dp
                do i = 0, n
                    ss = ss + &
                        f(i)*cos((pi/(n+1))*j*(i+0.5_dp))
                end do
                a(j) = 2._dp*ss/(n+1)
            end do
            a(0) = 0.5_dp*a(0)
            
        end subroutine chebyex
      
        subroutine chebyderiv(a,n,der,ya,yb) ! 

            integer :: n
            real(dp) :: ya, yb
            complex(dp) :: a(0:maxch), der(0:maxch)
            integer :: j
    
            der(n) = 0._dp; der(n-1) = 2*n*a(n)
            do j = n-1, 1, -1
                der(j-1) = der(j+1)+2*j*a(j)
            end do
            der(0) = der(0)/2
            der(0:n-1) = der(0:n-1)*2/(yb-ya)

        end subroutine chebyderiv
      
        function cheby(y,a,n,ya,yb) result(t) 
        ! func(y) =  sum_{i=0}^n  a_i T_i (x)

            implicit none
            integer :: n
            real(dp) :: ya, yb, aa, bb
            complex(dp) :: a(0:maxch)
            complex(dp) :: x, y, t, y0, y1
            integer :: k
    
            aa = (yb-ya)/2; bb = (yb+ya)/2
            x = (y-bb)/aa
            y1 = 0._dp; y0 = a(n)
            do k = n-1, 0, -1
                t = y1; y1 = y0
                y0 = a(k)+2*x*y1-t
            end do
            t = y0-x*y1

        end function cheby
                
        subroutine chebyzero(n,a,ya,yb,z0,iz0) 
        
            integer :: n, iz0
            real(dp) :: ya, yb
            complex(dp), dimension(0:maxch) :: a
            integer :: j
            complex(dp), dimension(maxch) :: w0,z0
            
            
            call boyd(n,a,w0)
            w0(1:n) = w0(1:n)*(yb-ya)/2+(yb+ya)/2

            
            iz0 = 0
            do j = 1, n
                if( ya <= dble(w0(j)) .and. dble(w0(j)) <= yb ) then 
                        iz0 = iz0+1;  z0(iz0) = w0(j)
                end if
            end do
                        
            contains
                
                subroutine boyd(n,a,w) !finds the zeros of a function using companion matrix
                
                    integer :: n, j, ie
                    complex(dp) :: a(0:maxch)
                    complex(dp) :: w(maxch)
                    integer, parameter :: lwork=4*maxch
                    complex(dp) :: aamat(maxch,maxch),  &
                        work(lwork), rwork(lwork), &
                        vl(1), vr(1)
                    
                    if (abs(a(n)) == 0._dp) stop 'a(n)=0'
                    aamat(1:n,1:n) = 0._dp
                    aamat(1,2) = 1._dp
                    do j = 2, n-1
                        aamat(j,j-1) = 0.5_dp
                        aamat(j,j+1) = 0.5_dp
                    end do
                    aamat(n,1:n) = -a(0:n-1)/(2*a(n))
                    aamat(n,n-1) = aamat(n,n-1) + 0.5_dp
                    
                    ie = 0
                    call zgeev('n','n',n,aamat,maxch,w,&
                        vl,1,vr,1,work,lwork,rwork,ie)                      
                    if( ie /= 0 ) stop ' boyd: ie /= 0 ' 
                            
                end subroutine boyd
                
        end subroutine chebyzero
        
        subroutine root_polish(func,zz,dz,eps,maxf)
        
            complex(dp), external :: func
            real(dp) :: eps
            complex(dp) :: zz, dz, z1, z2, z3, &
                f1, f2, f3, a12, a23, a31
            integer :: i, maxf
            
            z1 = zz+dz;   f1 = func(z1)
            z2 = zz-dz;   f2 = func(z2)
            z3 = zz;      f3 = func(z3)

            do i = 1,maxf
                a23 = (z2-z3)*f2*f3
                a31 = (z3-z1)*f1*f3
                a12 = (z1-z2)*f1*f2
                zz = (z1*a23+z2*a31+z3*a12)/(a23+a31+a12)
                if ( abs(zz-z3) < eps ) exit
                z1 = z2;  f1 = f2
                z2 = z3;  f2 = f3
                z3 = zz;  f3 = func(z3)
            end do
        
        end subroutine root_polish
      
end module chebyshev


