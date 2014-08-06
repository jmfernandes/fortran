
module cf_approx

	use numtype
	implicit none
    integer,    parameter   ::  n_basis=150, steps=50
	integer, parameter :: ncf = n_basis, ncf2 =steps
    real(dp), parameter :: tiny = 1.e-30_dp, eps = 1.e-13_dp
	
	contains
	
	    subroutine taylor_cfrac(f,n,e) !	f(n): Taylor coefficients
	    
            use numtype
            implicit none
            real(dp), dimension(0:ncf,0:ncf2) :: f, e, r, s
            integer :: n, i, nl, j, k, m
        
        do m=1,ncf
            e(m,0) = f(m,0)
            e(m,1) = f(m,1)
            e(m,2) = f(m,2)/(f(m,1)+tiny)
            do i = 3, n
                e(m,i) = f(m,i)/(f(m,i-1)+tiny)
                r(m,i) = e(m,i)-e(m,i-1)
            end do
            nl = 1
            do j = 4, n
                if(nl == 1) then
                    do k = j, n
                        s(m,k) = r(m,k)/r(m,k-1)*e(m,k-1)
                    end do
                else
                    do k = j, n
                        s(m,k) = r(m,k)-r(m,k-1)+e(m,k-1)
                    end do
                end if
                e(m,j-1:n) = r(m,j-1:n)
                r(m,j:n) = s(m,j:n)
                nl = -nl
            end do

        end do
                
        end subroutine taylor_cfrac
        
        function evalcf(e,n,x) result(g1)   
        
            use numtype
            implicit none
            real(dp), dimension(0:ncf,0:ncf2) :: e
            integer :: n, i, m
            real(dp) :: x
            real(dp), dimension(0:ncf) :: g0, g1, c0, c1, d0, d1, delta


        do m=1,ncf2    
            g0(m) = e(m,0) + tiny
            c1(m) = 1+e(m,1)*x/g0(m) + tiny
            g1(m) = g0(m)*c1(m)
            d0(m) = 1
            c0(m) = c1(m)
            g0(m) = g1(m)
            do i = 2, n
                d1(m) = 1-e(m,i)*x*d0(m) + tiny
                c1(m) = 1-e(m,i)*x/c0(m) + tiny
                d1(m) = 1/d1(m)
                delta(m) = c1(m)*d1(m)
                g1(m) = g0(m)*delta(m)
                !if (abs(delta-1) < eps) return
                d0(m) = d1(m)
                c0(m) = c1(m)
                g0(m) = g1(m)
            end do
        end do
    
        end function evalcf
              
end module cf_approx
