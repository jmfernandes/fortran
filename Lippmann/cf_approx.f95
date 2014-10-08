
module cf_approx

	use numtype
	implicit none
	integer, parameter :: ncf = 50
    real(dp), parameter :: tiny = 1.e-30_dp, eps = 1.e-13_dp
	
	contains
	
	    subroutine taylor_cfrac(f,n,e) !	f(n): Taylor coefficients
	    
            use numtype
            implicit none
            real(dp), dimension(0:ncf) :: f, e, r, s
            integer :: n, i, nl, j, k
        
            e(0) = f(0)
            e(1) = f(1)
            e(2) = f(2)/(f(1)+tiny)
            do i = 3, n
                e(i) = f(i)/(f(i-1)+tiny)
                r(i) = e(i)-e(i-1)
            end do
            nl = 1
            do j = 4, n
                if(nl == 1) then
                    do k = j, n
                        s(k) = r(k)/r(k-1)*e(k-1)
                    end do
                else
                    do k = j, n
                        s(k) = r(k)-r(k-1)+e(k-1)
                    end do
                end if
                e(j-1:n) = r(j-1:n)
                r(j:n) = s(j:n)
                nl = -nl
            end do
                
        end subroutine taylor_cfrac
        
        function evalcf(e,n,x) result(g1)   
        
            use numtype
            implicit none
            real(dp), dimension(0:ncf) :: e
            integer :: n, i
            real(dp) :: x, g0, g1, c0, c1, d0, d1, delta
        
            g0 = e(0) + tiny
            c1 = 1+e(1)*x/g0 + tiny
            g1 = g0*c1
            d0 = 1
            c0 = c1
            g0 = g1
            do i = 2, n
                d1 = 1-e(i)*x*d0 + tiny
                c1 = 1-e(i)*x/c0 + tiny
                d1 = 1/d1
                delta = c1*d1
                g1 = g0*delta
                !if (abs(delta-1) < eps) return
                d0 = d1
                c0 = c1
                g0 = g1
            end do
    
        end function evalcf
              
end module cf_approx
