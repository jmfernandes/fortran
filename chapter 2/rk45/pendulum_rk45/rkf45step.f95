
subroutine rkf45step(t,y,h)  ! 4-th order Runge-Kutta step
	
	use setup, only : dp, n_eq
	implicit none
	real(dp), intent(inout) :: t, h
	real(dp), dimension(n_eq), intent(inout) :: y 
	real(dp), dimension(n_eq) :: k1, k2, k3, k4, k5, k6, y1, y2
	real(dp), parameter :: epsilon = 1.e-6_dp, tiny = 1.e-20_dp
	real(dp) :: rr, delta
		
	call deriv(t,	     h,	 y,			                                                k1)
	call deriv(t+h/4,	 h,  y+ k1/4, 	                                                k2)
	call deriv(t+3*h/8,	 h,  y+ (3*k1+9*k2)/32,	                                        k3)
	call deriv(t+12*h/13,h,	 y+ (1932*k1-7200*k2+7296*k3)/2197 ,	                    k4)	
	call deriv(t+h,	     h,  y+ (439*k1/216-8*k2+3680*k3/513-845*k4/4104),	            k5)	
	call deriv(t+h/2,	 h,  y+ (-8*k1/27 +2*k2-3544*k3/2565 +1859*k4/4104 -11*k5/40),	k6)
	
    y1 = y + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104  - k5/5
    y2 = y + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55
    
    rr = sqrt(dot_product(y1-y2,y1-y2))/h + tiny
    
    if ( rr < epsilon ) then
        t = t + h
        y = y1
        delta = 0.92_dp * (epsilon/rr)**(0.2_dp)
        h = delta*h
        write (unit = 3,fmt='(3f20.10)') t, y(1)
    else
        delta = 0.92_dp * (epsilon/rr)**(0.25_dp)
        h = delta*h
    end if
	
    	
    contains
    
        subroutine deriv(t,h,y,k)   ! derivative

	        use setup, only : dp, n_eq, g, lenght
	        implicit none
	        real(dp), intent(in) :: t, h
	        real(dp), dimension(n_eq), intent(in) :: y
	        real(dp), dimension(n_eq) :: f, k

	        f(1)= y(2)
	        f(2)= -g/lenght*sin(y(1))
	        
	        k(1:n_eq) = h*f(1:n_eq)
	        
        end subroutine deriv
	
end subroutine rkf45step
