!mega important, the reason why papp likes fortran! ARRAYS!!!!


subroutine rk4step(t,y,dt)

    use setup
    implicit none
    real(dp), intent(in):: dt !t is changing in subroutine, dt is not
    real(dp), intent(inout):: t !goes in changes comes out
    real(dp), DIMENSION(n_eq), intent(inout):: y
    real(dp), dimension(n_eq) :: k1, k2, k3, k4, dy
    
    call deriv(t,   y, dt,  k1) !this is not for us
    call deriv(t+dt/2, y+k1/2, dt,  k2)
    call deriv(t+dt/2, y+k2/2, dt,  k3)
    call deriv(t+dt, y+k3, dt,  k4)
    
    dy = (k1 + 2*k2 + 2*k3 + k4)/6
    
    t = t+dt
    y = y + dy
    
    contains
        subroutine deriv(t,y,dt,k)
        
        use setup, only:dp,n_eq,hbar2,xm,energy,l

        real(dp), intent(in):: t, dt
        real(dp), dimension(n_eq), intent(in):: y
        real(dp), dimension(n_eq) :: k, f

      
        

        f(1) = y(2)
        f(2) = -2*xm/hbar2*(energy &
            -hbar2/(2*xm)*l*(l+1)/t**2 &
            -potential(t))*y(1)
        f(3) = y(1)**2
        
        k= dt*f
        
        end subroutine deriv

        function potential(x) result(v)
            use setup, only : dp
            implicit none
            real(dp), intent(in) :: x
            real(dp) :: v

            v = 120._dp*exp(-2*x) - 100._dp*exp(-x**2)
        end function potential
    
end subroutine rk4step