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
        
        real(dp), intent(in):: t, dt
        real(dp), dimension(n_eq), intent(in):: y
        real(dp), dimension(n_eq) :: k, f
        
        real(dp) :: v(2), v_rel(2), v_norm(2), vv, fdrag(2)
        
        v(1) = y(2)/mass; v(2) = y(4)/mass
        v_rel = v - wind
        vv =  sqrt(v_rel(1)**2 + v_rel(2)**2)
        v_norm = v_rel/vv
        fdrag = -drag * vv**2 * v_norm

        f(1) = v(1)
        f(2) = fdrag(1)
        f(3) = v(2)
        f(4) = -mass*g + fdrag(2)
        
        k= dt*f
        
        end subroutine deriv
    
end subroutine rk4step