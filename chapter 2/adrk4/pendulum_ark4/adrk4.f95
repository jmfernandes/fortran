!decides if the step is too big or small and adjusts

subroutine ark4(t,y,tau)

    use setup, only : dp, number_of_equations
    implicit none
    real(dp), intent(inout) :: t, tau
    real(dp), dimension(number_of_equations) :: y, y1, y2
    real(dp) :: tsave, halftau, erat, tau_old
    real(dp), parameter :: safe1 = 0.9_dp, safe2 = 4._dp, &
        eps = 1.e-16_dp, er = 1.e-10_dp
    integer :: maxtry=100, itry
   
    tsave = t
    do itry = 1, maxtry
        halftau = tau/2
        y2 = y
        t = tsave
        call rk4step(t,y2,halftau)
        call rk4step(t,y2,halftau)
        y1 = y
        t = tsave
        call rk4step(t,y1,tau)
        erat = maxval(abs(y1-y2)/(er*(abs(y1)+abs(y2))/2+eps))
        tau_old = tau
        tau = safe1*tau_old*erat**(-0.2_dp)
        tau = max(tau,tau_old/safe2)
        tau = min(tau,safe2*tau_old)
        if ( erat > 1 ) cycle
        y = y2
        exit
    end do
    
end subroutine ark4

