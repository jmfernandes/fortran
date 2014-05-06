

module setup

    use numtype
    implicit none
    real(dp) :: offset
    integer, parameter :: n_eq = 2
    real(dp), parameter :: g = 10._dp, mass = 1._dp, eps =5E-1
    real(dp), parameter :: x0 = 0._dp, y0 = 1.98_dp
    real(dp), parameter :: distance = 28._dp, height = 3._dp
    real(dp), parameter :: v0 =555._dp
    real(dp), parameter :: drag = 0.4_dp, wind = (3._dp, 0._dp)
    integer :: printvar

end module setup

program ball_shot

    use setup
    use chebyshev  !just chebyshev
    implicit none
    real(dp) :: ya, yb
    real(dp), external :: shoot  ! This should be changed back to real
    integer :: nch, i
    
    ya = 1
    yb = 10
    
    nch = 50

    offset = 0

    printvar = 0

    open(unit = 3, file='y.data', action= 'write', status ='replace')
    open(unit = 13, file='xy.data', action= 'write',status ='replace')
    
    call chebyex(shoot, nch, cheb, ya, yb)

    call chebyzero(nch, cheb, ya, yb, z0, iz0)
    
    print *, iz0, z0(1:iz0)*180._dp/pi


    printvar = 1

    do i=1,iz0
        offset = shoot(z0(i))
        print *, offset, z0(i)
    end do

    
end program ball_shot

function shoot(alpha) result(fail)

    use setup
    implicit none
    real(dp) :: y(n_eq), t, dt, fail, alpha, tmax
    
!     print *, "banans"
    t = 0._dp
    dt =0.01_dp
    tmax = 10._dp

    
    y(1) = 1
    y(2) = 1

    
    do while(t<= tmax)
!         print *, y(1)
!         if (printvar == 1) then
!             write(3, fmt='(2f20.10)') t, y(3)
!             write(13, fmt='(2f20.10)') y(1), y(3)
!         end if
        call rk4step(t,y,dt)
    
    end do
    
    fail = y(1)
    
end function shoot