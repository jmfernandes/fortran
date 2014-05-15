

module setup

    use numtype
    implicit none
    real(dp) :: offset
    integer, parameter :: n_eq = 2
    real(dp), parameter :: mass=1._dp, hbar=1._dp
    real(dp), parameter :: g = 10._dp, eps =5E-2
    real(dp) :: energy
    integer :: printvar
    real(dp), allocatable, dimension(:,:) :: wf
!     real(dp), parameter :: x0 = 0._dp, y0 = 0._dp
!     real(dp), parameter :: distance = 6._dp, height = 3._dp
!     real(dp), parameter :: v0 = 10._dp
!     real(dp), parameter :: drag = 0._dp, wind = (0._dp, 0._dp)

end module setup

program ball_shot

    use setup
    use chebyshev  !just chebyshev
    implicit none
    real(dp) :: ya, yb
    real(dp), external :: shoot  ! This should be changed back to real
    integer :: nch, i
    
    ya = -1
    yb = 0
    
    nch = 50

    printvar = 0

    open(unit = 3, file='y.data', action= 'write', status ='replace')
    open(unit = 13, file='xy.data', action= 'write',status ='replace')
    
    call chebyex(shoot, nch, cheb, ya, yb)

    call chebyzero(nch, cheb, ya, yb, z0, iz0)
    
    print *, iz0, z0(1:iz0)*180._dp/pi

    printvar = 1

    do i=1,iz0
        offset = shoot(z0(i))
        print *, offset
    end do

!     offset = shoot(1.0_dp)

    
end program ball_shot

function shoot(input) result(fail)

    use setup
    implicit none
    real(dp) :: y(n_eq), x, dx,xmax, xmin, fail, input
    integer :: imax, i
    
    dx =0.01_dp
    xmin = 0._dp
    xmax = 30
    imax = abs(xmax-xmin)/dx

    x = xmax
    i = imax

!     allocate(wf(0:imax,2))
    
    y(1) = 0.00000001_dp
    y(2) = -0.00000001_dp

    energy = input
    
    
    do while(x >= 0)
!         wf(i,1) = x
!         wf(i,2) = y(1)
!         i = i-1
!         if (printvar == 1) then
!             write(3, fmt='(2f20.10)') x, y(3)
!             write(13, fmt='(2f20.10)') y(1), y(3)
!         end if
        if (printvar == 1) then
            write(1,*) x,0._dp,y(1)
        end if
        call rk4step(x,y,-dx)
    
    end do

!     print *, y(3)

!     do i =1,imax
!     write(2,*) wf(i,1), wf(i,2)
!     end do

!     wf(0:imax,2) = wf(0:imax,2)/sqrt(2*y(3)) 

!     do i =1,imax
!     write(2,*) wf(i,1), 0._dp, wf(i,2)
!     end do
    
    fail =  y(1)
    
end function shoot