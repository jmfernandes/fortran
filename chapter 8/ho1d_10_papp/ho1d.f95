
module setup

    use NumType
    implicit none
    integer, parameter :: n_eq = 3
    real(dp), parameter :: hbar = 1._dp, mass = 1._dp, omega = 1._dp
    
    real(dp), allocatable, dimension(:,:) :: wf
    
    real(dp) :: energy
    integer :: parity

end module setup

program ho1d

    use setup
    implicit none
    real(dp) :: x, dx, xmax, y(n_eq)
    integer :: imax, i
    
    energy = 1.5
    xmax = 5._dp
    dx = 0.01_dp
    imax= abs(xmax/dx)
    
    allocate(wf(-imax:imax,2))
    
    x = xmax
    
    y(1) = 0.00000001
    y(2) = -0.0002
    y(3) = 0._dp
    
    
    i = imax
    do while ( x >= 0 )
    
        wf(i,1) = x
        wf(i,2) = y(1)
                
        i = i-1
        call rk4step(x,y,-dx)
        
    end do
    
    !print *,y(3),y(1),y(2)
    
    wf(0:imax,2) = wf(0:imax,2)/ sqrt(2*y(3))
    if ( abs(y(1)) > abs(y(2)) ) then
        parity = 1
    else if ( abs(y(2)) > abs(y(1)) ) then
        parity = -1
    end if
    print *,' E = ', energy,' P = ', parity
    
    do i = -imax, 0
        wf(i,1) = -wf(-i,1)
        wf(i,2) = parity * wf(-i,2)
    end do
    
    do i = -imax, imax
        write(1,*) wf(i,1),wf(i,2)
    end do

end program ho1d


