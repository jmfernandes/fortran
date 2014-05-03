
module setup

    use numtype
    implicit none

    integer, parameter :: n_eq = 3
    real(dp), parameter :: hbar = 1._dp, mass=1._dp, omega = 1._dp
    real(dp), allocatable, dimension(:,:) :: wf

    real(dp) :: energy

end module setup

program ho1d

    use setup
    implicit none

    real(dp) :: x, dx, xmax, xmin, y(n_eq)
    integer  :: imax, i

    energy = 2.5_dp
    xmin   = -5._dp
    xmax   = 5._dp
    dx     = 0.01_dp

    x    = xmax
    imax = abs(xmax-xmin)/dx

    allocate(wf(0:imax,2))

    y(1) = 0.000001_dp
    y(2) = -0.000002_dp
    y(3) = 0._dp

    i = imax

    do while (x >= xmin)

        wf(i,1) = x
        wf(i,2) = y(1)
        i = i-1
        call rk4step(x,y,-dx)

    end do

    print *, y(3)

    do i =1,imax
    write(2,*) wf(i,1), wf(i,2)
    end do

    wf(0:imax,2) = wf(0:imax,2)/sqrt(2*y(3)) 

    do i =1,imax
    write(1,*) wf(i,1), 0._dp, wf(i,2)
    end do

end program ho1d