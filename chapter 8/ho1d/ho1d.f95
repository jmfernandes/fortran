
module setup

    use numtype
    implicit none

    integer, parameter :: n_eq = 2
    real(dp), parameter :: hbar = 1._dp, mass=1._dp, omega = 1._dp
    
    real(dp) :: energy

end module setup

program ho1d

    use setup
    implicit none

    real(dp) :: x, dx, xmax, y(n_eq)

    energy = 4.5
    xmax   = 5._dp
    dx     = 0.01_dp

    x = 0._dp

    y(1) = 1
    y(2) = 0

    do while (x < xmax)

        write(1,*) x,0._dp,y(1)
        call rk4step(x,y,dx)

    end do

end program ho1d