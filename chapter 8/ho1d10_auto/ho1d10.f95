
module setup

    use NumType
    implicit none
    integer, parameter :: n_eq = 3
    real(dp), parameter :: hbar = 1._dp, hbar2 = hbar**2,&
                            mass = 1._dp, omega = 1._dp, &
                            x0 = sqrt(hbar/(mass*omega))
    real(dp) :: energy, xmax, dx, eps
    real(dp), allocatable, dimension(:,:) :: wf
    integer :: imax

end module setup

program ho1d

    use setup
    use chebyshev
    implicit none
    real(dp) :: eminx, emaxx, emin, emax, de, e0, de0
    real(dp), external :: psi0
    integer :: nch, iz, i, maxf, izz

    eminx = 0._dp
    emaxx = 10._dp

    xmax = 6.0_dp
    dx = 0.001_dp
    eps = 0.000001_dp
    de0 = 0.01_dp
    maxf = 10
    imax = abs(xmax/dx)
    allocate(wf(0:imax,2))
    de = 1._dp
    izz = 0
    nch = 5
    emin = eminx
    emax = emin+de

    do while (emin < emaxx)

        call chebyex(psi0, nch, cheb, emin, emax)
        call chebyzero(nch,cheb,emin,emax,z0,iz0)
        print *, z0(1:iz0)

        do iz=1,iz0
            e0=z0(iz)
            call root_polish(psi0,e0,de0,eps,maxf)
            izz = izz + 1
            print *, izz, e0
            call wavef(e0,izz)
        end do
        
        emin = emax
        emax = emin + de

    end do

end program ho1d

function psi0(eee) result(psi)

    use setup
    implicit none
    real(dp), intent(in) :: eee
    real(dp) :: x, psi
    real(dp), dimension(n_eq) :: y

    energy = eee
    x = xmax
    y(1) = exp(-0.5_dp*(x/x0)**2)
    y(2) = -(x/x0**2)*y(1)
    y(3) = 0._dp
    do while ( x > 0._dp )
        call rk4step(x,y,-dx)
    end do 
    psi = y(1)*y(2)
    print *, eee,psi

end function psi0

subroutine wavef(eee,iz)

    use setup
    implicit none
    real(dp), intent(in) :: eee
    real(dp) :: x, parity
    real(dp), dimension(n_eq) :: y
    integer :: iz, i, imin
    energy = eee
    x = xmax
    y(1) = exp(-0.5_dp*(x/x0)**2)
    y(2) = -(x/x0**2)*y(1)
    y(3) = 0._dp
    do while ( x > 0._dp)
        call rk4step(x,y,-dx)
    end do
   
    x = xmax
    i = imax + 1
    y(1) = exp(-0.5_dp*(x/x0)**2) / sqrt(2*y(3))
    y(2) = -(x/x0**2)*y(1)
    y(3) = 0._dp
    do while ( x > 0._dp)
        i = i-1
        wf(i,1) = x
        wf(i,2) = y(1)
        write(unit=20+iz, fmt='(2f15.5)') x,y(1)
        call rk4step(x,y,-dx)
    end do
    imin=i
    if( abs(y(1)) > abs(y(2)) ) then
        parity = 1
    else 
        parity = -1
    end if 

    do i = imin, imax
        write(unit=20+iz,fmt='(2f15.5)') - wf(i,1), &
            parity*wf(i,2)

    end do

end subroutine wavef

