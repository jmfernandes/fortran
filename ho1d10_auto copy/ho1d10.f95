
module setup

    use NumType
    implicit none
    integer, parameter :: n_eq = 3
    real(dp), parameter :: hbar = 1._dp, hbar2 = hbar**2,&
                            mass = 1._dp, omega = 1._dp, &
                            x0 = sqrt(hbar/(mass*omega))
    real(dp) :: energy, xmax, dx, eps, total,val1,val2
    real(dp), allocatable, dimension(:,:) :: wf, w1
    integer :: imax

end module setup

program ho1d

    use setup
    use chebyshev
    implicit none
    real(dp) :: eminx, emaxx, emin, emax, de, e0, de0
    real(dp), external :: psi0
    integer :: nch, iz, i, maxf, izz,m,n,j,k

    eminx = -1._dp
    emaxx = 0._dp

    xmax = 40.0_dp
    dx = 0.01_dp
    eps = 0.000001_dp
    de0 = 0.01_dp
    maxf = 10
    imax = abs(xmax/dx)
    allocate(wf(0:imax,2))
    allocate(w1(0:2*imax,10))
    de = 0.01_dp
    izz = 0
    nch = 5
    emin = eminx
    emax = emin+de

    print *, '=================Eigenvalues================='

    do while (emin < emaxx)

        call chebyex(psi0, nch, cheb, emin, emax)
        call chebyzero(nch,cheb,emin,emax,z0,iz0)

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

    print *, '============Orthogonality Test============'

    total = 0
    val1=0
    val2=0
    j=1
    k=1

    do j = 1,10
    do k = 1,10
    do i = 1,2*imax
        m=j
        n=k
        val1=w1(i,m)*w1(i,n)
        val2=w1(i+1,m)*w1(i+1,n)
        !Rombint/trapazoidal method of integration
        total = total + (0.01*(val1+1/2._dp*(val2-val1)) )
    end do

    print '(a,1f5.1,a5,i2,a5,i2)','total=', total,'i=',j,'j=',k

    total = 0

    end do
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
    y(1) = 0.00001
    y(2) = -0.00001
    y(3) = 0._dp
    do while ( x > 0._dp )
        call rk4step(x,y,-dx)
    end do 
    psi = y(1)*y(2)

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
    y(1) = 0.00001
    y(2) = -0.00001
    y(3) = 0._dp
    do while ( x > 0._dp)
        call rk4step(x,y,-dx)
    end do
   
    x = xmax
    i = imax + 1
    y(1) = 0.00001
    y(2) = -0.00001
    y(3) = 0._dp
    do while ( x > 0._dp)
        i = i-1
        wf(i,1) = x
        wf(i,2) = y(1)
        call rk4step(x,y,-dx)
    end do

    !establish parity so the negative side of graph cam be plotted
    imin=i
    if( abs(y(1)) > abs(y(2)) ) then
        parity = 1
    else 
        parity = -1
    end if 

    ! make a master list of wavefunctions called w1. Used for integration
    w1(0:imax,iz) = parity*wf(0:imax,2)/sqrt(2*y(3)) 
    w1(imax:2*imax,iz) = wf(0:imax,2)/sqrt(2*y(3)) 

    !normalize the wavefunction
    wf(0:imax,2) = wf(0:imax,2)/sqrt(2*y(3)) 
    

    !make the graphs for positive and negative sides
    do i = imax, imin, -1
        write(unit=20+iz,fmt='(2f15.5)') wf(i,1), &
            wf(i,2)
    end do

    do i = imin, imax
        write(unit=20+iz,fmt='(2f15.5)') -wf(i,1), &
            parity*wf(i,2)
    end do


end subroutine wavef

