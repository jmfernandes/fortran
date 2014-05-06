
module setup

    use NumType
    implicit none
    integer, parameter :: n_eq = 3
    real(dp), parameter :: hbar2=41.47_dp, &
            mass=4.0_dp, xm=mass*mass/(mass+mass)
    integer :: l

    real(dp) :: energy, xmax, dx, eps, xmid
    real(dp), allocatable, dimension(:,:) :: wf
    integer :: imax

end module setup

program ho1d

    use setup
    use chebyshev
    implicit none
    real(dp) :: dmin,dmax,de,delta0,psi
    real(dp), external :: psi0
    integer :: nch, iz, i, maxf

    l=0
    energy=1._dp

    xmax=20._dp
    dx=0.001_dp
    eps=0.0000001_dp
    maxf=20
    dmin=0._dp
    dmax=pi

    imax=nint(xmax/dx)+1
    allocate(wf(0:imax,2))

    nch=6

    call chebyex(psi0,nch,cheb,dmin,dmax)
    call chebyzero(nch,cheb,dmin,dmax,z0,iz0)
    print *, z0(i:iz0)

    de=0.01_dp

    do iz=1,iz0
        delta0=z0(iz)
        call root_polish(psi0,delta0,de,eps,maxf)
        psi=psi0(delta0)
        print *, ' Delta=', delta0
        call wavefunction(delta0)
    end do

end program ho1d

function psi0(delta) result(psi)

    use setup
    implicit none
    real(dp), intent(in) :: delta
    real(dp) :: x, psi, k
    real(dp), dimension(n_eq) :: y

    k = sqrt(2*xm/hbar2*energy)
    
    x=dx
    y(1)=x**(l+1)
    y(2)=(l+1)*x**l

    do while (x <= xmax .and. y(2) > 0._dp)
        call rk4step(x,y,dx)
    end do
    xmid=x
    if( xmid >= xmax) then
        stop 'xmax is too small'
    end if

    x=xmax
    y(1) = sin(k*x - l*pi/2 +delta)
    y(2) = k*cos(k*x -l*pi/2 +delta)
    do while ( x>xmid)
        call rk4step(x,y,-dx)
    end do
    psi=y(2)
    print *, delta,psi

end function psi0

subroutine wavefunction(delta)

    use setup
    implicit none
    real(dp), intent(in) :: delta
    real(dp) :: x, k, y12,y11
    real(dp), dimension(n_eq) :: y
    integer :: n, i
    
    k=sqrt(2*xm/hbar2*energy)

    x=xmax
    y(1)=sin(k*x-l*pi/2+delta)
    y(2)=k*cos(k*x-l*pi/2+delta)
    do while (x>xmid)
        n=nint(x/dx)
        wf(n,1)=x
        wf(n,2)=y(1)
        call rk4step(x,y,-dx)
    end do
    y12=y(1)

    wf(0,1)=0._dp
    wf(0,2)=0._dp
    x=dx
    y(1)=x**(l+1)
    y(2)=(l+1)*x**l
    do while (x<=xmid)
        n=nint(x/dx)
        wf(n,1)=x
        wf(n,2)=y(1)
        call rk4step(x,y,dx)
    end do
    y11=y(1)

    wf(0:n,2)=y12/y11*wf(0:n,2)
    n=nint(xmax/dx)
    do i = 0,n
        write(10,*) wf(i,1),wf(i,2)
    end do

end subroutine wavefunction

