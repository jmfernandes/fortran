
module setup

    use NumType
    implicit none
    integer, parameter :: n_eq = 3
    real(dp), parameter :: hbar2=1._dp, &
            mass=1.0_dp, xm=mass*mass/(mass+mass)
    integer :: l

    real(dp) :: energy, xmax, dx, eps, xmid
    real(dp), allocatable, dimension(:,:) :: wf
    integer :: imax

end module setup

program bound

    use setup
    use chebyshev
    implicit none
    real(dp) :: eminx,emin,emaxx,emax,deltae, &
        de,e0,psi
    real(dp), external :: psi0
    integer :: nch, izz, i, maxf, n, nstep

    l=0
    
    xmax=5.0_dp
    dx=0.001_dp
    eps=0.0001_dp
    maxf=20
    eminx=-100._dp
    emaxx=0._dp
    de=0.2_dp
    nstep=5
    nch=5

    imax=nint(xmax/dx)+1
    allocate(wf(0:imax,2))

    e0=eminx
    do
        psi=psi0(e0)
        if(psi /= xmax .or. e0 > emaxx) exit
        e0=e0+de
    end do

    eminx=e0
!     print *, eminx
    deltae=(emaxx-eminx)/nstep
!     print *, deltae

    izz = 0

    print *, 'The bound state energy are '
    do n=1,nstep

        emin = eminx+(n-1)*deltae
        emax = eminx+n*deltae
!         print *, emin, emax
        call chebyex(psi0,nch,cheb,emin,emax)
        call chebyzero(nch,cheb,emin,emax,z0,iz0)
!         print *, z0(1:iz0)

        de=0.1_dp
        do i=1,iz0
            e0=z0(i)
            call root_polish(psi0,e0,de,eps,maxf)
            psi=psi0(e0)
            izz=izz+1
            print *, izz,'E=',e0
            call wavef(e0,izz)
        end do

    end do

end program bound

function psi0(eee) result(psi)

    use setup
    implicit none
    real(dp), intent(in) :: eee
    real(dp) :: x, psi, k
    real(dp), dimension(n_eq) :: y

    energy = eee
    k = sqrt(2*xm/hbar2*(-energy))

    x=dx
    y(1)=x**(l+1)
    y(2)=(l+1)*x**l
    y(3)=0._dp
    do while (x <= xmax .and. y(2) > 0._dp)
        call rk4step(x,y,dx)
    end do

    xmid=x
    if ( xmid >= xmax ) then
        psi=xmax
        return
    end if

    x=xmax
    y(1) = exp(-k*x)
    y(2) = -k*y(1)
    y(3) = 0._dp
    do while (x > xmid)
        call rk4step(x,y,-dx)
    end do
    psi=y(2)

!     print *, e,psi

end function psi0

subroutine wavef(eee,iz)

    use setup
    implicit none
    real(dp), intent(in) :: eee
    real(dp) :: x, psi,k,y12,y32,y11,y31,yy
    real(dp), dimension(n_eq) :: y
    integer :: iz, n, i
    
    energy = eee
    k=sqrt(2*xm/hbar2*(-energy))

    x=xmax
    y(1)=exp(-k*x)
    y(2)=-k*y(1)
    y(3) = 0._dp
    do while (x>xmid)
        n=nint(x/dx)
        wf(n,1)=x
        wf(n,2)=y(1)
        call rk4step(x,y,-dx)
    end do
    y12=y(1)
    y32=-y(3)

    wf(0,1)=0._dp
    wf(0,2)=0._dp
    x=dx
    y(1)=x**(l+1)
    y(2)=(l+1)*x**l
    y(3)=0._dp
    do while ( x <= xmid )
        n=nint(x/dx)
        wf(n,1)=x
        wf(n,2)=y(1)
        call rk4step(x,y,dx)
    end do
    y11=y(1)
    y31=y(3)

    wf(0:n,2)=y12/y11*wf(0:n,2)
    y31=(y12/y11)**2* y31
    yy=y31+y32
    n=nint(xmax/dx)
    wf(0:n,2)=wf(0:n,2)/sqrt(yy)
    do i = 0,n
        write(30*(l+1)+iz,*) wf(i,1),-wf(i,2)
    end do

end subroutine wavef

