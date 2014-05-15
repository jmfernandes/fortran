
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
    allocate(w1(0:imax,10))
    de = 0.01_dp
    izz = 0
    nch = 5
    emin = eminx
    emax = emin+de

    print *, '=================Eigenvalues================='

    do while (emin < emaxx)

        call chebyex(psi0, nch, cheb, emin, emax)
        call chebyzero(nch,cheb,emin,emax,z0,iz0)
!         print *, z0(1:iz0)

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

!     do i =1,imax
!         write(99,*) wf(i,1), w1(i,1)*w1(i,2)
!     end do

    total = 0
    val1=0
    val2=0
    j=1
    k=1


    print *, '============Orthogonality Test============'

    do j = 1,10
    do k = 1,10
    do i = 1,imax
        m=j
        n=k
        val1=w1(i,m)*w1(i,n)
        val2=w1(i+1,m)*w1(i+1,n)
        total = total + (0.01*(val1+1/2._dp*(val2-val1)) )
    end do

    print '(a,i3,a5,i2,a5,i2)','total=', nint(total+.01),'i=',j,'j=',k

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
!     print *, eee,psi

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
!         write(unit=20+iz, fmt='(2f15.5)') x,y(1)
        call rk4step(x,y,-dx)
    end do


    imin=i
    if( abs(y(1)) > abs(y(2)) ) then
        parity = 1
    else 
        parity = -1
    end if 

!     do i =1,imax
!     write(2,*) wf(i,1), wf(i,2)
!     end do

!     print *, y(3), 'babnans'
!     if (iz == 1) then
!         w1(0:imax,1) = wf(0:imax,2)/sqrt(2*y(3)) 
!     else if (iz == 2) then
!         w1(0:imax,2) = wf(0:imax,2)/sqrt(2*y(3)) 
!     end if

!     print *, iz , 'iz'

    w1(0:imax,iz) = wf(0:imax,2)/sqrt(2*y(3)) 

    wf(0:imax,2) = wf(0:imax,2)/sqrt(2*y(3)) 
    

!     do i =1,imax
!     write(99,*) w1(i,1), w1(i,2)
!     end do

!     do i =1,imax
!     write(1,*) wf(i,1), wf(i,2)
!     end do

!     print *, parity, 'parity is', iz
    do i = imax, imin, -1
        write(unit=20+iz,fmt='(2f15.5)') wf(i,1), &
            wf(i,2)
    end do

    do i = imin, imax
        write(unit=20+iz,fmt='(2f15.5)') -wf(i,1), &
            parity*wf(i,2)
    end do


end subroutine wavef

