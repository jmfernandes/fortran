
module setup

    use NumType
    implicit none
    integer, parameter :: ndim = 50, &
        npint=200,nlag=100,ncleg=50,n_eq=3
    real(dp), parameter :: hbar2=1, &
        mass = 1.0_dp, &
        xm = mass*mass/(mass+mass)/hbar2
    integer :: l, nm, det0

    real(dp) :: b,bp,z12,v(ndim,ndim)
    real(dp) :: xcleg(npint),wcleg(npint)
    real(dp) :: eeenergy

end module setup

program bound

    use setup
    use chebyshev
    implicit none
    real(dp) :: eminx,emin,emaxx,emax,eps, &
        de,deltae,e0,energy(10),dr,rmax,x,xmax,dx
    real(dp), external ::cfd
    integer :: nch,izz,i,maxf,n,nstep,ifail
    real(dp), dimension(ndim,ndim) :: vp,w,ww
    real(dp) :: y(n_eq)

    call d01bcf(0,0._dp,2*pi,0._dp,0._dp, &
        ncleg,wcleg,xcleg,ifail)
    if(ifail /= 0) stop ' d01bcf '

    l=0
    z12=0
    b=2.5_dp
    bp=2._dp
    nm=20

    eps=1.e-8_dp
    maxf=10
    eminx=-35._dp
    emaxx=-0.001_dp
    nstep=50
    nch=5
    dr=0.1_dp
    rmax=10._dp
    deltae=(emaxx-eminx)/nstep
    do nm = 1,48

        call csme(b,bp,l,nm,vp,w)
        call dgei(nm+1,w)
        ww(1:nm+1,1:nm+1)=matmul(w(1:nm+1,1:nm+1), &
            vp(1:nm+1,1:nm+1))
        v(1:nm+1,1:nm+1)=matmul(ww(1:nm+1,1:nm+1), &
            transpose(w(1:nm+1,1:nm+1)))

        izz=0
        do n=1,nstep

            emin = eminx+(n-1)*deltae
            emax = eminx+n*deltae
!             print *, ' [,emin,emax,]'
            det0=0
            call chebyex(cfd,nch,cheb,emin,emax)
            call chebyzero(nch,cheb,emin,emax,z0,iz0)
!             print *, '-------------->>>', z0(1:iz0)

            de=0.01_dp
            do i=1,iz0
                e0=z0(i)
                if (e0 < emin .or. e0 > emax) cycle
                call root_polish(cfd,e0,de,eps,maxf)
                izz=izz+1
                energy(izz)=e0
                write(10*(l+1)+izz,*) nm,e0 
!                 print *, '                  E   ', izz,e0
            end do

        end do

    end do

    do i=1,izz
        e0=energy(i)
        call wavef(i,e0,dr,rmax)
    end do

    x=0._dp
    dx=0.01_dp
    xmax = 5._dp

    y(1) = 0
    y(2) = 1
    y(3) = 0

    eeenergy = 0.639

    do while (x < xmax)

        write(1,*) x, y(1)
        call rk4step(x,y,dx)
    end do

end program bound

function cfd(eee) result(fdet) !fredholm determinant
    use setup
    implicit none
    real(dp), intent(in) :: eee
    complex(dp) :: z,gz(ndim,ndim)
    real(dp) :: gr(ndim,ndim), det(2), fdet

    z=eee
    call cgri(z,b,xm,z12,l,nm,gz)
    gr(1:nm+1,1:nm+1)=gz(1:nm+1,1:nm+1)-v(1:nm+1,1:nm+1)

    call dsydet(nm+1,gr,det)

    if( det0 == 0 ) then
        fdet = det(1); det0 = det(2)
    else
        fdet = det(1)*10._dp**int(det(2)-det0)
    endif
!     print '(2x,2d20.10)',eee,fdet

end function cfd

subroutine wavef(iz,eee,dr,rmax)

    use setup
    implicit none
    real(dp) :: eee,dr,rmax
    complex(dp), dimension(ndim,ndim) :: g
    real(dp) :: gp(ndim,ndim), phi(ndim), st(ndim)
    real(dp) :: t,r,x,r0,wf
    complex(dp) :: z,zz
    integer :: iz,i,imax

    print *, iz,eee
    gp(1:nm+1,1:nm+1) = 0._dp
    r0 =0.1_dp
    do i=1,ncleg
        t = xcleg(i)
        z = r0*exp(iic*t)
        zz = eee+z
        call green(zz,g)
        gp(1:nm+1,1:nm+1) = gp(1:nm+1,1:nm+1) + &
            wcleg(i)*z/(2*pi)*g(1:nm+1,1:nm+1)
    end do
    phi(1:nm+1)=gp(1:nm+1,1)/sqrt(gp(1,1))

    imax=nint(rmax/dr)
    do i =0,imax
        r=i*dr
        x=2*b*r
        call csturmx(x,l,nm,st)
        wf=exp(-b*r)*sum(phi(1:nm+1)*st(1:nm+1))
        write(100+10*(l+1)+iz,*) r,wf
    end do

end subroutine wavef

subroutine green(z,g)

    use setup
    implicit none
    complex(dp) :: z, g(ndim,ndim)

    call cgri(z,b,xm,z12,l,nm,g)
    g(1:nm+1,1:nm+1)=g(1:nm+1,1:nm+1)-v(1:nm+1,1:nm+1)
    call zgei(nm+1,g)

end subroutine green

subroutine csme(b,bp,l,nm,v,w)

    use setup, only : dp, ndim, npint, nlag
    implicit none
    real(dp), intent(in) :: b,bp
    integer, intent(in) :: l,nm
    real(dp), dimension(ndim,ndim) :: v,w
    real(dp), dimension(npint) :: xlag,wlag
    real(dp), dimension(ndim,npint) :: st, stp
    real(dp), dimension(npint) :: pf
    integer :: info,i,j
    real(dp) :: x

    info=0
    call d01bcf(3,0._dp,1._dp,0._dp,0._dp, &
            nlag,wlag,xlag,info)

    do i =1,nlag
        x = xlag(i)
        call csturmx(x,l,nm,st(1,i))
        pf(i) = potential(x/(2*bp))/(2*bp)*wlag(i)
    end do
    forall (i=1:nm+1,j=1:nm+1) &
        v(i,j)=sum(st(i,1:nlag)*st(j,1:nlag)*pf(1:nlag))

    do i = 1,nlag
        x = xlag(i) 
        call csturmx(2*bp/(b+bp)*x,l,nm,st(1,i))
        call csturmx(2*b/(b+bp)*x,l,nm,stp(1,i))
        pf(i)=wlag(i)/x
    end do
    forall (i=1:nm+1,j=1:nm+1) &
        w(i,j)=sum(st(i,1:nlag)*stp(j,1:nlag)*pf(1:nlag)) 

    contains

        function potential(r) result(v) !potential

            use setup, only : dp
            implicit none
            real(dp), intent(in) :: r
            real(dp) :: v

            v = -1/sqrt(1+r**2)

        end function potential

end subroutine csme

subroutine csturmx(x,l,nm,st)

    use setup, only : dp, ndim
    implicit none
    integer, intent(in) :: l, nm
    real(dp), intent(in) :: x
    real(dp), dimension(ndim) :: st
    integer :: i
    real(dp) :: st0, df

    st0 = 0._dp
    st(1) = x**(l+1)/sqrt(gamma(2*l+2._dp))
    do i = 1,nm
        df = sqrt(dble(i*(i+2*l+1)))
        st(i+1) = ((2*i+2*l-x)*st(i)-st0)/df
        st0 = df*st(i)
    end do 

end subroutine csturmx

subroutine dgei(n,matr)

    use setup , only : dp, ndim
    implicit none
    integer, intent(in) :: n
    real(dp) :: matr(ndim,ndim)
    integer, parameter :: lwork = 4*ndim
    real(dp), dimension(lwork) :: work
    integer :: info, ipvt(ndim), k

    info = 0
    call dgetrf(n,n,matr,ndim, ipvt,info)
    if(info /= 0 ) stop 'dgetrf'
    call dgetri(n,matr,ndim,ipvt,work,lwork,info) 
    if(info /= 0) stop 'dgetri'


end subroutine dgei

subroutine dsydet(n,matr,det)
        
            use setup , only : dp , ndim
            implicit none
            integer :: n
            real(dp) :: matr(ndim,ndim), det(2)
            integer, parameter :: lwork = 4*ndim
            real(dp), dimension(lwork) :: work
            integer :: info, ipvt(ndim), k
            real(dp) :: t, d
        
            info = 0
            call dsytrf('u',n,matr,ndim,ipvt,work,lwork,info)
        
            det(1) = 1.0_dp
            det(2) = 0.0_dp
            t = 0.0_dp
            do k = 1, n
                d = matr(k,k)
                if ( ipvt(k) <= 0 ) then
                    if ( t == 0.0_dp ) then
                        t = matr(k,k+1)
                        d = (d/t)*matr(k+1,k+1) - t
                    else
                        d = t
                        t = 0.0_dp
                    end if
                end if
                det(1) = d*det(1)
                if ( det(1) /= 0.0_dp ) then
                    do while ( abs(det(1)) <= 1 )
                        det(1) = 10*det(1)
                        det(2) = det(2) - 1
                    end do
                    do while ( abs(det(1)) >= 10 )
                        det(1) = det(1)/10
                        det(2) = det(2) + 1
                    end do
                end if
            end do
            
end subroutine dsydet


subroutine zgei(n,matr)

            use setup, only : dp,ndim
            implicit none
            integer, intent(in) :: n
            complex(dp) :: matr(ndim,ndim)
            integer, parameter :: lwork = 4*ndim
            complex(dp), dimension(lwork) :: work
            integer :: info, ipvt(ndim)
            
            info=0
            call zgetrf(n,n,matr,ndim,ipvt,info)
                if(info /= 0) stop ' dgetrf '
            call zgetri(n,matr,ndim,ipvt,work,lwork,info)
                if(info /= 0) stop ' dgetri '
        
end subroutine zgei

