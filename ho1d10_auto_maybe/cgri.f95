
subroutine cgri(eee,b,xm,z12,l,nm,g)

    use NumType
    implicit none
    integer, parameter :: ndim=50
    integer , intent(in) :: l, nm
    complex(dp), intent(in) :: eee
    complex(dp), dimension(ndim,ndim), intent(out) :: g 
    real(dp), intent(in) :: b, xm, z12
    integer :: n
    complex(dp) :: zk, zki, gammai, w0, w1, znm

    zk = sqrt(2*xm*eee)
    if(dble(eee) < 0._dp .AND. dimag(zk) < 0._dp) zk = -zk
        w0 = (zk*zk-b*b)/(2*xm*b) 
        w1 = -(zk*zk+b*b)/(4*xm*b) 
        g(1:nm+1,1:nm+1) = 0._dp 
        do n=0,nm-1
            g(n+1,n+1) = w0*(n+l+1) - z12
            g(n+1,n+2) = w1*sqrt((n+1)*(n+2*l+2._dp)) 
            g(n+2,n+1) = w1*sqrt((n+1)*(n+2*l+2._dp))
        end do
    g(nm+1,nm+1) = w0*(nm+l+1) - z12



    znm = nm +1
    zki = (0.0_dp,1.0_dp)*zk
    gammai = (0.0_dp,1.0_dp)* z12*xm/zk 
    g(nm+1,nm+1) = g(nm+1,nm+1) &
            + w1*w1*znm*(znm+2*l+1) &
            * 4*xm*b/((znm+l+1+gammai)*(b-zki)**2) &
            * f21(-l+gammai,znm,znm+l+1+gammai, &
                ((b+zki)/(b-zki))**2)
    contains

        function f21(a,b,c,z) result(f1)

            implicit none
            complex(dp), intent(in) :: a,b,c,z
            complex(dp) :: f0,f1,c0,c1,d0,d1,de,r
            real(dp), parameter :: eps =1.e-14_dp, &
                tiny=1.e-50_dp
            integer, parameter :: nn = 800
            integer :: n

            if ( abs(a) == 0._dp .AND. &
                    abs(b) == 0._dp ) then
                f1 = (1._dp,0._dp)
            else
                f0 =1
                c0 =f0
                d0 = 0

                do n =0,nn
                    r= -(a+n)*(c-b+n)/ & 
                        ((c+2*n)*(c+2*n+1))*z
                    d1 = 1 + r*d0 +tiny 
                    c1 = 1 + r/c0 +tiny 
                    d1 = 1/d1
                    de = c1*d1
                    f1 = f0*de 
                    c0 = c1
                    d0 = d1
                    f0 = f1
                    r = -(b+n+1)*(c-a+n+1)/ & 
                        ((c+2*n+1)*(c+2*n+2))*z
                    d1 = 1 + r*d0 +tiny 
                    c1 = 1 + r/c0 +tiny 
                    d1 = 1/d1
                    de = c1*d1
                    f1 = f0*de 
                    c0 = c1
                    d0 = d1
                    f0 = f1

                    if ( abs(1-de) < eps) exit

                end do
                f1=1/f1
                if (n >= nn) &
                    print *,abs(1-de)
                end if
            end function f21



end subroutine cgri