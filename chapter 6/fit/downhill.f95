
subroutine downhill(n,func,xstart,fstart, &
                        stepi,epsf,itmin,iter)
!
!   n           dimension of the problem
!   func        function
!   xstart      starting values
!   fstart      conrespoding function value
!   stepi       relative stepsize for initial simplex
!   epsf        esilon for termination
!   itmin       termination istested if itmin < it
!   iter        maximum number of iterations
!

    use NumType
    implicit none
    integer :: n, iter, itmin
    real(dp), external :: func
    real(dp) :: xstart(1:n), fstart, stepi, epsf
    real(dp), parameter :: alph=1._dp, gamm=2._dp, &
                            rho=0.5_dp, sig=0.5_dp
    real(dp) :: xi(1:n,1:n+1), x(1:n,1:n+1), &
        fi(1:n+1), f(1:n+1),  &
        x0(1:n), xr(1:n), xe(1:n), xc(1:n), &
        fxr, fxe, fxc, deltaf
    integer :: i, j, ii, it
    
    xi(1:n,1) = xstart(1:n);    fi(1) = fstart
    do i = 2, n+1
        xi(1:n,i)=xi(1:n,1)
        xi(i-1,i)=xi(i-1,i)*(1+stepi)
        fi(i)=func(xi(1:n,i))
    end do
       
    do it = 1, iter 
        
        do i = 1, n+1                           ! ordering
            ii = minloc(fi(1:n+1),dim=1)
            x(1:n,i) = xi(1:n,ii);  f(i) = fi(ii)
            fi(ii) = huge(0._dp)
        end do
        xi(1:n,1:n+1) = x(1:n,1:n+1)
        fi(1:n+1) = f(1:n+1)
                
        x0(1:n) = sum(x(1:n,1:n),dim=2)/n   ! central
        
        if ( itmin < it ) then              ! condition for exit
            deltaf = (f(n)-f(1))
            write(7,*) it,deltaf
            if(deltaf < epsf ) exit
        end if
                
        xr(1:n) = x0(1:n)+alph*(x0(1:n)-x(1:n,n+1))
        fxr = func(xr)
        if( fxr < f(n) .and. &              ! reflection
                f(1) <= fxr ) then 
            xi(1:n,n+1) = xr(1:n);  fi(n+1) = fxr
            cycle
        
        else if ( fxr < f(1) ) then         ! expansion
            xe(1:n) = x0(1:n)+gamm*(x0(1:n)-x(1:n,n+1))
            fxe = func(xe)
            if( fxe < fxr ) then
                xi(1:n,n+1) = xe(1:n);  fi(n+1) = fxe
                cycle
            else
                xi(1:n,n+1) = xr(1:n);  fi(n+1) = fxr
                cycle
            end if

        else if ( fxr >= f(n) ) then        ! contraction
            xc(1:n) = x(1:n,n+1)+rho*(x0(1:n)-x(1:n,n+1))
            fxc = func(xc)
            if( fxc <= f(n+1) ) then
                xi(1:n,n+1) = xc(1:n);  fi(n+1) = fxc
                cycle
            end if
            
        else                                 ! reduction
            do i = 2, n+1
                xi(1:n,i) = x(1:n,1)+sig*(x(1:n,i)-x(1:n,1))
                fi(i) = func(xi)
            end do
            cycle
            
        end if
        
    end do
    
    xstart(1:n)=xi(1:n,1); fstart = fi(1)

end subroutine downhill


