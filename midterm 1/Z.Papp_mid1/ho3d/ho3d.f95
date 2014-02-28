
program ho_test

    use NumType
    implicit none
    real(dp) :: x, x_min, x_max, dx
    integer :: imax, i, l
    
    x_min = 0._dp
    x_max = 12._dp
    imax = 300
    
    l = 0
    
    !print *, psi_ho3d(3,5,3.1415_dp)
    
    !stop
    
    
    dx = (x_max-x_min)/imax
    
    do i = 0, imax
    
        x = x_min + i * dx
!        print *, x , psi_ho3d(0,l,x), psi_ho3d(1,l,x), psi_ho3d(2,l,x), &
!            psi_ho3d(3,l,x), psi_ho3d(4,l,x)
        print *, x , psi_ho3d(20,l,x), psi_ho3d(21,l,x)
        
    end do
        
    contains
    
        recursive function psi_ho3d(n,l,r) result(psi)
        
            implicit none
            real(dp) :: r, psi
            integer :: n, l  
            real(dp) :: alpha 
            
            alpha = l + 1._dp/2
        
            if ( n < 0 ) then
                psi =  0._dp
            else if ( n == 0 ) then
                psi = sqrt(2._dp/lp12f(l)) * exp(-r**2/2) * r**(l+1)
            else 
                psi = -(2*n-1+alpha -r**2)/sqrt(n*(n+alpha)) *  psi_ho3d(n-1,l,r)  - &
                    sqrt((n-1)*(n-1+alpha)/(n*(n+alpha))) * psi_ho3d(n-2,l,r)
            end if
            
        end function psi_ho3d
        
        
        recursive function lp12f(l) result(ss)   ! (l+1/2)!
        
            implicit none
            integer :: l
            real(dp) :: ss
            
            if ( l == 0 ) then  
                ss = sqrt(pi)/2
            else
                ss = (l + 0.5_dp ) * lp12f(l-1)
            end if
        
        end function lp12f
    
end program ho_test



