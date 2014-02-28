
program exp_cf_rec

    use NumType
    implicit none
    real(dp) :: a, x
    real(dp) :: xmin, xmax, dx
    integer :: imax, i

    xmin = -3
    xmax = 3
    imax = 50
    
    dx = (xmax - xmin)/imax
    
    a = 0.5_dp
    
    do i = 0, imax-1
    
        x = xmin + (i+0.5) * dx 
        
        print *, x, fun(x) , erfc(x)
        
    end do

    
    contains
    
    
    
        function fun (x) result(f)
        
            real(dp) :: a, x, f
            
            a = 0.5_dp
            
            if ( x >= 0 ) then 
                f = igamma(a,x*x)/sqrt(pi)
            else if ( x < 0 ) then
                f = 2 - igamma(a,x*x)/sqrt(pi)
            end if
        
        end function fun
    
        function igamma(a,z) result(ss)
        
            implicit none
            real(dp) :: a, z, ss
            
            ss = exp(-z)*z**a * 1/(z+contf(1,a,z))         
        
        end function igamma
    
        recursive function contf (i,a,z)  result(scf)
        
            implicit none
            integer :: i, imax
            real(dp) :: a, z, scf
            
            imax = 2000
            
            if ( i >= imax ) then
                scf = 0._dp
            else
                scf = (i-a)/(1+ i/(z + contf(i+1,a,z) )  )
            end if
        
        end function contf

end program exp_cf_rec

