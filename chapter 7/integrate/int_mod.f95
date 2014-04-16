
module integr

    use NumType
    integer, parameter :: maxint = 300
    
    contains

        subroutine cc11(n,a,b,x,w)  ! Clenshaw-Curtis
        !   n    number of points
        !   a, b interval
        !   x   abscissas
        !   w   weights 

            implicit none
            integer :: n, j, m
            real(dp) :: a, b, tj, ss
            real(dp), dimension(maxint) :: x, w
    
            do j=1, n
                tj=(j*pi)/(n+1)
                x(j)=(b-a)/2*cos(tj)+(b+a)/2
                ss=0._dp
                do m=1,n
                    ss=ss+sin(m*tj)*(1-cos(m*pi))/m
                end do
                w(j)=(b-a)*sin(tj)/(n+1)*ss
            end do
        
        end subroutine cc11


        subroutine cc0inf(n,scale,x,w) ! Clenshaw-Curtis
        !  0 to infinity
        !  n number of points
        !  scale to scale the function
        !  x  abscissas
        !  w  weights

            implicit none
            integer :: n, j, m
            real(dp) :: tj, ss, scale
            real(dp), dimension(maxint) :: x, w
    
            do j=1, n
                tj=(j*pi)/(n+1)
                x(j)=scale*(cot(tj/2))**2
                ss=0._dp
                do m=1,n
                    ss=ss+sin(m*tj)*(1-cos(m*pi))/m
                end do
                w(j)=4*scale*sin(tj)/((1-cos(tj))**2*(n+1))*ss
            end do
    
            contains 
    
                function cot(x)        
                    use NumType
                    real(dp) :: x, cot
            
                        cot=1/tan(x)  
                
                end function cot
        
        end subroutine cc0inf

    
      
        subroutine rombint( a, b, func, res, n, eps )

        !   res = int_a^b f(x) dx
        !   eps  required accuracy
        !   n    n_max in approx (input)
        !   n    accuaracy reached after n steps

            integer :: n 
            real(dp) :: a, b, eps, res
            real(dp), external :: func
            integer ::  np, i, j, k, m
            real(dp) :: h, sumt, r(maxint,maxint)

            h = b - a 
            np = 1 
            r(1,1) = h/2 * ( func(a) + func(b) )
            res = r(1,1)

            do i=2,n 
                h = h/2
                np = 2*np 
                sumt = 0.0_dp
                do k=1,(np-1),2
                    sumt = sumt + func( a + k*h)
                end do
                r(i,1) = 0.5_dp * r(i-1,1) + h * sumt
                m = 1
                do j=2,i 
                    m = 4*m
                    r(i,j) = r(i,j-1) + (r(i,j-1)-r(i-1,j-1))/(m-1)
                end do
                if ( abs(res-r(i,i)) < eps ) then
                    n = i
                    res = r(i,i)
                    return
                end if
                res = r(i,i)
            end do
            print *,' romint :',eps,r(i-1,i-1),res

        end subroutine rombint


end module integr

