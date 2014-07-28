
program exp_cf_rec

    use NumType
    use cf_approx
    implicit none
    real(dp) :: x,result, sign
    real(dp), dimension(0:2) :: aa,bb,cc,dd
    integer :: imax, i, j, n, info

    integer,    parameter   ::  n_basis=2,lwork=2*n_basis+1,steps=4
    integer :: ipiv(n_basis)
    real(dp),   parameter   ::  mass=1.0_dp, hbar=1.0_dp,&
                                omega_h=1._dp,omega_b=1/2._dp
    complex(dp) ::  x_mat(0:n_basis,0:n_basis+1), p_mat(0:n_basis,0:n_basis+1), &
                    x2_mat(0:n_basis,0:n_basis), p2_mat(0:n_basis,0:n_basis), &
                    h_mat(0:n_basis,0:n_basis),v_mat(0:n_basis,0:n_basis), &
                    work(lwork), g_mat(0:n_basis,0:n_basis),&
                    mult_mat(0:n_basis,0:n_basis), &
                    mult2_mat(0:n_basis), &
                    h0(0:n_basis,0:n_basis), &
                    unit_mat(0:n_basis,0:n_basis), &
                    result_mat(0:n_basis,0:n_basis), &
                    e_mat(0:n_basis,0:n_basis), &
                    inv_mat(0:n_basis,0:n_basis), &
                    f(0:n_basis,0:n_basis), &
                    c(0:n_basis,0:n_basis), &
                    d(0:n_basis,0:n_basis), &
                    psi_mat(0:n_basis),&
                    ss_mat(0:n_basis) 


    real(dp), dimension(0:n_basis,0:steps) :: taylor, cf                    
                    

    real(dp) :: rwork(3*(n_basis-2)), w_eigen(n_basis+1)

    x=4
    result = fun1(x)
    print *, 'exp(4) =',result, 'true value =', exp(x)

    x=3
    result = fun1_alt(x)
    print *, 'exp(4) =',result, 'true value =', exp(x)

    x=1
    result = fun2(x)
    print *, '1/(1+x) =',result, 'true value =', 1/(1._dp+x)

    x=1
    result = fun3(x)
    print *, 'log(1+x) =',result, ' true vale =', log(1+x)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    x_mat = 0._dp
    p_mat = 0._dp
    g_mat = 0._dp
    mult_mat = 0._dp
    unit_mat = 0._dp
    psi_mat = 0._dp
    result_mat = 0._dp


    do n=0,n_basis-1
        x_mat(n,n+1) = sqrt(hbar/(2*mass*omega_b))*sqrt(n+1._dp)
        x_mat(n+1,n) = x_mat(n,n+1)     
    end do
    x_mat(n_basis,n_basis+1) = sqrt(hbar/(2*mass*omega_b))*sqrt(n_basis+1._dp) !add last point


    do n=0,n_basis-1
        p_mat(n,n+1) = -iic*sqrt(mass*hbar*omega_b/2)*sqrt(n+1._dp)
        p_mat(n+1,n) = -p_mat(n,n+1)     
    end do
    p_mat(n_basis,n_basis+1) = -sqrt(mass*hbar*omega_b/2)*sqrt(n_basis+1._dp) !add last point

!     print *, 'x matrix-------------------------'

!     do n=0,n_basis
!         print '(12f10.3)', x_mat(n,0:n_basis+1)
!     end do

!     v_mat(0:n_basis,0:n_basis) = (1/2._dp)*mass*omega_h**2*&
!         matmul(x_mat(0:n_basis,0:n_basis+1),transpose(x_mat(0:n_basis,0:n_basis+1)))

!     print *, 'V matrix-------------------------'
!     do n=0,n_basis
!         print '(12f10.2)', v_mat(n,0:n_basis)
!     end do

    !!!!!!!!!!!!!!!!!!!

!     do n=0,n_basis
!         g_mat(n,n) = (n+1)/2._dp    
!     end do

    !!!!!!!!!!!!!!!!!!

!     do n=0,n_basis
!         unit_mat(n,n) = 1._dp   
!     end do

!     do n=0,n_basis
!         psi_mat(n,n) = 1._dp   
!     end do

!     print *, 'g matrix-------------------------'
!     do n=0,n_basis
!         print '(12f12.2)', g_mat(n,0:n_basis)
!     end do

    h_mat(0:n_basis,0:n_basis) = 1/(2*mass) * &
        matmul(p_mat(0:n_basis,0:n_basis+1),conjg(transpose(p_mat(0:n_basis,0:n_basis+1)))) + &
        mass*omega_h**2/2 * &
        matmul(x_mat(0:n_basis,0:n_basis+1),transpose(x_mat(0:n_basis,0:n_basis+1)))

!     h_mat(0:n_basis,0:n_basis) = 2._dp
    h0(0:n_basis,0:n_basis) = h_mat(0:n_basis,0:n_basis)
    v_mat(0:n_basis,0:n_basis) = h_mat(0:n_basis,0:n_basis)

    do i=0,n_basis
        v_mat(i,i) = 0
    end do

    do i=0,n_basis
        do j=0,n_basis
            h0(i,j)=h_mat(i,j)-v_mat(i,j)
        end do    
    end do


    e_mat(0:n_basis,0:n_basis) = 1._dp

    g_mat(0:n_basis,0:n_basis) = e_mat(0:n_basis,0:n_basis)-h0(0:n_basis,0:n_basis)

    inv_mat(0:n_basis,0:n_basis) = inv(g_mat(0:n_basis,0:n_basis))


!         print *, 'h matrix-------------------------'

!     h_mat=matmul(h_mat,h_mat)

    print *, 'h matrix--------------------'
    do n=0,n_basis
        print '(12f12.2)', dble(h_mat(n,0:n_basis))
    end do
    print *, 'h0 matrix--------------------'
    do n=0,n_basis
        print '(12f12.2)', dble(h0(n,0:n_basis))
    end do
    print *, 'v matrix--------------------'
    do n=0,n_basis
        print '(12f12.2)', dble(v_mat(n,0:n_basis))
    end do
    print *, '--------------------'

    !=======================================================


    x = 1._dp

    mult_mat(0:n_basis,0:n_basis) = matmul(inv_mat(0:n_basis,0:n_basis),&
        v_mat(0:n_basis,0:n_basis))

    psi_mat(0:n_basis) = 1._dp 

    mult2_mat(0:n_basis) = matmul(mult_mat(0:n_basis,0:n_basis),&
        psi_mat(0:n_basis)) 

    print *, 'mult matrix--------------------'
    do n=0,n_basis
        print '(12f12.2)', dble(mult_mat(n,0:n_basis))
    end do

!     aa = (/1,2,3/)
!     bb = (/5,6,7/)

!     aa = bb/aa

!     print *, aa

!   Taylor coefficients of logarithmic function

!     taylor(1:2,1:2) = reshape((/ 1._dp,2._dp,3._dp,4._dp /), &
!         (/2,2/))

!     taylor(0) = 0._dp
!     sign = 1._dp

!     ss_mat = multiply_mat(3,mult_mat,psi_mat)

    print *, 'psi matrix--------------------'
    do n=0,n_basis
        print '(12f12.2)', dble(psi_mat(n))
    end do


    do i = 0, 2
        ss_mat = multiply_mat(i,mult_mat,psi_mat)
        do j =0,steps
            if (i==0) then
                taylor(j,i) = psi_mat(j)
            else
                taylor(j,i+1) = ss_mat(j)
            end if 
        end do
    end do

    print *, 'taylor matrix--------------------'
    do n=0,n_basis
        print '(12f12.2)', dble(taylor(n,0:3))
    end do
    
!   print *,'  taylor sum ',x,'=', horner(taylor,n,x)
    
!   call taylor_cfrac(taylor,n,cf)  

!   print *,'          cf ',x,' =',evalcf(cf,n,x)
!     print *,'         log ',x,'=',log(1+x)    

!     print *,'--------------------------'



    !========================================================


!     e_mat=h_mat

!     call zheev('v','u',n_basis+1,e_mat,n_basis+1,w_eigen,work,lwork,rwork,info)

!     print *, '----------'
!     print *, info

!     print *, '------- Eigen Values ------- '

!     print *, w_eigen(1:n_basis)

!     print *, '------- Eigen Vectors ------- '

!     do i = 1,n_basis
!         print '(a,200f6.2)', 'vector', dble(e_mat(1:n_basis,i)) 
!     end do

!     f(1:n_basis,1:n_basis) = matmul(e_mat(1:n_basis,1:n_basis),&
!         transpose(h_mat(1:n_basis,1:n_basis)))

!     print *, 'f matrix-------------------------'
!     do n=1,n_basis
!         print '(12f10.2)', dble(f(n,1:n_basis))
!     end do


!     print *, '------- Orthogonality ------- '

!     do j = 1,n_basis
!             do i = 1,n_basis
!                 print '(i5,i5,2f5.0)', i, j, dot_product(e_mat(1:n_basis,i) ,e_mat(1:n_basis,j) ) !determines if matrix is orthogonal
            
!             end do
!         end do

!     mult_mat(0:n_basis,0:n_basis) = matmul(g_mat(0:n_basis,0:n_basis), &
!         v_mat(0:n_basis,0:n_basis))


!     print *, '------- completeness & spectral decomp------- '
        
!         do j = 1,n_basis
!             do i = 1,n_basis
!                 f(i,j)= dot_product(e_mat(i,1:n_basis) ,&
!                 e_mat(j,1:n_basis)*w_eigen(1:n_basis))  !determines if matrix is complete
!             end do
!         end do 
        
!         Do i= 1,n_basis
!             print '(10f7.2)', f(i,1:n_basis)
!         end do 
        
!         print *, '+++++++++++++++++++++++++++++++++++'

!         c(1:n_basis,1:n_basis) = conjg(transpose(e_mat(1:n_basis,1:n_basis)))
        
!         f(1:n_basis,1:n_basis) = matmul(c(1:n_basis,1:n_basis),e_mat(1:n_basis,1:n_basis))
!         f(1:n_basis,1:n_basis) = matmul(e_mat(1:n_basis,1:n_basis),c(1:n_basis,1:n_basis))
        
!         Do i= 1,n_basis
!             print '(10f7.2)', f(i,1:n_basis)
!         end do 
        
!         print *, '++++++++++++++++++++++++++++++++++++++++++++'

!         f(1:n_basis,1:n_basis) = matmul(h_mat(1:n_basis,1:n_basis),e_mat(1:n_basis,1:n_basis))
!         d(1:n_basis,1:n_basis) = matmul(c(1:n_basis,1:n_basis),f(1:n_basis,1:n_basis))
        
!         Do i= 1,n_basis
!             print '(10f7.2)', d(i,1:n_basis)
!         end do 


!     print *, 'multiplied matrix-------------------------'
!     do n=0,n_basis
!         print '(12f12.3)', mult_mat(n,0:n_basis)
!     end do


!     print *, 'result matrix-------------------------'
!     x=1
!     result_mat = fun4(x,psi_mat)
! !     print *, 'result =',result
!     do n=0,n_basis
!         print '(12f12.3)', result_mat(n,0:n_basis)
!     end do


    !!!!!!!!!!!!!!!!!

    contains

        ! Returns the inverse of a matrix calculated by finding the LU
        ! decomposition.  Depends on LAPACK.

        function inv(A) result(Ainv)
          complex(dp), dimension(:,:), intent(in) :: A
          complex(dp), dimension(size(A,1),size(A,2)) :: Ainv

          real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
          integer, dimension(size(A,1)) :: ipiv   ! pivot indices
          integer :: n, info

          ! External procedures defined in LAPACK
          external ZGETRF
          external ZGETRI

          ! Store A in Ainv to prevent it from being overwritten by LAPACK
          Ainv = A
          n = size(A,1)

          ! DGETRF computes an LU factorization of a general M-by-N matrix A
          ! using partial pivoting with row interchanges.
          call ZGETRF(n, n, Ainv, n, ipiv, info)

          if (info /= 0) then
             stop 'Matrix is numerically singular!'
          end if

          ! DGETRI computes the inverse of a matrix using the LU factorization
          ! computed by DGETRF.
          call ZGETRI(n, Ainv, n, ipiv, work, n, info)

          if (info /= 0) then
             stop 'Matrix inversion failed!'
          end if
        end function inv

        function multiply_mat(n,mult_mat,psi_mat) result(ss)
        
            implicit none
            integer :: n,i
            complex(dp) ::  mult_mat(0:n_basis,0:n_basis),&
                            res(0:n_basis,0:n_basis),&
                            ss(0:n_basis),&
                            psi_mat(0:n_basis)
            
            res(0:n_basis,0:n_basis) = mult_mat(0:n_basis,0:n_basis)

            do i=1,n
            res(0:n_basis,0:n_basis) = matmul(res(0:n_basis,0:n_basis),&
                mult_mat(0:n_basis,0:n_basis))
            end do         
        
            ss(0:n_basis)=matmul(res(0:n_basis,0:n_basis),psi_mat(0:n_basis))

        end function multiply_mat


    
        function fun1(x) result(ss)
        
            implicit none
            real(dp) :: x, ss
            
            ss = 1/(1-contf1(1,x))         
        
        end function fun1
    
        recursive function contf1(i,x)  result(scf)
        
            implicit none
            integer :: i, imax
            real(dp) :: scf, x

            imax = 2000
            
            if (i >= imax) then
                scf = 0._dp
            else
                scf = (x/i)/(1+(x/i)-contf1(i+1,x))
            end if 
        
        end function contf1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        function fun1_alt(x) result(ss)
        
            implicit none
            real(dp) :: x, ss
            
            ss = 1/(1-(x/contf1_alt(1,x)))         
        
        end function fun1_alt
    
        recursive function contf1_alt(i,x)  result(scf)
        
            implicit none
            integer :: i, imax
            real(dp) :: scf, x

            imax = 2000
            
            if (i >= imax) then
                scf = 0._dp
            else
                scf = i+x-((i*x)/contf1_alt(i+1,x))
            end if 
        
        end function contf1_alt

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        function fun2(x) result(ss)
        
            implicit none
            real(dp) :: x, ss
            
            ss = 1+contf2(1,x)         
        
        end function fun2
    
        recursive function contf2(i,x)  result(scf)
        
            implicit none
            integer :: i, imax
            real(dp) :: scf, x

            imax = 500
            
            if (i >= imax) then
                scf = 0._dp
            else
                scf = ((-1)**i)*x/(1+contf2(i+1,x))
            end if 
        
        end function contf2

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        function fun3(x) result(ss)
        
            implicit none
            real(dp) :: x, ss
            
            ss = x/(1+contf3(1,x))         
        
        end function fun3
    
        recursive function contf3(i,x)  result(scf)
        
            implicit none
            integer :: i, imax
            real(dp) :: scf, x

            imax = 500
            
            if (i >= imax) then
                scf = 0._dp
            else
                scf = ((i**2)*x)/((i+1)-(i*x)+contf3(i+1,x))
            end if 
        
        end function contf3

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        function fun4(x,psi_mat) result(ss)
        
            implicit none
            real(dp) :: x
            complex(dp) ::  psi_mat(0:n_basis,0:n_basis),&
                            ss(0:n_basis,0:n_basis)
            
            ss = psi_mat/(unit_mat-contf4(1,x))         
        
        end function fun4
    
        recursive function contf4(i,x)  result(scf)
        
            implicit none
            integer :: i, imax
            real(dp) :: x
            complex(dp) :: scf(0:n_basis,0:n_basis)

            imax = 500
            
            if (i >= imax) then
                scf = 0._dp
            else
                scf = (mult_mat**i)/(unit_mat+(mult_mat**i)-contf4(i+1,x))
            end if 
        
        end function contf4

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end program exp_cf_rec

