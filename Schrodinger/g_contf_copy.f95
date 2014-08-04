
program exp_cf_rec

    use NumType
    use cf_approx
    implicit none
    real(dp) :: x,result, sign, num1,num2,num3
    real(dp), dimension(0:2) :: aa,bb,cc,dd
    integer  :: imax, i, j, n, info, ii, int1, int2
    integer, parameter :: entries =10
    real(dp) :: yy1(entries),yy2(entries),yy(entries)

    integer,    parameter   :: lwork=2*n_basis+1
    integer  :: ipiv(n_basis)
    real(dp),   parameter   ::  mass=1.0_dp, hbar=1.0_dp,&
                                omega_h=1._dp,omega_b=1/2._dp
    complex(dp) ::  x_mat(0:n_basis,0:n_basis+1),       &
                    p_mat(0:n_basis,0:n_basis+1),       &
                    h_mat(0:n_basis,0:n_basis),         &
                    v_mat(0:n_basis,0:n_basis),         &
                    work(lwork),                        &
                    sub_mat(0:n_basis,0:n_basis),       &
                    h0_mat(0:n_basis,0:n_basis),        &
                    result_mat(0:n_basis,0:n_basis),    &
                    e_mat(0:n_basis,0:n_basis),         &
                    mult_mat(0:n_basis,0:n_basis),      &
                    g0_mat(0:n_basis,0:n_basis),       &
                    psi_mat(0:n_basis),                 &
                    ss_mat(0:n_basis),                  &
                    pade_mat(0:n_basis)


    real(dp), dimension(0:n_basis,0:steps) :: taylor, cf                    
                    

    real(dp) :: rwork(3*(n_basis-2)), w_eigen(n_basis+1)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    x_mat       = 0._dp
    p_mat       = 0._dp
    h_mat       = 0._dp
    v_mat       = 0._dp
    sub_mat     = 0._dp
    h0_mat      = 0._dp
    result_mat  = 0._dp
    e_mat       = 0._dp
    mult_mat    = 0._dp
    g0_mat     = 0._dp
    psi_mat     = 0._dp
    ss_mat      = 0._dp


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


    !build the hamiltonian
    h_mat(0:n_basis,0:n_basis) = 1/(2*mass) * &
        matmul(p_mat(0:n_basis,0:n_basis+1),conjg(transpose(p_mat(0:n_basis,0:n_basis+1)))) + &
        mass*omega_h**2/2 * &
        matmul(x_mat(0:n_basis,0:n_basis+1),transpose(x_mat(0:n_basis,0:n_basis+1)))

!============================================================================
    
    !Build the h0 and V matrix
    h0_mat(0:n_basis,0:n_basis) = h_mat(0:n_basis,0:n_basis)
    v_mat(0:n_basis,0:n_basis) = h_mat(0:n_basis,0:n_basis)

    !Need to change this for bigger matrices. Only works for 3x3
    do i=0,n_basis
        v_mat(i,i) = 0
    end do

    !this is fine for bigger matrices
    do i=0,n_basis
        do j=0,n_basis
            h0_mat(i,j)=h_mat(i,j)-v_mat(i,j)
        end do    
    end do

!===============================================================================

    open(unit=2,file='pfaff_Hmatrix_N8_l6.5.dat')
    read (2,100) int1, int2, num3
    100 format(1I1,1x,1I1,1x,1F1.5)
!         150 format (I3, 1x, I3, 1x, F7.2)
!     end do
    close(2)


!============================================================================

    !set up matrices for energy and psi - Need to get from PETERSON
    e_mat(0:n_basis,0:n_basis) = 1._dp

    psi_mat(0:n_basis) = 10._dp 

    !set up g0
    sub_mat(0:n_basis,0:n_basis) = e_mat(0:n_basis,0:n_basis)-h0_mat(0:n_basis,0:n_basis)

    g0_mat(0:n_basis,0:n_basis) = inv(sub_mat(0:n_basis,0:n_basis))


!===============================================================================

    !initialize x
    x = 1._dp

    !this is the matrix that represents (G0*V)
    mult_mat(0:n_basis,0:n_basis) = matmul(g0_mat(0:n_basis,0:n_basis),&
        v_mat(0:n_basis,0:n_basis))

    !build the coefficients for the power series
    do i = 0, n_basis+1
        ss_mat = multiply_mat(i,mult_mat,psi_mat)
        do j =0,steps
            if (i==0) then
                taylor(j,i) = psi_mat(j)
            else
                taylor(j,i+1) = ss_mat(j)
            end if 
        end do
    end do

!==================================================================================

    print *, 'hamiltonian matrix-------------------------'
    do n=0,n_basis
        print '(12f10.2)', dble(h_mat(n,0:n_basis))
    end do

    print *, 'h0 matrix-------------------------'
    do n=0,n_basis
        print '(12f10.2)', dble(h0_mat(n,0:n_basis))
    end do

    print *, 'potential matrix-------------------------'
    do n=0,n_basis
        print '(12f10.2)', dble(v_mat(n,0:n_basis))
    end do

    print *, 'G0 matrix-------------------------'
    do n=0,n_basis
        print '(12f10.2)', dble(g0_mat(n,0:n_basis))
    end do

    print *, 'taylor matrix-------STOPS at 4 steps instead of 50'
    do n=0,n_basis
        print '(12f12.2)', dble(taylor(n,0:steps))
    end do
    

!====================================================================

    !I'm not sure if this is right
    print *, 'pade approx--------------------------------'
    n = 3
    call taylor_cfrac(taylor,n,cf)  

    pade_mat = evalcf(cf,n,x)
    print *, 'result array='
    do i=0,n_basis
        print *, dble(pade_mat(i))
    end do

!=======================================================================


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

!     mult_mat(0:n_basis,0:n_basis) = matmul(sub_mat(0:n_basis,0:n_basis), &
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


end program exp_cf_rec

