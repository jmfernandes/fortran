
program exp_cf_rec

	use NumType
	use cf_approx
	implicit none

	!=========================== CONSTANTS ===========================

	integer,	parameter	::	n_basis		=	8,					&
								bs			=	3,					&
								nb 			=	(n_basis+1)/bs-1,	&
								steps 		=	ncf

	real(dp),	parameter	::	mass		=	1.0_dp,				& 
								hbar		=	1.0_dp,				&
								omega_h		=	1.0_dp,				&
								omega_b		=	2._dp,				&
								x			=	1.0_dp

	!=========================== MATRICES ===========================

	real(dp)				::	x_mat(0:n_basis,0:n_basis),			&
								p_mat(0:n_basis,0:n_basis),			&
								h_mat(0:n_basis,0:n_basis),			&
								eigen_mat(0:n_basis,0:n_basis),		&
								h0_mat(0:bs-1,0:bs-1,0:nb),			&
								v_mat(0:n_basis,0:n_basis),			&
								e_mat(0:bs-1,0:bs-1),				&
								g0_mat(0:bs-1,0:bs-1,0:nb),			&
								mult_mat(0:n_basis,0:n_basis),		&
								ss(0:n_basis,0:steps),			&
								psi_mat(0:n_basis),					&
								taylor(0:ncf),						&
								cf(0:ncf),							&
								gff(0:n_basis,0:n_basis,0:steps),	&
								gcc(0:n_basis,0:steps),				&
								g0_mat2(0:bs-1,0:bs-1,0:nb),		&
								g0_mat3(0:bs-1,0:bs-1,0:nb)

	!=========================== PARAMETERS ==========================

	integer,	parameter	::  lwork=(n_basis+2)*n_basis
	integer					::  n,i,j, block_size, info
	real(dp)				::  w_eigen(n_basis+1),work(lwork)



	!=================================================================

	x_mat	=	0._dp
	p_mat	=	0._dp
	h_mat	=	0._dp

	!build X
	do n=0,n_basis
		x_mat(n,n) = (hbar/(2*mass*omega_b))*(2*n+1)
		if (n>=2) then
			x_mat(n,n-2) = hbar/(2*mass*omega_b)*sqrt(real(n*(n-1)))
		end if
		if (n<=(n_basis-2)) then
			x_mat(n,n+2) = hbar/(2*mass*omega_b)*sqrt(real((n+1)*(n+2)))
		end if
	end do

	!build P
	do n=0,n_basis
		p_mat(n,n) = (hbar*mass*omega_b)/2*(2*n+1)
		if (n>=2) then
			p_mat(n,n-2) = -(hbar*mass*omega_b)/2*sqrt(real(n*(n-1)))
		end if
		if (n<=(n_basis-2)) then
			p_mat(n,n+2) = -(hbar*mass*omega_b)/2*sqrt(real((n+1)*(n+2)))
		end if
	end do

	h_mat(0:n_basis,0:n_basis) = p_mat(0:n_basis,0:n_basis)/(2*mass)+ &
							mass*omega_h**2*x_mat(0:n_basis,0:n_basis)/2

	print *, 'hamiltonian matrix-------------------------'
	do n=0,n_basis
		print '(20f10.2)', h_mat(n,0:n_basis)
	end do

	!============================================================================

	!   compute eigenvalues
	eigen_mat(0:n_basis,0:n_basis) = h_mat(0:n_basis,0:n_basis)

	call dsyev('V','U',n_basis+1,eigen_mat,n_basis+1,w_eigen,work,lwork,info)

! 	print *, 'eigenvalues-------------------------'
! 	print '(20f10.2)', w_eigen

	!============================================================================

	!build V matrix

	v_mat(0:n_basis,0:n_basis) = h_mat(0:n_basis,0:n_basis)

	do n=0,nb
		h0_mat(0:bs-1,0:bs-1,n)=h_mat(bs*n:bs*n+bs-1,bs*n:bs*n+bs-1)
		v_mat(bs*n:bs*n+bs-1,bs*n:bs*n+bs-1) = &
		h_mat(bs*n:bs*n+bs-1,bs*n:bs*n+bs-1) - h0_mat(0:bs-1,0:bs-1,n)
	end do


	print *, 'h0 matrix-------------------------'
	do n=0,nb
		print *, 'the n= ', n,' matrix'
		do i=0,bs-1
			print '(12f10.2)', h0_mat(i,0:bs-1,n)
		end do
	end do


	print *, "h' matrix-------------------------"
	do n=0,n_basis
		print '(20f10.2)', v_mat(n,0:n_basis)
	end do

	!============================================================================

	do n=0,bs-1
		e_mat(n,n)	= 0.5_dp
	end do


	do n=0,nb
		g0_mat(0:bs-1,0:bs-1,n)	=e_mat(0:bs-1,0:bs-1)-h0_mat(0:bs-1,0:bs-1,n)
		g0_mat3(0:bs-1,0:bs-1,n) = g0_mat(0:bs-1,0:bs-1,n)
		call inverse(g0_mat(0:bs-1,0:bs-1,n),g0_mat2(0:bs-1,0:bs-1,n),3)
		g0_mat3(0:bs-1,0:bs-1,n)=matmul(g0_mat3(0:bs-1,0:bs-1,n),g0_mat2(0:bs-1,0:bs-1,n))
	end do

! 	do n=0,nb
! 		g0_mat2(0:bs-1,0:bs-1,n) = inv(g0_mat(0:bs-1,0:bs-1,n))
! 	end do

	    print *, 'g0 matrix-------------------------'
	    do n=0,nb
	    	print *, 'the n= ', n,' matrix'
	    	do i=0,bs-1
	        	print '(12f10.2)', g0_mat(i,0:bs-1,n)
	    	end do
	    end do

! 	    call inverse(g0_mat(0:bs-1,0:bs-1,0),g0_mat2(0:bs-1,0:bs-1,0),3)

	    print *, 'g02 matrix-------------------------'
	    do n=0,nb
	    	print *, 'the n= ', n,' matrix'
	    	do i=0,bs-1
	        	print '(12f10.2)', g0_mat2(i,0:bs-1,n)
	    	end do
	    end do

	 	print *, 'g03 matrix-------------------------'
	    do n=0,nb
	    	print *, 'the n= ', n,' matrix'
	    	do i=0,bs-1
	        	print '(12f10.2)', g0_mat3(i,0:bs-1,n)
	    	end do
	    end do



	do n=0,nb
			do i=0,nb
				mult_mat(bs*n:bs*n+bs-1,bs*i:bs*i+bs-1)	= &
				matmul(g0_mat(0:bs-1,0:bs-1,n), v_mat(bs*n:bs*n+bs-1,bs*i:bs*i+bs-1))
			end do 
	end do

! 		print *, 'mult mat----------'
! 		do n=0,n_basis
! 	        	print '(12f10.2)', mult_mat(n,0:n_basis)
! 	    end do


	!======================================================================

	psi_mat(0:n_basis) = 1.0_dp/n_basis

	!doing it this way takes too long
! 	do n=0,steps
! 		ss(0:n_basis,n) = multiply_mat(n,mult_mat,psi_mat)
! 	end do

	gff(0:n_basis,0:n_basis,0) = mult_mat(0:n_basis,0:n_basis)
	do n=1,steps
		gff(0:n_basis,0:n_basis,n) = matmul(gff(0:n_basis,0:n_basis,n-1),mult_mat(0:n_basis,0:n_basis))
	end do

	gcc(0:n_basis,0) = psi_mat(0:n_basis)
	do n=1,steps
		gcc(0:n_basis,n) = matmul(gff(0:n_basis,0:n_basis,n-1),psi_mat(0:n_basis))
		taylor(n)=dot_product(gcc(0:n_basis,0),gcc(0:n_basis,n))
! 		print *, taylor(n)
	end do

	print *, '--------'

! 	do n=0,ncf
! 		taylor(n)=dot_product(multiply_mat(0,mult_mat,psi_mat),multiply_mat(n,mult_mat,psi_mat))
! 		print *, taylor(n)
! 	end do

	call taylor_cfrac(taylor,steps,cf)

	print *, 'cf=', 1/evalcf(cf,steps,x)

	print *, 'done'
	!=======================================================================

	contains

		! Returns the inverse of a matrix calculated by finding the LU
		! decomposition.  Depends on LAPACK.

		function inv(A) result(Ainv)
		  real(dp), dimension(0:bs-1,0:bs-1), intent(in) :: A
		  real(dp), dimension(0:bs-1,0:bs-1) :: Ainv

		  real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
		  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
		  integer :: n, info


		  ! External procedures defined in LAPACK
		  external DGETRF
		  external DGETRI

		  ! Store A in Ainv to prevent it from being overwritten by LAPACK
		  Ainv = A
		  n = size(A,1)

		  ! DGETRF computes an LU factorization of a general M-by-N matrix A
		  ! using partial pivoting with row interchanges.
		  call DGETRF(n, n, Ainv, n, ipiv, info)

		  if (info /= 0) then
			 stop 'Matrix is numerically singular!'
		  end if

		  ! DGETRI computes the inverse of a matrix using the LU factorization
		  ! computed by DGETRF.
		  call DGETRI(n, Ainv, n, ipiv, work, n, info)

		  if (info /= 0) then
			 stop 'Matrix inversion failed!'
		  end if
		end function inv

! 		function multiply_mat(n,mult_mat,psi_mat) result(ss)
		
! 			implicit none
! 			integer :: n,i
! 			real(dp) ::  mult_mat(0:n_basis,0:n_basis),&
! 							res(0:n_basis,0:n_basis),&
! 							ss(0:n_basis),&
! 							psi_mat(0:n_basis)
			
! 			!create placeholder matrix that will overwrite itself
! 			res(0:n_basis,0:n_basis) = mult_mat(0:n_basis,0:n_basis)

! 			!multiply (G0*V) with itself for n times (n is input)
! 			do i=1,n
! 			res(0:n_basis,0:n_basis) = matmul(res(0:n_basis,0:n_basis),&
! 				mult_mat(0:n_basis,0:n_basis))
! 			end do         
		
! 			!return the result as {(G0*V)^n}*psi
! 			if (n == 0) then
! 				ss(0:n_basis) = psi_mat(0:n_basis)
! 			else
! 				ss(0:n_basis)=matmul(res(0:n_basis,0:n_basis),psi_mat(0:n_basis))
! 			end if

! 		end function multiply_mat

end program exp_cf_rec

  subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse
