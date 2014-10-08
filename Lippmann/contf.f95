
program exp_cf_rec

	use NumType
	use cf_approx
	implicit none

	!=========================== CONSTANTS ===========================

	integer,	parameter	::	n_basis		=	99,					&
								bs			=	10,					&
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
								cf(0:ncf)

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

! 	print *, 'hamiltonian matrix-------------------------'
! 	do n=0,n_basis
! 		print '(20f10.2)', h_mat(n,0:n_basis)
! 	end do

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


! 	print *, 'h0 matrix-------------------------'
! 	do n=0,nb
! 		print *, 'the n= ', n,' matrix'
! 		do i=0,bs-1
! 			print '(12f10.2)', h0_mat(i,0:bs-1,n)
! 		end do
! 	end do


! 	print *, "h' matrix-------------------------"
! 	do n=0,n_basis
! 		print '(20f10.2)', v_mat(n,0:n_basis)
! 	end do

	!============================================================================

	do n=1,bs
		e_mat(n,n)	= 7._dp
	end do


	do n=0,nb
		g0_mat(0:bs-1,0:bs-1,n)	=inv(e_mat(0:bs-1,0:bs-1)-h0_mat(0:bs-1,0:bs-1,n))
	end do

	!     print *, 'g0 matrix-------------------------'
	!     do n=0,nb
	!     	print *, 'the n= ', n,' matrix'
	!     	do i=0,bs-1
	!         	print '(12f10.2)', g0_mat(i,0:bs-1,n)
	!     	end do
	!     end do

	do n=0,nb
			do i=0,nb
				mult_mat(bs*n:bs*n+bs-1,bs*i:bs*i+bs-1)	= &
				matmul(g0_mat(0:bs-1,0:bs-1,n), v_mat(bs*n:bs*n+bs-1,bs*i:bs*i+bs-1))
			end do 
	end do

	! 	print *, 'mult mat----------'
	! 	do n=0,n_basis
	!         	print '(12f10.2)', mult_mat(n,0:n_basis)
	!     end do


	!======================================================================

	psi_mat(0:n_basis) = 1.0_dp/n_basis

	!doing it this way takes too long
! 	do n=0,steps
! 		ss(0:n_basis,n) = multiply_mat(n,mult_mat,psi_mat)
! 	end do

	do n=0,ncf
		taylor(n)=dot_product(multiply_mat(0,mult_mat,psi_mat),multiply_mat(n,mult_mat,psi_mat))
		print *, taylor(n)
	end do

	call taylor_cfrac(taylor,steps,cf)

	print *, 'cf=', evalcf(cf,steps,x)

	!=======================================================================

	contains

		! Returns the inverse of a matrix calculated by finding the LU
		! decomposition.  Depends on LAPACK.

		function inv(A) result(Ainv)
		  real(dp), dimension(1:bs,1:bs), intent(in) :: A
		  real(dp), dimension(1:bs,1:bs) :: Ainv

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

		function multiply_mat(n,mult_mat,psi_mat) result(ss)
		
			implicit none
			integer :: n,i
			real(dp) ::  mult_mat(0:n_basis,0:n_basis),&
							res(0:n_basis,0:n_basis),&
							ss(0:n_basis),&
							psi_mat(0:n_basis)
			
			!create placeholder matrix that will overwrite itself
			res(0:n_basis,0:n_basis) = mult_mat(0:n_basis,0:n_basis)

			!multiply (G0*V) with itself for n times (n is input)
			do i=1,n
			res(0:n_basis,0:n_basis) = matmul(res(0:n_basis,0:n_basis),&
				mult_mat(0:n_basis,0:n_basis))
			end do         
		
			!return the result as {(G0*V)^n}*psi
			if (n == 0) then
				ss(0:n_basis) = psi_mat(0:n_basis)
			else
				ss(0:n_basis)=matmul(res(0:n_basis,0:n_basis),psi_mat(0:n_basis))
			end if

		end function multiply_mat

end program exp_cf_rec

