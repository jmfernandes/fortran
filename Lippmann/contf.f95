
program exp_cf_rec

	use NumType
	use cf_approx
	use chebyshev
	implicit none

	!=========================== CONSTANTS ===========================

	integer,	parameter	::	n_basis		=	99,					&
								bs			=	10,					&
								nb 			=	(n_basis+1)/bs-1,	&
								steps 		=	ncf,				&
								nch			=	4

	real(dp),	parameter	::	mass		=	1.0_dp,				& 
								hbar		=	1.0_dp,				&
								omega_h		=	1.0_dp,				&
								omega_b		=	2.0_dp,				&
								x			=	1.0_dp

	!=========================== MATRICES ===========================

	real(dp)				::	eigen_mat(0:n_basis,0:n_basis)

	!=========================== PARAMETERS ==========================

	integer,	parameter	::  lwork=(n_basis+2)*n_basis
	integer					::  n,i,j, block_size, info
	real(dp)				::  w_eigen(n_basis+1),work(lwork), offset
	real(dp)				::	ya,yb, smallthing, dz, zz, r, num


	!=================================================================

	!   compute eigenvalues using LAPACK
! 	eigen_mat(0:n_basis,0:n_basis) = h_mat(0:n_basis,0:n_basis)

! 	call dsyev('V','U',n_basis+1,eigen_mat,n_basis+1,w_eigen,work,lwork,info)

! 	print *, 'LAPACK eigenvalues-------------------------'
! 	print '(10f10.2)', w_eigen

	! 	do n = 0,n_basis
	! 		print '(A10,100f10.2)', '  vector  ', eigen_mat(0:n_basis,n) 
	! 	end do


	!============================================================================


	!build X
! 	do n=0,n_basis
! 		x_mat(n,n) = (hbar/(2*mass*omega_b))*(2*n+1)
! 		if (n>=2) then
! 			x_mat(n,n-2) = hbar/(2*mass*omega_b)*sqrt(real(n*(n-1)))
! 		end if
! 		if (n<=(n_basis-2)) then
! 			x_mat(n,n+2) = hbar/(2*mass*omega_b)*sqrt(real((n+1)*(n+2)))
! 		end if
! 	end do

! 	!build P
! 	do n=0,n_basis
! 		p_mat(n,n) = (hbar*mass*omega_b)/2*(2*n+1)
! 		if (n>=2) then
! 			p_mat(n,n-2) = -(hbar*mass*omega_b)/2*sqrt(real(n*(n-1)))
! 		end if
! 		if (n<=(n_basis-2)) then
! 			p_mat(n,n+2) = -(hbar*mass*omega_b)/2*sqrt(real((n+1)*(n+2)))
! 		end if
! 	end do

! 	h_mat(0:n_basis,0:n_basis) = p_mat(0:n_basis,0:n_basis)/(2*mass)+ &
! 							mass*omega_h**2*x_mat(0:n_basis,0:n_basis)/2

! 	!============================================================================

! 	!build V matrix

! 	v_mat(0:n_basis,0:n_basis) = h_mat(0:n_basis,0:n_basis)

! 	do n=0,nb
! 		h0_mat(0:bs-1,0:bs-1,n)=h_mat(bs*n:bs*n+bs-1,bs*n:bs*n+bs-1)
! 		v_mat(bs*n:bs*n+bs-1,bs*n:bs*n+bs-1) = &
! 		h_mat(bs*n:bs*n+bs-1,bs*n:bs*n+bs-1) - h0_mat(0:bs-1,0:bs-1,n)
! 	end do

! 	!build psi matrix
! 	psi_mat(0:n_basis) = 1.0_dp/sqrt(real(n_basis))

	!============================================================================

	!calculate the eigenvalues using chebyshev
	print *, 'calculated eigenvalues-------------------------'
	do j=0,9
		ya=j
		yb=j+1
		call chebyex(mcalc, nch, cheb, ya, yb)
		call chebyzero(nch, cheb, ya, yb, z0, iz0)

! 		print '(A6,1X,I2.1,A4,I2.1,5X,5f10.5)', 'range=',j,' to ', j+1, z0(1:iz0)
		print *, 'range=',ya,' to ', yb, z0(1:iz0)
	end do


! 	print *, floor(100*rand()),floor(100*rand()),floor(100*rand()),floor(100*rand())

!	smallthing = 1.e-15_dp
!    dz = 0.01_dp

!     do i = 1,1
!     	zz = 5.0
!     	print *, zz
!     call root_polish(mcalc,zz,dz,smallthing,100) 
!     print *, zz
! 	end do

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

		function mcalc(energy) result(res_sum)

			real(dp) :: energy, res_sum, cf_num, res
			real(dp) :: g0_mat(0:bs-1,0:bs-1,0:nb),			&
						mult_mat(0:n_basis,0:n_basis),		&
						taylor(0:ncf),						&
						cf(0:ncf),							&
						gff(0:n_basis,0:n_basis,0:steps),	&
						gcc(0:n_basis,0:steps),				&
						e_mat(0:bs-1,0:bs-1), &
						x_mat(0:n_basis,0:n_basis),			&
						p_mat(0:n_basis,0:n_basis),			&
						h_mat(0:n_basis,0:n_basis),			&
						h0_mat(0:bs-1,0:bs-1,0:nb),			&
						v_mat(0:n_basis,0:n_basis),			&
						psi_mat(0:n_basis)


			g0_mat(0:bs-1,0:bs-1,0:nb)= 0._dp
			mult_mat(0:n_basis,0:n_basis)=0._dp
			taylor(0:ncf)=0._dp
			cf(0:ncf)=0._dp
			gff(0:n_basis,0:n_basis,0:steps)=0._dp
			gcc(0:n_basis,0:steps)=0._dp
			e_mat(0:bs-1,0:bs-1)=0._dp
			x_mat(0:n_basis,0:n_basis)= 0._dp
			p_mat(0:n_basis,0:n_basis)= 0._dp
			h_mat(0:n_basis,0:n_basis)= 0._dp
			h0_mat(0:bs-1,0:bs-1,0:nb)= 0._dp
			v_mat(0:n_basis,0:n_basis)= 0._dp
			psi_mat(0:n_basis)=0._dp

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

			!============================================================================

			!build V matrix

			v_mat(0:n_basis,0:n_basis) = h_mat(0:n_basis,0:n_basis)

			do n=0,nb
				h0_mat(0:bs-1,0:bs-1,n)=h_mat(bs*n:bs*n+bs-1,bs*n:bs*n+bs-1)
				v_mat(bs*n:bs*n+bs-1,bs*n:bs*n+bs-1) = &
				h_mat(bs*n:bs*n+bs-1,bs*n:bs*n+bs-1) - h0_mat(0:bs-1,0:bs-1,n)
			end do

			!build psi matrix
			psi_mat(0:n_basis) = 1.0_dp/sqrt(real(n_basis))


			do n=0,bs-1
				e_mat(n,n)	=  energy
			end do

			do n=0,nb
				g0_mat(0:bs-1,0:bs-1,n)	=inv(e_mat(0:bs-1,0:bs-1)-h0_mat(0:bs-1,0:bs-1,n))
			end do


			do n=0,nb
					do i=0,nb
						mult_mat(bs*n:bs*n+bs-1,bs*i:bs*i+bs-1)	= &
						matmul(g0_mat(0:bs-1,0:bs-1,n), v_mat(bs*n:bs*n+bs-1,bs*i:bs*i+bs-1))
					end do 
			end do

			gff(0:n_basis,0:n_basis,0) = 1._dp
			gff(0:n_basis,0:n_basis,1) = mult_mat(0:n_basis,0:n_basis)
			do n=2,steps
				gff(0:n_basis,0:n_basis,n) = matmul(gff(0:n_basis,0:n_basis,n-1),mult_mat(0:n_basis,0:n_basis))
			end do
			res = 0._dp
			gcc(0:n_basis,0) = psi_mat(0:n_basis)
			do n=0,steps
				gcc(0:n_basis,n) = matmul(gff(0:n_basis,0:n_basis,n),psi_mat(0:n_basis))
				taylor(n) = dot_product(gcc(0:n_basis,0),gcc(0:n_basis,n))
! 				print *, taylor(n), 'taylor'
				res = res + taylor(n)
				!add small variation
! 				taylor(n) = taylor(n)+(rand()*(taylor(n)*1e-12))
			end do
! 			print *, energy, '-', res, 'res'

			call taylor_cfrac(taylor,steps,cf)

			cf_num = evalcf(cf,steps,x)

! 			print *, energy, cf_num, 'final energy'

			res_sum = abs(1/cf_num)

		end function mcalc



end program exp_cf_rec

