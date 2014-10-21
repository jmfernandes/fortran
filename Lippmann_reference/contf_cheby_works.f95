
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
								nch			=	10

	real(dp),	parameter	::	mass		=	1.0_dp,				& 
								hbar		=	1.0_dp,				&
								omega_h		=	1.0_dp,				&
								omega_b		=	1.1_dp,				&
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
								gcc(0:n_basis,0:steps)

	!=========================== PARAMETERS ==========================

	integer,	parameter	::  lwork=(n_basis+2)*n_basis
	integer					::  n,i,j, block_size, info
	real(dp)				::  w_eigen(n_basis+1),work(lwork), offset
	real(dp)				::	ya,yb


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

	print *, 'eigenvalues-------------------------'
	print '(20f10.2)', w_eigen

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

	do n=0,bs-1
		e_mat(n,n)	=  0.5_dp
	end do


	do n=0,nb
		g0_mat(0:bs-1,0:bs-1,n)	=inv(e_mat(0:bs-1,0:bs-1)-h0_mat(0:bs-1,0:bs-1,n))
	end do

! 	    print *, 'g0 matrix-------------------------'
! 	    do n=0,nb
! 	    	print *, 'the n= ', n,' matrix'
! 	    	do i=0,bs-1
! 	        	print '(12f10.2)', g0_mat(i,0:bs-1,n)
! 	    	end do
! 	    end do

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

	psi_mat(0:n_basis) = 1.0_dp/sqrt(real(n_basis))

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
! ! 		print *, taylor(n)
! 	end do

! 	call taylor_cfrac(taylor,steps,cf)

! 	print *, 'cf=', abs(evalcf(cf,steps,x))

	do j=0,20
		ya=j
		yb=j+1
		call chebyex(mcalc, nch, cheb, ya, yb)
		call chebyzero(nch, cheb, ya, yb, z0, iz0)

		print *, 'j=',j, z0(1:iz0)
	end do

! 	do i=1,iz0
!         offset = mcalc(z0(i))
!         print *, offset, z0(i), 'banana'
!     end do
! 	call chebyex(1/evalcf(cf,steps,x), nch, cheb, ya, yb)

!     call chebyzero(nch, cheb, ya, yb, z0, iz0)
    
!     print *, iz0, z0(1:iz0)*180._dp/pi


!     do i=1,iz0
!         offset = shoot(z0(i))
!         print *, offset
!     end do

! 	print *, 'done'
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

		function mcalc(energy) result(mine)

			real(dp) :: e, de, mine, cf_num, energy
			real(dp) :: g0_mat(0:bs-1,0:bs-1,0:nb),			&
						mult_mat(0:n_basis,0:n_basis),		&
						ss(0:n_basis,0:steps),			&
						psi_mat(0:n_basis),					&
						taylor(0:ncf),						&
						cf(0:ncf),							&
						gff(0:n_basis,0:n_basis,0:steps),	&
						gcc(0:n_basis,0:steps),				&
						e_mat(0:bs-1,0:bs-1)

			e = energy
! 			de = 0.1_dp

			do n=0,bs-1
				e_mat(n,n)	=  energy
			end do

			

! 			do while (e <= 10)
! 				do n=0,bs-1
! 					e_mat(n,n)	= energy(n,n)
! 				end do

				do n=0,nb
					g0_mat(0:bs-1,0:bs-1,n)	=inv(e_mat(0:bs-1,0:bs-1)-h0_mat(0:bs-1,0:bs-1,n))
				end do


				do n=0,nb
						do i=0,nb
							mult_mat(bs*n:bs*n+bs-1,bs*i:bs*i+bs-1)	= &
							matmul(g0_mat(0:bs-1,0:bs-1,n), v_mat(bs*n:bs*n+bs-1,bs*i:bs*i+bs-1))
						end do 
				end do

				psi_mat(0:n_basis) = 1.0_dp/sqrt(real(n_basis))


				gff(0:n_basis,0:n_basis,0) = mult_mat(0:n_basis,0:n_basis)
				do n=1,steps
					gff(0:n_basis,0:n_basis,n) = matmul(gff(0:n_basis,0:n_basis,n-1),mult_mat(0:n_basis,0:n_basis))
				end do

				gcc(0:n_basis,0) = psi_mat(0:n_basis)
				do n=1,steps
					gcc(0:n_basis,n) = matmul(gff(0:n_basis,0:n_basis,n-1),psi_mat(0:n_basis))
					taylor(n)=dot_product(gcc(0:n_basis,0),gcc(0:n_basis,n))
				end do



				call taylor_cfrac(taylor,steps,cf)

				cf_num = abs(evalcf(cf,steps,x))

! 				e = e + de
! 			end do

! 			print *, e, 'e', energy(0,0), 'emat', cf_num, 'num'

! 		print *, 'e mat----------'
! 		do n=0,bs-1
! 	        	print '(12f10.2)', e_mat(n,0:bs-1)
! 	    end do

			mine = 1/cf_num

		end function mcalc

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

