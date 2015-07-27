program exp_cf_rec

	use NumType
	use cf_approx
	use chebyshev
	implicit none

	!=========================== CONSTANTS ===========================

	integer,	parameter	::	dimensions	=	151,					&
								n_basis		=	dimensions-1,		&
								bs			=	10,					&
								nb 			=	(n_basis+1)/bs-1,	&
								steps 		=	ncf,				&
								nch			=	50

	real(dp),	parameter	::	mass		=	1.0_dp,				& 
								hbar		=	1.0_dp,				&
								omega_h		=	1.0_dp,				&
								omega_b		=	1.1_dp,				&
								x			=	1.0_dp

	!=========================== MATRICES ===========================

	real(dp)				::	eigen_mat(0:n_basis,0:n_basis)

	!=========================== PARAMETERS ==========================

	real, parameter :: set_dx = 0.001
	integer,	parameter	::  lwork=(n_basis+2)*n_basis,	&
								nmin=0,	&
								nmax=2,	&
								esteps = (nmax-nmin)/set_dx
	integer					::  n,i,j, block_size, info,num_points,m,int_num
	real(dp)				::  w_eigen(n_basis+1),work(lwork), offset
	real(dp)				::	ya,yb, smallthing, dz, zz, r, num, set_x,eps_n,ya2,deriv(0:esteps)
	integer					::  ind1(0:11759),ind2(0:11759)
	real(dp)				::  ind3(0:11759), &
								hmat_load(0:n_basis,0:n_basis)

	!=================================================================

!   compute eigenvalues using LAPACK
! 	eigen_mat(0:n_basis,0:n_basis) = 1.0_dp
! 			eigen_mat(0,0) = 1.0_dp
! 			eigen_mat(1,1) = 2.0_dp
! 			eigen_mat(2,2) = 3.0_dp
! 			eigen_mat(3,3) = 4.0_dp
! 			eigen_mat(4,4) = 5.0_dp
! 			eigen_mat(5,5) = 6.0_dp

	open(unit=9, file='pfaff_Hmatrix_N8_l6.5.dat')

    do i =1, 11759
        read(9,*) ind1(i),ind2(i),ind3(i)
        hmat_load(ind1(i)-1,ind2(i)-1) = ind3(i)
    end do
    close(9)

    eigen_mat(0:n_basis,0:n_basis) = hmat_load(0:n_basis,0:n_basis)


	call dsyev('V','U',n_basis+1,eigen_mat,n_basis+1,w_eigen,work,lwork,info)

	print *, 'LAPACK eigenvalues-------------------------'
	print '(10f10.2)', w_eigen
	print *, ' ====================================='

! 	do n = 0,n_basis
! 		print '(A10,100f10.2)', '  vector  ', eigen_mat(0:n_basis,n) 
! 	end do

!  stop

	!calculate the eigenvalues using chebyshev
	print *, 'calculated eigenvalues-------------------------'
! 	nmin=0
! 	nmax=6
! 	set_dx = 0.01
! 	esteps = (nmax-nmin)/set_dx
	do j=nmin,esteps,1
		int_num = 0
		ya=int_num+j*set_dx
		ya2=int_num+(j+1)*set_dx
! 		yb=j+1
! 		ya = 1
! 		yb = 10
! 		call chebyex(mcalc, nch, cheb, ya, yb)

! 		num_points = 100
! 		set_dx = (yb-ya)/num_points 
! 		do i = 1, num_points+1
! 			set_x = ya+(i-1)*set_dx
! 			write(1,*) set_x, cheby(set_x,cheb,nch,ya,yb)
! 		end do

! 		call chebyzero(nch, cheb, ya, yb, z0, iz0)

! 		print *, 'range=',ya,' to ', yb, z0(1:iz0)
		yb = mcalc(ya)
		write(3,*) ya, yb
! 		print *, ya, 'and the value is ', yb
		deriv(j) = mcalc(ya2)-yb
		write(4,*) ya, deriv(j)
		if (deriv(j) <= 0 .and. deriv(j-1) >= 0 .and. j > 0) then
			print '(10f10.2)', ya
		end if 
	end do

    
!     print *, hmat_load(1,1)
!     print *, hmat_load(1,2)
!     print *, hmat_load(1,3)
!     print *, hmat_load(1,4)
!     print *, hmat_load(1,5)
!     print *, hmat_load(2,2)
!     print *, hmat_load(3,3)

! 	ya = 0
! 	yb = 2
! 	set_dx = (yb-ya)/nch
! ! 	print *, set_dx
! 	call chebyex(mcalc,nch,cheb,ya,yb)
! ! 	print *, cheb
! 	do i=0,nch
! 	set_x = ya+i*set_dx
! ! 	print *, set_x
! 	write(4,*) set_x, cheb(i)
! 	end do



! 	ya = mcalc(2.53_dp)
! 	print *, '2.53 ', 'and the value is ', ya
! 	ya = mcalc(2.54_dp)
! 	print *,  '2.54 ', 'and the value is ', ya
! 	ya = mcalc(2.55_dp)
! 	print *,  '2.55 ', 'and the value is ', ya
! 	ya = mcalc(2.3_dp)
! 	print *,  '2.3 ', 'and the value is ', ya
! 	ya = mcalc(2.39_dp)
! 	print *,  '2.39 ', 'and the value is ', ya
! 	ya = mcalc(2.40_dp)
! 	print *,  '2.40 ', 'and the value is ', ya
! 	ya = mcalc(2.49_dp)
! 	print *,  '2.49 ', 'and the value is ', ya
! 	ya = mcalc(2.5_dp)
! 	print *,  '2.5 ', 'and the value is ', ya
! 	ya = mcalc(2.51_dp)
! 	print *,  '2.51 ', 'and the value is ', ya
! 	ya = mcalc(2.52_dp)
! 	print *,  '2.52 ', 'and the value is ', ya



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

		!========================================================================

		function mcalc(energy) result(res_sum)

			real(dp) :: energy, res_sum, cf_num, res
			real(dp) :: g0_mat(0:bs-1,0:bs-1,0:nb),			&
						mult_mat(0:n_basis,0:n_basis),		&
						taylor(0:ncf),						&
						cf(0:ncf),							&
						gff(0:n_basis,0:n_basis,0:steps),	&
						gcc(0:n_basis,0:steps),				&
						e_mat(0:bs-1,0:bs-1),				&
						x_mat(0:n_basis,0:n_basis),			&
						p_mat(0:n_basis,0:n_basis),			&
						h_mat(0:n_basis,0:n_basis),			&
						h0_mat(0:bs-1,0:bs-1,0:nb),			&
						v_mat(0:n_basis,0:n_basis),			&
						psi_mat(0:n_basis)


			g0_mat(0:bs-1,0:bs-1,0:nb)				=	0._dp
			mult_mat(0:n_basis,0:n_basis)			=	0._dp
			taylor(0:ncf)							=	0._dp
			cf(0:ncf)								=	0._dp
			gff(0:n_basis,0:n_basis,0:steps)		=	0._dp
			gcc(0:n_basis,0:steps)					=	0._dp
			e_mat(0:bs-1,0:bs-1)					=	0._dp
			x_mat(0:n_basis,0:n_basis)				=	0._dp
			p_mat(0:n_basis,0:n_basis)				=	0._dp
			h_mat(0:n_basis,0:n_basis)				=	0._dp
			h0_mat(0:bs-1,0:bs-1,0:nb)				=	0._dp
			v_mat(0:n_basis,0:n_basis)				=	0._dp
			psi_mat(0:n_basis)						=	0._dp

			!============================================================================

			!build X**2
! 			do n=0,n_basis
! 				x_mat(n,n) = (hbar/(2._dp*mass*omega_b))*(2._dp*n+1._dp)
! 				if (n>=2) then
! 					x_mat(n,n-2) = hbar/(2._dp*mass*omega_b)*sqrt(real(n*(n-1._dp)))
! 				end if
! 				if (n<=(n_basis-2)) then
! 					x_mat(n,n+2) = hbar/(2._dp*mass*omega_b)*sqrt(real((n+1._dp)*(n+2._dp)))
! 				end if
! 			end do

! 			!build P**2
! 			do n=0,n_basis
! 				p_mat(n,n) = (hbar*mass*omega_b)/2._dp*(2._dp*n+1._dp)
! 				if (n>=2) then
! 					p_mat(n,n-2) = -(hbar*mass*omega_b)/2._dp*sqrt(real(n*(n-1._dp)))
! 				end if
! 				if (n<=(n_basis-2)) then
! 					p_mat(n,n+2) = -(hbar*mass*omega_b)/2._dp*sqrt(real((n+1._dp)*(n+2._dp)))
! 				end if
! 			end do

! 			h_mat(0:n_basis,0:n_basis) = p_mat(0:n_basis,0:n_basis)/(2._dp*mass) + &
! 									mass*omega_h**2._dp*x_mat(0:n_basis,0:n_basis)/2._dp


! 			h_mat(0:n_basis,0:n_basis) = 1.0_dp
! 			h_mat(0,0) = 1.0_dp
! 			h_mat(1,1) = 2.0_dp
! 			h_mat(2,2) = 3.0_dp
! 			h_mat(3,3) = 4.0_dp
! 			h_mat(4,4) = 5.0_dp
! 			h_mat(5,5) = 6.0_dp

			h_mat(0:n_basis,0:n_basis) = hmat_load(0:n_basis,0:n_basis)

! 			print *, 'MAIN MATRIX'
! 			do n=0,n_basis
! 					print *, h_mat(n,0:n_basis)
! 			end do

			

			

			!============================================================================

			!build V & H_0 matrix

			v_mat(0:n_basis,0:n_basis) = h_mat(0:n_basis,0:n_basis)

			do n=0,nb
				h0_mat(0:bs-1,0:bs-1,n)=h_mat(bs*n:bs*n+bs-1,bs*n:bs*n+bs-1)
				v_mat(bs*n:bs*n+bs-1,bs*n:bs*n+bs-1) = &
				h_mat(bs*n:bs*n+bs-1,bs*n:bs*n+bs-1) - h0_mat(0:bs-1,0:bs-1,n)
			end do

! 			print *, 'V MATRIX'
! 			do n=0,n_basis
! 					print *, v_mat(n,0:n_basis)
! 			end do

! 			print *, 'H0 MATRIX'
! 			do j=0,nb
! 				print *, 'AND J IS ', j
! 			do n=0,bs-1
! 					print *, h0_mat(n,0:bs-1,j)
! 				end do
! 			end do

			!============================================================================

			!Define phi and energy

			psi_mat(0:n_basis) = 1.0_dp/sqrt(real(n_basis))

			do n=0,bs-1
				e_mat(n,n)	=  energy
			end do

! 			print *, 'E MATRIX'
! 			do n=0,bs-1
! 					print *, e_mat(n,0:bs-1)
! 			end do

			!============================================================================

			!Define G_0
			do n=0,nb
				g0_mat(0:bs-1,0:bs-1,n)	=e_mat(0:bs-1,0:bs-1)-h0_mat(0:bs-1,0:bs-1,n)
			end do


! 			print *, 'G0 MATRIX BEFORE!!!'
! 			do j=0,nb
! 				print *, 'AND J IS ', j
! 			do n=0,bs-1
! 					print *, g0_mat(n,0:bs-1,j)
! 				end do
! 			end do


			do n=0,nb
				g0_mat(0:bs-1,0:bs-1,n)	=inv(e_mat(0:bs-1,0:bs-1)-h0_mat(0:bs-1,0:bs-1,n))
			end do


			do n=0,nb
					do i=0,nb
						mult_mat(bs*n:bs*n+bs-1,bs*i:bs*i+bs-1)	= &
						matmul(g0_mat(0:bs-1,0:bs-1,n), v_mat(bs*n:bs*n+bs-1,bs*i:bs*i+bs-1))
					end do 
			end do

! 			print *, 'G0 MATRIX'
! 			do j=0,nb
! 				print *, 'AND J IS ', j
! 			do n=0,bs-1
! 					print *, g0_mat(n,0:bs-1,j)
! 				end do
! 			end do

! 			print *, '============================'


! 			print *, 'MULT MATRIX'
! 			do n=0,n_basis
! 					print *, mult_mat(n,0:n_basis)
! 			end do

			!============================================================================

			!Build Taylor series

			gff(0:n_basis,0:n_basis,0) = 1._dp
			gff(0:n_basis,0:n_basis,1) = mult_mat(0:n_basis,0:n_basis)
			do n=2,steps
				gff(0:n_basis,0:n_basis,n) = matmul(gff(0:n_basis,0:n_basis,n-1),mult_mat(0:n_basis,0:n_basis))
			end do

! 			print *, '1st MATRIX'
! 			do n=0,n_basis
! 					print *, gff(n,0:n_basis,0)
! 			end do

! 			print *, '2nd MATRIX'
! 			do n=0,n_basis
! 					print *, gff(n,0:n_basis,1)
! 			end do

! 			print *, '3rd MATRIX'
! 			do n=0,n_basis
! 					print *, gff(n,0:n_basis,2)
! 			end do

			gcc(0:n_basis,0) = psi_mat(0:n_basis)

			do n=1,steps
				gcc(0:n_basis,n) = matmul(gff(0:n_basis,0:n_basis,n),psi_mat(0:n_basis))
! 				 taylor(n) = dot_product(gcc(0:n_basis,0),gcc(0:n_basis,n))
! 				print *, ' n is ', n, ' value is ', taylor(n)
			end do

			res = 0.0_dp
			do n=0,steps
				taylor(n) = dot_product(gcc(0:n_basis,0),gcc(0:n_basis,n))
				res = res + taylor(n)
				call init_random_seed()
				call random_number(r)
! 				print *, 1+r*1e-15
				taylor(n) = taylor(n)*(1+r*1e-10)
! 				print *, ' n is ', n, ' value is ', taylor(n)
			end do

! 			print *, '1st MATRIX GCC'
! 			do n=0,n_basis
! 					print *, gcc(n,0)
! 			end do

! 			print *, '2nd MATRIX GCC'
! 			do n=0,n_basis
! 					print *, gcc(n,1)
! 			end do

! 			print *, '3rd MATRIX GCC'
! 			do n=0,n_basis
! 					print *, gcc(n,2)
! 			end do

! 			print *, '4th MATRIX GCC'
! 			do n=0,n_basis
! 					print *, gcc(n,3)
! 			end do

			!============================================================================

! 			do n=0,n_basis
! 				print '(12f10.2)', gcc(n,0)
! 			end do
! 			print *, energy, '-', res, 'res'

			call taylor_cfrac(taylor,steps,cf)

			cf_num = evalcf(cf,steps,x)

! 			print *, energy, cf_num, '< - final energy taylor', res, ' <- final energy res'

! 			res_sum = abs(1._dp/cf_num)

			res_sum = abs(cf_num)

		end function mcalc

		function horner(f,n,x) result(y)
			implicit none

			real(dp), dimension(0:ncf) :: f 
			integer :: n, i
			real(dp) :: x, y

			y = f(n)
			do i = n-1, 0, -1
				y = f(i) + x*y
			end do

		end function horner

		function func(x) result(f)

		use numtype , only : dp 
		implicit none

		real(dp) :: x, f

		f = -x*exp(-(x-2)**2/(2*x**2))

		end function func


end program exp_cf_rec
