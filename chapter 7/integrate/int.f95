
program integral

	use numtype
	use integr
	implicit none

	real(dp) :: a,b,c,d,res,eps,ifail
	real(dp), dimension(maxint) :: w_legendre,x_legendre
	integer :: nint, itype, i,n_legendre

	print *, 'int_-5^5 1/(1+x^2) dx'

	a = -5._dp
	b = 5._dp
	eps = 1.e-10_dp
	nint = 20

	call rombint(a,b,runge,res,nint,eps)

	print *, ' romberg result :', res, nint,2**nint, 2*atan(b)

	itype = 0
	a = -5._dp
	b = 5._dp
	c = 0._dp
	d = 0._dp
	ifail = 0._dp
	n_legendre = 50._dp

	call d01bcf(itype,a,b,c,d,n_legendre,w_legendre,x_legendre,ifail)

	print *,'ifail', ifail

! 	do i = 1,n_legendre
! 		print *, i, x_legendre(i), w_legendre(i)
! 	end do

	print *, '==========result========='
	res = 0._dp
	do i = 1,n_legendre
		res = res + runge(x_legendre(i))*w_legendre(i)
! 		print *, i, x_legendre(i), w_legendre(i)
	end do
	print *, 'gauss', res, 2*atan(b)

	contains

		function runge(x) result(fx)

			real(dp) :: x,fx


			fx = 1/(1+x**2)

		end function runge
	

end program integral