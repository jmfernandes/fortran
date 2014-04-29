
program integral

	use numtype
	use integr
	implicit none

	real(dp) :: a,b,c,d,res,eps,ifail,scale
	real(dp), dimension(maxint) :: w_legendre,x_legendre, w_cc, x_cc, &
								   w_lag, x_lag
	integer :: nint, itype, i,n_legendre, n_cc, n_lag

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


	!GAUSSIAN IS BETTER THAN ROMBINT

	print *, '==========result========='
	res = 0._dp
	do i = 1,n_legendre
		res = res + runge(x_legendre(i))*w_legendre(i)
! 		print *, i, x_legendre(i), w_legendre(i)
	end do
	print *, 'gauss', res, 2*atan(b)


	n_cc = 100

	!CC IS ALSO GOOD

	call cc11(n_cc,a,b,x_cc,w_cc) 

	res = 0._dp
	do i = 1,n_cc
		res = res + runge(x_cc(i))*w_cc(i)
	end do

	print *, 'CC', res, 2*atan(b)

	print *, '-------------------NEW FUNCTION---------------------------'

	! function is (x*2exp(-x))/(sqrt(1+x)) incorporate the x**2exp(-x) into weight

	print *, 'Gauss normal weights'

	itype = 3
	a = 0._dp
	b = 1._dp
	c = 2._dp
	d = 0._dp 
	ifail = 0._dp
	n_lag = 50._dp

	call d01bcf(itype,a,b,c,d,n_lag,w_lag,x_lag,ifail)

	res = 0._dp
	do i = 1,n_lag
		res = res + sq(x_lag(i))*w_lag(i)
	end do
	print *, 'gauss', res

	print *, '================================================================'
	print *, 'Gauss adjusted weights'

	itype = -3
	a = 0._dp
	b = 1._dp
	c = 2._dp
	d = 0._dp 
	ifail = 0._dp
	n_lag = 50._dp

	call d01bcf(itype,a,b,c,d,n_lag,w_lag,x_lag,ifail)

	res = 0._dp
	do i = 1,n_lag
		res = res + sqa(x_lag(i))*w_lag(i)
	end do
	print *, 'gauss', res

	print *, '================================================================'
	print *, 'NEW STUFF'

	n_cc = 120
	scale = 1._dp

	call cc0inf(n_cc,scale,x_cc,w_cc)

	res = 0._dp
	do i = 1,n_cc
		res = res + sqa(x_cc(i))*w_cc(i)
	end do

	print *, ' CC', res


	contains

		function runge(x) result(fx)

			real(dp) :: x,fx


			fx = 1/(1+x**2)

		end function runge

		function sq(x) result(fx)

			real(dp) :: x,fx


			fx = 1/sqrt(1+x)

		end function sq

		function sqa(x) result(fx)

			real(dp) :: x,fx


			fx = (x**2 * exp(-x))/sqrt(1+x)

		end function sqa
	

end program integral