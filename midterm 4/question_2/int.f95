
program integral

	use numtype
	use integr
	implicit none

	real(dp) :: a,b,c,d,res,eps,ifail,scale
	real(dp), dimension(maxint) :: w_legendre,x_legendre, w_cc, x_cc, &
								   w_lag, x_lag, w_new, x_new
	integer :: nint, itype, i,n_legendre, n_cc, n_lag, n_new

! 	print *, 'int_-5^5 1/(1+x^2) dx'

! 	a = -5._dp
! 	b = 5._dp
! 	eps = 1.e-10_dp
! 	nint = 20

! 	call rombint(a,b,runge,res,nint,eps)

! 	print *, ' romberg result :', res, nint,2**nint, 2*atan(b)

! 	itype = 0
! 	a = -5._dp
! 	b = 5._dp
! 	c = 0._dp
! 	d = 0._dp
! 	ifail = 0._dp
! 	n_legendre = 50._dp

! 	call d01bcf(itype,a,b,c,d,n_legendre,w_legendre,x_legendre,ifail)

! 	print *,'ifail', ifail

! ! 	do i = 1,n_legendre
! ! 		print *, i, x_legendre(i), w_legendre(i)
! ! 	end do


! 	!GAUSSIAN IS BETTER THAN ROMBINT

! 	print *, '==========result========='
! 	res = 0._dp
! 	do i = 1,n_legendre
! 		res = res + runge(x_legendre(i))*w_legendre(i)
! ! 		print *, i, x_legendre(i), w_legendre(i)
! 	end do
! 	print *, 'gauss', res, 2*atan(b)


! 	n_cc = 100

! 	!CC IS ALSO GOOD

! 	call cc11(n_cc,a,b,x_cc,w_cc) 

! 	res = 0._dp
! 	do i = 1,n_cc
! 		res = res + runge(x_cc(i))*w_cc(i)
! 	end do

! 	print *, 'CC', res, 2*atan(b)

! 	print *, '-------------------NEW FUNCTION---------------------------'

! 	! function is (x*2exp(-x))/(sqrt(1+x)) incorporate the x**2exp(-x) into weight

! 	print *, 'Gauss normal weights'

! 	itype = 3
! 	a = 0._dp
! 	b = 1._dp
! 	c = 2._dp
! 	d = 0._dp 
! 	ifail = 0._dp
! 	n_lag = 50._dp

! 	call d01bcf(itype,a,b,c,d,n_lag,w_lag,x_lag,ifail)

! 	res = 0._dp
! 	do i = 1,n_lag
! 		res = res + sq(x_lag(i))*w_lag(i)
! 		write(24,*) x_lag(i), sq(x_lag(i))
! 	end do
! 	print *, 'gauss', res

! 	print *, '================================================================'
! 	print *, 'Gauss adjusted weights'

! 	itype = -3
! 	a = 0._dp
! 	b = 1._dp
! 	c = 2._dp
! 	d = 0._dp 
! 	ifail = 0._dp
! 	n_lag = 50._dp

! 	call d01bcf(itype,a,b,c,d,n_lag,w_lag,x_lag,ifail)

! 	res = 0._dp
! 	do i = 1,n_lag
! 		res = res + sqa(x_lag(i))*w_lag(i)
! 		write(34,*) x_lag(i), sqa(x_lag(i))*w_lag(i)
! 	end do
! 	print *, 'gauss', res

! 	print *, '================================================================'
! 	print *, 'NEW STUFF'

! 	n_cc = 120
! 	scale = 1._dp

! 	call cc0inf(n_cc,scale,x_cc,w_cc)

! 	res = 0._dp
! 	do i = 1,n_cc
! 		res = res + sqa(x_cc(i))*w_cc(i)
! 	end do

! 	print *, ' CC', res

! 	print *, '================================================================'
! 	print *, 'My Stuff'

! 	itype = 4
! 	a = 0._dp
! 	b = 1._dp
! 	c = 1._dp
! 	d = 0._dp 
! 	ifail = 0._dp
! 	n_new = 50._dp

! 	call d01bcf(itype,a,b,c,d,n_new,w_new,x_new,ifail)

! 	res = 0._dp
! 	do i = 1,n_new
! 		res = res + sqas(x_new(i))*w_new(i)
! 		write(14,*) x_new(i), sqas(x_new(i))
! 	end do
! 	print *, 'gauss', res

! 	print *, '================================================================'
! 	print *, 'My Stuff2'


! 	itype = 0
! 	a = 1._dp
! 	b = 5._dp
! 	c = 3._dp
! 	d = 0._dp 
! 	ifail = 0._dp
! 	n_new = 50._dp

! 	call d01bcf(itype,a,b,c,d,n_new,w_new,x_new,ifail)

! 	res = 0._dp
! 	do i = 1,n_new
! 		res = res + sqa1(x_new(i))*w_new(i)
! 	end do
! 	print *, 'gauss', res

! 	print *, '================================================================'
! 	print *, 'My Stuff3'


! 	itype = 3
! 	a = 0._dp
! 	b = 2._dp
! 	c = 3._dp
! 	d = 0._dp 
! 	ifail = 0._dp
! 	n_new = 50._dp

! 	call d01bcf(itype,a,b,c,d,n_new,w_new,x_new,ifail)

! 	res = 0._dp
! 	do i = 1,n_new
! 		res = res + sqa2(x_new(i))*w_new(i)
! 	end do
! 	print *, 'gauss', res


! 	n_cc = 120
! 	scale = 1._dp

! 	call cc0inf(n_cc,scale,x_cc,w_cc)

! 	res = 0._dp
! 	do i = 1,n_cc
! 		res = res + sqa2a(x_cc(i))*w_cc(i)
! 	end do

! 	print *, ' CC', res

! 	print *, '================================================================'
! 	print *, 'My Stuff4'


! 	itype = 0
! 	a = 0.1_dp
! 	b = 100._dp
! 	c = 3._dp
! 	d = 0._dp 
! 	ifail = 0._dp
! 	n_new = 50._dp

! 	call d01bcf(itype,a,b,c,d,n_new,w_new,x_new,ifail)

! 	res = 0._dp
! 	do i = 1,n_new
! 		res = res + sqa3(x_new(i))*w_new(i)
! 		write(44,*) x_new(i),sqa3(x_new(i))*w_new(i)
! 	end do
! 	print *, 'gauss', res


! 	n_cc = 120
! 	scale = 1._dp

! 	call cc0inf(n_cc,scale,x_cc,w_cc)

! 	res = 0._dp
! 	do i = 1,n_cc
! 		res = res + sqa3(x_cc(i))*w_cc(i)
! 	end do
! 	print *, 'CC', res


! 	eps = 1.e-10_dp
! 	nint = 20

! 	res = 0._dp
! 	call rombint(a,b,sqa3,res,nint,eps)

! 	print *, ' romberg result :', res

	print *, '================================================================'
	print *, 'Midterm4'

	print *, 'first integration'

	itype = 1
	a = 0._dp
	b = 1._dp
	c = 3/2._dp
	d = 1/2._dp 
	ifail = 0._dp
	n_new = 50._dp

	call d01bcf(itype,a,b,c,d,n_new,w_new,x_new,ifail)

	res = 0._dp
	do i = 1,n_new
		res = res + mid1(x_new(i))*w_new(i)
! 		write(44,*) x_new(i),sqa3(x_new(i))*w_new(i)
	end do
	print *, 'gauss my  value =', res
	print *, 'true value=pi/16=', pi/16
	print *, 'accurate to 15 digits!!!!'

	print *, '================================================================'
	print *, 'Midterm4'

	print *, 'second integration'

	itype = 3
	a = 0._dp
	b = 1._dp
	c = -1/2._dp
	d = 0._dp 
	ifail = 0._dp
	n_new = 50._dp

	call d01bcf(itype,a,b,c,d,n_new,w_new,x_new,ifail)

	res = 0._dp
	do i = 1,n_new
		res = res + mid1(x_new(i))*w_new(i)
! 		write(44,*) x_new(i),sqa3(x_new(i))*w_new(i)
	end do

	print *, '  gauss  my  value =', res
	print *, 'true value=sqrt(pi)=', sqrt(pi)
	print *, 'accurate to 14 digits!!!!'
	


	contains

! 		function runge(x) result(fx)

! 			real(dp) :: x,fx


! 			fx = 1/(1+x**2)

! 		end function runge

! 		function sq(x) result(fx)

! 			real(dp) :: x,fx


! 			fx = 1/sqrt(1+x)

! 		end function sq

! 		function sqa(x) result(fx)

! 			real(dp) :: x,fx


! 			fx = (x**2 * exp(-x))/sqrt(1+x)

! 		end function sqa

! 		function sqas(x) result(fx)

! 			real(dp) :: x,fx


! 			fx = x**2

! 		end function sqas

! 		function sqa1(x) result(fx)

! 			real(dp) :: x,fx


! 			fx = x**2+(1/(x**2-x+5))+cos(x)

! 		end function sqa1

! 		function sqa2(x) result(fx)

! 			real(dp) :: x,fx


! 			fx = 1/x

! 		end function sqa2

! 		function sqa2a(x) result(fx)

! 			real(dp) :: x,fx


! 			fx = abs(x)**3*exp(-2*x)*(1/x)

! 		end function sqa2a

! 		function sqa3(x) result(fx)

! 			real(dp) :: x,fx


! 			fx = x**3/(exp(x)-1)

! 		end function sqa3

		function mid1(x) result(fx)

			real(dp) :: x,fx


			fx = 1

		end function mid1
	

end program integral