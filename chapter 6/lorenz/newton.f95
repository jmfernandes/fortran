!this program solves a set of linear equations
module setup

	use numtype
	use matrix
	implicit none
	integer, parameter 	:: nv = 3, np = 3
	real(dp) 			:: a(1:np)

end module setup


program newton

	use setup
	implicit none

	real(dp) 			:: x(nv), dx(nv), f(nv), jacobi(nv,nv), ff
	real(dp), parameter :: eps = 1.e-10_dp
	integer				:: maxstep, i

	a(1:np) = (/ 10._dp, 28._dp, 8._dp/3 /)
	x(1:nv) = (/ 40._dp, 60._dp, 50._dp /)

	maxstep = 15

	do i = 1, maxstep

		call func(nv,x,f,jacobi,ff)
		print '(3f12.4, 3x, 3e12.4, 3x, e12.4)', &
				x(1:nv), f(1:nv), ff

		if (ff <= eps ) exit

		call dgei(nv,jacobi,nv) 

		dx(1:nv) = matmul( jacobi(1:nv,1:nv), -f(1:nv))

		x(1:nv) = x(1:nv)+dx(1:nv)

	end do


end program newton

subroutine func(n,x,f,jmat,ff)

	use setup
	implicit none
	real(dp) :: x(n), f(n), jmat(n,n), ff
	integer :: n

	f(1) = a(1)*(x(2)-x(1))
	f(2) = x(1)*(a(2)-x(3))-x(2)
	f(3) = x(1)*x(2) - a(3)*x(3)

	ff = sqrt(dot_product(f(1:n),f(1:n)))

	jmat(1,1) = -a(1)
	jmat(2,1) = a(2) - x(3)
	jmat(3,1) = x(2)

	jmat(1,2) = a(1)
	jmat(2,2) = -1
	jmat(3,2) = x(1)

	jmat(1,3) = 0
	jmat(2,3) = -x(1)
	jmat(3,3) = -a(3)

end subroutine func