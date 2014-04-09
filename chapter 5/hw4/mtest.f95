
module setup

	use NumType
	implicit none
	integer,	parameter	::	ndim=10,lwork=5*ndim
	real(dp),	parameter	::	mass=1.0_dp, hbar=1.0_dp,&
							  	omega1=0.5_dp,omega2=1._dp

end module setup

program matrix

	use setup
	implicit none

	complex(dp), dimension(ndim,ndim)	:: A,B,C,D,E,H1,H2,H3,H4,H5
	real(dp),	 dimension(ndim)		:: w
	integer 							:: i, nn, info
	complex(dp) 						:: work(lwork)
	real(dp) 							:: rwork(lwork)

	A(1:10,1:10) = reshape((/	zero, sqrt(1*one), zero, zero, zero,		&
							zero, zero, zero, zero, zero,					&
							sqrt(1*one), zero, sqrt(2*one), zero, zero,		&
							zero, zero, zero, zero, zero,					&
							zero, sqrt(2*one), zero, sqrt(3*one), zero,		&
							zero, zero, zero, zero, zero,					&
							zero, zero, sqrt(3*one), zero, sqrt(4*one),		&
							zero, zero, zero, zero, zero,					&
							zero, zero, zero, sqrt(4*one), zero,			&
							sqrt(5*one), zero, zero, zero, zero,			&
							zero, zero, zero, zero, sqrt(5*one),			&
							zero, sqrt(6*one), zero, zero, zero,			&
							zero, zero, zero, zero, zero,					&
							sqrt(6*one), zero, sqrt(7*one), zero, zero,		&
							zero, zero, zero, zero, zero,					&
							zero, sqrt(7*one), zero, sqrt(8*one), zero,		&
							zero, zero, zero, zero, zero,					&
							zero, zero, sqrt(8*one), zero, sqrt(9*one),		&
							zero, zero, zero, zero, zero,					&
							zero, zero, zero, sqrt(9*one), zero 			&
						/), 												&
						(/10,10/))

	C(1:10,1:10) = reshape((/	zero, -sqrt(1*one), zero, zero, zero,		&
							zero, zero, zero, zero, zero,					&
							sqrt(1*one), zero, -sqrt(2*one), zero, zero,	&
							zero, zero, zero, zero, zero,					&
							zero, sqrt(2*one), zero, -sqrt(3*one), zero,	&
							zero, zero, zero, zero, zero,					&
							zero, zero, sqrt(3*one), zero, -sqrt(4*one),	&
							zero, zero, zero, zero, zero,					&
							zero, zero, zero, sqrt(4*one), zero,			&
							-sqrt(5*one), zero, zero, zero, zero,			&
							zero, zero, zero, zero, sqrt(5*one),			&
							zero, -sqrt(6*one), zero, zero, zero,			&
							zero, zero, zero, zero, zero,					&
							sqrt(6*one), zero, -sqrt(7*one), zero, zero,	&
							zero, zero, zero, zero, zero,					&
							zero, sqrt(7*one), zero, -sqrt(8*one), zero,	&
							zero, zero, zero, zero, zero,					&
							zero, zero, sqrt(8*one), zero, -sqrt(9*one),	&
							zero, zero, zero, zero, zero,					&
							zero, zero, zero, sqrt(9*one), zero 			&
						/), 												&
						(/10,10/))

	H1(1:10,1:10) = reshape((/ 1/2._dp*one, zero, zero, zero, zero,			&
							zero, zero, zero, zero, zero,					&
							zero, 3/2._dp*one, zero, zero, zero,			&
							zero, zero, zero, zero, zero,					&
							zero, zero, 5/2._dp*one, zero, zero,			&
							zero, zero, zero, zero, zero,					&
							zero, zero, zero, 7/2._dp*one, zero,			&
							zero, zero, zero, zero, zero,					&
							zero, zero, zero, zero, 9/2._dp*one,			&
							zero, zero, zero, zero, zero,					&
							zero, zero, zero, zero, zero,					&
							11/2._dp*one, zero, zero, zero, zero,			&
							zero, zero, zero, zero, zero,					&
							zero, 13/2._dp*one, zero, zero, zero,			&
							zero, zero, zero, zero, zero,					&
							zero, zero, 15/2._dp*one, zero, zero,			&
							zero, zero, zero, zero, zero,					&
							zero, zero, zero, 17/2._dp*one, zero,			&
							zero, zero, zero, zero, zero,					&
							zero, zero, zero, zero, 19/2._dp*one 			&
						/), 												&
						(/10,10/))

	H2(1:10,1:10) = reshape((/ 1/4._dp*one, zero, zero, zero, zero,			&
							zero, zero, zero, zero, zero,					&
							zero, 3/4._dp*one, zero, zero, zero,			&
							zero, zero, zero, zero, zero,					&
							zero, zero, 5/4._dp*one, zero, zero,			&
							zero, zero, zero, zero, zero,					&
							zero, zero, zero, 7/4._dp*one, zero,			&
							zero, zero, zero, zero, zero,					&
							zero, zero, zero, zero, 9/4._dp*one,			&
							zero, zero, zero, zero, zero,					&
							zero, zero, zero, zero, zero,					&
							11/4._dp*one, zero, zero, zero, zero,			&
							zero, zero, zero, zero, zero,					&
							zero, 13/4._dp*one, zero, zero, zero,			&
							zero, zero, zero, zero, zero,					&
							zero, zero, 15/4._dp*one, zero, zero,			&
							zero, zero, zero, zero, zero,					&
							zero, zero, zero, 17/4._dp*one, zero,			&
							zero, zero, zero, zero, zero,					&
							zero, zero, zero, zero, 19/4._dp*one 			&
						/), 												&
						(/10,10/))

	A(1:10,1:10)  = sqrt(hbar/(2*mass*omega1))*A(1:10,1:10)
	B(1:10,1:10)  = sqrt(hbar/(2*mass*omega2))*A(1:10,1:10)

	C(1:10,1:10)  = iic*sqrt((mass*omega1*hbar)/2)*C(1:10,1:10) 
	D(1:10,1:10)  = iic*sqrt((mass*omega2*hbar)/2)*C(1:10,1:10) 

	H3(1:10,1:10) = H1 + H2 + 1/2._dp*mass*(omega2-omega1)**2*(B-A)**2

	H4(1:10,1:10) = (C*C)/(2*mass) + mass*omega1**2*(A*A)/2
	H5(1:10,1:10) = (D*D)/(2*mass) + mass*omega2**2*(B*B)/2

	nn = 10
	info = 0
	E(1:nn,1:nn) = H3(1:nn,1:nn)

	print '(20f6.2)', A(1:10,1:10)
	print *, '__________'
	print '(20f6.2)', C(1:10,1:10)

	call zheev('v','u',nn,E,ndim,w,work,lwork,rwork,info)

	print *, 'info=', info

	do i = 1,10
		print '(a,f15.8,a,20f6.0)','eigenvalues',w(i),'  vector  ', dble(e(1:nn,i)) 
	end do

	print *, "wa-bam"

end program matrix