

module setup

	use numtype
	implicit none
	integer, parameter :: nsp = 1025, npar = 8
	integer :: yy(nsp), minsp, maxsp, iprint, ifcal
	
		
end module setup



program gamma_spectrum

	use setup
	implicit none
	integer :: i, ii, itmin, itmax
	real(dp) :: fstart, xstart(npar), stepi, epsf
	real(dp), external :: least2
	

	open(unit=2, file='na22.dat')
	open(unit=3, file='na22.d')
	do i = 1, nsp
		read(2,*) ii
		write(3,*) i, ii
		yy(i) = ii
	end do
	close(2)
	close(3)


	minsp = 230
	maxsp = 900
	xstart(1:npar) = (/ -0.01, 120.0, 3000.0, 300.0, 50.0, 400.0, 750.0, 40.0  /)
	
	ifcal = 0
	iprint = 1
	fstart = least2(xstart(1:npar))
	
	itmin = 100
	itmax = 1000
	epsf = 0.001_dp
	stepi = 0.1_dp
	iprint = 0
	
	call downhill(npar, least2, xstart, fstart, stepi, epsf, itmin, itmax)

	iprint = 2
	fstart = least2(xstart(1:npar))


end program gamma_spectrum


function least2(par) result(ss)


	use setup
	implicit none
	real(dp) :: par(npar), a, b, y1, x1, sig1, y2, x2, sig2, ss, fi
	integer :: i


	ifcal = ifcal +1
	a = par(1); b = par(2);
	y1 = par(3);	x1 = par(4);	sig1 = par(5)
	y2 = par(6);	x2 = par(7);	sig2 = par(8)

	ss = 0._dp
	
	do i = minsp, maxsp
		fi = a*i + b + y1*exp(-((i-x1)/sig1)**2) + y2*exp(-((i-x2)/sig2)**2)
		ss = ss + (fi-yy(i))**2/sqrt(yy(i)+1.0)
	end do
	ss = ss/(maxsp-minsp)
	print '(i5,8f12.4, f20.5)', ifcal, par(1:npar), ss
	write(8,*) ifcal, ss
	
	if ( iprint /= 0 ) then
		do i = minsp, maxsp
			fi = a*i + b + y1*exp(-((i-x1)/sig1)**2) + y2*exp(-((i-x2)/sig2)**2)
			write(unit=iprint,fmt='(i4,i10,f15.2)') i,yy(i),fi
		end do

	end if



end function least2

