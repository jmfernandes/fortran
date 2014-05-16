
module setup

	use NumType
	implicit none
	integer,	parameter	::	n_basis=50,lwork=2*n_basis+1,n_eq=3
	real(dp),	parameter	::	mass=1.0_dp, hbar=1.0_dp,&
							  	omega_h=10/10._dp,omega_b=4/5._dp, &
                                hbar2 = hbar**2
    real(dp), allocatable, dimension(:,:) :: wf
    real(dp) :: energy

end module setup

program matrix

	use setup
	implicit none

	complex(dp) ::	x_mat(0:n_basis,0:n_basis+1), p_mat(0:n_basis,0:n_basis+1), &
					x2_mat(0:n_basis,0:n_basis), p2_mat(0:n_basis,0:n_basis), &
					h_mat(0:n_basis,0:n_basis), &
                    work(lwork),e(0:n_basis,0:n_basis)
    integer :: n,m, info, i, j,imax, parity, imin, iz
    real(dp) :: rwork(3*(n_basis-2)), w_eigen(n_basis+1)
    real(dp) :: x,xmax,dx, y(n_eq)


	x_mat = 0._dp
    p_mat = 0._dp

	do n=0,n_basis-1
		x_mat(n,n+1) = sqrt(hbar/(2*mass*omega_b))*sqrt(n+1._dp)
		x_mat(n+1,n) = x_mat(n,n+1)		
	end do
	x_mat(n_basis,n_basis+1) = sqrt(hbar/(2*mass*omega_b))*sqrt(n_basis+1._dp) !add last point


    do n=0,n_basis-1
        p_mat(n,n+1) = -iic*sqrt(mass*hbar*omega_b/2)*sqrt(n+1._dp)
        p_mat(n+1,n) = -p_mat(n,n+1)     
    end do
    p_mat(n_basis,n_basis+1) = -sqrt(mass*hbar*omega_b/2)*sqrt(n_basis+1._dp) !add last point



    h_mat(0:n_basis,0:n_basis) = 1/(2*mass) * &
        matmul(p_mat(0:n_basis,0:n_basis+1),conjg(transpose(p_mat(0:n_basis,0:n_basis+1)))) + &
        1/sqrt(1 + &
        matmul(x_mat(0:n_basis,0:n_basis+1),transpose(x_mat(0:n_basis,0:n_basis+1))))


    e(1:n_basis,1:n_basis) = h_mat(1:n_basis,1:n_basis)/10

    call zheev('n','u',n_basis+1,e,n_basis+1,w_eigen,work,lwork,rwork,info)


    print *, '------- Orthogonality ------- '
        
        do j = 1,50
            do i = 1,50
                print *, i, j, dot_product(e(1:50,i) ,e(1:50,j) ) !determines if matrix is orthogonal
            end do
        end do


    print *, '===First Ten Energy States========'
    do i=1,10
    print *, w_eigen(i)
    end do

    print *, '===First ten entries for the Ten Eigenvectors========='
    do i=1,10
    print '(10f8.2)', dble(e(1:10,i)) 
    end do

    !Plot Results

    xmax=40
    dx=0.01
    imax = abs(xmax/dx)

    allocate(wf(0:imax,2))

    !only plot the first 10 intead of all 50
!     iz = size(w_eigen)
    iz = 10

    !plot the wavefunctions for the different eigenvalues
    do j=1,iz

    energy = w_eigen(j)

    x=xmax

    y(1) = 0.00001
    y(2) = -0.00001
    y(3) = 0

    i = imax+1
    do while(x > 0) 
        i = i -1
        wf(i,1) = x
        wf(i,2) = y(1)
        call rk4step(x,y,-dx)
    end do

    imin=i
    if( abs(y(1)) > abs(y(2)) ) then
        parity = 1
    else 
        parity = -1
    end if

    wf(0:imax,2) = wf(0:imax,2)/sqrt(2*y(3)) 

    
        do i = imax, imin, -1
            write(unit=20+j,fmt='(2f15.5)') wf(i,1), &
                wf(i,2)
        end do

        do i = imin, imax
            write(unit=20+j,fmt='(2f15.5)') -wf(i,1), &
                parity*wf(i,2)
        end do
    end do


end program matrix