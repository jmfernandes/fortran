
module setup

	use NumType
	implicit none
	integer,	parameter	::	n_basis=10,lwork=2*n_basis+1
	real(dp),	parameter	::	mass=1.0_dp, hbar=1.0_dp,&
							  	omega_h=0.5_dp,omega_b=1._dp

end module setup

program matrix

	use setup
	implicit none

	complex(dp) ::	x_mat(0:n_basis,0:n_basis+1), p_mat(0:n_basis,0:n_basis+1), &
					x2_mat(0:n_basis,0:n_basis), p2_mat(0:n_basis,0:n_basis), &
					h_mat(0:n_basis,0:n_basis), &
                    work(lwork)
    integer :: n,m, info
    real(dp) :: rwork(3*(n_basis-2)), w_eigen(n_basis+1)


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


    do n=0,n_basis
        print '(12f10.5)', x_mat(n,0:n_basis+1)
    end do

    print *, '----------'

    do n=0,n_basis
        print '(12f10.5)', p_mat(n,0:n_basis+1)
    end do


     print *, '==========================================================='


    h_mat(0:n_basis,0:n_basis) = 1/(2*mass) * &
        matmul(p_mat(0:n_basis,0:n_basis+1),conjg(transpose(p_mat(0:n_basis,0:n_basis+1)))) + &
        mass*omega_h**2/2 * &
        matmul(x_mat(0:n_basis,0:n_basis+1),transpose(x_mat(0:n_basis,0:n_basis+1)))

    do n=0,n_basis
        print '(12f10.5)', h_mat(n,0:n_basis)
    end do

    call zheev('n','u',n_basis+1,h_mat,n_basis+1,w_eigen,work,lwork,rwork,info)


    print *, '----------'
    print *, info

    print *, w_eigen(1:n_basis+1)

end program matrix