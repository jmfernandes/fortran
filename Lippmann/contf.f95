
program exp_cf_rec

    use NumType
    implicit none
    integer 					::	n
    integer,    parameter       ::  n_basis=5
    real(dp),   parameter		::  mass=1.0_dp, hbar=2.0_dp,&
                                	omega_h=1._dp,omega_b=1._dp
    real(dp) ::     x_mat(0:n_basis,0:n_basis),     &
                    p_mat(0:n_basis,0:n_basis),     &
                    h_mat(0:n_basis,0:n_basis)



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    x_mat   = 0._dp
    p_mat   = 0._dp
    h_mat   = 0._dp



    do n=0,n_basis
        x_mat(n,n) = (hbar/(2*mass*omega_b))*(2*n+1)
        if (n>=2) then
            x_mat(n,n-2) = hbar/(2*mass*omega_b)*sqrt(real(n*(n-1)))
        end if
        if (n<=(n_basis-2)) then
            x_mat(n,n+2) = hbar/(2*mass*omega_b)*sqrt(real((n+1)*(n+2)))
        end if
    end do

    do n=0,n_basis
        p_mat(n,n) = (hbar*mass*omega_b)/2*(2*n+1)
        if (n>=2) then
            p_mat(n,n-2) = -(hbar*mass*omega_b)/2*sqrt(real(n*(n-1)))
        end if
        if (n<=(n_basis-2)) then
            p_mat(n,n+2) = -(hbar*mass*omega_b)/2*sqrt(real((n+1)*(n+2)))
        end if
    end do

    h_mat(0:n_basis,0:n_basis) = p_mat(0:n_basis,0:n_basis)/(2*mass)+mass*omega_h**2*x_mat(0:n_basis,0:n_basis)/2



    print *, 'x^2 matrix-------------------------'
    do n=0,n_basis
        print '(12f10.2)', x_mat(n,0:n_basis)
    end do

    print *, 'p^2 matrix-------------------------'
    do n=0,n_basis
        print '(12f10.2)', p_mat(n,0:n_basis)
    end do

    print *, 'hamiltonian matrix-------------------------'
    do n=0,n_basis
        print '(12f10.2)', h_mat(n,0:n_basis)
    end do

end program exp_cf_rec

