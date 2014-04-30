
module setup

	use NumType
	implicit none
	integer,	parameter	::	n_basis=5,lwork=2*n_basis+1
	real(dp),	parameter	::	mass=1.0_dp, hbar=1.0_dp,&
							  	omega_h=2._dp,omega_b=1._dp

end module setup

program matrix

	use setup
	implicit none

	complex(dp) ::	x_mat(0:n_basis,0:n_basis+1), p_mat(0:n_basis,0:n_basis+1), &
					x2_mat(0:n_basis,0:n_basis), p2_mat(0:n_basis,0:n_basis), &
					h_mat(0:n_basis,0:n_basis),f(0:n_basis,0:n_basis), &
                    work(lwork),c(0:n_basis,0:n_basis),d(0:n_basis,0:n_basis), &
                    a(0:n_basis,0:n_basis)
    integer :: n,m, info, i, j,ipiv(n_basis)
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


   ! do n=0,n_basis
!        print '(12f10.5)', x_mat(n,0:n_basis+1)
!    end do
!
!    print *, '----------'
!
!    do n=0,n_basis
!        print '(12f10.5)', p_mat(n,0:n_basis+1)
!    end do


     print *, '==========================================================='


    h_mat(0:n_basis,0:n_basis) = 1/(2*mass) * &
        matmul(p_mat(0:n_basis,0:n_basis+1),conjg(transpose(p_mat(0:n_basis,0:n_basis+1)))) + &
        mass*omega_h**2/2 * &
        matmul(x_mat(0:n_basis,0:n_basis+1),transpose(x_mat(0:n_basis,0:n_basis+1)))
        
    h_mat(0,5) = h_mat(0,5) + (hbar*omega_h/2) !add contribution
    h_mat(5,0) = h_mat(5,0) + (hbar*omega_h/2)
    
    print *, '---------------hamiltonian matrix---------------'
    do n=0,n_basis
        print '(12f10.5)', h_mat(n,0:n_basis)
    end do

    a(1:n_basis,1:n_basis) = h_mat(1:n_basis,1:n_basis)

    call zheev('v','u',n_basis+1,h_mat,n_basis+1,w_eigen,work,lwork,rwork,info)


    print *, '----------'
    print *, info

    print *, 'eigenvalue', w_eigen(1:n_basis+1)
    
    do i = 1,n_basis
		print '(a,20f10.2)', ' vector  ', dble(h_mat(0:n_basis,i)) 
	end do
	
	!do n=1,n_basis
!        print '(12f10.5)', h_mat(1:n_basis,n)
!    end do
    
    
    print *, '------- Orthogonality ------- '
    
 	do j = 1,n_basis
            do i = 1,n_basis
                print '(i5,i5,20f10.1)', i, j, &
                dot_product(h_mat(1:n_basis,i) ,h_mat(1:n_basis,j) ) !determines if matrix is orthogonal
            
            end do
        end do

        
        print *, '------- completeness & spectral decomp------- '
        
        do j = 1,n_basis
            do i = 1,n_basis
                f(i,j)= dot_product(h_mat(i,1:n_basis) ,h_mat(j,1:n_basis)&
                *w_eigen(1:n_basis))  !determines if matrix is complete
            end do
        end do 
        
 
        Do i= 1,n_basis
            print '(10f7.2)', f(i,1:n_basis)
        end do 
        
        
        print *, '+++++++++++++++++++++++++++++++++++'
        
        c(1:n_basis,1:n_basis) = conjg(transpose(h_mat(1:n_basis,1:n_basis)))
        
        f(1:n_basis,1:n_basis) = matmul(c(1:n_basis,1:n_basis),h_mat(1:n_basis,1:n_basis))
        f(1:n_basis,1:n_basis) = matmul(h_mat(1:n_basis,1:n_basis),c(1:n_basis,1:n_basis))
        
        Do i= 1,n_basis
            print '(10f7.2)', f(i,1:n_basis)
        end do 
        
        print *, '++++++++++++++++++++++++++++++++++++++++++++'
        
        f(1:n_basis,1:n_basis)=matmul(a(1:n_basis,1:n_basis),h_mat(1:n_basis,1:n_basis))
        d(1:n_basis,1:n_basis) = matmul(c(1:n_basis,1:n_basis),f(1:n_basis,1:n_basis))
        
        Do i= 1,n_basis
            print '(10f7.2)', d(i,1:n_basis)
        end do 
        
        
        print *, '++++++++ inversions of e++++++++++++++ '
        info = 0
        call zgetrf(n_basis,n_basis,h_mat, n_basis, ipiv,info)
        
        if(info/= 0) stop 'stop'
        
        call zgetri(n_basis,h_mat,n_basis,ipiv,work, lwork,info)
        
        if(info/= 0) stop 'stop'
        
        Do i= 1,n_basis
            print '(6f7.2,5x,6f7.2)', c(i,1:n_basis)-h_mat(i,1:n_basis)
        end do 
            
 
     

end program matrix