
module thiele_approx

	use numtype
	implicit none
	integer, parameter :: maxpt = 50
	
	contains
        
        subroutine thiele_coef( nn, zn, fn, an )
        
            use numtype
            implicit none
            real(dp), dimension(maxpt) :: zn, fn, an
            real(dp), dimension(maxpt,maxpt) :: gn
            integer :: nn, n, nz
            
            gn(1,1:nn) = fn(1:nn)           
            do n = 2, nn
                do nz = n, nn
                    gn(n,nz) = ( gn(n-1,n-1) - gn(n-1,nz) )/ &
                                            ( (zn(nz)-zn(n-1) ) * gn(n-1,nz) )
                end do
            end do
            forall ( n = 1:nn ) an(n) = gn(n,n)
                    
        end subroutine thiele_coef
        
        function thiele_cf (z, nn, zn, an)  result(cfrac)
        
            use numtype
            implicit none
            real(dp) :: z
            real(dp), dimension(maxpt) :: zn, an
            integer :: nn, n
            real(dp) :: cf0(2), cf1(2), cf(2), cfrac
            
            cf0(1) = 0._dp; cf0(2) = 1._dp
            cf1(1) = an(1); cf1(2) = 1._dp
            do n = 1, nn-1
                cf = cf1 + (z - zn(n)) * an(n+1) * cf0
                cf0 = cf1;  cf1 = cf
            end do
            cfrac = cf(1)/cf(2)
            
        end function thiele_cf
              
end module thiele_approx

