
program money_prog

	use numtype
    use chebyshev
	implicit none
    integer :: n, min, max,i,j
    real(dp) :: m,dm,T

    integer, parameter :: npart = 5, nmin = 1, nmax=5
    real(dp), parameter :: Total_M = 100._dp, Total_N=10._dp

    integer :: B(0:npart), B_old(0:npart)
    
    !=================================================================


    T = real(Total_M)/real(Total_N)




    !=================================================================

    min = 0
    max = 100

    m = 0._dp
    n = 100
    dm = abs(max-min)/real(n)

    !=================================================================

    do while (m <= max)
        write(1,*) m, func(m,T)
        m = m + dm
    end do

    !=================================================================


    B = (/ 1,2,3,3,4,5 /)
    B_old = B

    i = 1
    do while (i<=1)
        do j = 0,npart
            B(j) = B(j) -1
            if (B(j) < 1) then 
                B(j) = B_old(j)
            else if (B(j) > 5) then
                B(j) = 5
            else 
                B(j) = B(j)
            end if 
        end do
        B_old = B
        i = i + 1
    end do

    print *, B

contains 


    function func(x,T) result(f)

    use numtype , only : dp 
    implicit none

    real(dp) :: x,f
    real(dp) :: T

    f = x**B_func()*exp(-x/T)

    end function func

    function B_func() result(f)

    use numtype , only : dp 
    implicit none

    real(dp) :: f
    real(dp), parameter :: G = 1/3._dp

    f = -1._dp-log(2._dp)/log(1._dp-G)

    end function B_func



end program money_prog


