program montecarlo

    use numtype
    implicit none
    integer :: i, count, seed
    real(dp) :: r2(2), R3(3), SUM, PI_value
    
    count = 10
    call random_seed(seed)
    
    do i = 1, count
    
        !print *, i
        call random_number(r2)
        print *, r2
    
    end do 
    
    
    count = 1000
    sum = 0._dp
    
    do i = 1, count
    
    call random_number(r2)
        if(r2(1)**2 + r2(2)**2 <1) then
            sum = sum +1
            pi_value = 4*sum/i
            write(1, *) r2
        else 
            write(2,*) r2
        end if
        
        if(mod(i,100) ==0) print*, i, pi_value
    
    end do 
    
    
    Print *, '----------------ZOMG----------------'
    
    count = 100000
    sum = 0._dp
    
    do i = 1, count
    
    call random_number(r3)
        if(r3(1)**2 + r3(2)**2 +r3(3)**2 <1) then
            sum = sum +1
            pi_value = 6*sum/i
        end if
        
    end do
        
        if(mod(i,1000) < 100000) print *, i, pi_value

end program montecarlo