!this is a program to return the change for a money transaction
program change
	integer :: quarters, dimes, nickels, pennies, ierr
	real :: q, d, n , p , money, price, diff, total
	character(1) :: yn
	q = 0.25
	d = 0.10
	n = 0.05
	p = 0.01
	interactive_loop: do
		write(*,*) 'Enter the price of item'
		read(*,*,iostat=ierr) price
		!stops user from using non-numbers
		if (ierr /= 0) then 
      		write(*,*) 'Error, invalid input.'
      		cycle interactive_loop
    	end if
		write(*,*) 'Enter the amount of money given'
		read(*,*,iostat=ierr) money
		!stops user from using non-numbers
		if (ierr /= 0) then 
      		write(*,*) 'Error, invalid input.'
      		cycle interactive_loop
    	end if
		diff = money - price
		total = diff
		!two lines below correct a float error that causes the number of quarters to be reduced by 1.
		ii = ANINT(diff*100.0) !Multiply by 100 and round to int
        diff  = ii / 100.0     !Floating point divide by 100 
		quarters = diff / q
		diff = diff - (q * quarters)
		!two lines below correct a float error that causes the number of dimes to be reduced by 1.
		ii = ANINT(diff*100.0) !Multiply by 100 and round to int
        diff  = ii / 100.0     !Floating point divide by 100 
		dimes = diff / d
		diff = diff - (d * dimes)
		!two lines below correct a float error that causes the number of nickels to be reduced by 1.
		ii = ANINT(diff*100.0) !Multiply by 100 and round to int
        diff  = ii / 100.0     !Floating point divide by 100 
		nickels = diff / n
		diff = diff - (n * nickels)
		!two lines below correct a float error that causes the number of pennies to be reduced by 1.
		ii = ANINT(diff*100.0) !Multiply by 100 and round to int
        diff  = ii / 100.0     !Floating point divide by 100 
		pennies = diff / p
		diff = diff - (p * pennies)
		write(*,*)'Customer recieves' , quarters , 'quarters,' , dimes , 'dimes,' , nickels , 'nickels,' , 'and ' , pennies , 'pennies.'
		yn_loop: do
      		write(*,*) 'Perform another calculation? y[n]'
      		read(*,'(a1)') yn
      		if (yn=='y' .or. yn=='Y') exit yn_loop
      		if (yn=='n' .or. yn=='N' .or. yn==' ') exit interactive_loop
   		end do yn_loop
	end do interactive_loop
end program change