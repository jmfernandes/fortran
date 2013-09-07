!this is a program to return the change for a money transaction
program change
	integer :: quarters, dimes, nickels, pennies
	real :: q, d, n , p , money, price, diff, total
	character(1) :: yn
	q = 0.25
	d = 0.10
	n = 0.05
	p = 0.01
	interactive_loop: do
		print *, 'Enter the price of item'
		read *, price
		print *, 'Enter the amount of money given'
		read *, money
		diff = money - price
		total = diff
		write(*,110) 'total=', diff 
  		110 format (A,F5.3)
  		print *, diff
		!do (while diff > 0)
		quarters = diff / q
		diff = diff - (q * quarters)
		dimes = diff / d
		diff = diff - (d * dimes)
		nickels = diff / n
		diff = diff - (n * nickels)
		pennies = diff / p
		diff = diff - (p * pennies)
		write(*,*)'Customer recieves' , quarters , 'quarters,' , dimes , 'dimes,' , nickels , 'nickels,' , 'and ' , pennies , 'pennies.'
		!end do
		yn_loop: do
      		write(*,*) 'Perform another calculation? y[n]'
      		read(*,'(a1)') yn
      		if (yn=='y' .or. yn=='Y') exit yn_loop
      		if (yn=='n' .or. yn=='N' .or. yn==' ') exit interactive_loop
   		end do yn_loop
	end do interactive_loop
end program change