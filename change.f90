!this is a program to return the change for a money transaction
program change
	integer :: quarters, dimes, nickels, pennies
	real :: q, d, n , p , money, price, diff
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
		!do (while diff > 0)
		quarters = diff / q
		diff = diff - (q * quarters)
		print *, diff, quarters
		!end do
		yn_loop: do
      		write(*,*) 'Perform another calculation? y[n]'
      		read(*,'(a1)') yn
      		if (yn=='y' .or. yn=='Y') exit yn_loop
      		if (yn=='n' .or. yn=='N' .or. yn==' ') exit interactive_loop
   		end do yn_loop
	end do interactive_loop
end program change