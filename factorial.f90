!program to calulate the factorial of a number
program factorial
	!declare my constants
	integer :: number, result
	character(1) :: yn
	!start a master do loop so that I can do multiple functions within it
	interactive_loop: do
		print *, 'Enter the Number to Factorize'
		read *, number
		result = 1
		do i=1, number, 1
			result = result*i
		end do
		print *, 'The answer is' , result
		yn = ' '
   		yn_loop: do
      		write(*,*) 'Perform another calculation? y[n]'
      		read(*,'(a1)') yn
      		if (yn=='y' .or. yn=='Y') exit yn_loop
      		if (yn=='n' .or. yn=='N' .or. yn==' ') exit interactive_loop
   		end do yn_loop
   	end do interactive_loop
end program factorial