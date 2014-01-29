!this code multiplies two complex numbers by i
program complex1
	implicit none
  	! Define variables and constants
  	complex, parameter :: i = (0, 1)
  	complex :: x, y
  	x = (1, 1); y = (1, -1)
  	write(*,*) i * x * y
end program complex1