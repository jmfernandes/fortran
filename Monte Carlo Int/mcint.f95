
program mc_demo
  implicit none
  integer, parameter :: n=10000

  call random_seed()  ! needed to initialize the random number generator used in MC eval
  call MC_integration(n,3.0)

  contains

    subroutine MC_integration(n,end_val)
      implicit none
      integer :: n
      real :: end_val
      real :: x, integral, integral_err
      real (kind=8) :: f, f2 
      integer :: i

      integral = 0.0

      f= 0.0d0
      f2 = 0.0d0

      do i=1,n
         call random_number(x)
         x=x*end_val  ! random_number only returns uniformly distributed from [0.0, 1.0]

         f = f+integrand(x)
         f2 = f2+ (integrand(x)**2)
      end do

      f=f/n
      f2=f2/n

      integral = (end_val-0.0)*f
      integral_err =  (end_val-0.0)*SQRT((f2 - f**2)/n)
      write (*,*) "# MC integration = ",integral,"+/-",integral_err

    end subroutine MC_integration

    function integrand(x) result (value)
      implicit none
      real :: x
      real :: value

      if (x .lt. 0.00001) then
         x = 0.00001
      end if

      value = (x**4)*EXP(X)/((EXP(X)-1.0)**2)
    end function integrand

  end program mc_demo