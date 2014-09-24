

module NumType

    save
    integer, parameter :: dp = kind(1.d0)       ! double precision
    real(dp), parameter :: pi = 4*atan(1._dp)   ! Pi= 3.1415926
    complex(dp), parameter :: iic = (0._dp,1._dp) ! complex unit


end module NumType


