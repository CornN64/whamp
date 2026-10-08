program probe19
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer, parameter :: d2p = 8
  real(kind=d2p) :: nanv, pinf, ninf, big, tinyv, z
  nanv  = ieee_value(1.d0, ieee_quiet_nan)
  pinf  = ieee_value(1.d0, ieee_positive_inf)
  ninf  = ieee_value(1.d0, ieee_negative_inf)
  big   = 1.23456789d20
  tinyv = 1.23456789d-20
  z = 0.d0
  write(*,'("|",e9.3,"|",e9.3,"|",e9.3,"|")') nanv, pinf, ninf
  write(*,'("|",e9.3,"|",e9.3,"|",e9.3,"|")') big, tinyv, z
  write(*,'("|",e8.2,"|",e8.2,"|",e8.2,"|")') nanv, pinf, 0.0d0
  write(*,'("|",e9.3,"|",f9.3,"|",g9.3,"|")') ninf, ninf, ninf
  write(*,'("|",e12.5,"|",e12.5,"|")') 9.999999d0, 0.0000001d0
  write(*,'("|",e9.3,"|")') 123456789012345678901234567890.0d0
  write(*,'("|",e9.3,"|")') -1.2345d20
  write(*,'("|",e9.3,"|")') 0.12345d20
  write(*,'("|",e9.3,"|")') 0.012345d20
  write(*,'("|",e9.3,"|")') 1.0d0
  write(*,'("|",e9.3,"|")') -1.0d0
  write(*,'("|",e9.3,"|")') 999999.9d0
  write(*,'("|",e9.3,"|")') 9999999.0d0
  write(*,'("|",e20.10,"|")') nanv
  write(*,'("|",e30.20,"|")') pinf
  write(*,'("|",e9.3,"|",e9.3,"|")') 0.0d0, -0.0d0
end program probe19
