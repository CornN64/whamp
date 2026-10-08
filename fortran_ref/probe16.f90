program probe16
  implicit none
  integer, parameter :: d2p = 8
  complex(kind=d2p) :: z1, z2
  real(kind=d2p) :: r1
  z1 = (0.1d0, 0.0d0)
  z2 = (1.23d-11, -4.5d0)
  r1 = 3.5d0
  write (*, *) 'A', z1, 'B', z2
  write (*, *) z1, z2
  write (*, *) z1, r1
  write (*, *) r1, z1
  write (*, *) 'START:', '. X=', z1, 'D=', z2
  write (*, '(I2,A ,I2 ,A,2E16.8,A,2E16.8)') 1, '.', 1, '. X=', z1, ' CX=', z2
end program
