program probe17
  implicit none
  integer, parameter :: d2p = 8
  real(kind=d2p) :: a, b, c
  real(kind=d2p) :: plg
  a = 1836.1
  b = 1836.1_d2p
  write (*, '(e26.17,e26.17)') a, b
  write (*, *) 'diff:', a - b
  c = 0.01
  write (*, *) 'c=', c
  plg = 0.5d0
  c = 10.**plg
  write (*, *) '10.**plg=', c
  c = 10.d0**plg
  write (*, *) '10.d0**plg=', c
end program
