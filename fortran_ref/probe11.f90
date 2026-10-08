program probe11
   implicit none
   integer, parameter :: d2p = 8
   real(kind=d2p) :: rarr(3)
   complex(kind=d2p) :: carr(4), c
   integer :: i
   rarr = [1.d0, 2.d0, 3.d0]
   c = (1.5d0,-2.25d0)
   carr = [(1.d0,2.d0), (0.d0,0.d0), (-1.23456789d10, 3.1d-11), (123456.789d0, 0.0001d0)]
   i = 42
   write (*, *) 'PM=', rarr
   write (*, *) 'cycleZfirst=', i
   write (*, *) 'fOUT=', carr
   write (*, *) ' TOO HEAVILY DAMPED!'
   write (*, *) '   '
   write (*)
   write (*, *) c
   write (*, *) 'START:', '. X=', c, 'D=', c, 'DX=', c
end program
