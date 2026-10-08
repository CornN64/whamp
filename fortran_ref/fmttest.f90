program fmttest
   implicit none
   integer, parameter :: d2p = 8
   real(kind=d2p) :: r, rarr(3)
   complex(kind=d2p) :: c, carr(4)
   integer :: i
   character(len=20) :: s
   logical :: l
   r = 1.23456789012345d0
   c = (1.5d0, -2.25d0)
   i = 42
   l = .true.
   s = 'abc'
   rarr = [1.d0, 2.d0, 3.d0]
   carr = [(1.d0,2.d0), (0.d0,0.d0), (-1.23456789d10, 3.1d-11), (123456.789d0, 0.0001d0)]
   write (*, *) 'PM=', rarr
   write (*, *) 'cycleZfirst=', i
   write (*, *) 'fOUT=', carr
   write (*, *) 'Input parameter:', s
   write (*, *) 'Reading file: ', s
   write (*, *) 'Max iterations: ', i
   write (*, *) 'r=', r, ' c=', c, ' i=', i, ' l=', l
   write (*, *) ' TOO HEAVILY DAMPED!'
   write (*, *) '   '
   write (*)
   write (*, *) c
   write (*, *) 0.d0, -0.d0, 1.d0, 1.d-5, 1.d20, 1.d-20
   write (*, *) carr(1)
   write (*, *) 'START:', '. X=', c, 'D=', c, 'DX=', c
end program
