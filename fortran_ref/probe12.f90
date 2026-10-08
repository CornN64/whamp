program probe12
   implicit none
   integer, parameter :: d2p = 8
   complex(kind=d2p) :: c
   integer :: i
   real(kind=d2p) :: r
   c = (1.5d0,-2.25d0)
   r = 2.5d0
   i = 7
   write (*, *) 'START:', '. X=', c, 'D=', c, 'DX=', c
   print *
   write (*, *) 'A', 'B'
   write (*, *) 'A', r
   write (*, *) r, 'A'
   write (*, *) -1.d0, -2.d0
   write (*, *) .true., .false.
   write (*, *) 'Max iterations: ', i
   write (*, *) 'Option ', 'x', ' is unknown'
end program
