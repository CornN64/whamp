program probe
   implicit none
   integer, parameter :: d2p = 8
   real(kind=d2p) :: r
   complex(kind=d2p) :: c
   integer :: i
   i = 42
   do
      ! measure integer field
      write (*, '(a,i0,a)') '[', i, ']'
      exit
   end do
   ! integer list-directed single
   write (*, '(a,$)') '<'
   write (*, *) i
   write (*, '(a,$)') '>'
   ! real list-directed single
   write (*, '(a,$)') '<'
   write (*, *) 1.0d0
   write (*, '(a,$)') '>'
   ! real various
   r = 1.23456789012345d0
   write (*, '(a,$)') '<'
   write (*, *) r
   write (*, '(a,$)') '>'
   write (*, '(a,$)') '<'
   write (*, *) 1.d-5
   write (*, '(a,$)') '>'
   write (*, '(a,$)') '<'
   write (*, *) 1.d20
   write (*, '(a,$)') '>'
   write (*, '(a,$)') '<'
   write (*, *) 0.d0
   write (*, '(a,$)') '>'
   ! complex single
   c = (1.5d0,-2.25d0)
   write (*, '(a,$)') '<'
   write (*, *) c
   write (*, '(a,$)') '>'
   ! string single
   write (*, '(a,$)') '<'
   write (*, *) 'abc'
   write (*, '(a,$)') '>'
   ! logical single
   write (*, '(a,$)') '<'
   write (*, *) .true.
   write (*, '(a,$)') '>'
end program
