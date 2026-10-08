program probe5
   implicit none
   integer, parameter :: d2p = 8
   real(kind=d2p), parameter :: vals(*) = [ &
     1.0d0, 0.1d0, 1.d-4, 10.d0, 1.d16, 1.d17, -1.0d0, 123456.789d0 ]
   integer :: k
   do k = 1, size(vals)
      write (*, '(a,g26.17,a)') '<', vals(k), '>'
   end do
   write (*, '(a)') '--- multiple items list ---'
   write (*, *) 1.0d0, 2.0d0, 3.0d0
   write (*, '(a,$)') '<'
   write (*, *) 1.0d0, 2.0d0
   write (*, '(a,$)') '>'
   write (*, '(a,$)') '<'
   write (*, *) 1.5d0, (2.5d0,-1.d0), 'str', 7
   write (*, '(a,$)') '>'
end program
