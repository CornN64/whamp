program probe4
   implicit none
   integer, parameter :: d2p = 8
   integer :: k
   real(kind=d2p), parameter :: vals(*) = [ &
     1.0d0, 0.1d0, 0.01d0, 1.d-4, 1.d-5, 1.d-6, 1.d-10, &
     10.d0, 100.d0, 1.d10, 1.d16, 1.d17, 1.d18, 1.d20, &
     1.5d0, 3.141592653589793d0, 123456.789d0, 0.0001d0, &
     1.0d-16, 1.2345678901234567d0, -1.0d0, 123456789012345678.d0 ]
   do k = 1, size(vals)
      write (*, '(a,$)') '<'
      write (*, '(g0.17)') vals(k)
      write (*, '(a,$)') '>'
      write (*, '(a,$)') '['
      write (*, *) vals(k)
      write (*, '(a,$)') ']'
   end do
end program
