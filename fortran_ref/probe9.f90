program probe9
   implicit none
   integer, parameter :: d2p = 8
   real(kind=d2p), parameter :: vals(*) = [ &
     0.5d0, 0.05d0, 0.0d0, 1.0d0, 12.34d0, 123.4d0, -0.5d0, -12.34d0 ]
   integer :: k
   do k = 1, size(vals)
      write (*, '("<",f4.2,">")') vals(k)
      write (*, '("<",f5.2,">")') vals(k)
      write (*, '("<",f6.3,">")') vals(k)
      write (*, '("<",f11.4,">")') vals(k)
   end do
end program
