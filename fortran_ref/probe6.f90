program probe6
   implicit none
   integer, parameter :: d2p = 8
   real(kind=d2p) :: v
   integer :: k
   real(kind=d2p), parameter :: vals(*) = [ &
     8.9786d0, 2.8d0, 1.0d6, 0.01d0, 0.001d0, 1.0d0, 5.0d0, 0.0d0, 1.0640799d-1 ]
   do k = 1, size(vals)
      v = vals(k)
      write (*, '("<",f11.4,">")') v
      write (*, '("<",f10.4,">")') v
      write (*, '("<",e11.5,">")') v
      write (*, '("<",e12.5,">")') v
      write (*, '("<",f9.5,">")') v
      write (*, '("<",f4.2,">")') v
      write (*, '("<",f5.2,">")') v
      write (*, '("<",f10.7,">")') v
      write (*, '("<",e12.2,">")') v
      write (*, '("<",e13.5,">")') v
   end do
end program
