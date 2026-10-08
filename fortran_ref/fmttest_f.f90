program fmttest_f
   implicit none
   integer, parameter :: d2p = 8
   real(kind=d2p), parameter :: vals(*) = [ &
     8.9786d0, 2.8d0, 1.0e6_d2p, 0.01d0, 0.001d0, 1.0d0, 5.0d0, 0.0d0, 1.0640799d-1, &
     0.999d0, 9.999d0, 3.34d-11, -0.10640799d0, -1.2345d-100, 1.5d0, 3.141592653589793d0, &
     123456.789d0, 1.d-4, 1.d-16, 1.2345678901234567d0, -1.0d0, 1.2345678901234567d17, &
     1.d16, 1.d17, 1.d18, 1.d20, 1.d-5, 1.d-6, 1.d-10, 0.5d0, 0.05d0, 12.34d0, 123.4d0, -12.34d0, &
     0.0022d0, 1.0640799d-01, 3.34d-11, 1.0d-3, 2.0d0, 3.0d0 ]
   integer :: k
   character(len=26) :: buf
   do k = 1, size(vals)
      write (*, '("<",e12.2,">")') vals(k)
      write (*, '("<",e13.5,">")') vals(k)
      write (*, '("<",e12.3,">")') vals(k)
      write (*, '("<",1pe11.5,">")') vals(k)
      write (*, '("<",1pe12.5,">")') vals(k)
      write (*, '("<",e8.2,">")') vals(k)
      write (*, '("<",e9.3,">")') vals(k)
      write (*, '("<",1pe12.3,">")') vals(k)
      write (*, '("<",f11.4,">")') vals(k)
      write (*, '("<",f10.4,">")') vals(k)
      write (*, '("<",f9.5,">")') vals(k)
      write (*, '("<",f4.2,">")') vals(k)
      write (*, '("<",f5.2,">")') vals(k)
      write (*, '("<",f6.3,">")') vals(k)
      write (*, '("<",f10.7,">")') vals(k)
      write (buf, *) vals(k)
      write (*, '("<",a,">")') buf
   end do
end program
