program probe8
   implicit none
   integer, parameter :: d2p = 8
   real(kind=d2p) :: rarr(3)
   rarr = [1.d0, 2.d0, 3.d0]
   write (*, *) 'PM=', rarr
   write (*, '(a,f11.4,a)') '# PLASMA FREQ.:', 8.9786d0, 'KHZ'
   write (*, *)
   write (*, '(a)') 'X'
end program
