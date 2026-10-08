program p
  implicit none
  real(kind(1.0d0)) :: PM(3)
  complex(kind(1.0d0)), allocatable :: fOUT(:)
  integer :: cycleZFirst
  PM = [0.0d0, 0.0d0, 10.0d0]
  cycleZFirst = 1
  allocate(fOUT(4))
  fOUT = [cmplx(1.0d0,2.0d0), cmplx(0.0d0,0.0d0), cmplx(-3.5d0,0.0d0), cmplx(0.0d0,1.0d0)]
  write (*, *) 'PM=', PM
  write (*, *) 'cycleZfirst=', cycleZFirst
  write (*, *) 'fOUT=', fOUT
  deallocate(fOUT)
  write (*, *) 'fOUT=', fOUT
end program
