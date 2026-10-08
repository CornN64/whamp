program p
  implicit none
  real(kind(1.0d0)) :: a,b,c,d,r,i
  integer :: k
  do k = 1, 8
     a = 1.0d0 + k*0.13d0; b = 2.0d0 - k*0.07d0
     c = 3.0d0 + k*0.29d0; d = 0.5d0 - k*0.11d0
     r = a*c - b*d; i = a*d + b*c
     write(*,'(Z16)') transfer(r, a)
     write(*,'(Z16)') transfer(i, a)
     ! complex divide
     r = (a*c + b*d)/(c*c + d*d)
     i = (b*c - a*d)/(c*c + d*d)
     write(*,'(Z16)') transfer(r, a)
     write(*,'(Z16)') transfer(i, a)
  end do
end program
