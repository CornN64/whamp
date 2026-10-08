program p
  implicit none
  integer(8) :: k
  real(kind(1.0d0)) :: a,b,c,d
  complex(kind(1.0d0)) :: z,w,m,q
  do k=1,8
     a=1.0d0+ k*0.13d0; b=2.0d0- k*0.07d0
     c=3.0d0+ k*0.29d0; d=0.5d0- k*0.11d0
     z = cmplx(a,b,kind=kind(1.0d0))
     w = cmplx(c,d,kind=kind(1.0d0))
     m = z*w
     write(*,'(Z16)') transfer(real(m), a)
     write(*,'(Z16)') transfer(aimag(m), a)
     q = z/w
     write(*,'(Z16)') transfer(real(q), a)
     write(*,'(Z16)') transfer(aimag(q), a)
  end do
end program
