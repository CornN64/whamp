program p
  implicit none
  real(kind(1.0d0)) :: x
  integer :: i
  do i = 1, 10
     x = 0.1d0 * i + 0.37d0
     write(*,'(Z16)') transfer(exp(x), x)
     write(*,'(Z16)') transfer(log(x), x)
     write(*,'(Z16)') transfer(sin(x), x)
     write(*,'(Z16)') transfer(cos(x), x)
  end do
end program
