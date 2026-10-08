program probe20
  implicit none
  integer, parameter :: d2p = kind(1.0d0)
  complex(kind=d2p) :: E(6,4), XX(10)
  real(kind=d2p) :: PP(10), ZZ(10)
  integer :: I, J, N, M
  E = (1.0d0, 2.0d0)
  XX = (3.0d0, 4.0d0)
  PP = 5.0d0
  ZZ = 6.0d0
  DO 30 J = 1, 6
     N = 1 + J/4 + J/6
     M = J - J/4*2 - J/6
     WRITE (*, 29) (N, M, E(J, I), I=1, 4)
29   FORMAT(' E', 2I1, '=', 2(1PE10.2), '  EX', 2I1, '=', 2(1PE10.2),&
            & 'EZ', 2I1, '=', 2(1PE10.2), '  EP', 2I1, '=', 2(1PE10.2),/)
30 end do
  PRINT 6
  PRINT 31, XX, PP, ZZ
31 FORMAT(1P, ' XX=', 12E12.3/' PP=', 6E12.3/' ZZ=', 6E12.3/)
  PRINT *, '=== end ==='
6 FORMAT('  ')
end program probe20
