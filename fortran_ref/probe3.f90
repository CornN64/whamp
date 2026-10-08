program probe3
   implicit none
   integer, parameter :: d2p = 8
   real(kind=d2p) :: r
   complex(kind=d2p) :: c
   integer :: i
   character(len=500) :: buf
   integer :: n
   ! integer list-directed width
   write (buf, '(i0)') 42
   write (*, '(a,i0)') 'INT i0 width=', len_trim(buf)
   ! Write list directed to internal buffer, then examine
   write (buf, *) 42
   n = len(buf)
   write (*, '(a,i0,a,i0)') 'int ld: rawlen=', n, ' '
   write (*, '(a,a,a)') 'int ld: [', buf, ']'
   write (buf, *) 1.0d0
   write (*, '(a,a,a)') 'real ld: [', buf, ']'
   write (buf, *) 1.d-5
   write (*, '(a,a,a)') 'reale xp:[', buf, ']'
   write (buf, *) (1.5d0,-2.25d0)
   write (*, '(a,a,a)') 'cplx ld: [', buf, ']'
   write (buf, *) 'abc'
   write (*, '(a,a,a)') 'str  ld:[', buf, ']'
   write (buf, *) .true.
   write (*, '(a,a,a)') 'log  ld:[', buf, ']'
end program
