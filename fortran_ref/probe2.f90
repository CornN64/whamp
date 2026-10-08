program probe2
   implicit none
   integer, parameter :: d2p = 8
   real(kind=d2p) :: r
   complex(kind=d2p) :: c
   integer :: i
   character(len=200) :: line
   i = 42
   write (line, '(i0)') i
   write (*, '(a,i0)') 'int_width=', len_trim(line)
   write (line, '(g0)') 1.0d0
   write (*, '(a,a,a,i0)') 'real_nonexp=[', trim(line), '] len=', len_trim(line)
   ! measure raw list-directed by writing to internal then checking
   write (line, '(a)') 'x'
   ! Use a big buffer and INQUIRE not possible; instead write with delimiters
   write (*, '("<",a,">")') ''
end program
