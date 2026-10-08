program probe15
  character(len=20) :: s
  character(len=5)  :: t
  s = 'abc'
  t = 'xy'
  write (*, *) "Input parameter:", s
  write (*, *) "Option '", trim(s), "' is unknown"
  write (*, *) "short:", t
  write (*, *) "two:", t, t
  write (*, *) trim(t), "end"
  write (*, *) "Max iterations: ", 42
  print *, '   '
  print *
  write (*, *) "Reading file: ", "Ex1"
end program
