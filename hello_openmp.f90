program hello

use OMP_LIB

implicit none

integer :: thread_id = 0

!$OMP PARALLEL

print '("Hello World from thread ", i0)', OMP_GET_THREAD_NUM()

!$OMP END PARALLEL

end program