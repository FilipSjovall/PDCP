module memory_util

! last modified
! M. Ristinmaa 2011-04-13
! M. Ristinmaa 2011-05-07
!  - Initial version that allow change of columns in allocation
! -----------------------------------------------------------------------------

!
! these routines can be used both for increasing the array
! as well as decreasing the array	

implicit none

! debug purposes
integer    :: debug=0
private debug


interface reallocate
   module procedure reallocate_dp_array
   module procedure reallocate_dp_array_col
   module procedure reallocate_dp_vector
   module procedure reallocate_int_array
   module procedure reallocate_int_array_col
   module procedure reallocate_int_vector
end interface


!------------------------------------------------------------------------------

contains


subroutine reallocate_dp_array_col(a,n,m)
  implicit none
  double precision, allocatable    :: a(:,:)
  integer                          :: n,m

  double precision, allocatable    :: b(:,:)
  integer                          :: rows, cols, ierr

  if (debug.eq.1) write(*,*)'---------------entering reallocate_dp_array_col'

  rows=size(a,1)
  cols=size(a,2)
  allocate(b(n,m), stat=ierr)
  b=0d0
  write(*,*)'rows cols ', rows,cols
  write(*,*)'n m', n,m
 
  if ((n.gt.rows).and.(m.gt.cols)) then
write(*,*)'1'
     b(1:rows,1:cols)=a(1:rows,1:cols)
  elseif ((n.le.rows).and.(m.le.cols)) then
write(*,*)'2'
     b=a(1:n,1:m) 
  elseif (n.le.rows) then
write(*,*)'3', n,cols
     b=a(1:n,1:cols)
  elseif (m.le.cols) then
write(*,*)'4'
     b=a(1:rows,1:m)
  end if
write(*,*)'5'

  deallocate(a, stat=ierr)
  allocate(a(n,m), stat=ierr)
  
  a=b

  deallocate(b, stat=ierr)

  if (debug.eq.1) write(*,*)'---------------leaving  reallocate_dp_array_col'
  return
end subroutine reallocate_dp_array_col

subroutine reallocate_dp_array(a,n)
  implicit none
  double precision, allocatable    :: a(:,:)
  integer                          :: n

  double precision, allocatable    :: b(:,:)
  integer                          :: rows, cols, ierr

  if (debug.eq.1) write(*,*)'---------------entering reallocate_dp_array'

  rows=size(a,1)
  cols=size(a,2)
  allocate(b(n,cols), stat=ierr)
  
  if (n.gt.rows) then
     b=0d0
     b(1:rows,1:cols)=a(1:rows,1:cols)
  else
     b=a(1:n,1:cols)
  end if

  deallocate(a, stat=ierr)
  allocate(a(n,cols), stat=ierr)
  
  a=b

  deallocate(b, stat=ierr)

  if (debug.eq.1) write(*,*)'---------------leaving  reallocate_dp_array'
  return
end subroutine reallocate_dp_array


subroutine reallocate_dp_vector(a,n)
  implicit none
  double precision, allocatable    :: a(:)
  integer                          :: n

  double precision, allocatable    :: b(:)
  integer                          :: rows, ierr

  	
  if (debug.eq.1) write(*,*)'---------------entering reallocate_dp_vector'
  
  rows=size(a,1)
  allocate(b(n), stat=ierr)
!  write(*,*)n,rows
  if (n.gt.rows) then
     b=0d0
     b(1:rows)=a(1:rows)
  else
     b=a(1:n)
  end if

  deallocate(a, stat=ierr)
  allocate(a(n), stat=ierr)
  
  a=b
!write(*,*)'a ',a
!write(*,*)'b ',b

  deallocate(b, stat=ierr)

  if (debug.eq.1) write(*,*)'---------------leaving  reallocate_dp_vector'
  return
end subroutine reallocate_dp_vector

subroutine reallocate_int_array_col(a,n,m)
  implicit none
  integer, allocatable             :: a(:,:)
  integer                          :: n,m

  double precision, allocatable    :: b(:,:)
  integer                          :: rows, cols, ierr

  if (debug.eq.1) write(*,*)'---------------entering reallocate_int_array_col'

  rows=size(a,1)
  cols=size(a,2)
  allocate(b(n,m), stat=ierr)
  b=0 
  
  write(*,*)'rows cols ', rows,cols
  write(*,*)'n m', n,m

  if ((n.gt.rows).and.(m.gt.cols)) then
write(*,*)'1'
     b(1:rows,1:cols)=a(1:rows,1:cols)
  elseif ((n.le.rows).and.(m.le.cols)) then
write(*,*)'2'
     b=a(1:n,1:m)
  elseif (n.le.rows) then
write(*,*)'3'
     b=a(1:n,1:cols)
  elseif (m.le.cols) then
write(*,*)'4'
     b=a(1:rows,1:m)
  end if

  deallocate(a, stat=ierr)
  allocate(a(n,m), stat=ierr)
  
  a=b

  deallocate(b, stat=ierr)

  if (debug.eq.1) write(*,*)'---------------leaving  reallocate_int_array_col'
  return
end subroutine reallocate_int_array_col


subroutine reallocate_int_array(a,n)
  implicit none
  integer, allocatable             :: a(:,:)
  integer                          :: n

  integer, allocatable             :: b(:,:)
  integer                          :: rows, cols, ierr

  	
  if (debug.eq.1) write(*,*)'---------------entering reallocate_int_array'

  rows=size(a,1)
  cols=size(a,2)
  allocate(b(n,cols), stat=ierr)
  
  if (n.gt.rows) then
     b=0
     b(1:rows,1:cols)=a(1:rows,1:cols)
  else
     b=a(1:n,1:cols)
  end if

  deallocate(a, stat=ierr)
  allocate(a(n,cols), stat=ierr)
  
  a=b

  deallocate(b, stat=ierr)

  if (debug.eq.1) write(*,*)'---------------leaving  reallocate_int_array'
  return
end subroutine reallocate_int_array


subroutine reallocate_int_vector(a,n)
  implicit none
  integer, allocatable    :: a(:)
  integer                          :: n

  integer, allocatable    :: b(:)
  integer                          :: rows, ierr

  	
  if (debug.eq.1) write(*,*)'---------------entering reallocate_int_vector'

  rows=size(a,1)
  allocate(b(n), stat=ierr)
  if (n.gt.rows) then
     b=0
     b(1:rows)=a(1:rows)
  else
     b=a(1:n)
  end if
  deallocate(a, stat=ierr)
  allocate(a(n), stat=ierr)
  
  a=b

  deallocate(b, stat=ierr)

  if (debug.eq.1) write(*,*)'---------------leaving  reallocate_int_vector'
  return
end subroutine reallocate_int_vector

end module memory_util
