program Init
!
use mesh
!
implicit none
!
character(4)                :: fname="mesh"
real(kind=8), allocatable   :: coord(:,:)
integer, allocatable        :: enod(:,:)
integer                     :: nelm, nnod, el
!
type(sparse)                :: K
integer                             :: ii, node
integer, allocatable                :: cellSym(:,:)
double precision, allocatable       :: valSym(:) 
!
call read_mesh(fname,coord,enod)
!
print * , shape(coord)
!
nelm = size(enod,1)
nnod = size(coord,1)
!
print * , nelm, nnod
!
!
! Init stiffness matrix K



end program Init