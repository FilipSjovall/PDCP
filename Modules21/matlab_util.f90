module matlab_util

! Last rev.:
! H. Hallberg, 2010-04-14
! M. Ristinmaa, 2011-02-02
! M. Ristinmaa 2011-04-17
!  - included matdrawmesh, matinit
! M. Ristinmaa 2011-04-20
!  - added matplot, matclose
! M. Ristinmaa 2011-04-25
!  - moved f2mat into this module, renamed to matWrt2f
!  - renamed the module from rw2mat to matlab_util
! M. Ristinmaa 2011-04-25
!  - new call matdispmesh
!  - 3d mesh implemented in matdrawmesh 
!  - 3d plot implemented in matplot
! M. Ristinmaa 2011-05-03
!  - bugg in matfilldefmesh 'tr3' fixed
!  - change coordinate, column is node coordninates
!  - change enod, column is nodes in element
! M. Ristinmaa 2011-05-1
!  - small change in matput_r_vector
! M. Ristinmaa 2011-09-14
!  - renamed coordxtr to coordxtr
! M. Ristinmaa 2011-11-10
!  - introduced matcommand
! M. Ristinmaa 2011-11-26
!  - change in matfilldisp to account for problems with dof>size(coord,1) 
!    per node, it is then assumed that the the first dofs are the displacements
!    dofs
!------------------------------------------------------------------------------

!
!   Module matlab_functions contains subroutines for
!   executing various Matlab commands by calls
!   from Fortran
!

use fem_util
use fem_system

implicit none

interface matplot
   module procedure matplot2
   module procedure matplot3
end interface
private matplot2, matplot3

interface matdispmesh
   module procedure matdispmesh1
   module procedure matdispmesh2
end interface
private matdispmesh1,matdispmesh2

interface matfilldispmesh
   module procedure matfilldispmesh1
   module procedure matfilldispmesh2
end interface
private matfilldispmesh1,matfilldispmesh2


interface matgetval
   module procedure matval_r_matrix
   module procedure matval_r_vector
   module procedure matval_r_scalar
!   module procedure matval_i_matrix
!   module procedure matval_i_vector
!   module procedure matval_i_scalar
end interface 
private matval_r_matrix, matval_r_vector, matval_r_scalar
 
 
interface matputval
   module procedure matput_r_matrix
   module procedure matput_r_vector
   module procedure matput_r_scalar
!   module procedure matput_i_matrix
!   module procedure matput_i_vector
!   module procedure matput_i_scalar
end interface
private matput_r_matrix, matput_r_vector, matput_r_scalar
 

interface matwrt2f
   module procedure F2MAT_REAL2
   module procedure F2MAT_INT2
   module procedure F2MAT_REAL_V2
   module procedure F2MAT_INT_V2
end interface
private F2MAT_REAL2,F2MAT_INT2,F2MAT_REAL_V2,F2MAT_INT_V2


!------------------------------------------------------------------------------
	
contains

subroutine matcommand(matptr,string)
  integer*8                         :: matptr
  character                         :: string*(*)

  integer*8                       :: status, engEvalString

write(*,*)string
! test string  string='clf;fill(ex,ey,ed);axis(''equal'')'
  status = engEvalString(matptr,string)

end subroutine matcommand


subroutine matclose(matptr)
   integer*8                         :: matptr

   call engClose(matptr)

end subroutine matclose


subroutine matinit(matptr)
   integer*8                         :: engOpen, matptr

! Start the Matlab engine
 !  matptr = engOpen(' matlab -nodesktop -nosplash -nodisplay  ')
   matptr = engOpen(' matlab -nodesktop -nosplash  ')
   if (matptr .eq. 0) then
      write(*,*) 'Can''t start Matlab engine'
      stop
   end if

end subroutine matinit

subroutine matplot2(ep,x,y)
  implicit none
  integer*8                       :: ep
  double precision                :: x(:), y(:)

  integer*8                       :: status, engEvalString

  call matputval(ep,'x',x)
  call matputval(ep,'y',y)

  status = engEvalString(ep, 'clf;plot(x,y)')
  
end subroutine matplot2

subroutine matplot3(ep,x,y,z)
  implicit none
  integer*8                       :: ep
  double precision                :: x(:), y(:), z(:)

  integer*8                       :: status, engEvalString

  call matputval(ep,'x',x)
  call matputval(ep,'y',y)
  call matputval(ep,'z',z)

  status = engEvalString(ep, 'clf;plot3(x,y,z)')
  
end subroutine matplot3

subroutine matfilldispmesh1(ep,eltype,var,coord,a,enod)
  implicit none
  integer*8                       :: ep
  double precision                :: coord(:,:), var(:), a(:)
  integer                         :: enod(:,:)
  character                       :: eltype*(*)

  integer                         :: ierr
  double precision, allocatable   :: ar(:,:)
  integer                         :: ndme, dfnd, dime

  ndme=maxval(enod)
  dime=size(coord,1)
  dfnd=size(a)/ndme

  allocate(ar(dfnd,ndme),stat=ierr)
  ar=reshape(a,[dfnd,ndme],order=[1,2])
  call matfillmesh(ep,eltype,var,coord+ar(1:dime,:),enod)


  return
end subroutine matfilldispmesh1


subroutine matfilldispmesh2(ep,eltype,var,coord,a,enod,dof)
  implicit none
  integer*8                       :: ep
  double precision                :: coord(:,:), var(:), a(:)
  integer                         :: enod(:,:), dof(:,:)
  character                       :: eltype*(*)

  integer                         :: ierr
  integer                         :: id, nnod,edim
  double precision, allocatable   :: axy(:,:)
 
  nnod=size(coord,2)
  edim=size(coord,1)

  allocate(axy(edim,nnod),stat=ierr)
  axy=0d0
  do id=1,nnod
    axy(:,id)=a(dof(1:edim,id))
  end do

!  call matdrawmesh(ep,eltype,coord+axy,enod)
  call matfillmesh(ep,eltype,var,coord+axy,enod)

  deallocate(axy)


  return
end subroutine matfilldispmesh2


subroutine matfillmesh(ep,eltype,var,coord,enod)
  implicit none
  integer*8                       :: ep
  double precision                :: coord(:,:), var(:)
  integer                         :: enod(:,:)
  character                       :: eltype*(*)

  double precision, allocatable   :: ex(:,:)
  double precision, allocatable   :: ey(:,:)
  double precision, allocatable   :: ez(:,:)
  double precision, allocatable   :: ed(:,:)
  integer*8                       :: status, engEvalString
  character                       :: string*140
  integer                         :: ierr

  if (eltype.eq.'tr3') then
    allocate(ex(3,size(enod,2)),stat=ierr)
    allocate(ey(3,size(enod,2)),stat=ierr)
    allocate(ed(3,size(enod,2)),stat=ierr)
    call coordxtr(ex,ey,coord,enod)
    call extract(ed,var,enod)
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    call matputval(ep,'ed',ed)
    string='clf;fill(ex,ey,ed);axis(''equal'')'
    status = engEvalString(ep,string)
  elseif (eltype.eq.'qu4') then
    allocate(ex(4,size(enod,2)),stat=ierr)
    allocate(ey(4,size(enod,2)),stat=ierr)
    allocate(ed(4,size(enod,2)),stat=ierr)
    call coordxtr(ex,ey,coord,enod)
    call extract(ed,var,enod)
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    call matputval(ep,'ed',ed)
    string='clf;fill(ex,ey,ed);axis(''equal'')'
    status = engEvalString(ep,string)
  elseif (eltype.eq.'qu8') then
    allocate(ex(8,size(enod,2)),stat=ierr)
    allocate(ey(8,size(enod,2)),stat=ierr)
    allocate(ed(8,size(enod,2)),stat=ierr)
    call coordxtr(ex,ey,coord,enod)
    call extract(ed,var,enod)
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    call matputval(ep,'ed',ed)
    string='clf;fill(ex([1 5 2 6 3 7 4 8],:),ey([1 5 2 6 3 7 4 8],:),ed([1 5 2 6 3 7 4 8],:));axis(''equal'')'
    status = engEvalString(ep,string)
  elseif (eltype.eq.'brick8') then
    allocate(ex(4,size(enod,2)),stat=ierr)
    allocate(ey(4,size(enod,2)),stat=ierr)
    allocate(ez(4,size(enod,2)),stat=ierr)
    allocate(ed(4,size(enod,2)),stat=ierr)

    call coordxtr(ex,ey,ez,coord,enod([1,2,3,4],:))
    call extract(ed,var,enod([1,2,3,4],:))
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    call matputval(ep,'ez',ez)
    call matputval(ep,'ed',ed)
!    string='clf;fill3(ex,ey,ez,ed);axis(''equal'');hold on'
    string='fill3(ex,ey,ez,ed);axis(''equal'');hold on'
    status = engEvalString(ep, string)

    call coordxtr(ex,ey,ez,coord,enod([5,6,7,8],:))
    call extract(ed,var,enod([5,6,7,8],:))
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    call matputval(ep,'ez',ez)
    call matputval(ep,'ed',ed)
    string='fill3(ex,ey,ez,ed);axis(''equal'')'
    status = engEvalString(ep, string)

    call coordxtr(ex,ey,ez,coord,enod([2,6,7,3],:))
    call extract(ed,var,enod([2,6,7,3],:))
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    call matputval(ep,'ez',ez)
    call matputval(ep,'ed',ed)
    string='fill3(ex,ey,ez,ed);axis(''equal'')'
    status = engEvalString(ep, string)

    call coordxtr(ex,ey,ez,coord,enod([1,5,8,4],:))
    call extract(ed,var,enod([1,5,8,4],:))
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    call matputval(ep,'ez',ez)
    call matputval(ep,'ed',ed)
    string='fill3(ex,ey,ez,ed);axis(''equal'')'
    status = engEvalString(ep, string)

    call coordxtr(ex,ey,ez,coord,enod([1,5,6,2],:))
    call extract(ed,var,enod([1,5,6,2],:))
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    call matputval(ep,'ez',ez)
    call matputval(ep,'ed',ed)
    string='fill3(ex,ey,ez,ed);axis(''equal'')'
    status = engEvalString(ep, string)

    call coordxtr(ex,ey,ez,coord,enod([4,8,7,3],:))
    call extract(ed,var,enod([4,8,7,3],:))
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    call matputval(ep,'ez',ez)
    call matputval(ep,'ed',ed)
    string='fill3(ex,ey,ez,ed);axis(''equal'')'
    status = engEvalString(ep, string)

  else
    stop 'matdrawmesh/matdispmesh element type is not implemeneted'
  endif

  return
end subroutine matfillmesh




subroutine matdrawmesh(ep,eltype,coord,enod)
  implicit none
  integer*8                       :: ep
  double precision                :: coord(:,:)
  integer                         :: enod(:,:)
  character                       :: eltype*(*)

  double precision, allocatable   :: ex(:,:)
  double precision, allocatable   :: ey(:,:)
  double precision, allocatable   :: ez(:,:)
  integer*8                       :: status, engEvalString
  character                       :: string*140
  integer                         :: ierr

  if (eltype.eq.'tr3') then
    allocate(ex(3,size(enod,2)),stat=ierr)
    allocate(ey(3,size(enod,2)),stat=ierr)
    call coordxtr(ex,ey,coord,enod)
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    string='clf;plot([ex; ex(1,:)],[ey; ey(1,:)],''b'');axis off;axis(''equal'')'
    status = engEvalString(ep,string)
  elseif (eltype.eq.'qu4') then
    allocate(ex(4,size(enod,2)),stat=ierr)
    allocate(ey(4,size(enod,2)),stat=ierr)
    call coordxtr(ex,ey,coord,enod)
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    string='clf;plot([ex; ex(1,:)],[ey; ey(1,:)],''b'');axis(''equal'')'
    status = engEvalString(ep,string)
  elseif (eltype.eq.'qu8') then
    allocate(ex(8,size(enod,2)),stat=ierr)
    allocate(ey(8,size(enod,2)),stat=ierr)
    call coordxtr(ex,ey,coord,enod)
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    string='clf;plot([ex([1 5 2 6 3 7 4 8],:); ex(1,:)], &
                  [ey([1 5 2 6 3 7 4 8],:); ey(1,:)],''b'');axis(''equal'')'
    status = engEvalString(ep,string)
  elseif (eltype.eq.'tet4') then
    allocate(ex(3,size(enod,2)),stat=ierr)
    allocate(ey(3,size(enod,2)),stat=ierr)
    allocate(ez(3,size(enod,2)),stat=ierr)
    call coordxtr(ex,ey,ez,coord,enod([1,2,3],:))
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    call matputval(ep,'ez',ez)
    string='clf;plot3([ex; ex(1,:)],[ey; ey(1,:)],[ez; ez(1,:)],''b'');axis(''equal'');hold on'
    status = engEvalString(ep, string)

    call coordxtr(ex,ey,ez,coord,enod([2,3,4],:))
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    call matputval(ep,'ez',ez)
    string='plot3([ex; ex(1,:)],[ey; ey(1,:)],[ez; ez(1,:)],''b'');axis(''equal'')'
    status = engEvalString(ep,string)

    call coordxtr(ex,ey,ez,coord,enod([3,1,4],:))
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    call matputval(ep,'ez',ez)
    string='plot3([ex; ex(1,:)],[ey; ey(1,:)],[ez; ez(1,:)],''b'');axis(''equal'')'
    status = engEvalString(ep,string)

  elseif (eltype.eq.'brick8') then
    allocate(ex(4,size(enod,2)),stat=ierr)
    allocate(ey(4,size(enod,2)),stat=ierr)
    allocate(ez(4,size(enod,2)),stat=ierr)
    call coordxtr(ex,ey,ez,coord,enod([1,2,3,4],:))
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    call matputval(ep,'ez',ez)
    string='clf;plot3([ex; ex(1,:)],[ey; ey(1,:)],[ez; ez(1,:)],''b'');axis(''equal'');hold on'
    status = engEvalString(ep, string)

    call coordxtr(ex,ey,ez,coord,enod([5,6,7,8],:))
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    call matputval(ep,'ez',ez)
    string='plot3([ex; ex(1,:)],[ey; ey(1,:)],[ez; ez(1,:)],''b'');axis(''equal'')'
    status = engEvalString(ep,string)

    call coordxtr(ex,ey,ez,coord,enod([2,6,7,3],:))
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    call matputval(ep,'ez',ez)
    status = engEvalString(ep,string)

    call coordxtr(ex,ey,ez,coord,enod([1,5,8,4],:))
    call matputval(ep,'ex',ex)
    call matputval(ep,'ey',ey)
    call matputval(ep,'ez',ez)
    status = engEvalString(ep,string)
  else
    stop 'matdrawmesh/matdispmesh element type is not implemeneted'
  endif
  
  if (allocated(ex)) deallocate(ex)
  if (allocated(ey)) deallocate(ey)
  if (allocated(ez)) deallocate(ez)

  return
end subroutine matdrawmesh


subroutine matdispmesh1(ep,eltype,coord,a,enod)
  implicit none
  integer*8                       :: ep
  double precision                :: coord(:,:), a(:)
  integer                         :: enod(:,:)
  character                       :: eltype*(*)

  integer                         :: ierr
  double precision, allocatable   :: ar(:,:)
  integer                         :: ndme, dfnd

  ndme=maxval(enod)
  dfnd=size(a)/ndme

  allocate(ar(dfnd,ndme),stat=ierr)
  ar=reshape(a,[dfnd,ndme],order=[1,2])
  call matdrawmesh(ep,eltype,coord+ar,enod)

  deallocate(ar)
  return
end subroutine matdispmesh1


subroutine matdispmesh2(ep,eltype,coord,a,enod,dof)
  implicit none
  integer*8                       :: ep
  double precision                :: coord(:,:), a(:)
  integer                         :: enod(:,:), dof(:,:)
  character                       :: eltype*(*)

  integer                         :: ierr, id, nnod, edim
  double precision, allocatable   :: axy(:,:)
 
  nnod=size(coord,2)
  edim=size(coord,1)

  allocate(axy(edim,nnod),stat=ierr)
  axy=0d0
  do id=1,nnod
    axy(:,id)=a(dof(1:edim,id))
  end do

  call matdrawmesh(ep,eltype,coord+axy,enod)

  deallocate(axy)
  return
end subroutine matdispmesh2


subroutine matput_r_matrix(ep,name,data)
   implicit none

   integer*8                                :: ep
   character*(*) name
   double precision, dimension(:,:)         :: data
   integer*8         :: rows, cols

   integer*8                     :: engPutVariable
   integer*8                     :: mxCreateDoubleMatrix, mxGetPr
   
   integer*8                     :: status, pt


   rows=size(data,1)
   cols=size(data,2)
   pt = mxCreateDoubleMatrix(rows, cols, 0)
   call mxCopyReal8ToPtr(data,mxGetPr(pt),rows*cols)
   status = engPutVariable(ep,name,pt)
   call mxDestroyArray(pt)

   return
   
end subroutine matput_r_matrix


subroutine matput_r_vector(ep,name,data)
   implicit none

   integer*8                                :: ep
   character*(*) name
   double precision, dimension(:)           :: data
   integer*8                                :: rows, cols

   integer*8                     :: engPutVariable
   integer*8                     :: mxCreateDoubleMatrix, mxGetPr
   
   integer*8                     :: status, pt

   rows=size(data)
   cols=1
   pt = mxCreateDoubleMatrix(rows, cols, 0)
   call mxCopyReal8ToPtr(data,mxGetPr(pt),rows)
   status = engPutVariable(ep,name,pt)
   call mxDestroyArray(pt)

   return
   
end subroutine matput_r_vector

subroutine matput_r_scalar(ep,name,data)
   implicit none

   integer*8                                :: ep
   character*(*) name
   double precision           :: data
   integer*8                                :: rows

   integer*8                     :: engPutVariable
   integer*8                     :: mxCreateDoubleMatrix, mxGetPr
   
   integer*8                     :: status, pt

   rows=1
   pt = mxCreateDoubleMatrix(rows, 1, 0)
   call mxCopyReal8ToPtr(data,mxGetPr(pt),rows)
   status = engPutVariable(ep,name,pt)
   call mxDestroyArray(pt)

   return
   
end subroutine matput_r_scalar



!subroutine matput_i_matrix(ep,name,data,rows,cols)
!   implicit none

!   integer*8                                :: ep
!   character*(*) name
!   integer*8, dimension(:,:),pointer 	    :: data
!   integer*8, optional, intent(in)          :: rows, cols

!   integer*8                     :: engPutVariable, mxGetM, mxGetN
!   integer*8                     :: mxCreateDoubleMatrix, mxGetPr
   
!   double precision, dimension(:),pointer :: datar
!   integer*8                     :: status, pt
!   integer*8                     :: irow, icol


!   pt = mxCreateDoubleMatrix(rows, cols, 0)

!   allocate(datar(irow,icol), STAT=status)
!   call mxCopyReal8ToPtr(datar,mxGetPr(pt),rows*cols)
!   status = engPutVariable(ep,name,pt)
!   call mxDestroyArray(pt)

!   allocate(data(irow,icol),STAT=status)
!   data=int(datar)
!   deallocate(datar)
  
! return
   
!end subroutine matput_i_matrix





subroutine matval_r_matrix(ep,name,data,rows,cols)
   implicit none

   integer*8                                :: ep
   character*(*) name
   double precision, dimension(:,:), allocatable :: data
   integer*8, optional, intent(out)         :: rows, cols

   integer*8                     :: engGetVariable, mxGetM, mxGetN
   integer*8                     :: mxGetPr, mxGetNumberOfElements
   
   integer*8                     :: status, pt
   integer                     :: irow, icol
    
! must deallocate data if it exists
   if (allocated(data)) deallocate(data)

!   write(*,*)'real', name
   pt = engGetVariable(ep,name)
   irow=mxGetM(pt)
   icol=mxGetN(pt)
!   write(*,*) irow, icol
   allocate(data(irow,icol), STAT=status)
   call mxCopyPtrToReal8(mxGetPr(pt), data, mxGetNumberOfElements(pt))
   call mxDestroyArray(pt)


   if (present(rows)) then
     rows=irow
   end if
   if (present(cols)) then
     cols=icol
   end if

   return
   
end subroutine matval_r_matrix


!subroutine matval_i_matrix(ep,name,data,rows,cols)
!   implicit none
!
!   integer*8                                :: ep
!   character*(*) name
!   integer, dimension(:,:),pointer :: data
!   integer*8, optional, intent(out)         :: rows, cols
!
!   integer*8                     :: engGetVariable, mxGetM, mxGetN
!   integer*8                     :: mxGetPr, mxGetNumberOfElements
!   
!   integer*8, dimension(:,:),pointer :: datar
!   integer*8                     :: status, pt
!   integer                     :: irow, icol
!    
!   write(*,*) name
!   pt = engGetVariable(ep,name)
!   irow=mxGetM(pt)
!   icol=mxGetN(pt)
!   write(*,*) irow, icol
!   allocate(datar(irow,icol), STAT=status)
!   call mxCopyPtrToInteger8(mxGetPr(pt), datar, mxGetNumberOfElements(pt))
!   call mxDestroyArray(pt)
!
!   allocate(data(irow,icol),STAT=status)
!   data=int(datar)
!   write(*,*)'data ir ', datar(1,:)
!   deallocate(datar)
!
!   write(*,*)'data i  ', data(1,:)
!   if (present(rows)) then
!     rows=irow
!   end if
!   if (present(cols)) then
!     cols=icol
!   end if
!
!   return
!
!end subroutine matval_i_matrix

subroutine matval_r_vector(ep,name,data,rows,cols)
   implicit none

   integer*8                                :: ep
   character*(*) name
   double precision, dimension(:), allocatable :: data
   integer*8, optional, intent(out)         :: rows, cols

   integer*8                     :: engGetVariable, mxGetM, mxGetN
   integer*8                     :: mxGetPr, mxGetNumberOfElements
   
   integer*8                     :: status, pt
   integer*8                     :: irow, icol
    
! must deallocate data if it exists
   if (allocated(data)) deallocate(data)

!   write(*,*) name
   pt = engGetVariable(ep,name)
   irow=mxGetM(pt)
   icol=mxGetN(pt)
!   write(*,*) irow, icol
   allocate(data(irow), STAT=status)
   call mxCopyPtrToReal8(mxGetPr(pt), data, mxGetNumberOfElements(pt))
   call mxDestroyArray(pt)

   if (present(rows)) then
     rows=irow
   end if
   if (present(cols)) then
     cols=icol
   end if

   return
   
end subroutine matval_r_vector


subroutine matval_r_scalar(ep,name,data)
   implicit none

   integer*8                                :: ep
   character*(*) name
   double precision :: data

   integer*8                     :: engGetVariable, mxGetM, mxGetN
   integer*8                     :: mxGetPr, mxGetNumberOfElements
   
   integer*8                     :: pt
   integer*8                     :: irow, icol
    
!   write(*,*) name
   pt = engGetVariable(ep,name)
   irow=mxGetM(pt)
   icol=mxGetN(pt)
!   write(*,*) irow, icol
!   allocate(data, STAT=status)
   call mxCopyPtrToReal8(mxGetPr(pt), data, mxGetNumberOfElements(pt))
   call mxDestroyArray(pt)

   return
   
end subroutine matval_r_scalar


subroutine f2mat_int2(filename,v_in,name,state)

   implicit none
   
   integer*8                :: nx, ny !, nn, fnn
   integer, intent(in)      :: v_in(:,:)
   integer*8                :: matOpen, matClose, matPutVariable
   integer*8                :: mxCreateNumericMatrix, mxClassIDFromClassName
   integer*8                :: status, mp, pa1, mxGetPr, p
   character(*), intent(in) :: filename, name, state     
      
   nx = size(v_in,1)
   ny = size(v_in,2)
      
   mp  = matOpen(filename, state)
   pa1 = mxCreateNumericMatrix(nx,ny,mxClassIDFromClassName('int32'),0)
   
   call mxCopyInteger4ToPtr(v_in,mxGetPr(pa1) , nx*ny)
   p = matPutVariable(mp,name, pa1)

   call mxDestroyArray(pa1)
   status = matClose(mp)
   
   return
   
end subroutine f2mat_int2


subroutine f2mat_int_v2(filename,v_in,name,state)

   implicit none
   
   integer*8                :: nx, ny !, nn, fnn
   integer, intent(in)      :: v_in(:)
   integer*8                :: matOpen, matPutVariable, matClose !, mxGetPi, mxCreateDoubleMatrix
   integer*8                :: mxCreateNumericMatrix, mxClassIDFromClassName
   integer*8                :: status, mp, pa1, p, mxGetPr
   character(*), intent(in) :: filename, name, state     
      
   nx = size(v_in,1)
   ny = 1
      
   mp  = matOpen(filename, state)
   !pa1 = mxCreateDoubleMatrix(nx,ny,0)
   pa1 = mxCreateNumericMatrix(nx,ny,mxClassIDFromClassName('int32'),0)
   call mxCopyInteger4ToPtr(v_in,mxGetPr(pa1) , nx*ny)
   !call mxCopyReal8ToPtr(dble(v_in),mxGetPr(pa1) , nx*ny)
   p = matPutVariable(mp,name, pa1)

   call mxDestroyArray(pa1)
   status = matClose(mp)

   return

end subroutine f2mat_int_v2


subroutine f2mat_real2(filename,v_in,name,state)

   implicit none

   integer*8                    :: nx, ny
   double precision, intent(in) :: v_in(:,:)
   integer*8                    :: matOpen, mxCreateDoubleMatrix, mxGetPr, matPutVariable, matClose
   integer*8                    :: status, mp, pa1, p
   character(*), intent(in)     :: filename, name, state     
   nx = size(v_in,1)
   ny = size(v_in,2)
   
   mp  = matOpen(filename, state)
   pa1 = mxCreateDoubleMatrix(nx,ny,0)
   call mxCopyReal8ToPtr(v_in, mxGetPr(pa1), nx*ny)
   p = matPutVariable(mp,name, pa1)
   call mxDestroyArray(pa1)
   status = matClose(mp)
   
   return
   
end subroutine f2mat_real2


subroutine f2mat_real_v2(filename,v_in,name,state)

   implicit none
   
   integer*8                    :: nx, ny
   double precision, intent(in) :: v_in(:)
   integer*8                    :: matOpen, mxCreateDoubleMatrix,mxGetPr, matPutVariable, matClose
   integer*8                    :: status, mp, pa1, p
   character(*), intent(in)     :: filename, name, state     
      
   nx = size(v_in,1)
   ny = 1

   mp  = matOpen(filename, state)
   pa1 = mxCreateDoubleMatrix(nx,ny,0)
   call mxCopyReal8ToPtr(v_in, mxGetPr(pa1), nx*ny)
   p   = matPutVariable(mp,name, pa1)
   
   call mxDestroyArray(pa1)
   status = matClose(mp)
   
   return
   
end subroutine f2mat_real_v2



end module matlab_util
