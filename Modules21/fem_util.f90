module fem_util
 
! last modified
! M. Ristinmaa 2011-04-12
! M. Ristinmaa 2011-04-17
!  - bugg in coordxtr assumed 3 nods
!  - 3D coordxtr added
! M. Ristinmaa 2011-04-28
!  - defdof rewritten
!  - neigbel introduced, previously in adaptive	
! M. Ristinmaa 2011-05-01
!  - change order of call variables in coordxtr
! M. Ristinmaa 2011-05-03
!  - implemented integ2nod averaging 'tr3' and 'qu4'
!  - change coordinate, column is node coordninates
!  - change enod, column is nodes in element
! M. Ristinmaa 2011-05-07
!  - small change in updcoord
! M. Ristinmaa 2011-05-07
!  - small change in integr2nod such that several variables 
!    can be used in call
! M. Ristinmaa 2011-05-11
!  - included inv3
! M. Ristinmaa 2011-09-14
!  - renamed coordxtr to coordel
!  - change order arguments inv3
! M. Ristinmaa 2011-11-26
!  - change in updcoord to account for problems with dof>size(coord,1) 
!    per node, it is then assumed that the the first dofs are the displacements
!    dofs
! M. Ristinmaa 2012-01-13
!  - changed in defedof such that different dofs at nodes in one element
!    can be used
! M. Ristinmaa 2012-06-07
!  - Bugg in defedof1 argument enod switched
! M. Ristinmaa 2012-09-05
!  - Bugg in defedof2  indx started with 1, giving first noddof=2
! ------------------------------------------------------------------------------
!
implicit none

interface integ2nod
  module procedure integ2nod1
  module procedure integ2nod1v
  module procedure integ2nod2
  module procedure integ2nod2v
end interface

! getEleCoord
interface coordxtr
  module procedure coordxtr2d
  module procedure coordxtr3d
end interface
private coordxtr2d, coordxtr3d

interface defedof
  module procedure defedof1
  module procedure defedof2
end interface

interface finddof
  module procedure bcfinddof1
  module procedure bcfinddof2
end interface

integer :: debug=0

private debug

!private inv3

!----------------------------------------------------------------------

contains


! THIS SHOULD BE REMOVED
subroutine inv3(res,Ca)
  implicit none
  double precision                :: Ca(:,:), res(:,:)

  double precision                :: t4, t6, t8, t10, t12, t14, t17

  t4 = Ca(1,1)*Ca(2,2)
  t6 = Ca(1,1)*Ca(2,3)
  t8 = Ca(1,2)*Ca(2,1)
  t10 = Ca(1,3)*Ca(2,1)
  t12 = Ca(1,2)*Ca(3,1)
  t14 = Ca(1,3)*Ca(3,1)
  t17 = 1d0/(t4*Ca(3,3)-t6*Ca(3,2)-t8*Ca(3,3)+t10*Ca(3,2) &
         +t12*Ca(2,3)-t14*Ca(2,2))

  res(1,1) = (Ca(2,2)*Ca(3,3)-Ca(2,3)*Ca(3,2))*t17
  res(1,2) = -(Ca(1,2)*Ca(3,3)-Ca(1,3)*Ca(3,2))*t17
  res(1,3) = (Ca(1,2)*Ca(2,3)-Ca(1,3)*Ca(2,2))*t17
  res(2,1) = -(Ca(2,1)*Ca(3,3)-Ca(2,3)*Ca(3,1))*t17
  res(2,2) = (Ca(1,1)*Ca(3,3)-t14)*t17
  res(2,3) = -(t6-t10)*t17
  res(3,1) = (Ca(2,1)*Ca(3,2)-Ca(2,2)*Ca(3,1))*t17
  res(3,2) = -(Ca(1,1)*Ca(3,2)-t12)*t17
  res(3,3) = (t4-t8)*t17

  return
end subroutine inv3


! Special for elements with one integration point
subroutine integ2nod1v(SVM_Avarege,svm,eltype,enod)
  implicit none
   
  integer                         :: nelm, elnr, ii
  double precision                :: SVM_Avarege(:), SVM_index(size(SVM_Avarege))
  double precision                :: SVM(:)
  integer                         :: enod(:,:)
  character(len=*)                :: eltype

  double precision                :: resM(1,size(SVM_Avarege)), svmM(1,size(svm))

  resM(1,:)=SVM_Avarege
  svmM(1,:)=svm

  call integ2nod1(resM,svmM,eltype,enod)

  SVM_Avarege=resM(1,:)
  svm=svmM(1,:)

  return
end subroutine integ2nod1v
 

subroutine integ2nod1(SVM_Avarege,svm,eltype,enod)
  implicit none
   
  integer                         :: nelm, elnr, ii
  double precision                :: SVM_Avarege(:,:)
  double precision                :: SVM_index(size(SVM_Avarege,1),size(SVM_Avarege,2))
  double precision                :: SVM(:,:)
  integer                         :: enod(:,:)
  character(len=*)                :: eltype

  SVM_Avarege=0D0
  SVM_index=0D0 

  nelm=size(enod,2)
  if ((eltype.eq.'tr3').or.(eltype.eq.'brick8')) then
    do ii=1,size(SVM_Avarege,1)
      do elnr=1,nelm
        SVM_Avarege(ii,enod(:,elnr))=SVM_Avarege(ii,enod(:,elnr))+SVM(ii,elnr)
        SVM_index(ii,enod(:,elnr))=SVM_index(ii,enod(:,elnr))+1D0
      enddo
    enddo
  else
    stop 'element not implemented in integ2nod'
  endif

  do ii=1,size(SVM_Avarege,1)
    SVM_Avarege(ii,:)=SVM_Avarege(ii,:)/SVM_index(ii,:)
  enddo
  return

end subroutine integ2nod1


subroutine integ2nod2v(SVM_Avarege,svm,eltype,enod)
  implicit none
   
  integer                         :: nelm, elnr
  double precision                :: SVM_Avarege(:), SVM_index(size(SVM_Avarege))
  double precision                :: SVM(:,:)
  integer                         :: enod(:,:)
  character(len=*)                :: eltype

  SVM_Avarege=0D0
  SVM_index=0D0
 
  nelm=size(enod,2)
  if (eltype.eq.'qu4') then
    do elnr=1,nelm
      SVM_Avarege(enod(1,elnr))=SVM_Avarege(enod(1,elnr))+SVM(1,elnr)
      SVM_Avarege(enod(2,elnr))=SVM_Avarege(enod(2,elnr))+SVM(2,elnr)
      SVM_Avarege(enod(3,elnr))=SVM_Avarege(enod(3,elnr))+SVM(4,elnr) !This is correct
      SVM_Avarege(enod(4,elnr))=SVM_Avarege(enod(4,elnr))+SVM(3,elnr)
   
      SVM_index(enod(1,elnr))=SVM_index(enod(1,elnr))+1D0
      SVM_index(enod(2,elnr))=SVM_index(enod(2,elnr))+1D0
      SVM_index(enod(3,elnr))=SVM_index(enod(3,elnr))+1D0 
      SVM_index(enod(4,elnr))=SVM_index(enod(4,elnr))+1D0
    enddo
  elseif (eltype.eq.'qu8') then
! check number of integration points 9
    if (size(svm,1).ne.9) stop 'Only implemented for 9 integration points'
    do elnr=1,nelm
      SVM_Avarege(enod(1,elnr))=SVM_Avarege(enod(1,elnr))+SVM(1,elnr)
      SVM_Avarege(enod(2,elnr))=SVM_Avarege(enod(2,elnr))+SVM(3,elnr)
      SVM_Avarege(enod(3,elnr))=SVM_Avarege(enod(3,elnr))+SVM(9,elnr) 
      SVM_Avarege(enod(4,elnr))=SVM_Avarege(enod(4,elnr))+SVM(7,elnr)
      SVM_Avarege(enod(5,elnr))=SVM_Avarege(enod(5,elnr))+SVM(2,elnr)
      SVM_Avarege(enod(6,elnr))=SVM_Avarege(enod(6,elnr))+SVM(6,elnr)
      SVM_Avarege(enod(7,elnr))=SVM_Avarege(enod(7,elnr))+SVM(8,elnr) 
      SVM_Avarege(enod(8,elnr))=SVM_Avarege(enod(8,elnr))+SVM(4,elnr)
   
      SVM_index(enod(:,elnr))=SVM_index(enod(:,elnr))+1D0
    enddo
  elseif (eltype.eq.'brick8') then
    do elnr=1,nelm
      SVM_Avarege(enod(1,elnr))=SVM_Avarege(enod(1,elnr))+SVM(1,elnr)
      SVM_Avarege(enod(2,elnr))=SVM_Avarege(enod(2,elnr))+SVM(2,elnr)
      SVM_Avarege(enod(3,elnr))=SVM_Avarege(enod(3,elnr))+SVM(4,elnr) !This is correct
      SVM_Avarege(enod(4,elnr))=SVM_Avarege(enod(4,elnr))+SVM(3,elnr) !This is correct
      SVM_Avarege(enod(5,elnr))=SVM_Avarege(enod(5,elnr))+SVM(5,elnr)
      SVM_Avarege(enod(6,elnr))=SVM_Avarege(enod(6,elnr))+SVM(6,elnr)
      SVM_Avarege(enod(7,elnr))=SVM_Avarege(enod(7,elnr))+SVM(8,elnr) 
      SVM_Avarege(enod(8,elnr))=SVM_Avarege(enod(8,elnr))+SVM(7,elnr)
   
      SVM_index(enod(:,elnr))=SVM_index(enod(:,elnr))+1D0
    enddo
  else
    stop 'element not implemented in integ2nod'
  endif
  
  SVM_Avarege=SVM_Avarege/SVM_index
  return
end subroutine integ2nod2v


subroutine integ2nod2(SVM_Avarege,svm,eltype,enod)
  implicit none
   
  integer                         :: nelm, elnr, ii
  double precision                :: SVM_Avarege(:,:)
  double precision                :: SVM_index(size(SVM_Avarege,1),size(SVM_Avarege,2))
  double precision                :: SVM(:,:,:)
  integer                         :: enod(:,:)
  character(len=*)                :: eltype

  SVM_Avarege=0D0
  SVM_index=0D0
 
  nelm=size(enod,2)
  if (eltype.eq.'qu4') then
  do ii=1,size(SVM_Avarege,1)
    do elnr=1,nelm
      SVM_Avarege(ii,enod(1,elnr))=SVM_Avarege(ii,enod(1,elnr))+SVM(ii,1,elnr)
      SVM_Avarege(ii,enod(2,elnr))=SVM_Avarege(ii,enod(2,elnr))+SVM(ii,2,elnr)
      SVM_Avarege(ii,enod(3,elnr))=SVM_Avarege(ii,enod(3,elnr))+SVM(ii,4,elnr) !This is correct
      SVM_Avarege(ii,enod(4,elnr))=SVM_Avarege(ii,enod(4,elnr))+SVM(ii,3,elnr)
   
      SVM_index(ii,enod(1,elnr))=SVM_index(ii,enod(1,elnr))+1D0
      SVM_index(ii,enod(2,elnr))=SVM_index(ii,enod(2,elnr))+1D0
      SVM_index(ii,enod(3,elnr))=SVM_index(ii,enod(3,elnr))+1D0 
      SVM_index(ii,enod(4,elnr))=SVM_index(ii,enod(4,elnr))+1D0
    enddo
  enddo
  elseif (eltype.eq.'brick8') then
  do ii=1,size(SVM_Avarege,1)
    do elnr=1,nelm
      SVM_Avarege(ii,enod(1,elnr))=SVM_Avarege(ii,enod(1,elnr))+SVM(ii,1,elnr)
      SVM_Avarege(ii,enod(2,elnr))=SVM_Avarege(ii,enod(2,elnr))+SVM(ii,2,elnr)
      SVM_Avarege(ii,enod(3,elnr))=SVM_Avarege(ii,enod(3,elnr))+SVM(ii,4,elnr) !This is correct  
      SVM_Avarege(ii,enod(4,elnr))=SVM_Avarege(ii,enod(4,elnr))+SVM(ii,3,elnr) !This is correct
      SVM_Avarege(ii,enod(5,elnr))=SVM_Avarege(ii,enod(5,elnr))+SVM(ii,5,elnr)
      SVM_Avarege(ii,enod(6,elnr))=SVM_Avarege(ii,enod(6,elnr))+SVM(ii,6,elnr)
      SVM_Avarege(ii,enod(7,elnr))=SVM_Avarege(ii,enod(7,elnr))+SVM(ii,8,elnr) 
      SVM_Avarege(ii,enod(8,elnr))=SVM_Avarege(ii,enod(8,elnr))+SVM(ii,7,elnr)
   
      SVM_index(ii,enod(:,elnr))=SVM_index(ii,enod(:,elnr))+1D0
    enddo
  enddo
  else
    stop 'element not implemented in integ2nod'
  endif
  
  SVM_Avarege=SVM_Avarege/SVM_index
  return
end subroutine integ2nod2


subroutine coordxtr2d(ex,ey,coord,enod)
implicit none
  double precision    :: ex(:,:), ey(:,:), coord(:,:)
  integer             :: enod(:,:)

  integer             :: ie, nele, nod(size(enod,1))

  nele=size(enod,2)

  do ie=1,nele
    nod=enod(:,ie)
    ex(:,ie)=coord(1,nod)
    ey(:,ie)=coord(2,nod)
  end do

end subroutine coordxtr2d

subroutine coordxtr3d(ex,ey,ez,coord,enod)
implicit none
  double precision    :: ex(:,:), ey(:,:), ez(:,:), coord(:,:)
  integer             :: enod(:,:)

  integer             :: ie, nele, nod(size(enod,1))

  nele=size(enod,2)

  do ie=1,nele
    nod=enod(:,ie)
    ex(:,ie)=coord(1,nod)
    ey(:,ie)=coord(2,nod)
    ez(:,ie)=coord(3,nod)
  end do

end subroutine coordxtr3d
  

!function norm(a)
!  implicit none
!  double precision   :: a(:), norm
!  
!  norm=dsqrt(dot_product(a,a))
!end function norm

subroutine updcoord(newcoord,coord,a)
  implicit none
  double precision                :: newcoord(:,:),coord(:,:), a(:)
  integer                         :: dofnod !enod(:,:)

  double precision, allocatable   :: ar(:,:)
  integer                         :: ndme, dfnd, dime, ierr

  ndme=size(coord,2) !maxval(enod)
  dime=size(coord,1)
  dfnd=size(a)/ndme
  
  allocate(ar(dfnd,ndme), stat=ierr)
  ar=reshape(a,[dfnd,ndme],order=[1,2])

  newcoord=coord+ar(1:dime,:)

  return
end subroutine updcoord


subroutine defedof1(edof,enod,dofnod)
  implicit none
  integer                         :: edof(:,:), enod(:,:), dofnod
  integer                         :: nelm, noel
  integer                         :: ie, ii, tmp1, tmp2

  nelm=size(enod,2)
  noel=size(enod,1)
  do ie=1,nelm
    do ii=1,noel
      tmp1=dofnod*(ii-1)+1
      tmp2=dofnod*(enod(ii,ie)-1)+1
      edof(tmp1:tmp1+dofnod-1,ie)=(/tmp2:(tmp2+dofnod-1)/)
    enddo
  enddo

  return
end subroutine defedof1

 
! Used when the number of dofs in the nodes differ in the element
subroutine defedof2(edof,noddof,enod,dofnod)
  implicit none
  integer                         :: edof(:,:), enod(:,:), dofnod(:)

  integer                         :: nrdof, nrnod, mxdof, dofs, nrelm, ndelm
  integer                         :: ierr, ind, iel, i, indx, nodnr

  integer, allocatable            :: dofsnod(:), noddof(:,:)

  nrnod=maxval(enod)
  nrdof=sum(dofnod)
  nrelm=size(enod,2)
  ndelm=size(enod,1)
  mxdof=maxval(dofnod)

! setup list node number and number of dofs exist at this node
  allocate(dofsnod(nrnod),stat=ierr)
  do iel=1,nrelm
    do ind=1,ndelm
      dofsnod(enod(ind,iel))=dofnod(ind)
    end do
  end do

! steup list global dofs beloning to a node
!  allocate(noddof(mxdof,nrnod),stat=ierr)
  noddof=0
!  indx=1
  indx=0
  do ind=1,nrnod
    dofs=dofsnod(ind)
    do i=1,dofs
        noddof(i,ind)=indx+i
    end do
    indx=indx+dofs
  end do

! defined the topology in terms of global dofs
  do iel=1,nrelm
    indx=0
    do ind=1,ndelm
      nodnr=enod(ind,iel)
      dofs=dofsnod(nodnr)
      do i=1,dofs
        edof(indx+i,iel)=noddof(i,nodnr)
      end do
      indx=indx+dofs
    end do
  end do

  deallocate(dofsnod)
!  deallocate(noddof)

  return
end subroutine defedof2


subroutine defenod(enod,edof,dofnod)
  implicit none
  integer                :: enod(:,:), edof(:,:), dofnod

  integer                :: ndel, nrelm, ie, ii, tmp


! expected number of dof in a node nrdof
! expected number of nodes in one element  
 
  ndel=size(edof,1)/dofnod
  nrelm=size(edof,2)

  do ie=1,nrelm
    do ii=1,ndel
      tmp=dofnod*(ii-1)+1
      enod(ii,ie)=(edof(tmp,ie)-1)/dofnod+1
    enddo
  enddo

end subroutine defenod

subroutine bcfinddof1(bcdof,bcnod,dofnod)
  implicit none
  integer                         :: bcdof(:), bcnod(:,:), dofnod
  integer                         :: ii

  do ii=1,size(bcdof)
     bcdof(ii)=dofnod*(bcnod(ii,1)-1)+bcnod(ii,2)
  enddo

  return
end subroutine bcfinddof1

subroutine bcfinddof2(bcdof,bcnode,dof)
  implicit none
  integer                :: bcnode(:,:), bcdof(:), dof(:,:)

  integer                :: i1

  do i1=1,size(bcnode,1)
    bcdof(i1)=dof(bcnode(i1,2),bcnode(i1,1))
  enddo

  return
end subroutine bcfinddof2


subroutine bcfindnode(bcnode,bcdof,dof)
  implicit none
  integer                :: bcnode(:,:), bcdof(:), dof(:,:)

  integer                :: i1, i2, i3

  do i2=1,size(bcdof)
    do i1=1,size(dof,2)
      do i3=1,size(dof,1)
        if (bcdof(i2).eq.dof(i3,i1)) then
          bcnode(i2,:)=(/i3, i1/)
        endif
      enddo
    enddo
  enddo

  return
end subroutine bcfindnode

subroutine neighbel(lele,side,enod,iel)
  implicit none
  integer               :: enod(:,:)
  integer               :: lele(3), side(3), iel
  
  integer               :: n(4), is, ie, current_nr_elements

  if (debug.eq.1) write(*,*)'---------------entering neighbel '
  
  if (size(enod,1).ne.3) stop 'neighbel only implemented for 3node ele'

! find the neigbours to element iel
! can be maximum 3 elements
! 2 is connected to external boundary
! this is indicated with a 0 in lele
!
! lele is the list of elements
! side is the side on which these elements
! are connected to element iel

 current_nr_elements=size(enod,2)  

  n(1:3)=enod(1:3,iel)
  n(4)=n(1)
  lele=0
  side=0   
 
! loop over element sides
  do is=1,3

! find element connected to nodes
! on one side

! SHOULD BE REWRITTEN FOR SPEED

    do ie=1,current_nr_elements

      if ((enod(1,ie).eq.n(is+1)).and.(enod(2,ie).eq.n(is))) then
        lele(is)=ie
        side(is)=1
      elseif ((enod(2,ie).eq.n(is+1)).and.(enod(3,ie).eq.n(is))) then
        lele(is)=ie
        side(is)=2
      elseif ((enod(3,ie).eq.n(is+1)).and.(enod(1,ie).eq.n(is))) then
        lele(is)=ie
        side(is)=3
      endif
     
    end do
  end do

  if (debug.eq.1) write(*,*)'---------------leaving  neighbel'
  return
end subroutine neighbel


end module fem_util
