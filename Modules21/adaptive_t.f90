module ADAPTIVE

! last modified
! M. Ristinmaa 2011-01-20
!   -Initial version
! M. Ristinmaa 2011-03-27
!   -Should work for general initial unstructed meshes
! M. Ristinmaa 2011-03-30
!   -Fixed bug in create_element, level incorrectly calculated in elorg
! M. Ristinmaa 2011-04-07
!   -Fixed bug in create-coord2, split can not be done at the called 
!    element side in a recursive call.  
! M. Ristinmaa 2011-04-09
!   -Added a subroutine for releasing internal memory
! M. Ristinmaa 2011-04-12
!   -Move not related subroutines to fem_util	
! M. Ristinmaa 2011-04-13
!   -Move not related subroutines to memory_util	
! M. Ristinmaa 2011-04-14
!   -Fixed boundary not divided, due to fix 2011-04-07
! M. Ristinmaa 2011-04-21
!   -Renamed adaptive_mesh to adapMesh and adaptive_close adapClose
! M. Ristinmaa 2011-04-28
!   -Moved neigbel to fem_util	
!----------------------------------------------------------------------

! adaptivity routines for 3-noded elements
!
! In the initial topology is used as the base this can not be 
! made more course. The initial topology can be refined as well
! as added elements can be remove.
!
! The module is nodal based.
!
! -External load vector is not changed due to distributed loading
!
! *profil should be done to speedup the routines
! *when splitting elements, distorted elements should be splitted first
!  will give a better mesh. 


use memory_util	
use fem_util
	
implicit none

integer, allocatable             :: elorg(:,:)
integer, allocatable             :: elloc(:)
integer                          :: current_nr_nodes, max_nr_nodes, alloc_calls
integer                          :: current_nr_elements, max_nr_elements
integer                          :: added_nr_elements
private elorg, elloc
private current_nr_nodes, max_nr_nodes, alloc_calls
private current_nr_elements, max_nr_elements, added_nr_elements

! these parameters are used for allocating memory for elements
! that might be created
integer, parameter               :: inc_nodes=1000
integer, parameter               :: inc_elements=2000
integer, parameter               :: max_alloc_calls=6
private inc_nodes, inc_elements, max_alloc_calls

! debug purposes
integer    :: debug=0
private debug


private refineel
private create_element
!private neighbel
private createcoord2
private adaptive_remove
private reduce_arrays
private removeele
private remove1

!------------------------------------------------------------------------------

contains


subroutine adapClose()

 if (allocated(elorg)) deallocate(elorg)
 if (allocated(elloc)) deallocate(elloc)
 
 return 
end subroutine adapClose


subroutine adapMesh(rmv,enodx,coordx,a,fext,bcnode,bcval,adaplevel)
  implicit none
  integer, allocatable           :: enodx(:,:)
  double precision, allocatable  :: coordx(:,:)

  integer, allocatable           :: enod(:,:), bcnode(:,:)
  double precision, allocatable  :: coord(:,:), rmv(:), a(:), fext(:), bcval(:)
  integer                        :: adaplevel

  double precision, allocatable  :: ar(:,:), fr(:,:)
  integer                        :: status, nele
  integer                        :: dfnd, ndme, ndof

  if (debug.eq.1) write(*,*)'---------------entering adaptive_mesh'

  nele=size(enod,1)
!  ndel=size(edof,2)
  ndme=maxval(enod)
  dfnd=size(a)/ndme

! SHOULD BE REMOVED LATER ON
  allocate(enod(size(enodx,2),size(enodx,1)),stat=status)
  allocate(coord(size(coordx,2),size(coordx,1)),stat=status)
  enod=transpose(enodx)
  coord=transpose(coordx)

  alloc_calls=0

  if (.not.(allocated(elorg))) then
    allocate(elorg(inc_elements+1,4),stat=status)
    elorg=0
    added_nr_elements=1  ! use position 1 for all initially existing elements
  end if

! allocate extra space such that the values belonging to new
! nodes can be stored
	
  current_nr_nodes=ndme
  max_nr_nodes=ndme+inc_nodes
	
  allocate(ar(ndme,dfnd),stat=status)
  allocate(fr(ndme,dfnd),stat=status)
  ar=reshape(a,[ndme,dfnd],order=[2,1])
  fr=reshape(fext,[ndme,dfnd],order=[2,1])

  call reallocate(ar,max_nr_nodes)
  call reallocate(fr,max_nr_nodes)
  call reallocate(coord,max_nr_nodes)
  
! allocate extra space such that new elements can be stored

  current_nr_elements=nele
  max_nr_elements=nele+inc_elements
  if (.not.(allocated(elloc))) then
    allocate(elloc(max_nr_elements),stat=status)
    elloc=0
  end if

  call reallocate(enod,max_nr_elements)
  call reallocate(elorg,added_nr_elements+inc_elements)
  call reallocate(rmv,max_nr_elements)
  call reallocate(elloc,max_nr_elements)
	
! refine mesh based on rmv

  call adaptive_add(coord,ar,enod,fr,bcnode,bcval,rmv,adaplevel)

! remove elements based on rmv

  call adaptive_remove(coord,ar,enod,fr,bcnode,bcval,rmv)

  ndof=current_nr_nodes*dfnd
  deallocate(a)
  allocate(a(ndof),stat=status)
  deallocate(fext)
  allocate(fext(ndof),stat=status)

  call reallocate(ar,current_nr_nodes)
  call reallocate(fr,current_nr_nodes)
  call reallocate(coord,current_nr_nodes)

! matrix storage, convert vector to matrix storage
  a=pack(transpose(ar),mask=.true.)
  fext=pack(transpose(fr), mask=.true.)

  call reallocate(enod,current_nr_elements)
  call reallocate(elorg,added_nr_elements)
  call reallocate(rmv,current_nr_elements)
  call reallocate(elloc,current_nr_elements)

!  nele=size(enod,1)
!  ndel=size(enod,2)
!  nbca=size(bcnode,1)

  deallocate(ar)
  deallocate(fr)   

! SHOULD BE REMOVED LATER
  deallocate(enodx,coordx)
  allocate(enodx(size(enod,2),size(enod,1)),stat=status)
  allocate(coordx(size(coord,2),size(coord,1)),stat=status)
  enodx=transpose(enod)
  coordx=transpose(coord)

  if (debug.eq.1) write(*,*)'---------------leaving  adaptive_mesh'
  return
end subroutine adapMesh


subroutine adaptive_add(coord,ar,enod,fr,bcnode,bcval,eta,adaplevel)
  implicit none
  integer, allocatable           :: enod(:,:), bcnode(:,:)
  double precision, allocatable  :: coord(:,:), ar(:,:), bcval(:), fr(:,:), eta(:)
  integer                        :: adaplevel
  
  integer                        ::  iel(1), iban

  if (debug.eq.1) write(*,*)'---------------entering adaptive_add'

  do while (maxval(eta).gt.0.5d0)
    iel=maxloc(eta)     
    if (elorg(iel(1),1).lt.adaplevel) then
      if (debug.eq.1) write(*,*)'refine element ',iel(1)
       iban=0
      call refineel(iel(1),enod,coord,eta,ar,fr,bcnode,bcval,iban)
    endif
    eta(iel)=0.25d0
  end do

  if (debug.eq.1) write(*,*)'---------------leaving  adaptive_add'
  return
end subroutine adaptive_add


recursive subroutine refineel(iel,enod,coord,eta,ar,fr,bcnode,bcval,iban)
  implicit none
  integer,allocatable               :: enod(:,:), bcnode(:,:)
  double precision, allocatable     :: coord(:,:), eta(:), ar(:,:), fr(:,:), bcval(:) 
  integer                           :: iel, iban

! For refinement it splits the current element on largest 
! element side

  integer  :: maxnod, iel2, isb, isb2, loc_iel, loc_iel2
  integer  :: lele(3), nside(3), split, tmp(1)

  if (debug.eq.1) write(*,*)'---------------entering refineel'
  if (debug.eq.1) write(*,*)'element ',iel
! define coordinates at the midpoint at the
! longest element side and define dofs at node
! bcold=bc;

  call createcoord2(isb,coord,enod,iel,ar,fr,bcval,bcnode,iban)
  maxnod=current_nr_nodes

! return new coordinates, isb is the longest side on the element

! find neigbouring element numbers
! and detect if element is part of 
! external boundary

  call neighbel(lele,nside,transpose(enod),iel)

! if side to split is part of external boundary
! we don't need to consider neighbour elements
! use that the minimum is zero in lele and find this
! when outer boundary is present

   split=1
   tmp=minloc(lele) 
   if ((tmp(1).eq.isb).and.(lele(tmp(1)).eq.0)) then
     split=-1
   end if

  if (split.eq.1) then

    iel2=lele(isb)

! check level of refinement for element
! must be of the same order
    loc_iel=elloc(iel)
    if (loc_iel.eq.0) loc_iel=1 
    loc_iel2=elloc(iel2)
    if (loc_iel2.eq.0) loc_iel2=1 

    if ((elorg(loc_iel,1)-(elorg(loc_iel2,1))).gt.0) then
      call refineel(iel2,enod,coord,eta,ar,fr,bcnode,bcval,iel)
    elseif ((elorg(loc_iel,1)-(elorg(loc_iel2,1))).gt.0) then
      write(*,*)' This should not happen'
      stop
    end if
 
! split neighbouring element

    call neighbel(lele,nside,transpose(enod),iel)
    iel2=lele(isb)
    isb2=nside(isb)

! create element

    call create_element(enod,eta,iel2,isb2,maxnod,iel,coord)
  else
    iel2=0
  end if

  call create_element(enod,eta,iel,isb,maxnod,iel2,coord)

  if (debug.eq.1) write(*,*)'element ',iel
  if (debug.eq.1) write(*,*)'---------------leaving  refineel'
  return
end subroutine refineel


subroutine create_element(enod,eta,ele,side,maxnod,iel2,coord)  
  implicit none
  integer, allocatable   :: enod(:,:)
  double precision, allocatable :: eta(:), coord(:,:)
  integer  :: ele,side, maxnod, iel2

  integer  :: n(5), ie, level, ele_min, i1, ie1
  double precision :: ex(3), ey(3), deta

  if (debug.eq.1) write(*,*)'---------------entering create_element', ele, iel2

  current_nr_elements=current_nr_elements+1
 
  if (current_nr_elements.eq.max_nr_elements) then
!    write(*,*)'allocating more space for elements', current_nr_elements
    max_nr_elements=current_nr_elements+inc_elements
    call reallocate(enod,max_nr_elements)
    call reallocate(elloc,max_nr_elements)
    call reallocate(elorg,added_nr_elements+inc_elements)
    call reallocate(eta,max_nr_elements)
    alloc_calls=alloc_calls+1
    if (alloc_calls.eq.max_alloc_calls) then
      write(*,*)'Speedup adaptivity by increasing inc_elemens and inc_nodes'
    endif
  end if

! nodes defining element

  n(1:3)=enod(ele,1:3)
  n(4)=n(1)
  n(5)=n(2)

  enod(ele,:)=(/n(side), maxnod, n(side+2) /)
  enod(current_nr_elements,:)=(/maxnod, n(side+1), n(side+2)/)

  ele_min=ele
  ex=(/coord(enod(ele_min,1),1),coord(enod(ele_min,2),1),coord(enod(ele_min,3),1) /) 
  ey=(/coord(enod(ele_min,1),2),coord(enod(ele_min,2),2),coord(enod(ele_min,3),2) /) 
  DetA=ex(2)*ey(3)+ex(1)*ey(2)+ey(1)*ex(3) &
         -ex(2)*ey(1)-ex(3)*ey(2)-ey(3)*ex(1)
  if (DetA.le.0d0) then
    write(*,*)' Area negative'
    write(*,*)'enod ',enod(ele_min,:)
    write(*,*)'coord1 ',coord(enod(ele_min,1),1:2)
    write(*,*)'coord2 ',coord(enod(ele_min,2),1:2)
    write(*,*)'coord3 ',coord(enod(ele_min,3),1:2)
    write(*,*)'n ',n
    do i1=1,3
      write(*,*)'n ',n(i1),' coord ',coord(n(i1),1:2)
    enddo
    write(*,*)'maxnod ',maxnod,' coord ',coord(maxnod,1:2)
!    enod(ele,:)=(/maxnod, n(side), n(side+2) /)
    stop
  endif

  ele_min=current_nr_elements
  ex=(/coord(enod(ele_min,1),1),coord(enod(ele_min,2),1),coord(enod(ele_min,3),1) /) 
  ey=(/coord(enod(ele_min,1),2),coord(enod(ele_min,2),2),coord(enod(ele_min,3),2) /) 
  DetA=ex(2)*ey(3)+ex(1)*ey(2)+ey(1)*ex(3) &
         -ex(2)*ey(1)-ex(3)*ey(2)-ey(3)*ex(1)
  if (DetA.le.0d0) then
    write(*,*)' Area negative'
    write(*,*)'enod ',enod(ele_min,:)
    write(*,*)'coord1 ',coord(enod(ele_min,1),:)
    write(*,*)'coord2 ',coord(enod(ele_min,2),:)
    write(*,*)'coord3 ',coord(enod(ele_min,3),:)
    write(*,*)'n ',n
    do i1=1,3
      write(*,*)'n ',n(i1),' coord ',coord(n(i1),1:2)
    enddo
    write(*,*)'maxnod ',maxnod,' coord ',coord(maxnod,1:2)
 !   enod(current_nr_elements,:)=(/n(side+1), maxnod, n(side+2)/)
    stop
  endif

  if (iel2.eq.0) then
    ie=added_nr_elements
    do while ((elorg(ie,2).ne.ele).and.(elorg(ie,3).ne.ele).and.(ie.ne.1))
      ie=ie-1
    end do
    elloc(ele)=added_nr_elements+1
    elloc(current_nr_elements)=added_nr_elements+1
    added_nr_elements=added_nr_elements+1

    if (ie.eq.1) then
      level=1
    else
      level=elorg(ie,1)+1
    end if

  else

    ie=added_nr_elements
    do while ((elorg(ie,2).ne.ele).and.(elorg(ie,3).ne.ele).and.(ie.ne.1))
      ie=ie-1
    end do
    ie1=added_nr_elements
    do while ((elorg(ie1,4).ne.ele).and.(ie1.ne.1))
      ie1=ie1-1
    end do
    ie=maxval((/ie,ie1/))

    elloc(ele)=added_nr_elements+1
    elloc(current_nr_elements)=added_nr_elements+1
    added_nr_elements=added_nr_elements+1

    if (ie.eq.1) then
      level=1
    else
      if (elorg(ie,2).eq.iel2) then
        level=elorg(ie,1)
      else
        level=elorg(ie,1)+1
      endif
    end if
  end if

  elorg(added_nr_elements,:)=(/ level, ele, current_nr_elements, iel2 /)
	
! change eta values for element such that they 
! are not active in the adaptivity process any more

  eta(ele)=0.25d0
  eta(current_nr_elements)=0.25d0

  if (debug.eq.1) write(*,*)'---------------leaving  create_element'
  return
end subroutine create_element



subroutine createcoord2(isb,coord,enod,iel,ar,fr,bcval,bcnode,iban)
  implicit none
  integer, allocatable          :: enod(:,:), bcnode(:,:)
  double precision, allocatable :: coord(:,:), bcval(:), ar(:,:), fr(:,:)
  integer                       :: isb,iel,iban

  integer      :: n(3), ierr, side(1), n1, n2, ic, i1,i2, nbc, lele(3), dummy(3)
  double precision :: ex(3), ey(3), l1, l2 ,l3, xc(3), yc(3), lv(3)
  double precision, allocatable :: ea(:,:), ec(:,:)

  if (debug.eq.1) write(*,*)'---------------entering createcoord2', iban
! Create coordinate at the midpoint on the
! longest element side

! nodes defining element
  n=enod(iel,:)

! coordinates for nodes in the element

  ex=(/coord(n(1),1), coord(n(2),1), coord(n(3),1)/)
  ey=(/coord(n(1),2), coord(n(2),2), coord(n(3),2)/)

! nodal field values for the three nodes

  allocate(ea(3,size(ar,2)),stat=ierr)
  ea(1,:)=ar(n(1),:)
  ea(2,:)=ar(n(2),:)
  ea(3,:)=ar(n(3),:)

! split boundary connecting the elements on the longest element side 

  l1=(ex(1)-ex(2))**2d0+(ey(1)-ey(2))**2d0
  l2=(ex(2)-ex(3))**2d0+(ey(2)-ey(3))**2d0
  l3=(ex(3)-ex(1))**2d0+(ey(3)-ey(1))**2d0
  lv=(/l1,l2,l3/)
  side=maxloc(lv)
  isb=side(1)

  call neighbel(lele,dummy,transpose(enod),iel)
  if ((lele(isb).eq.iban).and.(iban.ne.0)) then
    lv(isb)=-lv(isb)
    side=maxloc(lv)
    isb=side(1)
  endif

! create new coord/node

  xc=[ex(1)+ex(2),ex(2)+ex(3),ex(3)+ex(1)]*0.5d0
  yc=[ey(1)+ey(2),ey(2)+ey(3),ey(3)+ey(1)]*0.5d0

  current_nr_nodes=current_nr_nodes+1
!
! check if more space is needed
  if (current_nr_nodes.eq.max_nr_nodes) then
!    write(*,*)'allocating more space for nodes'
    max_nr_nodes=current_nr_nodes+inc_nodes
    call reallocate(ar,max_nr_nodes)
    call reallocate(fr,max_nr_nodes)
    call reallocate(coord,max_nr_nodes)
    alloc_calls=alloc_calls+1
    if (alloc_calls.eq.max_alloc_calls) then
      write(*,*)'Speedup adaptivity by increasing inc_elemens and inc_nodes'
    endif
  end if

  if (size(coord,2).eq.3) then
    coord(current_nr_nodes,:)=(/xc(side(1)), yc(side(1)), 0/)
  else
    coord(current_nr_nodes,:)=(/xc(side(1)), yc(side(1))/)
  end if
    
    
! boundary conditions  
!
 if (side(1).eq.1) then
   n1=n(1)
   n2=n(2)
 elseif( side(1).eq.2) then
   n1=n(2)
   n2=n(3)
 else
   n1=n(3)
   n2=n(1)
 endif	

!
! if both nodes are specified in bcnode
! then the created node should also be in
! bcnode	
!
! note a node can exist more than one time
! in bcnode since it appears for every dof that
! is prescribed
!	
 nbc=size(bcval)	
 ic=0	
 do i1=1,nbc
   if (n1.eq.bcnode(i1,1)) then
     do i2=1,nbc
       if ((n2.eq.bcnode(i2,1)).and.(bcnode(i1,2).eq.bcnode(i2,2))) then
         ic=ic+1
         call reallocate(bcnode,nbc+ic)
         call reallocate(bcval,nbc+ic)
         bcnode(nbc+ic,:)=(/current_nr_nodes, bcnode(i2,2)/)
         bcval(nbc+ic)=(bcval(i1)+bcval(i2))*0.5d0	
!         write(*,*)'new bc ',bcnode(nbc+ic,:),bcval(nbc+ic)
       endif
     enddo
   endif	
 enddo
!	
!	
! create new field variables and set the values 
! to the mean value between the two nodes

  allocate(ec(3,size(ar,2)),stat=ierr)
  ec(1,:)=(ea(1,:)+ea(2,:))*0.5d0
  ec(2,:)=(ea(2,:)+ea(3,:))*0.5d0
  ec(3,:)=(ea(3,:)+ea(1,:))*0.5d0

  ar(current_nr_nodes,:)=ec(side(1),:)
  fr(current_nr_nodes,:)=0

  deallocate(ea)
  deallocate(ec)
  
  if (debug.eq.1) write(*,*)'---------------leaving  createcoord2'
  return
end subroutine createcoord2


subroutine adaptive_remove(coord,ar,enod,fr,bcnode,bcval,eta)
  implicit none
  integer, allocatable           :: enod(:,:), bcnode(:,:)
  double precision, allocatable  :: coord(:,:), ar(:,:), bcval(:), fr(:,:), eta(:)
  
  integer                        ::  ie
	
  if (debug.eq.1) write(*,*)'---------------entering adaptive_remove'

!
! sort such that elements that are most refined is first
! check if they can be removed
!
  do ie=current_nr_elements,1,-1
    if (eta(ie).lt.0) then
      call removeele(ie,enod,eta,coord)
    end if
  end do
! remove elements and nods from array lists

  call reduce_arrays(coord,ar,enod,fr,bcnode,bcval,eta)

  if (debug.eq.1) write(*,*)'---------------leaving  adaptive_remove'
  return
end subroutine adaptive_remove


subroutine reduce_arrays(coord,ar,enod,fr,bcnode,bcval,eta)
  implicit none
  integer, allocatable           :: enod(:,:), bcnode(:,:)
  double precision, allocatable  :: coord(:,:), ar(:,:), bcval(:), fr(:,:), eta(:)

  integer :: nele, ie, inc, in, nnod, ierr, nods(3), tmp, nbc
  integer, allocatable :: nodlist(:), elelist(:), nodmap(:), elemap(:)

  if (debug.eq.1) write(*,*)'---------------entering reduce_arrays'
!
! element numbers will change
!  
 allocate(elelist(current_nr_elements),stat=ierr)
 nele=current_nr_elements	
  do ie=current_nr_elements,1,-1
    elelist(ie)=1
    if (enod(ie,1).eq.-1) then
      enod(ie:nele-1,:)=enod(ie+1:nele,:)
      nele=nele-1
      elelist(ie)=0
    end if
  end do

  allocate(elemap(current_nr_elements),stat=ierr)
  elemap=0
  inc=0
  do ie=1,current_nr_elements
    if (elelist(ie).ne.0) then  
      inc=inc+1
      elemap(ie)=inc
    end if
  end do

  current_nr_elements=nele	

  nele=added_nr_elements
  do ie=added_nr_elements,2,-1
    if (elorg(ie,1).eq.-1) then
      elorg(ie:nele-1,:)=elorg(ie+1:nele,:)
      eta(ie:nele-1)=eta(ie+1:nele)
      nele=nele-1
    end if
  end do
  added_nr_elements=nele

  do ie=2,added_nr_elements
    if (elorg(ie,4).eq.0) then
      tmp=0
    else
      tmp=elemap(elorg(ie,4))
    end if
    elorg(ie,2:4)=(/elemap(elorg(ie,2)),elemap(elorg(ie,3)),tmp /)
    elloc(elorg(ie,2))=ie
    elloc(elorg(ie,3))=ie
  end do 

! remove nodes from list
! note that enod must be modified
! first check which nodes are remaining
  allocate(nodlist(current_nr_nodes),stat=ierr)
  nodlist=0
  do ie=1,current_nr_elements
    nods=enod(ie,:)
    nodlist(nods)=1
  end do
!
! maybe best solution is to make a mapping vector
! and then loop over all elements at once
!write(*,*)'nodemap'
  allocate(nodmap(current_nr_nodes),stat=ierr)
  nodmap=0
  inc=0
  do in=1,current_nr_nodes
    if (nodlist(in).ne.0) then  
      inc=inc+1
      nodmap(in)=inc
    end if
  end do
!
! loop over elements and map nodes
!write(*,*)'nodemapping'

  do ie=1,current_nr_elements
    enod(ie,1)=nodmap(enod(ie,1))
    enod(ie,2)=nodmap(enod(ie,2))
    enod(ie,3)=nodmap(enod(ie,3))
  end do
!
! finally coord must be changed
  nnod=current_nr_nodes
  do in=current_nr_nodes,1,-1
    if (nodmap(in).eq.0) then
      coord(in:nnod-1,:)=coord(in+1:nnod,:)
      ar(in:nnod-1,:)=ar(in+1:nnod,:)
      fr(in:nnod-1,:)=fr(in+1:nnod,:)
      nnod=nnod-1
    end if
  end do
  current_nr_nodes=nnod

! change boundary conditions
!
  nbc=size(bcval)
  do in=1,nbc   
      bcnode(in,1)=nodmap(bcnode(in,1))
  enddo	
  nnod=nbc	
  do in=nbc,1,-1
    if (bcnode(in,1).eq.0) then
      if (in.ne.nbc) then
        bcnode(in:nnod-1,:)=bcnode(in+1:nnod,:)
        bcval(in:nnod-1)=bcval(in+1:nnod)
      endif	
      nnod=nnod-1	
    endif
  enddo	
  call reallocate(bcnode,nnod)
  call reallocate(bcval,nnod)	
!	
 if (debug.eq.1) then
!  do ie=1,current_nr_nodes
!    write(*,*)'ar ', ie,ar(ie,:)
!  end do
!  do ie=1,current_nr_nodes
!    write(*,*)'fr ', ie,fr(ie,:)
!  end do
  do ie=1,added_nr_elements
    write(*,*)'elorg ', ie,elorg(ie,:)
  end do
!
  do ie=1,current_nr_elements
    write(*,*)'eloc eta ', ie,elloc(ie), eta(ie)
  end do
!  
  do ie=1,current_nr_elements
    write(*,*)'enod ', ie,enod(ie,:)
  end do 
  end if
  
  deallocate(elelist)
  deallocate(elemap)
  deallocate(nodlist)
  deallocate(nodmap)
  
  if (debug.eq.1) write(*,*)'---------------leaving reduce_arrays'
  return
end subroutine reduce_arrays


subroutine removeele(ele,enod,eta,coord)
  implicit none
  integer, allocatable           :: enod(:,:)
  double precision, allocatable  :: eta(:), coord(:,:)
  integer                        :: ele
 
  integer  :: ele_n1, ele_n2, loc_n1, loc_n2, ie, ele1, ele2, iec

  if (debug.eq.1) write(*,*)'---------------entering removeel'

!
! find last entry of ele in elorg
! 
  ie=elloc(ele)

! not allowed to modify the original mesh
  if (ie.ne.0) then

! check if brother can be removed
! and that both elements are at the same refinement level

    ele1=elorg(ie,2)
    ele2=elorg(ie,3)

    if ((eta(ele2).le.0d0).and.(elloc(ele1).eq.elloc(ele2))) then
    
! check if the split is on a line which is part of outer boundary
      if (elorg(ie,4).eq.0) then
        call remove1(enod,eta,ele1,ele2,ie,coord)
        elloc=0
        do iec=2,added_nr_elements
          if (elorg(iec,1).ne.-1) then
            elloc(elorg(iec,2))=iec
            elloc(elorg(iec,3))=iec
          end if
        end do 
      else
      	      
! check if neighbouring elements are on the same refinement level  
! if not no removal is possible  	

        ele_n1=elorg(ie,4)
        if (elorg(ie+1,2).eq.ele_n1) then
          loc_n1=ie+1
          ele_n1=elorg(loc_n1,2)
          ele_n2=elorg(loc_n1,3)
        elseif (elorg(ie-1,2).eq.ele_n1) then
          loc_n1=ie-1
          ele_n1=elorg(loc_n1,2)
          ele_n2=elorg(loc_n1,3)
        else
          write(*,*)'something is wrong'
          write(*,*)'elorg 1 ',ie,elorg(ie,:)
          write(*,*)'elorg 2 ',loc_n1,elorg(loc_n1,:)
          stop
        endif
        loc_n2=maxval((/elloc(ele_n2),elloc(ele_n1)/))

! check that the elements are allowed to be removed before continuing
        if ((eta(ele_n1).le.0d0).and.(eta(ele_n2).le.0d0)) then
        	
! check that all elements are at the same level of refinement
	
          if ((elorg(ie,1).eq.elorg(loc_n1,1)).and.(loc_n1.eq.loc_n2)) then
            call remove1(enod,eta,ele1,ele2,ie,coord)
            call remove1(enod,eta,ele_n1,ele_n2,loc_n1,coord)
            elloc=0
            do iec=2,added_nr_elements
              if (elorg(iec,1).ne.-1) then
                elloc(elorg(iec,2))=iec
                elloc(elorg(iec,3))=iec
              end if
            end do 
          end if   
     	
        end if 
      end if
    end if
  end if
  eta(ele)=0d0

  if (debug.eq.1) write(*,*)'---------------leaving  removeel'
  return
end subroutine removeele


subroutine remove1(enod,eta,ele,ele_b,elr,coord)
  implicit none
  integer, allocatable          :: enod(:,:)
  double precision, allocatable :: eta(:), coord(:,:)
  integer                       :: ele, ele_b, elr

  integer                :: ele_min, ele_max, ie

  integer :: i1,i2,i3,n1,n2, cnod(2), pos_max(2), pos_min(2)
  double precision :: xm, ym, d1, d2, DetA, ex(3), ey(3)

  if (debug.eq.1) write(*,*)'---------------entering remove1'
 
! tag that the elements have been considered, note
! that a zero means that the mesh can be further coarsend
  eta(ele)=0.25d0
  eta(ele_b)=0.25d0

! mark for removal, only one element will be removed, one is overwritten
! store element with lowest element number
  ele_min=minval((/ele, ele_b/))
  ele_max=maxval((/ele, ele_b/))

  elorg(elr,:)=(/-1, -1, -1, -1 /)

! note that enod must be modified, sufficient to
! modify the element with lowest element number, the other will 
! be removed later on
! node on a straight should be removed
! we know that a node exist between two other nodes


  i3=0
  do i1=1,3
    do i2=1,3
      if (enod(ele_max,i1).eq.enod(ele_min,i2)) then
        i3=i3+1
        cnod(i3)=enod(ele_max,i1)
        pos_max(i3)=i1
        pos_min(i3)=i2
      endif
    enddo
  enddo


  enod(ele_max,pos_max)=-enod(ele_max,pos_max)
  n1=maxval(enod(ele_max,:))
  enod(ele_max,pos_max)=-enod(ele_max,pos_max)

  enod(ele_min,pos_min)=-enod(ele_min,pos_min)
  n2=maxval(enod(ele_min,:))
  enod(ele_min,pos_min)=-enod(ele_min,pos_min)

  xm=(coord(n1,1)+coord(n2,1))*0.5d0
  ym=(coord(n1,2)+coord(n2,2))*0.5d0

  d1=(xm-coord(cnod(1),1))**2d0+(ym-coord(cnod(1),2))**2d0
  d2=(xm-coord(cnod(2),1))**2d0+(ym-coord(cnod(2),2))**2d0

  if (d1.lt.d2) then
! cnod(1) removed   
    enod(ele_min,pos_min(1))=n1
  else
! cnod(2) removed   
    enod(ele_min,pos_min(2))=n1
  endif

! check that area is positive
  ex=(/coord(enod(ele_min,1),1),coord(enod(ele_min,2),1),coord(enod(ele_min,3),1) /) 
  ey=(/coord(enod(ele_min,1),2),coord(enod(ele_min,2),2),coord(enod(ele_min,3),2) /) 
  DetA=ex(2)*ey(3)+ex(1)*ey(2)+ey(1)*ex(3) &
         -ex(2)*ey(1)-ex(3)*ey(2)-ey(3)*ex(1)
  if (DetA.lt.0d0) then
    write(*,*)' shift area negative'
    n1=enod(ele_min,1)
    enod(ele_min,1)=enod(ele_min,2)
    enod(ele_min,2)=n1
  endif

  enod(ele_max,:)=(/ -1, -1, -1 /)


!
! find new locations of ele_min and ele_max
!
  if (debug.eq.1) write(*,*)' removed element ', ele_max, elloc(ele_max)

  elloc(ele_min)=0
  elloc(ele_max)=0
  ie=added_nr_elements
  do while ((elloc(ele_min).eq.0).and.(elloc(ele_max).eq.0).and.(ie.gt.0)) 
    if ((elorg(ie,2).eq.ele_min).or.(elorg(ie,3).eq.ele_min)) then
      elloc(ele_min)=ie
    end if
    if ((elorg(ie,2).eq.ele_max).or.(elorg(ie,3).eq.ele_max)) then
      elloc(ele_max)=ie
    end if
    ie=ie-1
  end do
  
  if (debug.eq.1) write(*,*)'---------------leaving  remove1'
  return
end subroutine remove1


end module adaptive

