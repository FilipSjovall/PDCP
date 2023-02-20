module abaqus_util
! Last modified
! M. Ristinmaa 2011-04-16
!  -initial version based on old f77 code
! M. Ristinmaa 2011-04-17
!  -included lists for nset and elset. 
!  -bcnod and bcval are now read from nset
! M. Ristinmaa 2011-04-21
!  -key for concentrated loads implemented, cload
! M. Ristinmaa 2011-05-05
!   -Transpose of coord and enod, at end, this is a fix
!    and should be resolved as it slows down the routine
! M. Ristinmaa 2011-06-01
!   -"ENCASTRE" implemented in *Boundary	
! M. Ristinmaa 2011-10-27
!   -check for redudant boundary conditions implemented	
!   -implemented boundary condition 0d0 same as empty in abaqus inp-file
!   -local dof >4 is set to dof 4 in *Boundary
! M. Ristinmaa 2011-10-28
!   -increased len=120 for reading inp-line
! M. Ristinmaa 2012-01-16
!   -debubbed nset, error in generate: calculation of nn
!   -debubbed elset
!   -included *Surface
!   -included *Dsload
! M. Ristinmaa 2012-01-17
!   - bc dof>4 3D dof set to 4
!        dof>3 2D dof set to 3
! M. Ristinmaa 2012-05-14
!   - include *Dsload for element C3D8 should be move to element routines
!-----------------------------------------------------------------------------

! TODO

! The format for calls should be changed such that a call is 
! done for each entry, i.e.
! call abaqus_get('coord',coord)
! call abaqus_get('enod',enod)
! call abaqus_get('dload',dload)
! etc

! Should introduce a constant load vector and one that changes
	 
! this version can read nods and elements
! and displacement boundary conditions
	

use memory_util

implicit none

! maxnodinc : increase of number nodes in reallocation
! maxeleinc : increase of number elements in reallocation
integer, parameter    :: maxnodinc=1000
integer, parameter    :: maxeleinc=1000
integer               :: maxnod, maxele

private maxnodinc, maxeleinc, maxnod, maxele

! maxnodele : maximum number of nodes in one element
! maxcomlin : maximum number of commands in one line
integer, parameter    :: maxnodele=24
integer, parameter    :: maxcomlin=40

private maxnodele, maxcomlin

type list
   character(len=80)              :: name
   integer, allocatable           :: da(:)
end type

private list

integer, parameter               :: nrsub_surf=10
type list_surf
   character(len=80)              :: name
   integer                        :: nr
   type(list)                     :: surf(nrsub_surf)
end type

private list_surf

! maxnsetinc : increase of number of nset in reallocation
! maxelstinc : increase of number of elset in reallocation
! maxsurfinc : increase of number of surfaces in reallocation
integer, parameter     :: maxnsetinc=10
integer, parameter     :: maxelstinc=10
integer, parameter     :: maxsurfinc=10
integer                :: maxnset, maxelst, maxsurf

private maxnsetinc, maxelstinc, maxsurfinc, maxnset, maxelst, maxsurf

! maxbcinc : increase of number of bc in reallocation
integer, parameter     :: maxbcinc=10
integer                :: maxbc

private maxbcinc, maxbc

! maxflinc : increase of number of fl in reallocation
! concentrated loads
integer, parameter     :: maxflinc=10
integer                :: maxfl

private maxflinc, maxfl

! maxclinc : increase of number of fl in reallocation
! distributed loads
integer, parameter     :: maxclinc=30
integer                :: maxcl

private maxclinc, maxcl

interface abareadinp
   module procedure abaread1
   module procedure abaread2
   module procedure abaread3
   module procedure abaread4
end interface


private nsetcom
private elstcom
private elemen
private parsline
private cross3

!-----------------------------------------------------------------------------
contains

subroutine abaread1(filename,coord,enod)
  implicit none
  double precision, allocatable     :: coord(:,:), bcval(:), flval(:)
  integer, allocatable              :: enod(:,:), bcnod(:,:), flnod(:,:)
  character(len=*)                  :: filename

  call abaread3(filename,coord,enod,bcnod,bcval,flnod,flval)
  
  return
end subroutine abaread1  
	
subroutine abaread2(filename,coord,enod,bcnod,bcval)
  implicit none
  double precision, allocatable     :: coord(:,:), bcval(:), flval(:)
  integer, allocatable              :: enod(:,:), bcnod(:,:), flnod(:,:)
  character(len=*)                  :: filename

  call abaread3(filename,coord,enod,bcnod,bcval,flnod,flval)
  
  return
end subroutine abaread2
	

subroutine abaread3(filename,coordx,enodx,bcnod,bcval,flnod,flval)
  implicit none
  double precision, allocatable     :: coordx(:,:)
  integer, allocatable              :: enodx(:,:)
  double precision, allocatable     :: bcval(:), flval(:), fcval(:)
  integer, allocatable              :: bcnod(:,:), flnod(:,:), fcnod(:,:)
  character(len=*)                  :: filename

  call abaread4(filename,coordx,enodx,bcnod,bcval,flnod,flval,fcnod,fcval)

  return
end subroutine abaread3
	

subroutine abaread4(filename,coordx,enodx,bcnod,bcval,flnod,flval,fcnod,fcval)
  implicit none
  double precision, allocatable     :: coordx(:,:)
  integer, allocatable              :: enodx(:,:)
  double precision, allocatable     :: coord(:,:), bcval(:), flval(:), fcval(:)
  integer, allocatable              :: enod(:,:), bcnod(:,:), flnod(:,:), fcnod(:,:)
  character(len=*)                  :: filename

  character(len=120)                 :: strin, nsetname, elstname, surfname, eltype
  character(len=120)                 :: command(maxcomlin), command2(maxcomlin)
  integer                           :: nrcom
  
  integer                           :: nelem, nodel
  integer                           :: nrcoor
  double precision                  :: coore(3)

  integer                           :: nndim, knode(maxnodele), ierr, nndof
  integer                           :: nrnod, nrele, nrnset, nrelst, nrsurf
  integer                           :: nrbc, nrfl, nrcl
  
  integer                           :: dof1, dof2
  double precision                  :: val, Q
  
  type(list),allocatable            :: nset(:), elset(:)
  type(list_surf), allocatable      :: surfs(:)
  integer                           :: igener, nstart, nend, ninc, nns, isurf
  
  integer                           :: i, j, nn, jf, i2
  integer                           :: node_bc, dof_bc

  if (allocated(coord)) deallocate(coord)
  if (allocated(enod))  deallocate(enod)
  if (allocated(bcnod)) deallocate(bcnod)
  if (allocated(bcval)) deallocate(bcval)
  if (allocated(flnod)) deallocate(flnod)
  if (allocated(flval)) deallocate(flval)


  maxnod=maxnodinc
  maxele=maxeleinc
  maxnset=maxnsetinc
  maxelst=maxelstinc
  maxsurf=maxsurfinc
  maxbc=maxbcinc
  maxfl=maxflinc
  maxcl=maxclinc

  nndim=0
  nrnod=0
  nrele=0
  nrnset=0
  nrelst=0
  nrsurf=0
  nrbc=0
  nrfl=0
  nrcl=0
  
  allocate(nset(maxnset),stat=ierr)
  allocate(elset(maxelst),stat=ierr)
  allocate(surfs(maxsurf),stat=ierr)
  allocate(bcnod(maxbc,2),stat=ierr)
  allocate(bcval(maxbc),stat=ierr)
  allocate(flnod(maxfl,2),stat=ierr)
  allocate(flval(maxfl),stat=ierr)
  allocate(fcnod(maxcl,2),stat=ierr)
  allocate(fcval(maxcl),stat=ierr)

  open(unit=21,file=filename,status='UNKNOWN')

10 READ(21,"(A)",END=12,ERR=9999)STRIN

! Resolve * line	
	
   IF (STRIN(1:1).EQ.'*') THEN
!         write(*,*)strin
     IF (STRIN(2:2).NE.'*') then
      call parsline(STRIN,command,nrcom)
!      write(*,*)'command ',strin
!     do j=1,nrcom
!        write(*,*)'command ',nrcom,command(j)
!     enddo
 
      IF (command(1).EQ.'*Element') THEN
        CALL ELEMEN(command,eltype,nodel,nndof)
      elseif (command(1).eq.'*Nset') then
        nrnset=nrnset+1
        call nsetcom(command,nset,igener,nsetname)
      elseif (command(1).eq.'*Elset') then
        nrelst=nrelst+1
        call elstcom(command,elset,igener,elstname)
!write(*,*)'nr elset ',nrelst,elstname
      elseif (command(1).eq.'*Surface') then
        nrsurf=nrsurf+1
        isurf=0
        call surfcom(command,surfs,igener,surfname)
      ENDIF
     endif
    !read(*,*)j
   ELSE
!
! Resolve data line
!
! Node
!
     IF (command(1).EQ.'*Node') THEN

       if (nndim.eq.0) then
         do j=1,79
           if (strin(j:j).eq.',') nndim=nndim+1      
         enddo
       endif

       READ(STRIN,*,ERR=9999)NRCOOR,(COORE(J),J=1,nndim)
       if (.not.(allocated(coord))) allocate(coord(maxnod,nndim),stat=ierr)
       if (nrcoor.gt.maxnod) then
         maxnod=maxnod+maxnodinc
         call reallocate(coord,maxnod)
       endif
       coord(nrcoor,:)=COORE(1:nndim)
       nrnod=nrnod+1
!
! Element
!
     ELSEIF (command(1).EQ.'*Element') THEN

       READ(STRIN,*)NELEM,(KNODE(J),J=1,nodel)
       if (.not.(allocated(enod))) allocate(enod(maxele,nodel),stat=ierr)
       if (nelem.gt.maxele) then
         maxele=maxele+maxeleinc
         call reallocate(enod,maxele)
       endif
       enod(nelem,:)=knode(1:nodel)
       nrele=nrele+1
!
! Cload
!
     ELSEIF (command(1).EQ.'*Cload') THEN
     	     
       call parsline(STRIN,command2,nrcom)
       if (isnumber(command2(1))) then 
         write(*,*)'Not implemented'
         stop
       else
         jf=-1
         do j=1,nrnset
           if (command2(1).eq.nset(j)%name) jf=j
         enddo
         if (jf.ne.-1) then
           nn=size(nset(jf)%da)
           do while ((nrfl+nn).gt.maxfl) 
             maxfl=maxfl+maxflinc
!             maxfl=maxbc+maxflinc
             call reallocate(flnod,maxfl)
             call reallocate(flval,maxfl)
           enddo
           flnod((nrfl+1):(nrfl+nn),1)=nset(jf)%da
           flval((nrfl+1):(nrfl+nn))=0d0
! we assume that there is an integer in the second position
           read(command2(2),*,ERR=9999)dof1
           read(command2(3),*,ERR=9999)val           
           flnod((nrfl+1):(nrfl+nn),2)=dof1
           flval((nrfl+1):(nrfl+nn))=val
           nrfl=nrfl+size(nset(jf)%da)
         else
           write(*,*)'ERROR Node set not found', command2(1)
           stop
         endif  
       endif
!
! Boundary
!
     ELSEIF (command(1).EQ.'*Boundary') THEN
 !      write(*,*)'command ',strin
    	     
       call parsline(STRIN,command2,nrcom)
       if (isnumber(command2(1))) then 
         write(*,*)'Not implemented'
         stop
       else
         jf=-1
         do j=1,nrnset
           if (command2(1).eq.nset(j)%name) jf=j
         enddo
         if (jf.ne.-1) then
           nn=size(nset(jf)%da)
           do while ((nrbc+nn).gt.maxbc) 
             maxbc=maxbc+maxbcinc
             call reallocate(bcnod,maxbc)
             call reallocate(bcval,maxbc)
           enddo
           bcnod((nrbc+1):(nrbc+nn),1)=nset(jf)%da
           bcval((nrbc+1):(nrbc+nn))=0d0
           if (command2(2).eq.'XSYMM') then
              bcnod((nrbc+1):(nrbc+nn),2)=1           
           elseif (command2(2).eq.'YSYMM') then
              bcnod((nrbc+1):(nrbc+nn),2)=2
           elseif (command2(2).eq.'ZSYMM') then
              bcnod((nrbc+1):(nrbc+nn),2)=3
           elseif ((command2(2).eq.'ENCASTRE').or.(command2(2).eq.'PINNED')) then
              bcnod((nrbc+1):(nrbc+nn),2)=1
              nrbc=nrbc+nn
              do while ((nrbc+nn).gt.maxbc) 
                maxbc=maxbc+maxbcinc
                call reallocate(bcnod,maxbc)
                call reallocate(bcval,maxbc)
              enddo
              bcnod((nrbc+1):(nrbc+nn),1)=nset(jf)%da
              bcnod((nrbc+1):(nrbc+nn),2)=2
              if (nndim.eq.3) then
                nrbc=nrbc+nn
                do while ((nrbc+nn).gt.maxbc) 
                  maxbc=maxbc+maxbcinc
                  call reallocate(bcnod,maxbc)
                  call reallocate(bcval,maxbc)
                enddo
                bcnod((nrbc+1):(nrbc+nn),1)=nset(jf)%da
                bcnod((nrbc+1):(nrbc+nn),2)=3
              endif
           else  !if (isnumber(command2(2))) then
             read(command2(2),*,ERR=9999)dof1
             read(command2(3),*,ERR=9999)dof2
             if (nrcom.eq.4) then
               read(command2(4),*,ERR=9999)val      
             else
               val=0d0
             end if     
             if (dof1.eq.dof2) then
               if ((dof1.gt.4).and.(nndim.eq.3)) dof1=4
               if ((dof1.gt.3).and.(nndim.eq.2)) dof1=3
               bcnod((nrbc+1):(nrbc+nn),2)=dof1
               bcval((nrbc+1):(nrbc+nn))=val
             else
               write(*,*)' Do not know what the below means', strin
               stop
             endif
   !        else
   !           write(*,*)'Not implemented ',command2(2)
   !           stop
           endif
           nrbc=nrbc+size(nset(jf)%da)
         else
           write(*,*)'ERROR Node set not found', command2(1)
           stop
         endif  
       endif
!
! Nset
!
     ELSEIF (command(1).EQ.'*Nset') THEN
       if (nrnset.ge.maxnset) then
         stop 'Increase maxnsetinc or rewrite for dynamic reallocation'
       endif

       if (igener.eq.-1) then
! note can continue over several lines
	call parsline(STRIN,command2,nn)
 	 if (allocated(nset(nrnset)%da)) then
 	   nns=size(nset(nrnset)%da)
 	   call reallocate(nset(nrnset)%da,nns+nn)
 	 else	 
           nns=0
           allocate(nset(nrnset)%da(nn),stat=ierr)
         endif
         READ(STRIN,*,ERR=9999)(nset(nrnset)%da(j+nns),j=1,nn)
         nset(nrnset)%name=nsetname
       else
         READ(STRIN,*,ERR=9999)nstart,nend,ninc
         nn=(nend-nstart)/ninc+1
         allocate(nset(nrnset)%da(nn),stat=ierr)
         nset(nrnset)%da=(/nstart:nend:ninc/)
         nset(nrnset)%name=nsetname
       endif

!
! Elset
!
     ELSEIF (command(1).EQ.'*Elset') THEN
 
       if (nrelst.ge.maxelst) then
         stop 'Increase maxelstinc or rewrite for dynamic reallocation'
       endif
!
       if (igener.eq.-1) then
! simple list of elements
         nn=1
         call parsline(STRIN,command2,nn)
 	 if (allocated(elset(nrelst)%da)) then
 	   nns=size(elset(nrelst)%da)
 	   call reallocate(elset(nrelst)%da,nns+nn)
 	 else	 
           nns=0
           allocate(elset(nrelst)%da(nn),stat=ierr)
         endif
         READ(STRIN,*,ERR=9999)(elset(nrelst)%da(j+nns),j=1,nn)
         elset(nrelst)%name=elstname
       else
! list is generated
! write(*,*)'generated elset ',strin,elstname
       	 READ(STRIN,*,ERR=9999)nstart,nend,ninc
         nn=(nend-nstart)/ninc+1
         allocate(elset(nrelst)%da(nn),stat=ierr)
         elset(nrelst)%da=(/nstart:nend:ninc/)
         elset(nrelst)%name=elstname
       endif
!
! Surface
!
     ELSEIF (command(1).EQ.'*Surface') THEN
       if (nrsurf.ge.maxsurf) then
         stop 'Increase maxelstinc or rewrite for dynamic reallocation'
       endif
!      write(*,*)'command ',strin
       call parsline(STRIN,command2,nrcom)
 !      write(*,*)'command2 ',command2(1)
       if (isnumber(command2(1))) then 
         write(*,*)'Not implemented'
         stop
       else
         isurf=isurf+1
         allocate(surfs(nrsurf)%surf(isurf)%da(1),stat=ierr)
         if (command2(2).eq.'S1') then
           surfs(nrsurf)%surf(isurf)%da(1)=1
         elseif (command2(2).eq.'S2') then
           surfs(nrsurf)%surf(isurf)%da(1)=2
         elseif (command2(2).eq.'S3') then
           surfs(nrsurf)%surf(isurf)%da(1)=3
         elseif (command2(2).eq.'S4') then
           surfs(nrsurf)%surf(isurf)%da(1)=4
         elseif (command2(2).eq.'S5') then
           surfs(nrsurf)%surf(isurf)%da(1)=5
         elseif (command2(2).eq.'S6') then
           surfs(nrsurf)%surf(isurf)%da(1)=6
         else
           write(*,*)'surf name ',surfname, command2(2)
           stop 'surface not implemented'
         end if
         surfs(nrsurf)%name=surfname
         surfs(nrsurf)%nr=isurf
         surfs(nrsurf)%surf(isurf)%name=command2(1)
       end if
!
! Dsload
!
     ELSEIF (command(1).EQ.'*Dsload') THEN
       call parsline(STRIN,command2,nrcom)
!       write(*,*)'command2 ',command2(1)
         jf=-1
         do j=1,nrsurf
!           write(*,*)'surf name ',surfs(j)%name
           if (command2(1).eq.surfs(j)%name) jf=j
         enddo
         if (jf.ne.-1) then
!
! Calculate load using shape functions for element
! Only pressure is implemented command2(2)='P'
!
!          write(*,*)'surface found ',jf,surfs(jf)%name
!          write(*,*)'number',surfs(jf)%nr
          if (command2(2).eq.'P') then
            READ(command2(3),*,ERR=9999) Q
            if (eltype.eq.'CPE8RP') then
              call dload_2d8(fcnod,fcval,nrcl,Q,surfs,jf,elset,coord,nrelst,enod)
            elseif (eltype.eq.'C3D8') then
              call dload_3d8(fcnod,fcval,nrcl,Q,surfs,jf,elset,coord,nrelst,enod)
            else
              write(*,*)'distributed load for element', eltype
              stop 'not implemented'
            end if
          else
            write(*,*)strin
            stop 'loading not implemented'
          end if
         else
           write(*,*)strin
           stop 'elset not found'
         end if
!         do j=1,nrcl
!           write(*,*)'node dof val ',fcnod(j,:),fcval(j)
!         end do         

!
! Dsflow
!
     ELSEIF (command(1).EQ.'*Dsflow') THEN
 !      write(*,*)'command ',strin
    	     
       call parsline(STRIN,command2,nrcom)
       if (isnumber(command2(1))) then 
         write(*,*)'Not implemented'
         stop
       else
         jf=-1
         do j=1,nrsurf
!           write(*,*)'surf name ',surfs(j)%name
           if (command2(1).eq.surfs(j)%name) jf=j
         enddo
         if (command2(2).eq.'S') then
           read(command2(3),*,ERR=9999)val      
           call dsflow_2d(bcnod,bcval,nrbc,val,surfs,jf,elset,nrelst,enod)
         else
           write(*,*)'ERROR Node set not found', command2(1)
           stop
         endif  
       endif

!
! Rest
!
     ELSEIF (command(1).EQ.'*Heading') THEN
     ELSEIF (command(1).EQ.'*Solid') THEN
     ELSEIF (command(1).EQ.'*Static') THEN
!
! END STEP
!
     ELSEIF (command(1).EQ.'*End') THEN
       GOTO 12
!
! W A R N I N G
!
     ELSE
       WRITE(*,*)'WARNING NO MATCH ',command(1),STRIN
     ENDIF
   ENDIF
   GOTO 10
!
! E R R O R   R E A D
!
9999 WRITE(*,*)'   ERROR READING HLFFEM INPUTFILE '
   STOP
!
12 CLOSE (21)
!

! Release allocated memory	
  do j=1,nrnset
    deallocate(nset(j)%da)
  enddo
  deallocate(nset)
!   
  do j=1,nrelst
    deallocate(elset(j)%da)
  enddo
  deallocate(elset)

! Get rid of redundant boundary conditions

!do i=1,nrbc
!  write(*,*)'bc_aba ',i,bcnod(i,1),bcnod(i,2)
!end do

  do i=1,nrbc-1
    node_bc=bcnod(i,1)
    dof_bc=bcnod(i,2)
    if (node_bc.ne.-1) then
    do j=i+1,nrbc
      if ((bcnod(j,1).eq.node_bc).and.(bcnod(j,2).eq.dof_bc)) then
!        write(*,*)'redudant node dof',j,bcnod(j,1),bcnod(j,2)
        bcnod(j,1)=-1
      end if
    end do
    end if
  end do

  j=0
  do i=1,nrbc
     do while (bcnod(i,1).eq.-1)
       if ((i+(j-1)).ne.nrbc) then
        j=j+1
        do i2=i,nrbc-1
          bcnod(i2,1)=bcnod(i2+1,1)
          bcnod(i2,2)=bcnod(i2+1,2)
          bcval(i2)=bcval(i2+1)
        end do
        bcnod(nrbc,1)=0
        if (i.eq.nrbc) bcnod(i,1)=0
       end if
      end do
  end do
  nrbc=nrbc-j
!  write(*,*)'number redundant bc', j
!do i=1,nrbc
!  write(*,*)'bc_aba ',i,bcnod(i,1),bcnod(i,2)
!end do
! pause  
! 

  call reallocate(coord,nrnod)
  call reallocate(enod,nrele)
  call reallocate(bcnod,nrbc)
  call reallocate(bcval,nrbc)
  call reallocate(flnod,nrfl)
  call reallocate(flval,nrfl)
  call reallocate(fcnod,nrcl)
  call reallocate(fcval,nrcl)

  allocate(coordx(size(coord,2),size(coord,1)),stat=ierr)
  allocate(enodx(size(enod,2),size(enod,1)),stat=ierr)
  enodx=transpose(enod)
  coordx=transpose(coord)

!  write(*,*)'Number of nods                :',NRNOD
!  write(*,*)'Number of dofs                :',NRNOD*NNDIM
!  write(*,*)'Number of elements            :',NRELE
!  write(*,*)'Number of boundary conditions :',NRBC
!  write(*,*)'Number of concentrated loads  :',NRFL

  return
end subroutine abaread4


subroutine dsflow_2d(bcnod,bcval,nrbc,Q,surfs,isurf,elset,nrelst,enod)
  implicit none
  integer, allocatable           :: bcnod(:,:)
  integer                        :: isurf, nrelst, nrbc, enod(:,:)
  double precision, allocatable  :: bcval(:)
  double precision               :: Q
!  type(list),allocatable            :: nset(:), elset(:)
!  type(list_surf), allocatable      :: surfs(:)
  type(list_surf)                :: surfs(:)
  type(list)                     :: elset(:)
    
  integer                        :: i, jf, side, j, iel, ele
  integer                        :: n1, n2, n3, n4
  double precision               :: x1, x2, x3, x4, y1, y2, y3, y4
  double precision               :: normal(2), v_tmp(2), L
  
! quadratic shape function on boundary
! calculate size and direction in the three nodes
! put the values in flnod, flval
! loop over all element sets, need also to identify the nodes on the side
! distribution (1/6,2/3,1/6)*Q*L

! loop over all surfaces for this loading
  do i=1,surfs(isurf)%nr

! find element set for this surface
    jf=-1
    do j=1,nrelst
      if (surfs(isurf)%surf(i)%name.eq.elset(j)%name) jf=j
    enddo
    if (jf.eq.-1) stop 'dsflow - element set not found'
!     write(*,*)'elset name',elset(jf)%name
! extract node numbers
    side=surfs(isurf)%surf(i)%da(1)
 
! loop over all elements in elset
    iel=1
     do iel=1,size(elset(jf)%da)
!      write(*,*)'element number ',iel,elset(jf)%da(iel), side
      ele=elset(jf)%da(iel)
!if (1.eq.0) then
      if (side.eq.1) then
        n1=enod(ele,1)
        n2=enod(ele,5)
        n3=enod(ele,2)
      elseif (side.eq.2) then
        n1=enod(ele,2)
        n2=enod(ele,6)
        n3=enod(ele,3)
      elseif (side.eq.3) then
        n1=enod(ele,3)
        n2=enod(ele,7)
        n3=enod(ele,4)
      elseif (side.eq.4) then
        n1=enod(ele,4)
        n2=enod(ele,8)
        n3=enod(ele,1)
      end if  
!      write(*,*)'nods ',n1,n2,n3

!
! check allocation
! only corner nodes are prescribed, i.e. 2 
! THIS IS A SPECIAL CASE WORKS ONLY WITH QU8
      do while ((nrbc+2).gt.maxbc) 
         maxbc=maxbc+maxbcinc
         call reallocate(bcnod,maxbc)
         call reallocate(bcval,maxbc)
      enddo

! loading, note that it is a pore pressure 
! that is applied and dof 3 is assume
      bcnod(nrbc+1,1)=n1; bcnod(nrbc+1,2)=3; bcval(nrbc+1)=Q
      bcnod(nrbc+2,1)=n3; bcnod(nrbc+2,2)=3; bcval(nrbc+2)=Q
      nrbc=nrbc+2

    end do

  end do

end subroutine dsflow_2d



subroutine dload_3d8(fcnod,fcval,nrcl,Q,surfs,isurf,elset,coord,nrelst,enod)
  implicit none
  integer, allocatable           :: fcnod(:,:)
  integer                        :: isurf, nrelst, nrcl, enod(:,:)
  double precision, allocatable  :: fcval(:)
  double precision               :: coord(:,:), Q
!  type(list),allocatable            :: nset(:), elset(:)
!  type(list_surf), allocatable      :: surfs(:)
  type(list_surf)                :: surfs(:)
  type(list)                     :: elset(:)

  integer :: i, j, jf, side, ele, iel, igp
  integer :: nodes(4)
  double precision :: x(4), y(4), z(4)
  double precision :: gp1, gp2, xsi(9), eta(9), etv, xsv
  double precision :: wp1, wp2, wpxsi(9), wpeta(9)

  double precision :: dx(3), dy(3), dNdeta(4), dNdxsi(4)
  double precision :: normal(3), Jac, traction(3)
  double precision :: Nt(12,3), N1, N2, N3, N4, fload(12)

  parameter        (gp1=0.0D0)
  parameter        (gp2=0.774596669241483D0)
  parameter        (xsi=(/-gp2, gp1, gp2,-gp2, gp1, gp2,-gp2, gp1, gp2/))
  parameter        (eta=(/-gp2,-gp2,-gp2, gp1, gp1, gp1, gp2, gp2, gp2/))

  parameter        (wp1=0.888888888888889D0)
  parameter        (wp2=0.555555555555556D0)
  parameter        (wpxsi=(/wp2,wp1,wp2,wp2,wp1,wp2,wp2,wp1,wp2/))
  parameter        (wpeta=(/wp2,wp2,wp2,wp1,wp1,wp1,wp2,wp2,wp2/))


! loop over all surfaces for this loading
  do i=1,surfs(isurf)%nr

! find element set for this surface
    jf=-1
    do j=1,nrelst
      if (surfs(isurf)%surf(i)%name.eq.elset(j)%name) jf=j
    enddo
    if (jf.eq.-1) stop 'dload - element set not found'
! extract node numbers
    side=surfs(isurf)%surf(i)%da(1)
!     write(*,*)'elset name and side',elset(jf)%name, side

! loop over all elements in elset
     do iel=1,size(elset(jf)%da)
!      write(*,*)'element number ',iel,elset(jf)%da(iel), side
      ele=elset(jf)%da(iel)
! faces from Abaqus user's manual, vol II, page 14.1.4-15
      if (side.eq.1) then
        nodes=enod(ele,[1,2,3,4])
      elseif (side.eq.2) then
        nodes=enod(ele,[5,8,7,6])
      elseif (side.eq.3) then
        nodes=enod(ele,[1,5,6,2])
      elseif (side.eq.4) then
        nodes=enod(ele,[2,6,7,3])
      elseif (side.eq.5) then
        nodes=enod(ele,[3,7,8,4])
      elseif (side.eq.6) then
        nodes=enod(ele,[4,8,5,1])
      end if  
!      write(*,*)'nods ',nodes

! coordinates of the nodes on the side

      x=coord(nodes,1)
      y=coord(nodes,2)
      z=coord(nodes,3)

! loop over gauss points

      fload=0d0
      do igp=1,9
        
        etv=eta(igp)
        xsv=xsi(igp)
        dNdxsi=0.25d0*(/(etv-1d0),-1d0*(etv-1d0),(etv+1d0),-1d0*(etv+1)/)
        dNdeta=0.25d0*(/(xsv-1d0),-1d0*(xsv+1d0),(xsv+1d0),-1d0*(xsv-1)/)

        dx=(/dot_product(dNdxsi,x),dot_product(dNdxsi,y),dot_product(dNdxsi,z)/)
        dy=(/dot_product(dNdeta,x),dot_product(dNdeta,y),dot_product(dNdeta,z)/)

! define normal to the surface and traction load
        normal=cross3(dx,dy)
        Jac=dsqrt(normal(1)**2d0+normal(2)**2d0+normal(3)**2d0)
        traction=normal/Jac
!        write(*,*)'traction ',traction

        N1= 0.25d0*(xsv-1d0)*(etv-1d0)
        N2=-0.25d0*(xsv+1d0)*(etv-1d0)
        N3= 0.25d0*(xsv+1d0)*(etv+1d0)
        N4=-0.25d0*(xsv-1d0)*(etv+1d0)

        Nt=0d0
        Nt(1,1) =N1; Nt(2,2) =N1; Nt(3,3) =N1
        Nt(4,1) =N2; Nt(5,2) =N2; Nt(6,3) =N2
        Nt(7,1) =N3; Nt(8,2) =N3; Nt(9,3) =N3
        Nt(10,1)=N4; Nt(11,2)=N4; Nt(12,3)=N4

        fload=fload+Q*matmul(Nt,traction)*Jac*wpxsi(igp)*wpeta(igp)        
 
      end do

! loading
      if (nrcl+12.gt.maxcl) then
        maxcl=maxcl+maxclinc
        call reallocate(fcnod,maxcl)
        call reallocate(fcval,maxcl)
      end if             

! loading, note that it is a pressure that is applied

      fcnod([nrcl+[1,2,3]],1)=nodes(1)
      fcnod([nrcl+[4,5,6]],1)=nodes(2)
      fcnod([nrcl+[7,8,9]],1)=nodes(3)
      fcnod([nrcl+[10,11,12]],1)=nodes(4)

      fcnod([nrcl+[1,4,7,10]],2)=1
      fcnod([nrcl+[2,5,8,11]],2)=2
      fcnod([nrcl+[3,6,9,12]],2)=3
      fcval([nrcl+1:nrcl+12])=fload

!      write(*,*)'fcnod ',fcnod([nrcl+1:nrcl+12],1)
      nrcl=nrcl+12
 
   end do
  end do

end subroutine dload_3d8

 
subroutine dload_2d8(fcnod,fcval,nrcl,Q,surfs,isurf,elset,coord,nrelst,enod)
  implicit none
  integer, allocatable           :: fcnod(:,:)
  integer                        :: isurf, nrelst, nrcl, enod(:,:)
  double precision, allocatable  :: fcval(:)
  double precision               :: coord(:,:), Q
!  type(list),allocatable            :: nset(:), elset(:)
!  type(list_surf), allocatable      :: surfs(:)
  type(list_surf)                :: surfs(:)
  type(list)                     :: elset(:)
    
  integer                        :: i, jf, side, j, iel, ele
  integer                        :: n1, n2, n3, n4
  double precision               :: x1, x2, x3, x4, y1, y2, y3, y4
  double precision               :: normal(2), v_tmp(2), L
  
! quadratic shape function on boundary
! calculate size and direction in the three nodes
! put the values in flnod, flval
! loop over all element sets, need also to identify the nodes on the side
! distribution (1/6,2/3,1/6)*Q*L

! loop over all surfaces for this loading
  do i=1,surfs(isurf)%nr

! find element set for this surface
    jf=-1
    do j=1,nrelst
      if (surfs(isurf)%surf(i)%name.eq.elset(j)%name) jf=j
    enddo
    if (jf.eq.-1) stop 'dload - element set not found'
!     write(*,*)'elset name',elset(jf)%name
! extract node numbers
    side=surfs(isurf)%surf(i)%da(1)

! loop over all elements in elset
    iel=1
     do iel=1,size(elset(jf)%da)
!      write(*,*)'element number ',iel,elset(jf)%da(iel), side
      ele=elset(jf)%da(iel)
!if (1.eq.0) then
      if (side.eq.1) then
        n1=enod(ele,1)
        n2=enod(ele,5)
        n3=enod(ele,2)
      elseif (side.eq.2) then
        n1=enod(ele,2)
        n2=enod(ele,6)
        n3=enod(ele,3)
      elseif (side.eq.3) then
        n1=enod(ele,3)
        n2=enod(ele,7)
        n3=enod(ele,4)
      elseif (side.eq.4) then
        n1=enod(ele,4)
        n2=enod(ele,8)
        n3=enod(ele,1)
      end if  
!      write(*,*)'nods ',n1,n2,n3

! coordinates and length of element side
! for simplicity it will be assumed that the side is straigth
      x1=coord(n1,1)
      y1=coord(n1,2)
      x2=coord(n2,1)
      y2=coord(n2,2)
      x3=coord(n3,1)
      y3=coord(n3,2)
    
! TODO
! Note that isoparametric mapping should be used for taking care of
! curved boundaries.

      L=dsqrt((x1-x3)**2d0+(y1-y3)**2d0)

! define normal pointing out from the element    
      normal=(/ -(y3-y1)/L,(x3-x1)/L /)

! define a vector using any other node, a scalar multiplication between these should
! then be positive, select mid node on opposite element side
      if (side.eq.1) n4=enod(ele,7)
      if (side.eq.2) n4=enod(ele,8)
      if (side.eq.3) n4=enod(ele,5)
      if (side.eq.4) n4=enod(ele,6)
      x4=coord(n4,1)
      y4=coord(n4,2)
      v_tmp=(/ (x4-x2)/L,(y4-y2)/L /)
     
      if (dot_product(normal,v_tmp).gt.0d0) normal=-normal 

! loading
      if (nrcl+6.gt.maxcl) then
        maxcl=maxcl+maxclinc
        call reallocate(fcnod,maxcl)
        call reallocate(fcval,maxcl)
      end if             

! loading, note that it is a pressure that is applied
      fcnod(nrcl+1,1)=n1; fcnod(nrcl+1,2)=1; fcval(nrcl+1)=-normal(1)*Q*L/6d0
      fcnod(nrcl+2,1)=n1; fcnod(nrcl+2,2)=2; fcval(nrcl+2)=-normal(2)*Q*L/6d0
    
      fcnod(nrcl+3,1)=n2; fcnod(nrcl+3,2)=1; fcval(nrcl+3)=-normal(1)*Q*L*2d0/3d0
      fcnod(nrcl+4,1)=n2; fcnod(nrcl+4,2)=2; fcval(nrcl+4)=-normal(2)*Q*L*2d0/3d0

      fcnod(nrcl+5,1)=n3; fcnod(nrcl+5,2)=1; fcval(nrcl+5)=-normal(1)*Q*L/6d0
      fcnod(nrcl+6,1)=n3; fcnod(nrcl+6,2)=2; fcval(nrcl+6)=-normal(2)*Q*L/6d0
        
      nrcl=nrcl+6
    end do

  end do

end subroutine dload_2d8

subroutine surfcom(command,elset,igener,elstname)
  implicit none
  character(len=*)                  :: command(maxcomlin), elstname
  type(list_surf), allocatable      :: elset(:)
  integer                           :: igener

  integer      :: j, jfound

! check if generate statement is present
  igener=-1
  do j=1,maxcomlin
    if (command(j).eq.'generate') igener=j
  enddo

! name of the nset
  do j=2,maxcomlin
    if (command(j).eq.'name') jfound=j
  enddo
  elstname=command(jfound+1)

  return


end subroutine surfcom

subroutine elstcom(command,elset,igener,elstname)
  implicit none
  character(len=*)                  :: command(maxcomlin), elstname
  type(list), allocatable           :: elset(:)
  integer                           :: igener

  integer      :: j, jfound

! check if generate statement is present
  igener=-1
  do j=1,maxcomlin
    if (command(j).eq.'generate') igener=j
  enddo

! name of the nset
  do j=2,maxcomlin
    if (command(j).eq.'elset') jfound=j
  enddo
  elstname=command(jfound+1)

  return


end subroutine elstcom


subroutine nsetcom(command,nset,igener,nsetname)
  implicit none
  character(len=*)                  :: command(maxcomlin), nsetname
  type(list), allocatable           :: nset(:)
  integer                           :: igener

  integer      :: j, jfound

! check if generate statement is present
  igener=-1
  do j=1,maxcomlin
    if (command(j).eq.'generate') igener=j
  enddo

! name of the nset
  do j=2,maxcomlin
    if (command(j).eq.'nset') jfound=j
  enddo
  nsetname=command(jfound+1)

  return
end subroutine nsetcom


SUBROUTINE ELEMEN(command,eltype,nodel,nndof)
  implicit none
  integer                           :: nodel, nndof
  character(len=*)                  :: command(maxcomlin), eltype
  
  integer      :: j, jfound

  do j=1,maxcomlin
    if (command(j).eq.'type') jfound=j
  enddo
  nodel=0
  select case(command(jfound+1))
  case('CPE3')
    nodel=3
    nndof=2
  case('CPS3') 
    nodel=3
    nndof=2
  case('CPE4')
    nodel=4
    nndof=2
  case('CPE8RP')  ! 8node biquadratic displ, bilinear pressure
    nodel=8
    nndof=2
  case('C3D8')
    nodel=8
    nndof=3
  case('C3D4')
    nodel=4
    nndof=3
  case('C3D8T')
    nodel=8
    nndof=4
  end select

  if (nodel.eq.0) write(*,*) command(jfound+1)
  if (nodel.eq.0) stop 'abaqus_util element not found'
!write(*,*)'jfound ',jfound, nodel

  eltype=command(jfound+1)
!
  return
end subroutine elemen  



subroutine parsline(line,sym,nsym)
  implicit none
!
!    Separate the individual words in the line for
!    later interpretation. 
!
!    INPUT
!     line   -   line to be parsed
!     OUTPUT
!     sym    -   array containing individual words in the line
!     nsym   -   number of words found in the line

  integer            :: istart, iend, ierr, lline, nsym, lsym, nfnd, iexp, i, j
  character          :: line*(*),sym(*)*(*), le*120


! search for , and = and space as separators
	
  character (len=3)    :: sep=',= '

!
  nfnd=maxcomlin
  
  le=adjustl(line)    
  lsym=len(sym(1))
  lline=len(le)

  nsym=0
  ierr=0
  iend=0
  istart=1

!  write(*,*)line
!
! Loop to find all words (symbols) in the line
!
  iend = scan(le,sep)
  do while (iend.ne.1)
     nsym=nsym+1
     sym(nsym)=adjustl(le(1:iend-1))
     le=adjustl(le((iend+1):lline))
     iend = scan(le,sep)
  end do
  sym((nsym+1):nfnd)=' '

 return
end subroutine parsline


logical function isnumber(str)
!     check if the string argument contain a number
  implicit none
  character        :: str*(*)   
!  logical          :: isnumber
  
  integer          :: i
  
  isnumber=.true.
	
  do i=1,len(str)
    if((str(i:i).lt.'0'.or.str(i:i).gt.'9').and.  &
        str(i:i).ne.'.'.and.str(i:i).ne.'-'.and.  &
        str(i:i).ne.'+'.and.str(i:i).eq.' ')then
       isnumber=.false.
       return
     end if
   end do

   return  
end function isnumber


FUNCTION cross3(a, b)
  double precision :: cross3(3)
  double precision, INTENT(IN) :: a(:), b(:)

  cross3(1) = a(2) * b(3) - a(3) * b(2)
  cross3(2) = a(3) * b(1) - a(1) * b(3)
  cross3(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION cross3
      
end module abaqus_util
