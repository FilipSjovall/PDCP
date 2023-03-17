!!!
module mesh_reader
! Program that reads from ASCII-mesh file
!!!
!
!
! The file hej.msh is in ASCII - 1 format
!
!
use fem_system
use matlab_util
use sparse_util
use elem_large_cont_2d6
use mater_hyperel
use elem_large_cont_3d
use elem_flow_3d
!type(sparse)                        :: KSYM, K
!integer, allocatable                :: cellSym(:,:)
!double precision, allocatable       :: valSym(:)
integer                             :: ndof   
integer                             :: ierr
!
!
contains 
!
!#########################################!
!      Read .gmsh (format = ASCII 1)      !
!#########################################!
subroutine read_mesh(fname, coord, enod)
    character(len=*) :: fname
    integer :: m,n, iostat, nstart, nend, elstart, eled, idum, idum2
    integer, dimension(1,4) :: dummy_int
    real(kind=8) :: x
    character(len=100) :: dummy_char
    real(kind=8), allocatable, intent(inout) :: coord(:,:)
    integer, allocatable, intent(inout) :: enod(:,:)

    print * , '------------------ READING MESH ------------------------'
    n = 0
    eled = 0
    open(45, file=fname, status="old",action="read")
    !
    !
    ! Loop to find the number of lines in the file
    do
        read(45 ,* ,IOSTAT=iostat) dummy_char
        n = n+1
        if (dummy_char(1:6) =="$NOD" ) then
            nstart = n
        elseif ( dummy_char(1:9)=="$ENDNOD" ) then
            nend   = n
        elseif ( dummy_char(1:9)=="$ELM" ) then
            elstart = n
        elseif ( dummy_char(1:13)=="$ENDELM" ) then
            eled   = n
        endif
        if(iostat.lt.0) then
            exit
        endif
    end do
    !
    print *, "Number of nodes ", nend-nstart-2
    !
    print *, "Number of elements ", eled-elstart-3
    !
    1 rewind(45)
    !
    !
    if(allocated(coord)) deallocate(coord)
    if(allocated(enod)) deallocate(enod)
    !
    allocate(coord(nend-nstart-2,2),stat=ierr)
    !
    allocate(enod(eled-elstart-3,7),stat=ierr)
    ! Loop to find the number of lines in the file
    do n=1, eled
        if (( n>nstart+1).and.(n.le.nend-1)) then
            read(45 ,* ,IOSTAT=iostat), idum2, coord(n-nstart-1,1), coord(n-nstart-1,2), idum
            !print *, n-nstart-1, coord(n-nstart-1,:)
        elseif ((n > elstart+1).and.(n<eled-1)) then
            read(45 ,* ,IOSTAT=iostat), enod(n-elstart-1,1), dummy_int, enod(n-elstart-1,2:7)
            !print *, enod(n-elstart-1,:), "\\"
        else
            read(45 ,* ,IOSTAT=iostat) dummy_char
            !print *, dummy_char
        endif
    end do
    !
end subroutine read_mesh
!
!
!#########################################!
!      Initialize sparse pattern          !
!!#########################################!
! ###################################################### !
! Set Newton parameters                                # !
! ###################################################### !

subroutine NewtonSolver(K,u,edofn,coord,mp,TOL,bcval,bcdof,nelm,ndofx,enod,t)   
type(sparse)                    :: K
double precision, intent(inout) :: u(:)
double precision                :: coord(:,:) 
integer, intent(in)             :: edofn(:,:), enod(:,:)
double precision, intent(in)    :: TOL
double precision, intent(in)    :: mp(2)
double precision                :: bcval(:)
integer, intent(in)             :: bcdof(:)
integer, intent(in)             :: nelm

double precision, allocatable   :: dg(:,:,:), ef(:,:), ed(:,:), ke(:,:), fe(:), Deltaa(:), F(:), bc0(:), res(:)
double precision                :: Dgp(3,3,3),esh(4,3), es(3,3), residual, t
integer                         :: el,ie, iter, dofnod
 
! # Allocation
allocate(dg(ngp,2,2), stat = ierr)
allocate(ef(4,ngp), stat=ierr)
allocate(ed(12,nelm), stat=ierr)

allocate(ke(12,12),stat=ierr)
allocate(fe(12), stat=ierr)

allocate(Deltaa(ndofx), stat=ierr)
allocate(F(ndofx), stat=ierr)
allocate(bc0(6),stat=ierr)

call neohooke_init(mp)   
bc0    = bcval
dofnod = 2
imax   = 25
! ###################################################### !
! For each load step n                                 # !
! ###################################################### !
do n = 1,10

   ! ###################################################### !
   ! Reset variables                                      # !
   ! ###################################################### !
   res     = 0d0 
   bcval   = bc0
   residual= 1d0
   iter    = 0
   print * , 'Starting equlibrium iter at loadstep:', n

   ! Newton-loop                                          # !
   ! ###################################################### !
   do while ((iter.lt.(imax)).and.(residual.gt.TOL).or.(iter.lt.2))
       iter=iter+1
       K%a=0d0
       F=0d0

      ! ###################################################### !
      ! Extract element displacments                         # !
      ! ###################################################### !
      !call extract(ed,a,transpose(edofn(:,2:13) ),dofnod)

      ! ###################################################### !
      ! Assembly                                             # !
      ! ###################################################### !
      do ie = 1, nelm
         ! ###################################################### !
         ! Compute the deformation gradient                     # !
         ! ###################################################### !
         call c2tl6_d(dg,u(edofn(ie,2:13)),coord(enod(ie,2:7),:))

         ! ###################################################### !
         ! Reshape deformation gradient          (needed?)      # !
         ! ###################################################### !
         do gp=1,3
            ef(:,gp) = (/ dg(gp,1,1), dg(gp,1,2), dg(gp,2,1), dg(gp,2,2) /)
            !ef(:,gp) = (/ dg(gp,1,1), dg(gp,2,1), dg(gp,1,2), dg(gp,2,2) /)
         end do	
         !print * , ef(:,1)

         ! ###################################################### !
         ! Compute the material stiffness [nelm x ngp x 3 x 3]  # !
         ! ###################################################### !
         call dneohooke('tl',Dgp,ef)

         ! ###################################################### !
         ! Compute the 2nd Piola Kirchhoff stress               # !
         ! ###################################################### !
         call neohooke('2ndPiola',esh,ef)        

         ! ###################################################### !
         ! Reshape stress tensor                                # !
         ! ###################################################### !
         do gp=1,ngp
            es(1,gp) = esh(1,gp)
            es(2,gp) = esh(2,gp)
            es(3,gp) = esh(4,gp)					 
         end do
         !print * , es

         ! ###################################################### !
         ! Compute element stiffness                            # !
         ! ###################################################### !
         call c2tl6_e(ke,coord(enod(ie,2:7),:),t,Dgp,u(edofn(ie,2:13)),es)

         ! ###################################################### !
         ! Compute internal forces                              # !
         ! ###################################################### !
         call c2tl6_f(fe,coord(enod(ie,2:7),:),t,es,u(edofn(ie,2:13)))

         ! ###################################################### !
         ! Printout                                             # !
         ! ###################################################### !
         !if(ie.eq.1) then
            !print * , coord(enod(ie,2:7),:) 
            !print * , dg(1,:,:)
            !print * , Dgp(:,:,1)
            !print * , esh(:,:)
            !print * , ke
         !endif 

         ! ###################################################### !
         ! Assemble                                             # !
         ! ###################################################### !
         !call assem(K,Ke,edofn(ie,2:13))
         
         call assem(K,Ke,enod(ie,2:7),dofnod)

         F(edofn(ie,2:13)) = F(edofn(ie,2:13)) + fe
      end do

      
      ! ###################################################### !
      ! Solve                                                # !
      ! ###################################################### !
      res    = F
      Deltaa = -res
      call solveq(K,Deltaa,bcdof,bcval)	
      u      = u + Deltaa
      
      ! Uppdatera F_int?
      F = 0d0
      do ie = 1, nelm
         ! ###################################################### !
         ! Compute the deformation gradient                     # !
         ! ###################################################### !
         call c2tl6_d(dg,u(edofn(ie,2:13)),coord(enod(ie,2:7),:))

!         ! ###################################################### !
         ! Reshape deformation gradient          (needed?)      # !
         ! ###################################################### !
         do gp=1,ngp
            ! Stämmer ordningen???
            ef(:,gp) = (/ dg(gp,1,1), dg(gp,1,2), dg(gp,2,1), dg(gp,2,2) /)
            !ef(:,gp) = (/ dg(gp,1,1), dg(gp,2,1), dg(gp,1,2), dg(gp,2,2) /)
         end do	
   
!         ! ###################################################### !
         ! Compute the material stiffness [nelm x ngp x 3 x 3]  # !
         ! ###################################################### !
         call dneohooke('tl',Dgp,ef)

!         ! ###################################################### !
         ! Compute the 2nd Piola Kirchhoff stress               # !
         ! ###################################################### !
         call neohooke('2ndPiola',esh,ef)     

!         ! ###################################################### !
         ! Reshape stress tensor                                # !
         ! ###################################################### !
         do gp=1,ngp
            es(1,gp) = esh(1,gp)
            es(2,gp) = esh(2,gp)
            es(3,gp) = esh(4,gp)					
            !print * , "Stress in element ", ie , "Gauss point" , gp , " = " , es(:,gp) 
         end do

!         ! ###################################################### !
         ! Compute internal forces                              # !
         ! ###################################################### !
         call c2tl6_f(fe,coord(enod(ie,2:7),:),t,es,u(edofn(ie,2:13)))

!         ! ###################################################### !
         ! Printout                     # !
         ! ###################################################### !
         !if(ie.eq.1) then
            !print * , coord(enod(ie,2:7),:) 
            !print * , dg(1,:,:)
            !print * , Dgp(:,:,1)
            !print * , esh(:,:)
            !print * , ke
         !endif 

!         ! ###################################################### !
         ! Assemble                                             # !
         ! ###################################################### !
         !call insert(F,fe,edofn(ie,2:13))
         F(edofn(ie,2:13)) = F(edofn(ie,2:13)) + fe
      end do
      ! ###################################################### !
      ! Efter systemet är löst                               # !
      ! ###################################################### !
      bcval = 0d0
      res = F 
      res(bcdof) = 0d0
      residual   = dsqrt(dot_product(res,res))
      print * ,  'Iteration:', iter, '; Residual:', residual
      if(iter.eq.imax) then
         print * , "----------------------- Maximum iteration count reached -----------------------" 
      end if

   end do
end do

bcval = bc0

end subroutine NewtonSolver

end module mesh_reader