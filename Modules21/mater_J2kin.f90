module mater_J2kin

! last modified
! M. Wallin 2011-05-24
!  - Project initiated
!  - J2kin_init implemented
!  - J2kin_accept implemented 
!  - J2kin        implemented
!  - dJ2kin       implemented
!  - J2kin_getVal implemented
!  - M. Wallin 2011-09-15.
!  - Several bugs removed and program tested.
!  - M. Wallin 2011-09-16
!  - syminv removed
!  - inv3, inv3s det3 moved to module matrix_util
!   Ref. Wallin et al. (2003). Eur. J. Mech. A/Solids Vol. 21, pp. 341-356   \\
!    Wallin and Ristinmaa (2005), Int. J. Plast. Vol. 21, pp. 2025-2050
!  2011-09-23 Mixed formulation introduced, not debugged
! M. Ristinmaa 2012-02-14
!  - J2kin_init arguments changed, stype is still present differs from J2iso

! TODO mater_J2kin is implemented differently than mater_J2iso
! stype is not given in calls


use mater_large    
use matrix_util
    
implicit none

                 

double precision                  :: SIGY0, MU, KAPPA, KSI1, GAMMA1,KSI2, GAMMA2
private SIGY0, MU, KAPPA, KSI1, GAMMA1,KSI2, GAMMA2



double precision, allocatable       :: Isv(:,:,:,:)
double precision, allocatable       :: dgn(:,:,:)
character(len=80)                   :: type_defgr
integer                             :: New, Old
double precision                    :: ep(7)
private Isv, dgn, New, Old, ep
	

interface J2kin_getVal
  module procedure J2kin_getValScalar
  module procedure J2kin_getValVector
end interface
private J2kin_getValScalar, J2kin_getValVector


interface J2kin_init
  module procedure init_J2_kin1
end interface
private init_J2_kin1

interface J2kin_accept
  module procedure accept_J2_kin1
  module procedure accept_J2_kin2
  module procedure accept_J2_kin3
end interface
private accept_J2_kin1

interface dJ2kin
  module procedure dJ2_kin0
  module procedure dJ2_kin1
  module procedure dJ2_kin2
end interface
private dJ2_kin0, dJ2_kin1, dJ2_kin2


interface J2kin
  module procedure J2_kin0
  module procedure J2_kin1
  module procedure J2_kin2
end interface
private J2_kin0, J2_kin1, J2_kin2


!-------------------------------------------------------------------------------
contains


subroutine J2kin_getValScalar(stype,out)
  implicit none
  integer                         :: ierr
  character(len=*)                :: stype
  double precision, allocatable   :: out(:,:)
 
  if (stype.eq.'DL') then
    out=Isv(28,:,:,new)
  elseif (stype.eq.'plastic') then
    out=0d0
    where(Isv(28,:,:,new).gt.0d0) out=1d0   !dl
  else
    stop "Variable does not exist"
  end if
  
  return
end subroutine J2kin_getValScalar


subroutine J2kin_getValVector(stype,out)
  implicit none
  integer                         :: ierr
  character(len=*)                :: stype
  double precision, allocatable   :: out(:,:,:)

  
  if (stype.eq.'FP') then
    out=Isv(1:9,:,:,new)
  elseif (stype.eq.'FK1') then
    out=Isv(10:18,:,:,new)
  elseif (stype.eq.'FK2') then
    out=Isv(19:27,:,:,new)
    stop "Variable does not exist"
  end if
  
  return
end subroutine J2kin_getValVector

subroutine init_J2_kin1(stype,mp,Nelm,Ngp)
  implicit none
  character(len=*)                :: stype
  integer                         :: Nisv, Ngp, Nelm, ierr
  double precision                :: mp(7)

  Nisv=28
  allocate(Isv(Nisv,Ngp,nelm,2), stat=ierr)   
  New=1
  Old=2
  
  SIGY0  =MP(1)
  MU     =MP(2)
  KAPPA  =MP(3)
  KSI1   =MP(4)
  GAMMA1 =MP(5)
  KSI2   =MP(6)
  GAMMA2 =MP(7)

  Isv=0D0

  Isv(1,:,:,:)=1D0
  Isv(5,:,:,:)=1D0
  Isv(9,:,:,:)=1D0
  Isv(10,:,:,:)=1D0
  Isv(14,:,:,:)=1D0
  Isv(18,:,:,:)=1D0
  Isv(19,:,:,:)=1D0
  Isv(23,:,:,:)=1D0
  Isv(27,:,:,:)=1D0
 
  if (stype.eq.'total') then
    allocate(dgn(9,Ngp,nelm), stat=ierr)
    dgn=0d0
    dgn((/1,5,9/),:,:)=1d0
    type_defgr='total'
  elseif (stype.eq.'relative') then
    allocate(dgn(9,Ngp,nelm), stat=ierr)
    dgn=0d0
    dgn((/1,5,9/),:,:)=1d0
    type_defgr='relativ'
  else
    stop "choice is not possible"
  end if

  return
end subroutine init_J2_kin1


subroutine accept_J2_kin1
  implicit none

  Isv(:,:,:,Old)=Isv(:,:,:,New)
  if (New.eq.1) then
    New=2
    Old=1
  else
    New=1
    Old=2
  endif

  return
end subroutine accept_J2_kin1

! one gauss point per element
subroutine accept_J2_kin2(dg)
  implicit none
  double precision                :: dg(:,:)

  double precision                :: ddg(9)
  integer                         :: ie, ig


  Isv(:,:,:,Old)=Isv(:,:,:,New)
  if (New.eq.1) then
    New=2
    Old=1
  else
    New=1
    Old=2
  endif

  if (type_defgr.eq.'total') then

    if (size(dg,1).eq.9) then
      dgn(:,1,:)=dg(:,:)
    else
      dgn(1,1,:)=dg(1,:)
      dgn(2,1,:)=dg(2,:)
      dgn(4,1,:)=dg(3,:)
      dgn(5,1,:)=dg(4,:)
    end if
  elseif (type_defgr.eq.'relative') then

     do ie=1,size(dg,2)
       if (size(dg,1).eq.9) then
         ddg=dg(:,ie)
       else
         ddg=0d0
         ddg(1)=dg(1,ie)
         ddg(2)=dg(2,ie)
         ddg(4)=dg(3,ie)
         ddg(5)=dg(4,ie)
         ddg(9)=1d0
       end if
       call upddefgr(dgn(:,1,ie),ddg,dgn(:,1,ie))
     enddo
  end if

  return
end subroutine accept_J2_kin2

! more than one gauss point per element
subroutine accept_J2_kin3(dg)
  implicit none
  double precision                :: dg(:,:,:)

  double precision                :: ddg(9)
  integer                         :: ie, ig

  Isv(:,:,:,Old)=Isv(:,:,:,New)

  if (New.eq.1) then
    New=2
    Old=1
  else
    New=1
    Old=2
  endif


  if (type_defgr.eq.'total') then
    if (size(dg,1).eq.9) then
      dgn=dg
    else
      dgn(1,:,:)=dg(1,:,:)
      dgn(2,:,:)=dg(2,:,:)
      dgn(4,:,:)=dg(3,:,:)
      dgn(5,:,:)=dg(4,:,:)
    end if
  elseif (type_defgr.eq.'relative') then
     do ie=1,size(dg,3)
       do ig=1,size(dg,2)
         if (size(dg,1).eq.9) then
           ddg=dg(:,ig,ie)
         else
           ddg=0d0
           ddg(1)=dg(1,ig,ie)
           ddg(2)=dg(2,ig,ie)
           ddg(4)=dg(3,ig,ie)
           ddg(5)=dg(4,ig,ie)
           ddg(9)=1d0
         end if
         call upddefgr(dgn(:,ig,ie),ddg,dgn(:,ig,ie))
       enddo
     enddo
  end if

  return
end subroutine accept_J2_kin3


subroutine J2_kin0(stype,stress,ef_New,ie)
  implicit none
  double precision                :: stress(:), ef_New(:)
  character(len=*)                :: stype
  integer                         :: gp, nr_gp,ie


  call J2_kin2(stype, stress(:), Isv(:,1,ie,old), Isv(:,1,ie,new), & 
                  ef_New(:),ie,1)

  return
end subroutine J2_kin0


subroutine J2_kin1(stype,stress,ef_New,ie)
  implicit none
  double precision                :: stress(:,:), ef_New(:,:)
  character(len=*)                :: stype
  integer                         :: gp, nr_gp,ie

  nr_gp=size(ef_New,2)

  do gp=1,nr_gp
    call J2_kin2(stype, stress(:,gp), Isv(:,gp,ie,old), Isv(:,gp,ie,new), & 
                  ef_New(:,gp),ie,gp)
  enddo

  return
end subroutine J2_kin1


subroutine J2_kin2(stype,stress,Isv_Old,Isv_New,ef_New,ie,gp)
  implicit none
  double precision                :: stress(:), Isv_Old(:), Isv_New(:), th
  double precision                :: ef_New(:)
  character(len=*)                :: stype  
  integer                         :: ie, gp, iter


       DOUBLE PRECISION                                 :: J, SIGEFF, TOL, NORM_DY
       DOUBLE PRECISION                                 :: DLAMB
       DOUBLE PRECISION, DIMENSION(3)                   :: DWE, D, GA, TMP1, GAMMASE, EIGENVALUES, DWK1, DWK2
       DOUBLE PRECISION, DIMENSION(19)                  :: Y, DY, RES, RES1
       DOUBLE PRECISION, DIMENSION(3,3)                 :: FP, FBAR1, FBAR2,IFP_O, F_NEW,F_NEW_HAT,FE_TRIAL_HAT
       DOUBLE PRECISION, DIMENSION(3,3)                 :: FP_OLD, FBAR1_OLD, FBAR2_OLD, B1_DEV, B2_DEV
       DOUBLE PRECISION, DIMENSION(3,3)                 :: CE_TRIAL_HAT,CBAR1_TRIAL,CBAR2_TRIAL,LOGCE_TRIAL_HAT,UWe,VWe
       DOUBLE PRECISION, DIMENSION(3,3)                 :: LOGCBAR1_TRIAL,LOGCBAR2_TRIAL, MANDEL_RED_DEV_TRIAL, CE_TR
       DOUBLE PRECISION, DIMENSION(3,3)                 :: KR,APINV, AP, AB1, AB2, S_BAR, MANDEL_DEV
       DOUBLE PRECISION, DIMENSION(3,3)                 :: CE_new_hat,CBAR1,CBAR2,LOGCE_HAT, U, V, NE,NK1, NK2
       DOUBLE PRECISION, DIMENSION(3,3)                 :: LOGCBAR1,LOGCBAR2,M_RED, BETA1, BETA2, LP, TMP33
       DOUBLE PRECISION, DIMENSION(3,3)                 :: FP_NEW, FBAR1_NEW, FBAR2_NEW, CE, ICE, MANDEL, IFP, C, IC, SE
       DOUBLE PRECISION, DIMENSION(19,19)               :: JAC, JAC_INV, JAC_TEST      
       double precision                                 :: S(3,3), C_NEW(3,3), C_NEW_HAT(3,3), CORRECTION

   TOL    =1D-14
   KR     =getI()

   if (type_defgr.eq.'total') then
     F_New=getF(ef_New)
   elseif (type_defgr.eq.'relative') then
     stop "Not implemented, err 82256"   
   end if   

  
   C_NEW      =MATMUL(TRANSPOSE(F_NEW), F_NEW) 
   CALL INV3S(IC, C_NEW)  

   J          =DSQRT(det3(C_NEW))  
   th=J ! Displacement based formulation
   C_NEW_HAT  =J**(-2D0/3D0)*C_NEW

   FP_OLD=   TRANSPOSE(RESHAPE((/Isv_Old(1:3),   Isv_Old(4:6),   Isv_Old(7:9)/),(/3,3/)))
   FBAR1_OLD=TRANSPOSE(RESHAPE((/Isv_Old(10:12), Isv_Old(13:15), Isv_Old(16:18)/),(/3,3/)))
   FBAR2_OLD=TRANSPOSE(RESHAPE((/Isv_Old(19:21), Isv_Old(22:24), Isv_Old(25:27)/),(/3,3/)))

   CALL INV3(IFP_O, FP_OLD)      
   CE_TRIAL_HAT  =MATMUL(MATMUL(TRANSPOSE(IFP_O), C_NEW_HAT),IFP_O)

   CBAR1_TRIAL   =MATMUL(TRANSPOSE(FBAR1_OLD),FBAR1_OLD)
   CBAR2_TRIAL   =MATMUL(TRANSPOSE(FBAR2_OLD),FBAR2_OLD)

   CALL PADELOG(CE_TRIAL_HAT,LOGCE_TRIAL_HAT )
   CALL PADELOG(CBAR1_TRIAL, LOGCBAR1_TRIAL)
   CALL PADELOG(CBAR2_TRIAL, LOGCBAR2_TRIAL)

   MANDEL_RED_DEV_TRIAL=MU*LOGCE_TRIAL_HAT-KSI1*LOGCBAR1_TRIAL-KSI2*LOGCBAR2_TRIAL
  
   TMP33=MATMUL(MANDEL_RED_DEV_TRIAL,MANDEL_RED_DEV_TRIAL)
   SIGEFF=DSQRT(3D0/2D0*(TMP33(1,1)+TMP33(2,2)+TMP33(3,3)))
   
   IF (SIGEFF>SIGY0) THEN  
	
      Y=(/CE_TRIAL_HAT(1,1), CE_TRIAL_HAT(2,2), CE_TRIAL_HAT(3,3), &
          CE_TRIAL_HAT(1,2), CE_TRIAL_HAT(1,3), CE_TRIAL_HAT(2,3), &
          CBAR1_TRIAL(1,1),  CBAR1_TRIAL(2,2),  CBAR1_TRIAL(3,3), &
          CBAR1_TRIAL(1,2),  CBAR1_TRIAL(1,3),  CBAR1_TRIAL(2,3), &
          CBAR2_TRIAL(1,1),  CBAR2_TRIAL(2,2),  CBAR2_TRIAL(3,3), &
          CBAR2_TRIAL(1,2),  CBAR2_TRIAL(1,3),  CBAR2_TRIAL(2,3), (sigeff-sigy0)*4.115D-6 /)

      NORM_DY=99D0
      ITER=0
      DO WHILE (NORM_DY>TOL)
          ITER=ITER+1
          CALL ANALYTICALJACOBIAN(RES,JAC,Y,CE_TRIAL_HAT,CBAR1_TRIAL,CBAR2_TRIAL,19,ITER)
 !         CALL NUMERICALJACOBIAN (RES,JAC,Y,CE_TRIAL_HAT,CBAR1_TRIAL,CBAR2_TRIAL,19,ITER)	                             
          CALL INVERSE(JAC,JAC_INV,19)
          DY=MATMUL(JAC_INV,RES)
          Y=Y-DY		                         
          CALL NORM(DY,NORM_DY,19) 
          if (iter.gt.100) then
             stop "Local system is diverging, ERRCODE 99"
          endif                                                                                           
      END DO  
                                
      CALL ANALYTICALJACOBIAN(RES,JAC,Y,CE_TRIAL_HAT,CBAR1_TRIAL,CBAR2_TRIAL,19,ITER)
!      CALL NUMERICALJACOBIAN(RES,JAC,Y,CE_TRIAL_HAT,CBAR1_TRIAL,CBAR2_TRIAL,19,ITER)	     
      CALL INVERSE(JAC,JAC_INV,19)           	     	     
      CE_NEW_HAT=RESHAPE((/Y(1), Y(4), Y(5), &
                           Y(4), Y(2), Y(6), &
                           Y(5), Y(6), Y(3)/),(/3, 3/))
   
      CBAR1=RESHAPE((/Y(7),  Y(10), Y(11), &
                      Y(10), Y(8),  Y(12), &
                      Y(11), Y(12), Y(9)/),(/3, 3/))
                
      CBAR2=RESHAPE((/Y(13), Y(16), Y(17), &
                      Y(16), Y(14), Y(18), &
                      Y(17), Y(18), Y(15)/),(/3, 3/))
	  
      DLAMB=Y(19)

      CALL  PADELOG(CBAR1,LOGCBAR1)
      CALL  PADELOG(CBAR2,LOGCBAR2)
      CALL  PADELOG(CE_NEW_HAT,LOGCE_HAT )
	 
      B1_DEV=KSI1*LOGCBAR1 
      B2_DEV=KSI2*LOGCBAR2
                  
      M_RED=MU*LOGCE_HAT-B1_DEV -B2_DEV
      LP=3D0/(2D0*SIGY0)*M_RED
      BETA1=LP-GAMMA1*B1_DEV
      BETA2=LP-GAMMA2*B2_DEV                
        
      CALL PADEEXP(LP*DLAMB, AP)
      CALL PADEEXP(BETA1*DLAMB,AB1)
      CALL PADEEXP(BETA2*DLAMB,AB2)     

      CALL INV3S(APINV,AP)

      FP_NEW   =MATMUL(AP,FP_OLD)
      FP_NEW   =((det3(FP_NEW))**(-1D0/3D0))*FP_NEW
      FBAR1_NEW=MATMUL(FBAR1_OLD,AB1)  
      FBAR1_NEW=((det3(FBAR1_NEW))**(-1D0/3D0))*FBAR1_NEW
      FBAR2_NEW=MATMUL(FBAR2_OLD,AB2)    
      FBAR2_NEW=((det3(FBAR2_NEW))**(-1D0/3D0))*FBAR2_NEW
                				             	
      Isv_New(1:3)  =FP_NEW(1,1:3)
      Isv_New(4:6)  =FP_NEW(2,1:3)
      Isv_New(7:9)  =FP_NEW(3,1:3)
      Isv_New(10:12)=FBAR1_NEW(1,1:3)
      Isv_New(13:15)=FBAR1_NEW(2,1:3)
      Isv_New(16:18)=FBAR1_NEW(3,1:3)
      Isv_New(19:21)=FBAR2_NEW(1,1:3)
      Isv_New(22:24)=FBAR2_NEW(2,1:3)
      Isv_New(25:27)=FBAR2_NEW(3,1:3)
      Isv_New(28)   =DLAMB
   ELSE             
      Isv_New(:)     =Isv_Old   
      Isv_New(28)    =0D0
      FP_NEW         =FP_OLD
      CE_NEW_HAT     =CE_TRIAL_HAT
      LOGCE_HAT      =LOGCE_TRIAL_HAT
   ENDIF
		  
    CE   =(J**(2D0/3D0))*CE_NEW_HAT
    CALL INV3S(ICE,CE)
!    MANDEL=MU*LOGCE_HAT+KAPPA*DLOG(J)*KR  
    MANDEL=MU*LOGCE_HAT+KAPPA*DLOG(th)*KR  
    CALL INV3(IFP,FP_NEW)
    S=MATMUL(MATMUL(IFP,ICE),MATMUL(MANDEL,TRANSPOSE(IFP))) 

    MANDEL_DEV=MU*LOGCE_HAT
    CALL INV3(IFP,FP_NEW)
    S_BAR=MATMUL(MATMUL(IFP,ICE),MATMUL(MANDEL_DEV,TRANSPOSE(IFP))) 
    S=S_BAR+KAPPA*DLOG(J)*IC   

    if (stype.eq.'Cauchy') then
      S=MATMUL(MATMUL(F_NEW,S),TRANSPOSE(F_NEW))/J
    elseif (stype.eq.'2ndPiola') then
    ! Do nothing, S already calculated
    elseif (stype.eq.'Kirchhoff') then
      S=MATMUL(MATMUL(F_NEW,S),TRANSPOSE(F_NEW))/J
    else
      stop 'stype not implemented'
    endif


    if (size(ef_New).eq.4) then
      stress=(/S(1,1),S(2,2),S(3,3),S(1,2)/)
    else
      stress=(/S(1,1),S(2,2),S(3,3),S(1,2),S(1,3),S(2,3)/)
    endif

end subroutine J2_kin2


subroutine dJ2_kin0(stype,D,ef_New,ie,iter)
  implicit none
  character(len=*)                :: stype
  double precision                :: D(:,:), ef_New(:)
  integer                         :: gp,nr_gp, ie, iter

   call dJ2_kin2(stype,D(:,:),ef_New(:), &
                  Isv(:,1,ie,New), Isv(:,1,ie,Old),ie,1)
 

  return
end subroutine dJ2_kin0


subroutine dJ2_kin1(stype,D,ef_New,ie,iter)
  implicit none
  character(len=*)                :: stype
  double precision                :: D(:,:,:), ef_New(:,:)
  integer                         :: gp,nr_gp, ie,iter

  nr_gp=size(D,3)
 
  do gp=1,nr_gp
    call dJ2_kin2(stype,D(:,:,gp),ef_New(:,gp), &
                  Isv(:,gp,ie,Old), Isv(:,gp,ie,New),ie,gp)
  enddo

  return
end subroutine dJ2_kin1


subroutine dJ2_kin2(stype,Dout,ef_New,Isv_Old,Isv_New,ie,gp)
  implicit none
  character(len=*)                :: stype
  double precision                :: ef_New(:)
  double precision                :: Isv_New(:), Isv_Old(:), Dout(:,:) 
  integer                         :: ie, gp, ITER, M, N


  DOUBLE PRECISION                                 :: J, SIGEFF, TOL, NORM_DY, PER, R0, R1, ACC, th
  DOUBLE PRECISION                                 :: DLAMB
  DOUBLE PRECISION, DIMENSION(19)                  :: Y, DY, RES, RES1
  DOUBLE PRECISION, DIMENSION(28)                  :: Y_OLD, Y_NEW
  DOUBLE PRECISION, DIMENSION(3,3)                 :: FP, FBAR1, FBAR2,IFP_O, F_NEW,F_NEW_HAT,FE_TRIAL_HAT
  DOUBLE PRECISION, DIMENSION(3,3)                 :: FP_OLD, FBAR1_OLD, FBAR2_OLD, FE_NEW_HAT
  DOUBLE PRECISION, DIMENSION(3,3)                 :: CE_TR_HAT,CBAR1_TRIAL,CBAR2_TRIAL,LOGCE_TRIAL_HAT
  DOUBLE PRECISION, DIMENSION(3,3)                 :: LOGCBAR1_TRIAL,LOGCBAR2_TRIAL, MANDEL_RED_DEV_TRIAL, CE_TR
  DOUBLE PRECISION, DIMENSION(3,3)                 :: APINV, AP, AB1, AB2
  DOUBLE PRECISION, DIMENSION(3,3)                 :: CE_new_hat, CBAR1, CBAR2,LOGCE_HAT, NE,NK1, NK2
  DOUBLE PRECISION, DIMENSION(3,3)                 :: LOGCBAR1, LOGCBAR2, B1_DEV,B2_DEV,M_RED, BETA1, BETA2, LP
  DOUBLE PRECISION, DIMENSION(3,3)                 :: FP_NEW, FBAR1_NEW, FBAR2_NEW, CE, ICE, MANDEL, IFP, C, IC, SE
  DOUBLE PRECISION, DIMENSION(19,19)               :: JAC, JAC_INV, JAC_TEST      
  DOUBLE PRECISION, DIMENSION(6,6)                 :: DIAG_UNIT, DRCEDCE_TRIAL_HAT, DCE_HAT_TRIAL_DC,DCE_DC
  DOUBLE PRECISION, DIMENSION(6,6)                 :: PUSH_FW, TMP66, STIFF, TEST
  DOUBLE PRECISION, DIMENSION(6,1)                 :: TMP_M, TMP_N,  CE_IM, CE_M, C_IM, CE_TRM

  DOUBLE PRECISION, DIMENSION(19,6)                :: DRDCE_TRIAL_HAT, DYDCETRIAL
  DOUBLE PRECISION, DIMENSION(9)                   :: TMP_9

  DOUBLE PRECISION, DIMENSION(6,6)                 :: DE_HYD, CEI_X, CEI_KR, DLOGCEHAT_DCEHAT, DCEHAT_DCE, DE_DEV, XFP1, XFP2, XFP
  DOUBLE PRECISION, DIMENSION(6,6)                 :: DLOGCBAR1_DCBAR1, DLOGCBAR2_DCBAR2, Conny, DEXPM, MW66, Part2, QNEW,FPISFPI
  DOUBLE PRECISION, DIMENSION(3,3)                 :: X, EXPM,SFPI, C_NEW_HAT, TMP33, KR, FE_NEW, S, SE_BAR
  DOUBLE PRECISION, DIMENSION(6,6)                 :: TMP66_1,TMP66_2, TMP66_3, PUSH, KR_DYAD_KR, D_HYD	
  DOUBLE PRECISION, DIMENSION(6,1)                 :: TMP_61
  DOUBLE PRECISION, DIMENSION(6)                   :: TMP_6


   Y_OLD=Isv_Old
   Y_NEW=Isv_New  

   if (type_defgr.eq.'total') then
     F_NEW=getF(ef_New)
   elseif (type_defgr.eq.'relative') then
     stop "Not implemented, err 82256"   
   end if   

   KR     =getI()
   CALL EYE(DIAG_UNIT,6)
   DIAG_UNIT(4,4)=.5D0
   DIAG_UNIT(5,5)=.5D0
   DIAG_UNIT(6,6)=.5D0

   KR_DYAD_KR(1:3,1:3)=1D0

   DRDCE_TRIAL_HAT=0D0

   C=MATMUL(TRANSPOSE(F_NEW),F_NEW)         
   CALL INV3S(IC,C)

   J=DSQRT(det3(C))   
   th=J ! Displacmenet based formulation 
   C_NEW_HAT=J**(-2D0/3D0)*C
 
   FP_OLD=TRANSPOSE(   RESHAPE((/Y_OLD(1:3),   Y_OLD(4:6),   Y_OLD(7:9)/),(/3,3/)))
   FBAR1_OLD=TRANSPOSE(RESHAPE((/Y_OLD(10:12), Y_OLD(13:15), Y_OLD(16:18)/),(/3,3/)))
   FBAR2_OLD=TRANSPOSE(RESHAPE((/Y_OLD(19:21), Y_OLD(22:24), Y_OLD(25:27)/),(/3,3/)))
   FP_NEW=TRANSPOSE(   RESHAPE((/Y_NEW(1:3),   Y_NEW(4:6),   Y_NEW(7:9)/),(/3,3/)))
   FBAR1_NEW=TRANSPOSE(RESHAPE((/Y_NEW(10:12), Y_NEW(13:15), Y_NEW(16:18)/),(/3,3/)))
   FBAR2_NEW=TRANSPOSE(RESHAPE((/Y_NEW(19:21), Y_NEW(22:24), Y_NEW(25:27)/),(/3,3/)))
   DLAMB=Y_NEW(28)

   CALL INV3(IFP,FP_NEW)  
   FE_NEW=MATMUL(F_NEW, IFP)

   CE    =MATMUL(TRANSPOSE(FE_NEW), FE_NEW)
   CALL INV3S(ICE,CE)

   TMP_9=(/F_NEW(1,1), F_NEW(1,2), F_NEW(1,3), &
           F_NEW(2,1), F_NEW(2,2), F_NEW(2,3), &
           F_NEW(3,1), F_NEW(3,2), F_NEW(3,3) /)

   FE_NEW_HAT=J**(-1D0/3D0)*FE_NEW
   CE_NEW_HAT  =MATMUL(TRANSPOSE(FE_NEW_HAT), FE_NEW_HAT)
  
   CBAR1       =MATMUL(TRANSPOSE(FBAR1_NEW), FBAR1_NEW)
   CBAR2       =MATMUL(TRANSPOSE(FBAR2_NEW), FBAR2_NEW)

   CALL INV3(IFP_O,FP_OLD)  
   CE_TR_HAT  =MATMUL(MATMUL(TRANSPOSE(IFP_O), C_NEW_HAT),IFP_O)

   CBAR1_TRIAL   =MATMUL(TRANSPOSE(FBAR1_OLD),FBAR1_OLD)
   CBAR2_TRIAL   =MATMUL(TRANSPOSE(FBAR2_OLD),FBAR2_OLD)

   CE_TRM=J**(2D0/3D0)*RESHAPE((/CE_TR_HAT(1,1), CE_TR_HAT(2,2), CE_TR_HAT(3,3),CE_TR_HAT(1,2), CE_TR_HAT(1,3), CE_TR_HAT(2,3)/),(/6,1/))

   TMP66=KAPPA*(0.5D0*KR_DYAD_KR-DIAG_UNIT*DLOG(th))
   call pullback(D_HYD, TMP66,TMP_9)

   Y=(/CE_NEW_HAT(1,1), CE_NEW_HAT(2,2), CE_NEW_HAT(3,3), &
       CE_NEW_HAT(1,2), CE_NEW_HAT(1,3), CE_NEW_HAT(2,3), &
       CBAR1(1,1),  CBAR1(2,2),  CBAR1(3,3), &
       CBAR1(1,2),  CBAR1(1,3),  CBAR1(2,3), &
       CBAR2(1,1),  CBAR2(2,2),  CBAR2(3,3), &
       CBAR2(1,2),  CBAR2(1,3),  CBAR2(2,3), DLAMB /)
  
    CE_IM  =RESHAPE((/ICE(1,1), ICE(2,2), ICE(3,3), ICE(1,2), ICE(1,3), ICE(2,3)/),(/6,1/))
    CE_M   =RESHAPE((/ CE(1,1),  CE(2,2),  CE(3,3),  CE(1,2),  CE(1,3),  CE(2,3)/),(/6,1/))
    C_IM   =RESHAPE((/ IC(1,1),  IC(2,2),  IC(3,3),  IC(1,2),  IC(1,3),  IC(2,3)/),(/6,1/))      
             
    CALL PADELOG(CE_NEW_HAT, LOGCE_HAT)

    X=MATMUL(ICE,LOGCE_HAT)        
    CALL BAR_DYAD(CEI_X,ICE,X)
    CEI_X(:,4:6)=0.5D0*CEI_X(:,4:6)        
            
    CALL BAR_DYAD(CEI_KR,ICE,KR)

    CALL DPADELOG(CE_NEW_HAT, DLOGCEHAT_DCEHAT)	
        	     
    TMP66=MATMUL(CEI_KR,DLOGCEHAT_DCEHAT) 

    DCEHAT_DCE=J**(-2D0/3D0)*(DIAG_UNIT-1D0/3D0*MATMUL(CE_M,TRANSPOSE(CE_IM)))
       
    DE_DEV=MU*(-CEI_X+MATMUL(TMP66,DCEHAT_DCE))
   
    CALL BAR_DYAD(PUSH_FW,IFP,IFP)

    IF (DLAMB>0D0) THEN
 !          CALL NUMERICALJACOBIAN (RES,JAC,Y,CE_TR_HAT,CBAR1_TRIAL,CBAR2_TRIAL,19,ITER)	      
            CALL ANALYTICALJACOBIAN(RES,JAC,Y,CE_TR_HAT,CBAR1_TRIAL,CBAR2_TRIAL,19,ITER) 
            CALL INVERSE(JAC,JAC_INV,19)    

            CALL PADELOG(CBAR1, LOGCBAR1)
            CALL PADELOG(CBAR2, LOGCBAR2)

            B1_DEV=KSI1*LOGCBAR1 
            B2_DEV=KSI2*LOGCBAR2
                 
            M_RED=MU*LOGCE_HAT-B1_DEV-B2_DEV
            LP=3D0/(2D0*SIGY0)*M_RED
                
            CALL PADEEXP(LP*DLAMB, AP)
            CALL INV3(APINV,AP)  

            CALL DPADELOG(CBAR1, DLOGCBAR1_DCBAR1)
            CALL DPADELOG(CBAR2, DLOGCBAR2_DCBAR2)
            CALL DPADEEXP(-LP*DLAMB,DEXPM)
            
            CALL BAR_DYAD(PUSH,APINV,APINV)
  
            DRCEDCE_TRIAL_HAT=-MATMUL(PUSH,DIAG_UNIT)                  
            DRDCE_TRIAL_HAT(1:6,:)=DRCEDCE_TRIAL_HAT	       
            DYDCETRIAL=-MATMUL(JAC_INV,DRDCE_TRIAL_HAT)  

            CALL BAR_DYAD(PUSH,TRANSPOSE(IFP_o),TRANSPOSE(IFP_o))
            DCE_HAT_TRIAL_DC=(J**(-2D0/3D0))*(MATMUL(PUSH,DIAG_UNIT)-1D0/3D0*MATMUL(CE_TRM,TRANSPOSE(C_IM)))

            DCE_HAT_TRIAL_DC(4:6,:)=2D0*DCE_HAT_TRIAL_DC(4:6,:)
      
            TMP66=DYDCETRIAL(1:6,1:6)
            DCE_DC=J**(2D0/3D0)*MATMUL(TMP66,DCE_HAT_TRIAL_DC)+1D0/3D0*MATMUL(CE_M,TRANSPOSE(C_IM)) 
	    DCE_DC(4:6,:)=2D0*DCE_DC(4:6,:)

            TMP_N=RESHAPE((/dYdCetrial(19,1), dYdCetrial(19,2), dYdCetrial(19,3), &
                            dYdCetrial(19,4), dYdCetrial(19,5), dYdCetrial(19,6)/),(/6,1/))			  

            TMP_61=reshape((/LP(1,1), LP(2,2), LP(3,3), LP(1,2), LP(3,1),  LP(3,2)/),(/6,1/))
                      
            TMP66_1 =MATMUL(DLOGCEHAT_DCEHAT, dYdCetrial(1:6,:))	                      
            TMP66_2 =MATMUL(DLOGCBAR1_DCBAR1,dYdCetrial(7:12,:))
            TMP66_3 =MATMUL(DLOGCBAR2_DCBAR2,dYdCetrial(13:18,:))
         
            Qnew=3D0*DLAMB/(2D0*SIGY0)*(KSI1*TMP66_2+KSI2*TMP66_3-MU*TMP66_1)-MATMUL(TMP_61,TRANSPOSE(TMP_N))
          
            SE_BAR=MU*MATMUL(ICE,LOGCE_HAT)        
            SFpI=MATMUL(IFP, SE_BAR) 
  
            CALL BAR_DYAD(XFP1,IFP_O,SFpi)
            CALL BAR_DYAD(XFP2,SFpi,IFP_O)
            XFP=(XFP1+XFP2)
	                           
            Part2=MATMUL(QNEW,dCe_hat_trial_dC)
            Part2=MATMUL(DEXPM,PART2)
            PART2=MATMUL(XFP,PART2)         
            STIFF=2D0*(MATMUL(MATMUL(PUSH_FW,DE_DEV),DCE_DC)+Part2+D_HYD)            
        ELSE
           STIFF=2D0*(MATMUL(MATMUL(PUSH_FW,DE_DEV),TRANSPOSE(PUSH_FW))+D_HYD)         
        ENDIF       

  if (size(ef_New).eq.4) then
    Dout(1,:)=(/STIFF(1,1), STIFF(2,2), STIFF(1,4)/)
    Dout(2,:)=(/STIFF(2,1), STIFF(2,2), STIFF(2,4)/)
    Dout(3,:)=(/STIFF(4,1), STIFF(4,2), STIFF(4,4)/)
  else
    Dout=STIFF
  endif 
  

  if (stype.eq.'ul') then
    call pushforward(Dout,Dout,ef_new)
  elseif (stype.eq.'tl') then
  else
  endif
 

end subroutine dJ2_kin2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UP
subroutine J2_kin0up(stype,stress,ef_New,ie,th)
  implicit none
  double precision                :: stress(:), ef_New(:), th
  character(len=*)                :: stype
  integer                         :: gp, nr_gp,ie


  call J2_kin2up(stype, stress(:), Isv(:,1,ie,old), Isv(:,1,ie,new), & 
                  ef_New(:),ie,1,th)

  return
end subroutine J2_kin0up


subroutine J2_kin1up(stype,stress,ef_New,ie,th)
  implicit none
  double precision                :: stress(:,:), ef_New(:,:), th(:)
  character(len=*)                :: stype
  integer                         :: gp, nr_gp,ie

  nr_gp=size(ef_New,2)

  do gp=1,nr_gp
    call J2_kin2up(stype, stress(:,gp), Isv(:,gp,ie,old), Isv(:,gp,ie,new), & 
                  ef_New(:,gp),ie,gp,th(gp))
  enddo

  return
end subroutine J2_kin1up


subroutine J2_kin2up(stype,stress,Isv_Old,Isv_New,ef_New,ie,gp,th)
  implicit none
  double precision                :: stress(:), Isv_Old(:), Isv_New(:), th
  double precision                :: ef_New(:)
  character(len=*)                :: stype  
  integer                         :: ie, gp, iter

       DOUBLE PRECISION                                 :: J, SIGEFF, TOL, NORM_DY
       DOUBLE PRECISION                                 :: DLAMB
       DOUBLE PRECISION, DIMENSION(3)                   :: DWE, D, GA, TMP1, GAMMASE, EIGENVALUES, DWK1, DWK2
       DOUBLE PRECISION, DIMENSION(19)                  :: Y, DY, RES, RES1
       DOUBLE PRECISION, DIMENSION(3,3)                 :: FP, FBAR1, FBAR2,IFP_O, F_NEW,F_NEW_HAT,FE_TRIAL_HAT
       DOUBLE PRECISION, DIMENSION(3,3)                 :: FP_OLD, FBAR1_OLD, FBAR2_OLD, B1_DEV, B2_DEV
       DOUBLE PRECISION, DIMENSION(3,3)                 :: CE_TRIAL_HAT,CBAR1_TRIAL,CBAR2_TRIAL,LOGCE_TRIAL_HAT,UWe,VWe
       DOUBLE PRECISION, DIMENSION(3,3)                 :: LOGCBAR1_TRIAL,LOGCBAR2_TRIAL, MANDEL_RED_DEV_TRIAL, CE_TR
       DOUBLE PRECISION, DIMENSION(3,3)                 :: KR,APINV, AP, AB1, AB2, S_BAR, MANDEL_DEV
       DOUBLE PRECISION, DIMENSION(3,3)                 :: CE_new_hat,CBAR1,CBAR2,LOGCE_HAT, U, V, NE,NK1, NK2
       DOUBLE PRECISION, DIMENSION(3,3)                 :: LOGCBAR1,LOGCBAR2,M_RED, BETA1, BETA2, LP, TMP33
       DOUBLE PRECISION, DIMENSION(3,3)                 :: FP_NEW, FBAR1_NEW, FBAR2_NEW, CE, ICE, MANDEL, IFP, C, IC, SE
       DOUBLE PRECISION, DIMENSION(19,19)               :: JAC, JAC_INV, JAC_TEST      
       double precision                                 :: S(3,3), C_NEW(3,3), C_NEW_HAT(3,3), CORRECTION

   TOL    =1D-14
   KR     =getI()

   if (type_defgr.eq.'total') then
     F_New=getF(ef_New)
   elseif (type_defgr.eq.'relative') then
     stop "Not implemented, err 82256"   
   end if   

  
   C_NEW      =MATMUL(TRANSPOSE(F_NEW), F_NEW) 
   CALL INV3S(IC, C_NEW)  

   J          =DSQRT(det3(C_NEW))  
 !  th=J ! Mixed formulation since th differs from J
   C_NEW_HAT  =J**(-2D0/3D0)*C_NEW

   FP_OLD=   TRANSPOSE(RESHAPE((/Isv_Old(1:3),   Isv_Old(4:6),   Isv_Old(7:9)/),(/3,3/)))
   FBAR1_OLD=TRANSPOSE(RESHAPE((/Isv_Old(10:12), Isv_Old(13:15), Isv_Old(16:18)/),(/3,3/)))
   FBAR2_OLD=TRANSPOSE(RESHAPE((/Isv_Old(19:21), Isv_Old(22:24), Isv_Old(25:27)/),(/3,3/)))

   CALL INV3(IFP_O, FP_OLD)      
   CE_TRIAL_HAT  =MATMUL(MATMUL(TRANSPOSE(IFP_O), C_NEW_HAT),IFP_O)

   CBAR1_TRIAL   =MATMUL(TRANSPOSE(FBAR1_OLD),FBAR1_OLD)
   CBAR2_TRIAL   =MATMUL(TRANSPOSE(FBAR2_OLD),FBAR2_OLD)

   CALL PADELOG(CE_TRIAL_HAT,LOGCE_TRIAL_HAT )
   CALL PADELOG(CBAR1_TRIAL, LOGCBAR1_TRIAL)
   CALL PADELOG(CBAR2_TRIAL, LOGCBAR2_TRIAL)

   MANDEL_RED_DEV_TRIAL=MU*LOGCE_TRIAL_HAT-KSI1*LOGCBAR1_TRIAL-KSI2*LOGCBAR2_TRIAL
  
   TMP33=MATMUL(MANDEL_RED_DEV_TRIAL,MANDEL_RED_DEV_TRIAL)
   SIGEFF=DSQRT(3D0/2D0*(TMP33(1,1)+TMP33(2,2)+TMP33(3,3)))
   
   IF (SIGEFF>SIGY0) THEN  
	
      Y=(/CE_TRIAL_HAT(1,1), CE_TRIAL_HAT(2,2), CE_TRIAL_HAT(3,3), &
          CE_TRIAL_HAT(1,2), CE_TRIAL_HAT(1,3), CE_TRIAL_HAT(2,3), &
          CBAR1_TRIAL(1,1),  CBAR1_TRIAL(2,2),  CBAR1_TRIAL(3,3), &
          CBAR1_TRIAL(1,2),  CBAR1_TRIAL(1,3),  CBAR1_TRIAL(2,3), &
          CBAR2_TRIAL(1,1),  CBAR2_TRIAL(2,2),  CBAR2_TRIAL(3,3), &
          CBAR2_TRIAL(1,2),  CBAR2_TRIAL(1,3),  CBAR2_TRIAL(2,3), (sigeff-sigy0)*4.115D-6 /)

      NORM_DY=99D0
      ITER=0
      DO WHILE (NORM_DY>TOL)
          ITER=ITER+1
          CALL ANALYTICALJACOBIAN(RES,JAC,Y,CE_TRIAL_HAT,CBAR1_TRIAL,CBAR2_TRIAL,19,ITER)
 !         CALL NUMERICALJACOBIAN (RES,JAC,Y,CE_TRIAL_HAT,CBAR1_TRIAL,CBAR2_TRIAL,19,ITER)	                             
          CALL INVERSE(JAC,JAC_INV,19)
          DY=MATMUL(JAC_INV,RES)
          Y=Y-DY		                         
          CALL NORM(DY,NORM_DY,19) 
          if (iter.gt.100) then
             stop "Local system is diverging, ERRCODE 99"
          endif                                                                                           
      END DO  
                                
      CALL ANALYTICALJACOBIAN(RES,JAC,Y,CE_TRIAL_HAT,CBAR1_TRIAL,CBAR2_TRIAL,19,ITER)
!      CALL NUMERICALJACOBIAN(RES,JAC,Y,CE_TRIAL_HAT,CBAR1_TRIAL,CBAR2_TRIAL,19,ITER)	     
      CALL INVERSE(JAC,JAC_INV,19)           	     	     
      CE_NEW_HAT=RESHAPE((/Y(1), Y(4), Y(5), &
                           Y(4), Y(2), Y(6), &
                           Y(5), Y(6), Y(3)/),(/3, 3/))
   
      CBAR1=RESHAPE((/Y(7),  Y(10), Y(11), &
                      Y(10), Y(8),  Y(12), &
                      Y(11), Y(12), Y(9)/),(/3, 3/))
                
      CBAR2=RESHAPE((/Y(13), Y(16), Y(17), &
                      Y(16), Y(14), Y(18), &
                      Y(17), Y(18), Y(15)/),(/3, 3/))
	  
      DLAMB=Y(19)

      CALL  PADELOG(CBAR1,LOGCBAR1)
      CALL  PADELOG(CBAR2,LOGCBAR2)
      CALL  PADELOG(CE_NEW_HAT,LOGCE_HAT )
	 
      B1_DEV=KSI1*LOGCBAR1 
      B2_DEV=KSI2*LOGCBAR2
                  
      M_RED=MU*LOGCE_HAT-B1_DEV -B2_DEV
      LP=3D0/(2D0*SIGY0)*M_RED
      BETA1=LP-GAMMA1*B1_DEV
      BETA2=LP-GAMMA2*B2_DEV                
        
      CALL PADEEXP(LP*DLAMB, AP)
      CALL PADEEXP(BETA1*DLAMB,AB1)
      CALL PADEEXP(BETA2*DLAMB,AB2)     

      CALL INV3S(APINV,AP)

      FP_NEW   =MATMUL(AP,FP_OLD)
      FP_NEW   =((det3(FP_NEW))**(-1D0/3D0))*FP_NEW
      FBAR1_NEW=MATMUL(FBAR1_OLD,AB1)  
      FBAR1_NEW=((det3(FBAR1_NEW))**(-1D0/3D0))*FBAR1_NEW
      FBAR2_NEW=MATMUL(FBAR2_OLD,AB2)    
      FBAR2_NEW=((det3(FBAR2_NEW))**(-1D0/3D0))*FBAR2_NEW
                				             	
      Isv_New(1:3)  =FP_NEW(1,1:3)
      Isv_New(4:6)  =FP_NEW(2,1:3)
      Isv_New(7:9)  =FP_NEW(3,1:3)
      Isv_New(10:12)=FBAR1_NEW(1,1:3)
      Isv_New(13:15)=FBAR1_NEW(2,1:3)
      Isv_New(16:18)=FBAR1_NEW(3,1:3)
      Isv_New(19:21)=FBAR2_NEW(1,1:3)
      Isv_New(22:24)=FBAR2_NEW(2,1:3)
      Isv_New(25:27)=FBAR2_NEW(3,1:3)
      Isv_New(28)   =DLAMB
   ELSE             
      Isv_New(:)     =Isv_Old   
      Isv_New(28)    =0D0
      FP_NEW         =FP_OLD
      CE_NEW_HAT     =CE_TRIAL_HAT
      LOGCE_HAT      =LOGCE_TRIAL_HAT
   ENDIF
		  
    CE   =(J**(2D0/3D0))*CE_NEW_HAT
    CALL INV3S(ICE,CE)
!    MANDEL=MU*LOGCE_HAT+KAPPA*DLOG(J)*KR  
    MANDEL=MU*LOGCE_HAT+KAPPA*DLOG(th)*KR  
    CALL INV3(IFP,FP_NEW)
    S=MATMUL(MATMUL(IFP,ICE),MATMUL(MANDEL,TRANSPOSE(IFP))) 

    MANDEL_DEV=MU*LOGCE_HAT
    CALL INV3(IFP,FP_NEW)
    S_BAR=MATMUL(MATMUL(IFP,ICE),MATMUL(MANDEL_DEV,TRANSPOSE(IFP))) 
    S=S_BAR+KAPPA*DLOG(th)*IC  ! Note that the stress is calculated based on th 

    if (stype.eq.'Cauchy') then
      S=MATMUL(MATMUL(F_NEW,S),TRANSPOSE(F_NEW))/J
    elseif (stype.eq.'2ndPiola') then
    ! Do nothing, S already calculated
    elseif (stype.eq.'Kirchhoff') then
      S=MATMUL(MATMUL(F_NEW,S),TRANSPOSE(F_NEW))/J
    else
      stop 'stype not implemented'
    endif


    if (size(ef_New).eq.4) then
      stress=(/S(1,1),S(2,2),S(3,3),S(1,2)/)
    else
      stress=(/S(1,1),S(2,2),S(3,3),S(1,2),S(1,3),S(2,3)/)
    endif

end subroutine J2_kin2up


subroutine dJ2_kin0up(stype,D_BAR,D_HYD,ef_New,ie,iter,th)
  implicit none
  character(len=*)                :: stype
  double precision                :: D_BAR(:,:), D_HYD, ef_New(:),th
  integer                         :: gp,nr_gp, ie, iter

   call dJ2_kin2up(stype,D_BAR(:,:),D_HYD,ef_New(:), &
                  Isv(:,1,ie,New), Isv(:,1,ie,Old),ie,1,th)
 

  return
end subroutine dJ2_kin0up


subroutine dJ2_kin1up(stype,D_BAR,D_HYD,ef_New,ie,iter,th)
  implicit none
  character(len=*)                :: stype
  double precision                :: D_BAR(:,:,:),D_HYD(:), ef_New(:,:), th(:)
  integer                         :: gp,nr_gp, ie,iter

  nr_gp=size(D_BAR,3)
 
  do gp=1,nr_gp
    call dJ2_kin2up(stype,D_BAR(:,:,gp),D_HYD(gp),ef_New(:,gp), &
                  Isv(:,gp,ie,Old), Isv(:,gp,ie,New),ie,gp,th(gp))
  enddo

  return
end subroutine dJ2_kin1up


subroutine dJ2_kin2up(stype,Dout_BAR,Dout_HYD,ef_New,Isv_Old,Isv_New,ie,gp,th)
  implicit none
  character(len=*)                :: stype
  double precision                :: ef_New(:)
  double precision                :: Isv_New(:), Isv_Old(:), Dout_BAR(:,:) , Dout_HYD 
  integer                         :: ie, gp, ITER, M, N


  DOUBLE PRECISION                                 :: J, SIGEFF, TOL, NORM_DY, PER, R0, R1, ACC, th
  DOUBLE PRECISION                                 :: DLAMB
  DOUBLE PRECISION, DIMENSION(19)                  :: Y, DY, RES, RES1
  DOUBLE PRECISION, DIMENSION(28)                  :: Y_OLD, Y_NEW
  DOUBLE PRECISION, DIMENSION(3,3)                 :: FP, FBAR1, FBAR2,IFP_O, F_NEW,F_NEW_HAT,FE_TRIAL_HAT
  DOUBLE PRECISION, DIMENSION(3,3)                 :: FP_OLD, FBAR1_OLD, FBAR2_OLD, FE_NEW_HAT
  DOUBLE PRECISION, DIMENSION(3,3)                 :: CE_TR_HAT,CBAR1_TRIAL,CBAR2_TRIAL,LOGCE_TRIAL_HAT
  DOUBLE PRECISION, DIMENSION(3,3)                 :: LOGCBAR1_TRIAL,LOGCBAR2_TRIAL, MANDEL_RED_DEV_TRIAL, CE_TR
  DOUBLE PRECISION, DIMENSION(3,3)                 :: APINV, AP, AB1, AB2
  DOUBLE PRECISION, DIMENSION(3,3)                 :: CE_new_hat, CBAR1, CBAR2,LOGCE_HAT, NE,NK1, NK2
  DOUBLE PRECISION, DIMENSION(3,3)                 :: LOGCBAR1, LOGCBAR2, B1_DEV,B2_DEV,M_RED, BETA1, BETA2, LP
  DOUBLE PRECISION, DIMENSION(3,3)                 :: FP_NEW, FBAR1_NEW, FBAR2_NEW, CE, ICE, MANDEL, IFP, C, IC, SE
  DOUBLE PRECISION, DIMENSION(19,19)               :: JAC, JAC_INV, JAC_TEST      
  DOUBLE PRECISION, DIMENSION(6,6)                 :: DIAG_UNIT, DRCEDCE_TRIAL_HAT, DCE_HAT_TRIAL_DC,DCE_DC
  DOUBLE PRECISION, DIMENSION(6,6)                 :: PUSH_FW, TMP66, STIFF, TEST
  DOUBLE PRECISION, DIMENSION(6,1)                 :: TMP_M, TMP_N,  CE_IM, CE_M, C_IM, CE_TRM

  DOUBLE PRECISION, DIMENSION(19,6)                :: DRDCE_TRIAL_HAT, DYDCETRIAL
  DOUBLE PRECISION, DIMENSION(9)                   :: TMP_9

  DOUBLE PRECISION, DIMENSION(6,6)                 :: DE_HYD, CEI_X, CEI_KR, DLOGCEHAT_DCEHAT, DCEHAT_DCE, DE_DEV, XFP1, XFP2, XFP
  DOUBLE PRECISION, DIMENSION(6,6)                 :: DLOGCBAR1_DCBAR1, DLOGCBAR2_DCBAR2, Conny, DEXPM, MW66, Part2, QNEW,FPISFPI
  DOUBLE PRECISION, DIMENSION(3,3)                 :: X, EXPM,SFPI, C_NEW_HAT, TMP33, KR, FE_NEW, S, SE_BAR
  DOUBLE PRECISION, DIMENSION(6,6)                 :: TMP66_1,TMP66_2, TMP66_3, PUSH, KR_DYAD_KR, D_BAR	, D_HYD	
  DOUBLE PRECISION, DIMENSION(6,1)                 :: TMP_61
  DOUBLE PRECISION, DIMENSION(6)                   :: TMP_6


   Y_OLD=Isv_Old
   Y_NEW=Isv_New  

   if (type_defgr.eq.'total') then
     F_NEW=getF(ef_New)
   elseif (type_defgr.eq.'relative') then
     stop "Not implemented, err 82256"   
   end if   

   KR     =getI()
   CALL EYE(DIAG_UNIT,6)
   DIAG_UNIT(4,4)=.5D0
   DIAG_UNIT(5,5)=.5D0
   DIAG_UNIT(6,6)=.5D0

   KR_DYAD_KR(1:3,1:3)=1D0

   DRDCE_TRIAL_HAT=0D0

   C=MATMUL(TRANSPOSE(F_NEW),F_NEW)         
   CALL INV3S(IC,C)

   J=DSQRT(det3(C))   
!   th=J ! Mixed formulation, since th differs from J
   C_NEW_HAT=J**(-2D0/3D0)*C
 
   FP_OLD=TRANSPOSE(   RESHAPE((/Y_OLD(1:3),   Y_OLD(4:6),   Y_OLD(7:9)/),(/3,3/)))
   FBAR1_OLD=TRANSPOSE(RESHAPE((/Y_OLD(10:12), Y_OLD(13:15), Y_OLD(16:18)/),(/3,3/)))
   FBAR2_OLD=TRANSPOSE(RESHAPE((/Y_OLD(19:21), Y_OLD(22:24), Y_OLD(25:27)/),(/3,3/)))
   FP_NEW=TRANSPOSE(   RESHAPE((/Y_NEW(1:3),   Y_NEW(4:6),   Y_NEW(7:9)/),(/3,3/)))
   FBAR1_NEW=TRANSPOSE(RESHAPE((/Y_NEW(10:12), Y_NEW(13:15), Y_NEW(16:18)/),(/3,3/)))
   FBAR2_NEW=TRANSPOSE(RESHAPE((/Y_NEW(19:21), Y_NEW(22:24), Y_NEW(25:27)/),(/3,3/)))
   DLAMB=Y_NEW(28)

   CALL INV3(IFP,FP_NEW)  
   FE_NEW=MATMUL(F_NEW, IFP)

   CE    =MATMUL(TRANSPOSE(FE_NEW), FE_NEW)
   CALL INV3S(ICE,CE)

   TMP_9=(/F_NEW(1,1), F_NEW(1,2), F_NEW(1,3), &
           F_NEW(2,1), F_NEW(2,2), F_NEW(2,3), &
           F_NEW(3,1), F_NEW(3,2), F_NEW(3,3) /)

   FE_NEW_HAT=J**(-1D0/3D0)*FE_NEW
   CE_NEW_HAT  =MATMUL(TRANSPOSE(FE_NEW_HAT), FE_NEW_HAT)
  
   CBAR1       =MATMUL(TRANSPOSE(FBAR1_NEW), FBAR1_NEW)
   CBAR2       =MATMUL(TRANSPOSE(FBAR2_NEW), FBAR2_NEW)

   CALL INV3(IFP_O,FP_OLD)  
   CE_TR_HAT  =MATMUL(MATMUL(TRANSPOSE(IFP_O), C_NEW_HAT),IFP_O)

   CBAR1_TRIAL   =MATMUL(TRANSPOSE(FBAR1_OLD),FBAR1_OLD)
   CBAR2_TRIAL   =MATMUL(TRANSPOSE(FBAR2_OLD),FBAR2_OLD)

   CE_TRM=J**(2D0/3D0)*RESHAPE((/CE_TR_HAT(1,1), CE_TR_HAT(2,2), CE_TR_HAT(3,3),CE_TR_HAT(1,2), CE_TR_HAT(1,3), CE_TR_HAT(2,3)/),(/6,1/))

   TMP66=KAPPA*(0.5D0*KR_DYAD_KR-DIAG_UNIT*DLOG(th))
   call pullback(D_HYD, TMP66,TMP_9)
   D_HYD=2D0*D_HYD

   Y=(/CE_NEW_HAT(1,1), CE_NEW_HAT(2,2), CE_NEW_HAT(3,3), &
       CE_NEW_HAT(1,2), CE_NEW_HAT(1,3), CE_NEW_HAT(2,3), &
       CBAR1(1,1),  CBAR1(2,2),  CBAR1(3,3), &
       CBAR1(1,2),  CBAR1(1,3),  CBAR1(2,3), &
       CBAR2(1,1),  CBAR2(2,2),  CBAR2(3,3), &
       CBAR2(1,2),  CBAR2(1,3),  CBAR2(2,3), DLAMB /)
  
    CE_IM  =RESHAPE((/ICE(1,1), ICE(2,2), ICE(3,3), ICE(1,2), ICE(1,3), ICE(2,3)/),(/6,1/))
    CE_M   =RESHAPE((/ CE(1,1),  CE(2,2),  CE(3,3),  CE(1,2),  CE(1,3),  CE(2,3)/),(/6,1/))
    C_IM   =RESHAPE((/ IC(1,1),  IC(2,2),  IC(3,3),  IC(1,2),  IC(1,3),  IC(2,3)/),(/6,1/))      
             
    CALL PADELOG(CE_NEW_HAT, LOGCE_HAT)

    X=MATMUL(ICE,LOGCE_HAT)        
    CALL BAR_DYAD(CEI_X,ICE,X)
    CEI_X(:,4:6)=0.5D0*CEI_X(:,4:6)        
            
    CALL BAR_DYAD(CEI_KR,ICE,KR)

    CALL DPADELOG(CE_NEW_HAT, DLOGCEHAT_DCEHAT)	
        	     
    TMP66=MATMUL(CEI_KR,DLOGCEHAT_DCEHAT) 

    DCEHAT_DCE=J**(-2D0/3D0)*(DIAG_UNIT-1D0/3D0*MATMUL(CE_M,TRANSPOSE(CE_IM)))
       
    DE_DEV=MU*(-CEI_X+MATMUL(TMP66,DCEHAT_DCE))
   
    CALL BAR_DYAD(PUSH_FW,IFP,IFP)

    IF (DLAMB>0D0) THEN
 !          CALL NUMERICALJACOBIAN (RES,JAC,Y,CE_TR_HAT,CBAR1_TRIAL,CBAR2_TRIAL,19,ITER)	      
            CALL ANALYTICALJACOBIAN(RES,JAC,Y,CE_TR_HAT,CBAR1_TRIAL,CBAR2_TRIAL,19,ITER) 
            CALL INVERSE(JAC,JAC_INV,19)    

            CALL PADELOG(CBAR1, LOGCBAR1)
            CALL PADELOG(CBAR2, LOGCBAR2)

            B1_DEV=KSI1*LOGCBAR1 
            B2_DEV=KSI2*LOGCBAR2
                 
            M_RED=MU*LOGCE_HAT-B1_DEV-B2_DEV
            LP=3D0/(2D0*SIGY0)*M_RED
                
            CALL PADEEXP(LP*DLAMB, AP)
            CALL INV3(APINV,AP)  

            CALL DPADELOG(CBAR1, DLOGCBAR1_DCBAR1)
            CALL DPADELOG(CBAR2, DLOGCBAR2_DCBAR2)
            CALL DPADEEXP(-LP*DLAMB,DEXPM)
            
            CALL BAR_DYAD(PUSH,APINV,APINV)
  
            DRCEDCE_TRIAL_HAT=-MATMUL(PUSH,DIAG_UNIT)                  
            DRDCE_TRIAL_HAT(1:6,:)=DRCEDCE_TRIAL_HAT	       
            DYDCETRIAL=-MATMUL(JAC_INV,DRDCE_TRIAL_HAT)  

            CALL BAR_DYAD(PUSH,TRANSPOSE(IFP_o),TRANSPOSE(IFP_o))
            DCE_HAT_TRIAL_DC=(J**(-2D0/3D0))*(MATMUL(PUSH,DIAG_UNIT)-1D0/3D0*MATMUL(CE_TRM,TRANSPOSE(C_IM)))

            DCE_HAT_TRIAL_DC(4:6,:)=2D0*DCE_HAT_TRIAL_DC(4:6,:)
      
            TMP66=DYDCETRIAL(1:6,1:6)
            DCE_DC=J**(2D0/3D0)*MATMUL(TMP66,DCE_HAT_TRIAL_DC)+1D0/3D0*MATMUL(CE_M,TRANSPOSE(C_IM)) 
	    DCE_DC(4:6,:)=2D0*DCE_DC(4:6,:)

            TMP_N=RESHAPE((/dYdCetrial(19,1), dYdCetrial(19,2), dYdCetrial(19,3), &
                            dYdCetrial(19,4), dYdCetrial(19,5), dYdCetrial(19,6)/),(/6,1/))			  

            TMP_61=reshape((/LP(1,1), LP(2,2), LP(3,3), LP(1,2), LP(3,1),  LP(3,2)/),(/6,1/))
                      
            TMP66_1 =MATMUL(DLOGCEHAT_DCEHAT, dYdCetrial(1:6,:))	                      
            TMP66_2 =MATMUL(DLOGCBAR1_DCBAR1,dYdCetrial(7:12,:))
            TMP66_3 =MATMUL(DLOGCBAR2_DCBAR2,dYdCetrial(13:18,:))
         
            Qnew=3D0*DLAMB/(2D0*SIGY0)*(KSI1*TMP66_2+KSI2*TMP66_3-MU*TMP66_1)-MATMUL(TMP_61,TRANSPOSE(TMP_N))
          
            SE_BAR=MU*MATMUL(ICE,LOGCE_HAT)        
            SFpI=MATMUL(IFP, SE_BAR) 
  
            CALL BAR_DYAD(XFP1,IFP_O,SFpi)
            CALL BAR_DYAD(XFP2,SFpi,IFP_O)
            XFP=(XFP1+XFP2)
	                           
            Part2=MATMUL(QNEW,dCe_hat_trial_dC)
            Part2=MATMUL(DEXPM,PART2)
            PART2=MATMUL(XFP,PART2)         
            D_BAR=2D0*(MATMUL(MATMUL(PUSH_FW,DE_DEV),DCE_DC)+Part2)            
        ELSE
            D_BAR=2D0*(MATMUL(MATMUL(PUSH_FW,DE_DEV),TRANSPOSE(PUSH_FW)))         
        ENDIF       
                   

  if (size(ef_New).eq.4) then
    Dout_bar(1,:)=(/D_BAR(1,1), D_BAR(2,2), D_BAR(1,4)/)
    Dout_bar(2,:)=(/D_BAR(2,1), D_BAR(2,2), D_BAR(2,4)/)
    Dout_bar(3,:)=(/D_BAR(4,1), D_BAR(4,2), D_BAR(4,4)/)
  else
    Dout_bar=D_BAR
  endif 

  Dout_hyd=kappa/th


  if (stype.eq.'ul') then
    call pushforward(Dout_bar,Dout_bar,ef_new)
  elseif (stype.eq.'tl') then
  else
  endif
 

end subroutine dJ2_kin2up
 

SUBROUTINE    NUMERICALJACOBIAN(RES,JAC,Y,MP,CE_trial,CBAR1_trial,CBAR2_trial,N_Y,N_MP,ITER)
       IMPLICIT NONE
       DOUBLE PRECISION,    DIMENSION(N_Y)              :: Y, RES, FTY ,D 
       DOUBLE PRECISION,    DIMENSION(N_Y)              :: FTY_TMPB,FTY_TMPF,FTY_TMPBB,FTY_TMPFF
       DOUBLE PRECISION,    DIMENSION(N_Y,N_Y)          :: JAC,YDELF,YDELB,YDELFF,YDELBB
       DOUBLE PRECISION,    DIMENSION(N_MP)             :: MP
       DOUBLE PRECISION,    DIMENSION(3,3)              :: CE_trial, CBAR1_trial,CBAR2_trial,CE,CBAR1, CBAR2     
       INTEGER                                          :: N_Y, N_MP,ITER, I, J      
       DOUBLE PRECISION                                 :: DLAMB, DEL, DEL1
   
       CE   =RESHAPE((/Y(1), Y(4),  Y(5),  Y(4),  Y(2), Y(6),  Y(5),  Y(6),  Y(3)/),(/3,3/))
       CBAR1=RESHAPE((/Y(7), Y(10), Y(11), Y(10), Y(8), Y(12), Y(11), Y(12), Y(9)/),(/3,3/))      
       CBAR2=RESHAPE((/Y(13),Y(16), Y(17), Y(16), Y(14),Y(18), Y(17), Y(18), Y(15)/),(/3,3/))                    
       DLAMB=Y(19)
 
       DEL   = 2D-9
                 
       CALL YIELD(FTY,Y,CE_trial,CBAR1_trial,CBAR2_trial,N_Y)
	     
       DO I=1,N_Y
          YDELF(:,I) =Y  
          YDELB(:,I) =Y 
          YDELFF(:,I)=Y  
	       YDELBB(:,I)=Y 
	  	        
          DO J=1,N_Y
             IF (I==J) THEN 	       
	             D(I)= DEL*max(4D-2,abs(Y(I)))                
                YDELF(I,J)=YDELF(I,J)+D(I)
                YDELB(I,J)=YDELB(I,J)-D(I)
                YDELFF(I,J)=YDELFF(I,J)+2D0*D(I)
                YDELBB(I,J)=YDELBB(I,J)-2D0*D(I)				
	          ENDIF
	       END DO
       END DO       
      
       JAC = 0D0
       DO I=1,N_Y
          CALL YIELD(FTY_TMPF,YDELF(:,I)  ,CE_trial,CBAR1_trial,CBAR2_trial,N_Y)
          CALL YIELD(FTY_TMPB,YDELB(:,I)  ,CE_trial,CBAR1_trial,CBAR2_trial,N_Y)	
          CALL YIELD(FTY_TMPFF,YDELFF(:,I),CE_trial,CBAR1_trial,CBAR2_trial,N_Y)
          CALL YIELD(FTY_TMPBB,YDELBB(:,I),CE_trial,CBAR1_trial,CBAR2_trial,N_Y)	
	  
!          JAC(:,I)=(FTY_TMPF-FTY_TMPB)/(2D0*D(I))	  
!          JAC(:,I)=(FTY_TMP-FTY)/D(I)
          JAC(:,I)=(FTY_TMPBB-8D0*FTY_TMPB+8D0*FTY_TMPF-FTY_TMPFF)/(12D0*D(I))	  
       END DO            
       CALL YIELD(RES,Y,CE_trial,CBAR1_trial,CBAR2_trial,N_Y)

RETURN
END  subroutine NUMERICALJACOBIAN


SUBROUTINE    YIELD(RES,Y,CE_trial,CBAR1_trial,CBAR2_trial,N_Y)
       IMPLICIT NONE
       INTEGER                                          :: N_Y, INFO
       DOUBLE PRECISION                                 :: tol, dlamb, sigeff,fe, fk, ge, gk    
       DOUBLE PRECISION,    DIMENSION(N_Y,N_Y)          :: JAC 
       DOUBLE PRECISION,    DIMENSION(N_Y)              :: Y, RES
       DOUBLE PRECISION,    DIMENSION(3,3)              :: CE, CBAR1, CBAR2, LOGCE, LOGCBAR1, LOGCBAR2
       DOUBLE PRECISION,    DIMENSION(3,3)              :: CRAP1,CRAP2,CRAP3,CRAP4,CRAP5, B1_DEV, B2_DEV, MANDEL_REDE_DEV       
       DOUBLE PRECISION,    DIMENSION(3,3)              :: APINV, AB1, AB2, NE, NK1, NK2, CE_trial, CBAR1_trial, CBAR2_trial
       DOUBLE PRECISION,    DIMENSION(3)                :: DWe, DWk1,DWk2,DLe, DLk1, DLk2
            
                   
       CE   =RESHAPE((/Y(1), Y(4),  Y(5),  Y(4),  Y(2), Y(6),  Y(5),  Y(6),  Y(3)/),(/3,3/))
       CBAR1=RESHAPE((/Y(7), Y(10), Y(11), Y(10), Y(8), Y(12), Y(11), Y(12), Y(9)/),(/3,3/))      
       CBAR2=RESHAPE((/Y(13),Y(16), Y(17), Y(16), Y(14),Y(18), Y(17), Y(18), Y(15)/),(/3,3/))                    
       DLAMB=Y(19)

       CALL PADELOG(CE, LOGCE)
       CALL PADELOG(CBAR1, LOGCBAR1)
       CALL PADELOG(CBAR2, LOGCBAR2)

       B1_DEV=KSI1*LOGCBAR1
       B2_DEV=KSI2*LOGCBAR2
       MANDEL_REDE_DEV=MU*LOGCE-B1_DEV-B2_DEV
       CRAP1=MATMUL(MANDEL_REDE_DEV,MANDEL_REDE_DEV)
       SIGEFF=DSQRT(3D0/2D0*(CRAP1(1,1)+CRAP1(2,2)+CRAP1(3,3)))
              
       NE=3D0/(2D0*SIGY0)*MANDEL_REDE_DEV
       NK1=NE-GAMMA1*B1_DEV
       NK2=NE-GAMMA2*B2_DEV

       CALL PADEEXP(-NE*dlamb, APINV)
       CALL PADEEXP(NK1*dlamb,AB1 ) 
       CALL PADEEXP(NK2*dlamb,AB2 ) 

       CRAP1=CE-MATMUL(TRANSPOSE(APINV),MATMUL(CE_trial,APINV))
       RES(1:6)=(/CRAP1(1,1), CRAP1(2,2), CRAP1(3,3), CRAP1(1,2), CRAP1(1,3), CRAP1(2,3)/)
       CRAP1=CBAR1-MATMUL(TRANSPOSE(AB1),MATMUL(CBAR1_trial,AB1))
       RES(7:12)=(/CRAP1(1,1), CRAP1(2,2), CRAP1(3,3), CRAP1(1,2), CRAP1(1,3), CRAP1(2,3)/)
       CRAP1=CBAR2-MATMUL(TRANSPOSE(AB2),MATMUL(CBAR2_trial,AB2))
       RES(13:19)=(/CRAP1(1,1), CRAP1(2,2), CRAP1(3,3), CRAP1(1,2), CRAP1(1,3), CRAP1(2,3),sigeff-sigy0/)
     
RETURN
END subroutine Yield







SUBROUTINE    ANALYTICALJACOBIAN(RES,JAC,Y,CE_trial,CBAR1_trial,CBAR2_trial,N_Y,ITER)
       IMPLICIT NONE
       INTEGER                                       :: N_Y,   ITER
       DOUBLE PRECISION                              :: DLAMB, SIGEFF
       DOUBLE PRECISION,    DIMENSION(N_Y,N_Y)       :: JAC 
       DOUBLE PRECISION,    DIMENSION(N_Y)           :: Y, RES
       DOUBLE PRECISION,    DIMENSION(3,3)           :: CE_trial, CBAR1_trial, CBAR2_trial, CE, CBAR1, CBAR2
       DOUBLE PRECISION,    DIMENSION(6,6)           :: DEXPE, DEXPK1, DEXPK2, DLOGCBAR1_DCBAR1, DLOGCBAR2_DCBAR2, DLOGCEHAT_DCEHAT

       DOUBLE PRECISION,    DIMENSION(3,3)           :: LOGCE, LOGCBAR1, LOGCBAR2, TMP33, NE, NK1, NK2, APINV, AB1, AB2
       DOUBLE PRECISION,    DIMENSION(3,3)           :: KR, APICE, AB1CBAR1, AB2CBAR2, B1_DEV, B2_DEV, MANDEL_REDE_DEV
       DOUBLE PRECISION,    DIMENSION(6,6)           :: TMP66_1,TMP66_2, APICE_KR, AB1CBAR1_KR,  AB2CBAR2_KR
       DOUBLE PRECISION,    DIMENSION(6)             :: TMP6


       KR     =getI()      
      
       CE   =RESHAPE((/Y(1), Y(4),  Y(5),  Y(4),  Y(2), Y(6),  Y(5),  Y(6),  Y(3)/),(/3,3/))
       CBAR1=RESHAPE((/Y(7), Y(10), Y(11), Y(10), Y(8), Y(12), Y(11), Y(12), Y(9)/),(/3,3/))      
       CBAR2=RESHAPE((/Y(13),Y(16), Y(17), Y(16), Y(14),Y(18), Y(17), Y(18), Y(15)/),(/3,3/))                    
       DLAMB=Y(19)
     
       if (dlamb<0) then
          write(*,*) 'DLAMB<0 when solving local system'
       endif

       CALL PADELOG(CE,    LOGCE)
       CALL PADELOG(CBAR1, LOGCBAR1)
       CALL PADELOG(CBAR2, LOGCBAR2)
                 
       B1_DEV=KSI1*LOGCBAR1
       B2_DEV=KSI2*LOGCBAR2
              
       MANDEL_REDE_DEV=MU*LOGCE-B1_DEV-B2_DEV
       TMP33=MATMUL(MANDEL_REDE_DEV,MANDEL_REDE_DEV)

       SIGEFF=DSQRT(3D0/2D0*(TMP33(1,1)+TMP33(2,2)+TMP33(3,3)))
              
       NE=3D0/(2D0*SIGY0)*MANDEL_REDE_DEV
       NK1=NE-GAMMA1*B1_DEV
       NK2=NE-GAMMA2*B2_DEV       

       CALL PADEEXP(-NE*DLAMB, APINV)
       CALL PADEEXP(NK1*DLAMB,AB1)
       CALL PADEEXP(NK2*DLAMB,AB2)     

       CALL DPADEEXP(-NE*DLAMB, DEXPE)
       CALL DPADEEXP( NK1*DLAMB,DEXPK1)
       CALL DPADEEXP( NK2*DLAMB,DEXPK2)

       CALL DPADELOG(CBAR1, DLOGCBAR1_DCBAR1)
       CALL DPADELOG(CBAR2, DLOGCBAR2_DCBAR2)
       CALL DPADELOG(CE   , DLOGCEHAT_DCEHAT)

       APICE=MATMUL(APINV,CE_TRIAL)
       CALL BAR_DYAD(TMP66_1,APICE,KR)
       CALL BAR_DYAD(TMP66_2,KR,APICE)
       APICE_KR=TMP66_1+TMP66_2

       AB1CBAR1=MATMUL(AB1,CBAR1_TRIAL)
       CALL BAR_DYAD(TMP66_1,AB1CBAR1,KR)
       CALL BAR_DYAD(TMP66_2,KR,AB1CBAR1)
       AB1CBAR1_KR=TMP66_1+TMP66_2

       AB2CBAR2=MATMUL(AB2,CBAR2_TRIAL)
       CALL BAR_DYAD(TMP66_1,AB2CBAR2,KR)
       CALL BAR_DYAD(TMP66_2,KR,AB2CBAR2)
       AB2CBAR2_KR=TMP66_1+TMP66_2

       CALL EYE(JAC,19)
       JAC(19,19)=0D0
        

       JAC(1:6,1:6)  =JAC(1:6,1:6)+3D0*MU*DLAMB/(2D0*SIGY0)*MATMUL(APICE_KR,MATMUL(DEXPE, DLOGCEHAT_DCEHAT))
       JAC(1:6,7:12) =            -3D0*KSI1*DLAMB/(2D0*SIGY0)*MATMUL(APICE_KR,MATMUL(DEXPE, DLOGCBAR1_DCBAR1))
       JAC(1:6,13:18)=            -3D0*KSI2*DLAMB/(2D0*SIGY0)*MATMUL(APICE_KR,MATMUL(DEXPE, DLOGCBAR2_DCBAR2))
     
       JAC(7:12,1:6)  =           -MU*DLAMB*3D0/(2D0*SIGY0)*MATMUL(AB1CBAR1_KR,MATMUL(DEXPK1, DLOGCEHAT_DCEHAT))
       JAC(7:12,7:12) = JAC(7:12,7:12)+ KSI1*DLAMB*(3D0/(2D0*SIGY0)+GAMMA1)*MATMUL(AB1CBAR1_KR,MATMUL(DEXPK1, DLOGCBAR1_DCBAR1))
       JAC(7:12,13:18)=            KSI2*DLAMB*(3D0/(2D0*SIGY0))*MATMUL(AB1CBAR1_KR,MATMUL(DEXPK1, DLOGCBAR2_DCBAR2))

       JAC(13:18,1:6)  =           -MU*DLAMB*3D0/(2D0*SIGY0)*MATMUL(AB2CBAR2_KR,MATMUL(DEXPK2, DLOGCEHAT_DCEHAT))
       JAC(13:18,7:12) =            KSI1*DLAMB*(3D0/(2D0*SIGY0))*MATMUL(AB2CBAR2_KR,MATMUL(DEXPK2, DLOGCBAR1_DCBAR1))
       JAC(13:18,13:18)= JAC(13:18,13:18)+ KSI2*DLAMB*(3D0/(2D0*SIGY0)+GAMMA2)*MATMUL(AB2CBAR2_KR,MATMUL(DEXPK2, DLOGCBAR2_DCBAR2))
       

       TMP6=(/Ne(1,1), Ne(2,2), Ne(3,3), Ne(1,2), Ne(1,3), Ne(2,3)/)
       DLOGCEHAT_DCEHAT(4:6,:)=2D0*DLOGCEHAT_DCEHAT(4:6,:)  
       JAC(1:6,19)=MATMUL(MATMUL(APICE_KR,DEXPE),TMP6)    
       JAC(19,1:6)=MU*SIGY0/SIGEFF*MATMUL(DLOGCEHAT_DCEHAT,TMP6)

       DLOGCBAR1_DCBAR1(4:6,:)=2D0*DLOGCBAR1_DCBAR1(4:6,:)  
       JAC(19,7:12)=-KSI1*SIGY0/SIGEFF*MATMUL(DLOGCBAR1_DCBAR1,TMP6)

       DLOGCBAR2_DCBAR2(4:6,:)=2D0*DLOGCBAR2_DCBAR2(4:6,:)  
       JAC(19,13:18)=-KSI2*SIGY0/SIGEFF*MATMUL(DLOGCBAR2_DCBAR2,TMP6)


       TMP6=(/Nk1(1,1), Nk1(2,2), Nk1(3,3), Nk1(1,2), Nk1(1,3), Nk1(2,3)/)
       JAC(7:12,19)=-MATMUL(MATMUL(AB1CBAR1_KR,DEXPK1),TMP6)

       TMP6=(/Nk2(1,1), Nk2(2,2), Nk2(3,3), Nk2(1,2), Nk2(1,3), Nk2(2,3)/)
       JAC(13:18,19)=-MATMUL(MATMUL(AB2CBAR2_KR,DEXPK2),TMP6)


            
       TMP33=CE-MATMUL(TRANSPOSE(APINV),MATMUL(CE_trial,APINV))
       RES(1:6)=(/TMP33(1,1), TMP33(2,2), TMP33(3,3), TMP33(1,2), TMP33(1,3), TMP33(2,3)/)
       TMP33=CBAR1-MATMUL(TRANSPOSE(AB1),MATMUL(CBAR1_trial,AB1))
       RES(7:12)=(/TMP33(1,1), TMP33(2,2), TMP33(3,3), TMP33(1,2), TMP33(1,3), TMP33(2,3)/)
       TMP33=CBAR2-MATMUL(TRANSPOSE(AB2),MATMUL(CBAR2_trial,AB2))
       RES(13:19)=(/TMP33(1,1), TMP33(2,2), TMP33(3,3), TMP33(1,2), TMP33(1,3), TMP33(2,3),SIGEFF-SIGY0/)
     
RETURN
END SUBROUTINE ANALYTICALJACOBIAN












end module mater_J2kin
