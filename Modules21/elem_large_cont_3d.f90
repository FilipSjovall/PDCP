module elem_large_cont_3d

! Last rev.:
!  M. Ristinmaa 2011-05-04
!   - initial version 8-node brick element 
!  M. Ristinmaa 2011-05-06
!   - changed order of integration points to be consistent with Abaqus
!  M. Ristinmaa 2011-05-27
!   - changed location of nodes in c3dtl8
!  M. Ristinmaa 2011-05-27
!   - introducing 8 noded Q1P0 
!  M. Ristinmaa 2011-09-11
!   - debugged 8 noded Q1P0 
!  M. Ristinmaa 2011-10-10
!   - start to defined parameters
!  M. Ristinmaa 2011-10-12
!   - rewrote TL element to use parameters and data definitions
!     noticable spped increase
!   - changed such that matrix_util is used, for inv3 order of arguments
!     had to be changed
!  M. Ristinmaa 2011-10-13
!   - rewrote UL element to use paramaters and data definitions
!------------------------------------------------------------------------------
!  M. Wallin 2013-09-15
!   - elmvol added
!  M. Wallin 2013-09-18
!   - c3dtl8_fsens added
!------------------------------------------------------------------------------	
!  M. Wallin 2013-09-26
!   - c3dtl8_etilde added
!------------------------------------------------------------------------------	
!	A. Dalklint 2018-04-16
!	 - c3dtl8_mtilde added
! 	 - dc3dtl8_mtilde added


!   Module elem_large_cont_3d contains element subroutines 
!   large deformations.	

! 8 node brick element, trilinear displacement interpolation
! 2x2x2 gauss integration points (in total 8), the weight for every
! gauss point is 1 for this scheme, and it therefore omitted in the code

use matrix_util, only: inv3, det3

implicit none

! This is an attempt to predefine some vectors, to obtain speed
! Implemented for the TL and UL elements
double precision  xsi(8), eta(8), zet(8), G1, Ifm(9)

integer           index(3,8)
parameter        (G1=0.577350269189626D0)
parameter        (xsi=(/-1D0,  1D0, 1D0, -1D0, -1D0,  1D0, 1D0, -1D0/)*G1 )
parameter        (eta=(/-1D0, -1D0, 1D0,  1D0, -1D0, -1D0, 1D0,  1D0/)*G1 )
parameter        (zet=(/ 1D0,  1D0, 1D0,  1D0, -1D0, -1D0,-1D0, -1D0/)*G1 )
parameter        (Ifm=(/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/))
parameter        (index=[(/1,2,3/),(/4,5,6/),(/7,8,9/),(/10,11,12/), &
                         (/13,14,15/),(/16,17,18/),(/19,20,21/),(/22,23,24/)])
private xsi, eta, zet, G1, Ifm, index

double precision  DNR(24,8) ! Maybe one should use DNR(8,24) instead for speed
! derivate of shape functions with respect to xsi
data             (DNR(1:24:3,1)= (1D0+ETA)*(1D0-ZET)/8D0)
data             (DNR(1:24:3,2)=-(1D0+ETA)*(1D0-ZET)/8D0)
data             (DNR(1:24:3,3)=-(1D0-ETA)*(1D0-ZET)/8D0)
data             (DNR(1:24:3,4)= (1D0-ETA)*(1D0-ZET)/8D0)
data             (DNR(1:24:3,5)= (1D0+ETA)*(1D0+ZET)/8D0)
data             (DNR(1:24:3,6)=-(1D0+ETA)*(1D0+ZET)/8D0)
data             (DNR(1:24:3,7)=-(1D0-ETA)*(1D0+ZET)/8D0)
data             (DNR(1:24:3,8)= (1D0-ETA)*(1D0+ZET)/8D0)
! derivate of shape functions with respect to eta
data             (DNR(2:24:3,1)= (1D0+XSI)*(1D0-ZET)/8D0)
data             (DNR(2:24:3,2)= (1D0-XSI)*(1D0-ZET)/8D0)
data             (DNR(2:24:3,3)=-(1D0-XSI)*(1D0-ZET)/8D0)
data             (DNR(2:24:3,4)=-(1D0+XSI)*(1D0-ZET)/8D0)
data             (DNR(2:24:3,5)= (1D0+XSI)*(1D0+ZET)/8D0)
data             (DNR(2:24:3,6)= (1D0-XSI)*(1D0+ZET)/8D0)
data             (DNR(2:24:3,7)=-(1D0-XSI)*(1D0+ZET)/8D0)
data             (DNR(2:24:3,8)=-(1D0+XSI)*(1D0+ZET)/8D0)
! derivate of shape functions with respect to zet
data             (DNR(3:24:3,1)=-(1D0+XSI)*(1D0+ETA)/8D0)
data             (DNR(3:24:3,2)=-(1D0-XSI)*(1D0+ETA)/8D0)
data             (DNR(3:24:3,3)=-(1D0-XSI)*(1D0-ETA)/8D0)
data             (DNR(3:24:3,4)=-(1D0+XSI)*(1D0-ETA)/8D0)
data             (DNR(3:24:3,5)= (1D0+XSI)*(1D0+ETA)/8D0)
data             (DNR(3:24:3,6)= (1D0-XSI)*(1D0+ETA)/8D0)
data             (DNR(3:24:3,7)= (1D0-XSI)*(1D0-ETA)/8D0)
data             (DNR(3:24:3,8)= (1D0+XSI)*(1D0-ETA)/8D0) 
private DNR

double precision  NR(8,8) ! 8 Gauss points and 8 nodes
data             (NR(:,1)= (1D0+XSI)*(1D0+ETA)*(1D0-ZET)/8D0)
data             (NR(:,2)= (1D0-XSI)*(1D0+ETA)*(1D0-ZET)/8D0)
data             (NR(:,3)= (1D0-XSI)*(1D0-ETA)*(1D0-ZET)/8D0)
data             (NR(:,4)= (1D0+XSI)*(1D0-ETA)*(1D0-ZET)/8D0)
data             (NR(:,5)= (1D0+XSI)*(1D0+ETA)*(1D0+ZET)/8D0)
data             (NR(:,6)= (1D0-XSI)*(1D0+ETA)*(1D0+ZET)/8D0)
data             (NR(:,7)= (1D0-XSI)*(1D0-ETA)*(1D0+ZET)/8D0)
data             (NR(:,8)= (1D0+XSI)*(1D0-ETA)*(1D0+ZET)/8D0)
private NR

interface c3dul8p0_f
  module procedure c3dul8_f
end interface

interface c3dul8_d
  module procedure c3dtl8_d
end interface

interface c3dul8_emu
  module procedure c3dtl8_emu
  module procedure c3dtl8_emu2
end interface

interface c3dul8_emu2
  module procedure c3dtl8_emu2
end interface

interface c3dtl8_m
  module procedure c3dtl8_m
  module procedure dc3dtl8_m
end interface


!------------------------------------------------------------------------------
contains


! 8-node brick Q1P0 element based on updated Lagrangian formulation	
subroutine c3dul8p0_e(ke,coord,Diso,ddWdJJ,ed,es)
  implicit none
  double precision                :: ke(:,:), coord(:,:), Diso(:,:,:)
  double precision                :: ed(:), es(:,:), ddWdJJ(:)

  integer, parameter              :: NGP=8
  integer                         :: gp_nr

  double precision                :: JT(24,3), JT0(24,3), JTinv(3,3)
  double precision                :: DetJ, DetJ0
  double precision                :: DNX(3,8)
  double precision                :: B(6,24)
  double precision                :: Dgp(6,6)
  DOUBLE PRECISION                :: STRESS(3,3)
  DOUBLE PRECISION                :: HE(9,24)
  DOUBLE PRECISION                :: R(9,9)
  double precision                :: ucoord(3,8), Volume0, trS, DW
  double precision                :: PB(1,24)

  double precision, parameter     :: P(1,6)=(/1d0,1d0,1d0,0d0,0d0,0d0/)
  double precision, parameter     :: Q(6,6)=(/-1d0, 1d0, 1d0, 0d0, 0d0, 0d0, &
                                               1d0,-1d0, 1d0, 0d0, 0d0, 0d0, &
                                               1d0, 1d0,-1d0, 0d0, 0d0, 0d0, &
                                               0d0, 0d0, 0d0,-1d0, 0d0, 0d0, &
                                               0d0, 0d0, 0d0, 0d0,-1d0, 0d0, &
                                               0d0, 0d0, 0d0, 0d0, 0d0,-1d0/)

  ucoord(1,:)=coord(1,:)+ed([1,4,7,10,13,16,19,22])
  ucoord(2,:)=coord(2,:)+ed([2,5,8,11,14,17,20,23])
  ucoord(3,:)=coord(3,:)+ed([3,6,9,12,15,18,21,24])

  JT=MATMUL(DNR,transpose(ucoord))
  JT0=MATMUL(DNR,transpose(coord))

  KE=0D0
  PB=0d0
  Volume0=0d0
  DW=0d0
  DO GP_NR=1,NGP

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 
    DETJ0=det3(JT0(index(:,gp_nr),:)) 

    B=0D0
    B(1,1:24:3)=dNx(1,:)
    B(2,2:24:3)=dNx(2,:)
    B(3,3:24:3)=dNx(3,:)
    B(4,1:24:3)=dNx(2,:)
    B(4,2:24:3)=dNx(1,:)
    B(5,1:24:3)=dNx(3,:)
    B(5,3:24:3)=dNx(1,:)
    B(6,2:24:3)=dNx(3,:)
    B(6,3:24:3)=dNx(2,:)
  
    He=0d0
    He(1:3,1:24:3)=dNx(1:3,:)
    He(4:6,2:24:3)=dNx(1:3,:)
    He(7:9,3:24:3)=dNx(1:3,:)

    stress(1,:)=(/es(1,gp_nr), es(4,gp_nr), es(5,gp_nr)/)
    stress(2,:)=(/es(4,gp_nr), es(2,gp_nr), es(6,gp_nr)/)
    stress(3,:)=(/es(5,gp_nr), es(6,gp_nr), es(3,gp_nr)/)

    R=0D0
    R(1:3,1:3)=STRESS(:,:)
    R(4:6,4:6)=STRESS(:,:)
    R(7:9,7:9)=STRESS(:,:)

    trS=(es(1,gp_nr)+es(2,gp_nr)+es(3,gp_nr))/3d0
    Dgp=Diso(:,:,gp_nr)+trS*Q
  
    KE=KE+(MATMUL(TRANSPOSE(B),MATMUL(Dgp,B)) &
           +MATMUL(MATMUL(TRANSPOSE(He),R),He))*DETJ

!   P^T*B
    PB=PB+matmul(P,B)*DetJ
    Volume0=Volume0+DetJ0
!    DW=ddWdJJ(GP_NR)**DetJ0*WP(GP_NR)

  END DO

! Missing part
  KE=KE+matmul(transpose(PB),PB)*ddWdJJ(1)/Volume0
!  KE=KE+PB*PB*DW/Volume0/Volume0

  RETURN
END subroutine c3dul8p0_e


! calculates the total deformation gradient
subroutine c3dul8p0_d(dg,th,coord,ed)
  implicit none

  double precision                :: dg(:,:), th, coord(:,:), ed(:)

  integer, parameter              :: NGP=8
  integer                         :: gp_nr

  double precision                :: JT(24,3), JT0(24,3), JTinv(3,3)
  double precision                :: DNX(3,8), DetJ, DetJ0
  double precision                :: Volume, Volume0
  double precision                :: ucoord(3,8), He(9,24)

! current nodal values
  ucoord(1,:)=coord(1,:)+ed([1,4,7,10,13,16,19,22])
  ucoord(2,:)=coord(2,:)+ed([2,5,8,11,14,17,20,23])
  ucoord(3,:)=coord(3,:)+ed([3,6,9,12,15,18,21,24])

  Volume=0d0
  Volume0=0d0

  JT=MATMUL(DNR,transpose(uCOORD))
  JT0=MATMUL(DNR,transpose(coord))
   
  DO GP_NR=1,NGP

    CALL INV3(JTinv,JT0(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 
    DETJ0=det3(JT0(index(:,gp_nr),:)) 

    He=0D0
    He(1:3,1:24:3)=dNx(1:3,:)
    He(4:6,2:24:3)=dNx(1:3,:)
    He(7:9,3:24:3)=dNx(1:3,:)

    dg(:,gp_nr)=Ifm+MATMUL(HE,ED)
 
    Volume=Volume+DetJ
    Volume0=Volume0+DetJ0

  END DO
  
  th=Volume/Volume0

  RETURN
end subroutine c3dul8p0_d


! 8-node brick element  based on updated Lagrangian formulation	
subroutine c3dul8_e(ke,coord,D,ed,es)
  implicit none
  double precision                :: ke(:,:), coord(:,:), D(:,:,:)
  double precision                :: ed(:), es(:,:)

  integer, parameter              :: NGP=8
  integer                         :: gp_nr

  double precision                :: JT(24,3), DetJ, JTinv(3,3)
  double precision                :: DNX(3,8)
  double precision                :: BL0(6,24)
  double precision                :: Dgp(6,6)
  DOUBLE PRECISION                :: STRESS(3,3)
  DOUBLE PRECISION                :: HE(9,24)
  DOUBLE PRECISION                :: R(9,9)
  double precision                :: ucoord(3,8)

  ucoord(1,:)=coord(1,:)+ed([1,4,7,10,13,16,19,22])
  ucoord(2,:)=coord(2,:)+ed([2,5,8,11,14,17,20,23])
  ucoord(3,:)=coord(3,:)+ed([3,6,9,12,15,18,21,24])

  JT=MATMUL(DNR,transpose(ucoord))

  KE=0D0
  DO GP_NR=1,NGP

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:24:3)=dNx(1,:)
    BL0(2,2:24:3)=dNx(2,:)
    BL0(3,3:24:3)=dNx(3,:)
    BL0(4,1:24:3)=dNx(2,:)
    BL0(4,2:24:3)=dNx(1,:)
    BL0(5,1:24:3)=dNx(3,:)
    BL0(5,3:24:3)=dNx(1,:)
    BL0(6,2:24:3)=dNx(3,:)
    BL0(6,3:24:3)=dNx(2,:)
  
    He=0d0
    He(1:3,1:24:3)=dNx(1:3,:)
    He(4:6,2:24:3)=dNx(1:3,:)
    He(7:9,3:24:3)=dNx(1:3,:)

    stress(1,:)=(/es(1,gp_nr), es(4,gp_nr), es(5,gp_nr)/)
    stress(2,:)=(/es(4,gp_nr), es(2,gp_nr), es(6,gp_nr)/)
    stress(3,:)=(/es(5,gp_nr), es(6,gp_nr), es(3,gp_nr)/)
    R=0D0
    R(1:3,1:3)=STRESS(:,:)
    R(4:6,4:6)=STRESS(:,:)
    R(7:9,7:9)=STRESS(:,:)

    Dgp=D(:,:,gp_nr)
  
    KE=KE+(MATMUL(TRANSPOSE(BL0),MATMUL(Dgp,BL0)) &
           +MATMUL(MATMUL(TRANSPOSE(He),R),He))*DETJ

  END DO
  RETURN
END subroutine c3dul8_e


subroutine c3dul8_f(ef,coord,ed,es)
  implicit none

  double precision                :: eF(:), coord(:,:),  ed(:), es(:,:)

  integer, parameter              :: NGP=8
  integer                         :: gp_nr

  double precision                :: JT(24,3), JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8), DETJ
  double precision                :: BL0(6,24)
  double precision                :: ucoord(3,8)

  ucoord(1,:)=coord(1,:)+ed([1,4,7,10,13,16,19,22])
  ucoord(2,:)=coord(2,:)+ed([2,5,8,11,14,17,20,23])
  ucoord(3,:)=coord(3,:)+ed([3,6,9,12,15,18,21,24])

  JT=MATMUL(DNR,transpose(ucoord))

  ef=0D0
  DO GP_NR=1,NGP
    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:24:3) =dNx(1,:)
    BL0(2,2:24:3) =dNx(2,:)
    BL0(3,3:24:3) =dNx(3,:)
    BL0(4,1:24:3) =dNx(2,:)
    BL0(4,2:24:3) =dNx(1,:)
    BL0(5,1:24:3) =dNx(3,:)
    BL0(5,3:24:3) =dNx(1,:)
    BL0(6,2:24:3) =dNx(3,:)
    BL0(6,3:24:3) =dNx(2,:)

    EF=EF+MATMUL(TRANSPOSE(BL0),es(:,gp_nr))*detJ

  END DO
  RETURN
END subroutine c3dul8_f


subroutine c3dul8_eSens(fke,coord,D,ed,larg,rarg,es)
  implicit none
  double precision                :: fke(:), kegp, coord(:,:), D(:,:,:)
  double precision                :: ed(:), es(:,:), larg(:), rarg(:)

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8)
  double precision                :: BL0(6,24), BL(6,24)
  double precision                :: Dgp(6,6)
  DOUBLE PRECISION                :: STRESS(3,3)
  DOUBLE PRECISION                :: HE(9,24)
  DOUBLE PRECISION                :: A(6,9)
  DOUBLE PRECISION                :: dg(9)
  DOUBLE PRECISION                :: R(9,9), largTMP(24), rargTMP(24), V2(8)
  double precision                :: ucoord(3,8)

  ucoord(1,:)=coord(1,:)+ed([1,4,7,10,13,16,19,22])
  ucoord(2,:)=coord(2,:)+ed([2,5,8,11,14,17,20,23])
  ucoord(3,:)=coord(3,:)+ed([3,6,9,12,15,18,21,24])

  JT=MATMUL(DNR,transpose(ucoord))

  kegp = 0d0
  fke = 0d0
  DO GP_NR=1,NGP

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:24:3)=dNx(1,:)
    BL0(2,2:24:3)=dNx(2,:)
    BL0(3,3:24:3)=dNx(3,:)
    BL0(4,1:24:3)=dNx(2,:)
    BL0(4,2:24:3)=dNx(1,:)
    BL0(5,1:24:3)=dNx(3,:)
    BL0(5,3:24:3)=dNx(1,:)
    BL0(6,2:24:3)=dNx(3,:)
    BL0(6,3:24:3)=dNx(2,:)
  
    He=0d0
    He(1:3,1:24:3)=dNx(1:3,:)
    He(4:6,2:24:3)=dNx(1:3,:)
    He(7:9,3:24:3)=dNx(1:3,:)

    stress(1,:)=(/es(1,gp_nr), es(4,gp_nr), es(5,gp_nr)/)
    stress(2,:)=(/es(4,gp_nr), es(2,gp_nr), es(6,gp_nr)/)
    stress(3,:)=(/es(5,gp_nr), es(6,gp_nr), es(3,gp_nr)/)
    R=0D0
    R(1:3,1:3)=STRESS(:,:)
    R(4:6,4:6)=STRESS(:,:)
    R(7:9,7:9)=STRESS(:,:)

    V2(:)=NR(GP_NR,:)
    Dgp=D(:,:,gp_nr)
    largTMP = larg
    rargTMP = rarg
  
    KEGP=dot_product(largTMP,(matmul((MATMUL(TRANSPOSE(BL0),MATMUL(Dgp,BL0)) &
           +MATMUL(MATMUL(TRANSPOSE(He),R),He)),rargTMP)))*DETJ

    fke=fke+V2*KEGP

  END DO
  RETURN
END subroutine c3dul8_eSens


subroutine c3dul8_fSens(ffe,coord,ed,edvarphi,es)
  implicit none

  double precision                :: ffe(:), fegp, coord(:,:),  ed(:), es(:,:), edvarphi(:)

  integer, parameter              :: NGP=8
  integer                         :: gp_nr

  double precision                :: JT(24,3), JTinv(3,3)
  double precision                :: TMP(3,8), dNx(3,8), DETJ
  double precision                :: BL0(6,24), He(9,24)
  double precision                :: A(6,9), dg(9), varphiTMP(24), V2(8)
  double precision                :: ucoord(3,8)

  ucoord(1,:)=coord(1,:)+ed([1,4,7,10,13,16,19,22])
  ucoord(2,:)=coord(2,:)+ed([2,5,8,11,14,17,20,23])
  ucoord(3,:)=coord(3,:)+ed([3,6,9,12,15,18,21,24])

  JT=MATMUL(DNR,transpose(ucoord))
  fegp = 0d0
  ffe = 0d0
  DO GP_NR=1,NGP

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    dNx=MATMUL(JTINV,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:24:3) =dNx(1,:)
    BL0(2,2:24:3) =dNx(2,:)
    BL0(3,3:24:3) =dNx(3,:)
    BL0(4,1:24:3) =dNx(2,:)
    BL0(4,2:24:3) =dNx(1,:)
    BL0(5,1:24:3) =dNx(3,:)
    BL0(5,3:24:3) =dNx(1,:)
    BL0(6,2:24:3) =dNx(3,:)
    BL0(6,3:24:3) =dNx(2,:)

    V2(:)=NR(GP_NR,:)
    varphiTMP = edvarphi
 
    fegp=dot_product(varphiTMP,MATMUL(TRANSPOSE(BL0),es(:,gp_nr)))*detJ

    ffe=ffe+V2*fegp

  END DO
  RETURN
END subroutine c3dul8_fSens














! 8-node brick element  based on total Lagrangian formulation	
subroutine c3dtl8_e(ke,coord,D,ed,es)
  implicit none
  double precision                :: ke(:,:), coord(:,:), D(:,:,:)
  double precision                :: ed(:), es(:,:)

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8)
  double precision                :: BL0(6,24), BL(6,24)
  double precision                :: Dgp(6,6)
  DOUBLE PRECISION                :: STRESS(3,3)
  DOUBLE PRECISION                :: HE(9,24)
  DOUBLE PRECISION                :: A(6,9)
  DOUBLE PRECISION                :: dg(9)
  DOUBLE PRECISION                :: R(9,9)


  JT=MATMUL(DNR,transpose(COORD))

  KE=0D0
  DO GP_NR=1,NGP

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:24:3)=dNx(1,:)
    BL0(2,2:24:3)=dNx(2,:)
    BL0(3,3:24:3)=dNx(3,:)
    BL0(4,1:24:3)=dNx(2,:)
    BL0(4,2:24:3)=dNx(1,:)
    BL0(5,1:24:3)=dNx(3,:)
    BL0(5,3:24:3)=dNx(1,:)
    BL0(6,2:24:3)=dNx(3,:)
    BL0(6,3:24:3)=dNx(2,:)
  
    He=0d0
    He(1:3,1:24:3)=dNx(1:3,:)
    He(4:6,2:24:3)=dNx(1:3,:)
    He(7:9,3:24:3)=dNx(1:3,:)

! displacement gradient
    dg=matmul(He,ed)
    A=0d0
    A([1,4,5],1)=dg(1:3)
    A([2,4,6],2)=(/dg(2),dg(1),dg(3)/)
    A([3,5,6],3)=(/dg(3),dg(1),dg(2)/)
    A([1,4,5],4)=dg(4:6)
    A([2,4,6],5)=(/dg(5),dg(4),dg(6)/)
    A([3,5,6],6)=(/dg(6),dg(4),dg(5)/)
    A([1,4,5],7)=dg(7:9)
    A([2,4,6],8)=(/dg(8),dg(7),dg(9)/)
    A([3,5,6],9)=(/dg(9),dg(7),dg(8)/)

    BL=BL0+MATMUL(A,He)

    stress(1,:)=(/es(1,gp_nr), es(4,gp_nr), es(5,gp_nr)/)
    stress(2,:)=(/es(4,gp_nr), es(2,gp_nr), es(6,gp_nr)/)
    stress(3,:)=(/es(5,gp_nr), es(6,gp_nr), es(3,gp_nr)/)
    R=0D0
    R(1:3,1:3)=STRESS(:,:)
    R(4:6,4:6)=STRESS(:,:)
    R(7:9,7:9)=STRESS(:,:)

    Dgp=D(:,:,gp_nr)
  
    KE=KE+(MATMUL(TRANSPOSE(BL),MATMUL(Dgp,BL)) &
           +MATMUL(MATMUL(TRANSPOSE(He),R),He))*DETJ

  END DO
  RETURN
END subroutine c3dtl8_e


! MASS MATRIX rho-tilde
subroutine c3dtl8_m(Me,coord,val,rhogp)
  implicit none
  double precision                :: Me(:,:), coord(:,:)
  double precision                :: val

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR,ie

  double precision                :: rhogp(8), rhogpCurr
  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8)
  DOUBLE PRECISION                :: V(24,3)

  JT=MATMUL(DNR,transpose(COORD))
  Me=0D0
  DO GP_NR=1,NGP    

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 
   
    V = 0d0 
    V(1:24:3,1)=NR(GP_NR,:)
    V(2:24:3,2)=NR(GP_NR,:)
   V(3:24:3,3)=NR(GP_NR,:)
    rhogpCurr = rhogp(GP_NR)

    Me=Me+rhogpCurr*MATMUL(V,TRANSPOSE(V))*DETJ*val

  END DO
  RETURN
END subroutine c3dtl8_m


! Derivative of mass matrix with respect to rho-tilde
subroutine dc3dtl8_m(fme,coord,val,rhogp,larg,rarg)
  implicit none
  double precision                :: fme(:), coord(:,:), larg(:), rarg(:)
  double precision                :: val

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR,ie

  double precision                :: rhogp(8),rhogpCurr
  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8)
  DOUBLE PRECISION                :: V(24,3), V2(8), tmp2(24,24), tmp3

  JT=MATMUL(DNR,transpose(COORD))
  fme=0D0
  DO GP_NR=1,NGP    

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 
    
    V = 0d0
    V(1:24:3,1)=NR(GP_NR,:)
    V(2:24:3,2)=NR(GP_NR,:)
    V(3:24:3,3)=NR(GP_NR,:)
    rhogpCurr = rhogp(GP_NR)

    V2(:)=NR(GP_NR,:)
    
    tmp2 = rhogpCurr*MATMUL(V,TRANSPOSE(V))*DETJ*val
    tmp3 = dot_product(larg,MATMUL(tmp2,rarg))

    fme=fme+V2*tmp3

  END DO
  RETURN
END subroutine dc3dtl8_m













subroutine c3dtl8_egamma(ke,coord,D,ed,es,gammagp)
  implicit none
  double precision                :: ke(:,:), coord(:,:), D(:,:,:),gammagp(:)
  double precision                :: ed(:), es(:,:)

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8)
  double precision                :: BL0(6,24), BL(6,24)
  double precision                :: Dgp(6,6)
  DOUBLE PRECISION                :: STRESS(3,3)
  DOUBLE PRECISION                :: HE(9,24)
  DOUBLE PRECISION                :: A(6,9)
  DOUBLE PRECISION                :: dg(9)
  DOUBLE PRECISION                :: R(9,9)


  JT=MATMUL(DNR,transpose(COORD))

  KE=0D0
  DO GP_NR=1,NGP

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:24:3)=dNx(1,:)
    BL0(2,2:24:3)=dNx(2,:)
    BL0(3,3:24:3)=dNx(3,:)
    BL0(4,1:24:3)=dNx(2,:)
    BL0(4,2:24:3)=dNx(1,:)
    BL0(5,1:24:3)=dNx(3,:)
    BL0(5,3:24:3)=dNx(1,:)
    BL0(6,2:24:3)=dNx(3,:)
    BL0(6,3:24:3)=dNx(2,:)
  
    He=0d0
    He(1:3,1:24:3)=dNx(1:3,:)
    He(4:6,2:24:3)=dNx(1:3,:)
    He(7:9,3:24:3)=dNx(1:3,:)

! displacement gradient
    dg=matmul(He,ed)*gammagp(GP_NR)
    A=0d0
    A([1,4,5],1)=dg(1:3)
    A([2,4,6],2)=(/dg(2),dg(1),dg(3)/)
    A([3,5,6],3)=(/dg(3),dg(1),dg(2)/)
    A([1,4,5],4)=dg(4:6)
    A([2,4,6],5)=(/dg(5),dg(4),dg(6)/)
    A([3,5,6],6)=(/dg(6),dg(4),dg(5)/)
    A([1,4,5],7)=dg(7:9)
    A([2,4,6],8)=(/dg(8),dg(7),dg(9)/)
    A([3,5,6],9)=(/dg(9),dg(7),dg(8)/)

    BL=BL0+MATMUL(A,He)

    stress(1,:)=(/es(1,gp_nr), es(4,gp_nr), es(5,gp_nr)/)
    stress(2,:)=(/es(4,gp_nr), es(2,gp_nr), es(6,gp_nr)/)
    stress(3,:)=(/es(5,gp_nr), es(6,gp_nr), es(3,gp_nr)/)
    R=0D0
    R(1:3,1:3)=STRESS(:,:)
    R(4:6,4:6)=STRESS(:,:)
    R(7:9,7:9)=STRESS(:,:)

    Dgp=D(:,:,gp_nr)
  
    KE=KE+(MATMUL(TRANSPOSE(BL),MATMUL(Dgp,BL)) &
           +MATMUL(MATMUL(TRANSPOSE(He),R),He))*DETJ

  END DO
  RETURN
END subroutine c3dtl8_egamma






subroutine c3dtl8_egammaSens(fke,coord,D,ed,larg,rarg,es,gammagp)
  implicit none
  double precision                :: fke(:), kegp, coord(:,:), D(:,:,:),gammagp(:)
  double precision                :: ed(:), es(:,:), larg(:), rarg(:)

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8)
  double precision                :: BL0(6,24), BL(6,24)
  double precision                :: Dgp(6,6)
  DOUBLE PRECISION                :: STRESS(3,3)
  DOUBLE PRECISION                :: HE(9,24)
  DOUBLE PRECISION                :: A(6,9)
  DOUBLE PRECISION                :: dg(9)
  DOUBLE PRECISION                :: R(9,9), largTMP(24), rargTMP(24), V2(8)


  JT=MATMUL(DNR,transpose(COORD))
  kegp = 0d0
  fke = 0d0
  DO GP_NR=1,NGP

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:24:3)=dNx(1,:)
    BL0(2,2:24:3)=dNx(2,:)
    BL0(3,3:24:3)=dNx(3,:)
    BL0(4,1:24:3)=dNx(2,:)
    BL0(4,2:24:3)=dNx(1,:)
    BL0(5,1:24:3)=dNx(3,:)
    BL0(5,3:24:3)=dNx(1,:)
    BL0(6,2:24:3)=dNx(3,:)
    BL0(6,3:24:3)=dNx(2,:)
  
    He=0d0
    He(1:3,1:24:3)=dNx(1:3,:)
    He(4:6,2:24:3)=dNx(1:3,:)
    He(7:9,3:24:3)=dNx(1:3,:)

! displacement gradient
    dg=matmul(He,ed)*gammagp(GP_NR)
    A=0d0
    A([1,4,5],1)=dg(1:3)
    A([2,4,6],2)=(/dg(2),dg(1),dg(3)/)
    A([3,5,6],3)=(/dg(3),dg(1),dg(2)/)
    A([1,4,5],4)=dg(4:6)
    A([2,4,6],5)=(/dg(5),dg(4),dg(6)/)
    A([3,5,6],6)=(/dg(6),dg(4),dg(5)/)
    A([1,4,5],7)=dg(7:9)
    A([2,4,6],8)=(/dg(8),dg(7),dg(9)/)
    A([3,5,6],9)=(/dg(9),dg(7),dg(8)/)

    BL=BL0+MATMUL(A,He)

    stress(1,:)=(/es(1,gp_nr), es(4,gp_nr), es(5,gp_nr)/)
    stress(2,:)=(/es(4,gp_nr), es(2,gp_nr), es(6,gp_nr)/)
    stress(3,:)=(/es(5,gp_nr), es(6,gp_nr), es(3,gp_nr)/)
    R=0D0
    R(1:3,1:3)=STRESS(:,:)
    R(4:6,4:6)=STRESS(:,:)
    R(7:9,7:9)=STRESS(:,:)

    V2(:)=NR(GP_NR,:)
    Dgp=D(:,:,gp_nr)
    largTMP = larg!*gammagp(GP_NR)
    rargTMP = rarg!*gammagp(GP_NR)
  
    KEGP=dot_product(largTMP,(matmul((MATMUL(TRANSPOSE(BL),MATMUL(Dgp,BL)) &
           +MATMUL(MATMUL(TRANSPOSE(He),R),He)),rargTMP)))*DETJ

    fke=fke+V2*KEGP

  END DO
  RETURN
END subroutine c3dtl8_egammaSens













! calculates the total deformation gradient
subroutine c3dtl8_d(dg,coord,ed)
  implicit none

  double precision                :: dg(:,:), coord(:,:),  ed(:)

  integer                         :: gp_nr

  double precision                :: JT(24,3), JTinv(3,3)
  double precision                :: TMP(3,8), dNx(3,8)
  double precision                :: He(9,24)

  JT=MATMUL(DNR,transpose(COORD))
   
  do GP_NR=1,8

   CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
   dNx=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
 
   He=0d0
   He(1:3,1:24:3)=dNx(1:3,:)
   He(4:6,2:24:3)=dNx(1:3,:)
   He(7:9,3:24:3)=dNx(1:3,:)
     
   dg(:,gp_nr)=Ifm+MATMUL(HE,ED)

  end do

  return
end subroutine c3dtl8_d















subroutine c3dtl8_f(ef,coord,ed,es)
  implicit none

  double precision                :: eF(:), coord(:,:),  ed(:), es(:,:)

  integer, parameter              :: NGP=8
  integer                         :: gp_nr

  double precision                :: JT(24,3), JTinv(3,3)
  double precision                :: TMP(3,8), dNx(3,8), DETJ
  double precision                :: BL0(6,24), He(9,24)
  double precision                :: A(6,9), dg(9)

 
  JT=MATMUL(DNR,transpose(COORD))

  ef=0D0
  DO GP_NR=1,NGP

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    dNx=MATMUL(JTINV,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:24:3) =dNx(1,:)
    BL0(2,2:24:3) =dNx(2,:)
    BL0(3,3:24:3) =dNx(3,:)
    BL0(4,1:24:3) =dNx(2,:)
    BL0(4,2:24:3) =dNx(1,:)
    BL0(5,1:24:3) =dNx(3,:)
    BL0(5,3:24:3) =dNx(1,:)
    BL0(6,2:24:3) =dNx(3,:)
    BL0(6,3:24:3) =dNx(2,:)

    He=0d0
    He(1:3,1:24:3)=dNx(1:3,:)
    He(4:6,2:24:3)=dNx(1:3,:)
    He(7:9,3:24:3)=dNx(1:3,:)

! displacement gradient
    dg=matmul(He,ed)
    A=0d0
    A([1,4,5],1)=dg(1:3)
    A([2,4,6],2)=(/dg(2),dg(1),dg(3)/)
    A([3,5,6],3)=(/dg(3),dg(1),dg(2)/)
    A([1,4,5],4)=dg(4:6)
    A([2,4,6],5)=(/dg(5),dg(4),dg(6)/)
    A([3,5,6],6)=(/dg(6),dg(4),dg(5)/)
    A([1,4,5],7)=dg(7:9)
    A([2,4,6],8)=(/dg(8),dg(7),dg(9)/)
    A([3,5,6],9)=(/dg(9),dg(7),dg(8)/)
 
    EF=EF+MATMUL(TRANSPOSE(BL0+MATMUL(A,He)),es(:,gp_nr))*detJ

  END DO
  RETURN
END subroutine c3dtl8_f



subroutine c3dtl8_fgamma(ef,coord,ed,es,gammagp)
  implicit none

  double precision                :: eF(:), coord(:,:),  ed(:), es(:,:), gammagp(:)

  integer, parameter              :: NGP=8
  integer                         :: gp_nr

  double precision                :: JT(24,3), JTinv(3,3)
  double precision                :: TMP(3,8), dNx(3,8), DETJ
  double precision                :: BL0(6,24), He(9,24)
  double precision                :: A(6,9), dg(9)

 
  JT=MATMUL(DNR,transpose(COORD))

  ef=0D0
  DO GP_NR=1,NGP

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    dNx=MATMUL(JTINV,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:24:3) =dNx(1,:)
    BL0(2,2:24:3) =dNx(2,:)
    BL0(3,3:24:3) =dNx(3,:)
    BL0(4,1:24:3) =dNx(2,:)
    BL0(4,2:24:3) =dNx(1,:)
    BL0(5,1:24:3) =dNx(3,:)
    BL0(5,3:24:3) =dNx(1,:)
    BL0(6,2:24:3) =dNx(3,:)
    BL0(6,3:24:3) =dNx(2,:)

    He=0d0
    He(1:3,1:24:3)=dNx(1:3,:)
    He(4:6,2:24:3)=dNx(1:3,:)
    He(7:9,3:24:3)=dNx(1:3,:)

! displacement gradient
    dg=matmul(He,ed)*gammagp(GP_NR)
    A=0d0
    A([1,4,5],1)=dg(1:3)
    A([2,4,6],2)=(/dg(2),dg(1),dg(3)/)
    A([3,5,6],3)=(/dg(3),dg(1),dg(2)/)
    A([1,4,5],4)=dg(4:6)
    A([2,4,6],5)=(/dg(5),dg(4),dg(6)/)
    A([3,5,6],6)=(/dg(6),dg(4),dg(5)/)
    A([1,4,5],7)=dg(7:9)
    A([2,4,6],8)=(/dg(8),dg(7),dg(9)/)
    A([3,5,6],9)=(/dg(9),dg(7),dg(8)/)
 
    EF=EF+MATMUL(TRANSPOSE(BL0+MATMUL(A,He)),es(:,gp_nr))*detJ

  END DO
  RETURN
END subroutine c3dtl8_fgamma







subroutine c3dtl8_fgammaSens(ffe,coord,ed,edvarphi,es,gammagp)
  implicit none

  double precision                :: ffe(:), fegp, coord(:,:),  ed(:), es(:,:), gammagp(:), edvarphi(:)

  integer, parameter              :: NGP=8
  integer                         :: gp_nr

  double precision                :: JT(24,3), JTinv(3,3)
  double precision                :: TMP(3,8), dNx(3,8), DETJ
  double precision                :: BL0(6,24), He(9,24)
  double precision                :: A(6,9), dg(9), varphiTMP(24), V2(8)

 
  JT=MATMUL(DNR,transpose(COORD))
  fegp = 0d0
  ffe = 0d0
  DO GP_NR=1,NGP

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    dNx=MATMUL(JTINV,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:24:3) =dNx(1,:)
    BL0(2,2:24:3) =dNx(2,:)
    BL0(3,3:24:3) =dNx(3,:)
    BL0(4,1:24:3) =dNx(2,:)
    BL0(4,2:24:3) =dNx(1,:)
    BL0(5,1:24:3) =dNx(3,:)
    BL0(5,3:24:3) =dNx(1,:)
    BL0(6,2:24:3) =dNx(3,:)
    BL0(6,3:24:3) =dNx(2,:)

    He=0d0
    He(1:3,1:24:3)=dNx(1:3,:)
    He(4:6,2:24:3)=dNx(1:3,:)
    He(7:9,3:24:3)=dNx(1:3,:)

! displacement gradient
    dg=matmul(He,ed)*gammagp(GP_NR)
    A=0d0
    A([1,4,5],1)=dg(1:3)
    A([2,4,6],2)=(/dg(2),dg(1),dg(3)/)
    A([3,5,6],3)=(/dg(3),dg(1),dg(2)/)
    A([1,4,5],4)=dg(4:6)
    A([2,4,6],5)=(/dg(5),dg(4),dg(6)/)
    A([3,5,6],6)=(/dg(6),dg(4),dg(5)/)
    A([1,4,5],7)=dg(7:9)
    A([2,4,6],8)=(/dg(8),dg(7),dg(9)/)
    A([3,5,6],9)=(/dg(9),dg(7),dg(8)/)

    V2(:)=NR(GP_NR,:)
    varphiTMP = edvarphi!*gammagp(GP_NR)
 
    fegp=dot_product(varphiTMP,MATMUL(TRANSPOSE(BL0+MATMUL(A,He)),es(:,gp_nr)))*detJ

    ffe=ffe+V2*fegp

  END DO
  RETURN
END subroutine c3dtl8_fgammaSens










subroutine elmvol(evol,coord)
  implicit none
  double precision                :: evol, coord(:,:)

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: DNX(3,8)

  JT=MATMUL(DNR,transpose(COORD))

  evol=0D0
  DO GP_NR=1,NGP
    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 
    evol=evol+DETJ
  END DO
  RETURN
END subroutine elmvol






subroutine c3dtl8_emu2(coord,ed,edmutilde, Emutilde,dgmutilde, gammagp, ones)
! Input:  coord,ed,edmutilde, gammagp
! Outout: Emutilde,dgmutilde
  implicit none
  double precision                :: coord(:,:)
  double precision                :: ed(:), edmutilde(:), Emutilde(:,:), gammagp(:), ones(:)

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8)
  double precision                :: BL0(6,24), BL(6,24), Bphi(6,24)
  DOUBLE PRECISION                :: HE(9,24)
  DOUBLE PRECISION                :: A(6,9), Aphi(6,9)
  DOUBLE PRECISION                :: dg(9), dgmutilde(:,:)
  DOUBLE PRECISION                :: R(9,9),p,rho,sw, Rphi(9,9), C(3,3), id(3,3), iC(3,3)

  id=0
  id(1,1)=1D0
  id(2,2)=1D0
  id(3,3)=1D0

  JT=MATMUL(DNR,transpose(COORD))
  
  DO GP_NR=1,NGP

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:24:3)=dNx(1,:)
    BL0(2,2:24:3)=dNx(2,:)
    BL0(3,3:24:3)=dNx(3,:)
    BL0(4,1:24:3)=dNx(2,:)
    BL0(4,2:24:3)=dNx(1,:)
    BL0(5,1:24:3)=dNx(3,:)
    BL0(5,3:24:3)=dNx(1,:)
    BL0(6,2:24:3)=dNx(3,:)
    BL0(6,3:24:3)=dNx(2,:)

    He=0d0
    He(1:3,1:24:3)=dNx(1:3,:)
    He(4:6,2:24:3)=dNx(1:3,:)
    He(7:9,3:24:3)=dNx(1:3,:)

! displacement gradient
    dg=matmul(He,ed)* gammagp(gp_nr)


    A=0d0
    A([1,4,5],1)=dg(1:3)
    A([2,4,6],2)=(/dg(2),dg(1),dg(3)/)
    A([3,5,6],3)=(/dg(3),dg(1),dg(2)/)
    A([1,4,5],4)=dg(4:6)
    A([2,4,6],5)=(/dg(5),dg(4),dg(6)/)
    A([3,5,6],6)=(/dg(6),dg(4),dg(5)/)
    A([1,4,5],7)=dg(7:9)
    A([2,4,6],8)=(/dg(8),dg(7),dg(9)/)
    A([3,5,6],9)=(/dg(9),dg(7),dg(8)/)
    BL=BL0+MATMUL(A,He)


    dgmutilde(:,GP_Nr)=matmul(He,edmutilde)!*gammagp(gp_nr)
    Emutilde(:,GP_Nr) =matmul(BL,edmutilde)!*gammagp(gp_nr)! Lagrange strain-type tensor
  END DO
  RETURN
END subroutine c3dtl8_emu2









subroutine c3dtl8_emu(coord,ed,edmutilde, Emutilde,dgmutilde, gammagp)
! Input:  coord,ed,edmutilde, gammagp
! Outout: Emutilde,dgmutilde
  implicit none
  double precision                :: coord(:,:)
  double precision                :: ed(:), edmutilde(:), Emutilde(:,:), gammagp(:)

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8)
  double precision                :: BL0(6,24), BL(6,24), Bphi(6,24)
  DOUBLE PRECISION                :: HE(9,24)
  DOUBLE PRECISION                :: A(6,9), Aphi(6,9)
  DOUBLE PRECISION                :: dg(9), dgmutilde(:,:)
  DOUBLE PRECISION                :: R(9,9),p,rho,sw, Rphi(9,9), C(3,3), id(3,3), iC(3,3)

  id=0
  id(1,1)=1D0
  id(2,2)=1D0
  id(3,3)=1D0

  JT=MATMUL(DNR,transpose(COORD))
  
  DO GP_NR=1,NGP

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:24:3)=dNx(1,:)
    BL0(2,2:24:3)=dNx(2,:)
    BL0(3,3:24:3)=dNx(3,:)
    BL0(4,1:24:3)=dNx(2,:)
    BL0(4,2:24:3)=dNx(1,:)
    BL0(5,1:24:3)=dNx(3,:)
    BL0(5,3:24:3)=dNx(1,:)
    BL0(6,2:24:3)=dNx(3,:)
    BL0(6,3:24:3)=dNx(2,:)

    He=0d0
    He(1:3,1:24:3)=dNx(1:3,:)
    He(4:6,2:24:3)=dNx(1:3,:)
    He(7:9,3:24:3)=dNx(1:3,:)

! displacement gradient
    dg=matmul(He,ed)* gammagp(gp_nr)


    A=0d0
    A([1,4,5],1)=dg(1:3)
    A([2,4,6],2)=(/dg(2),dg(1),dg(3)/)
    A([3,5,6],3)=(/dg(3),dg(1),dg(2)/)
    A([1,4,5],4)=dg(4:6)
    A([2,4,6],5)=(/dg(5),dg(4),dg(6)/)
    A([3,5,6],6)=(/dg(6),dg(4),dg(5)/)
    A([1,4,5],7)=dg(7:9)
    A([2,4,6],8)=(/dg(8),dg(7),dg(9)/)
    A([3,5,6],9)=(/dg(9),dg(7),dg(8)/)
    BL=BL0+MATMUL(A,He)


    dgmutilde(:,GP_Nr)=matmul(He,edmutilde)*gammagp(gp_nr)
    Emutilde(:,GP_Nr) =matmul(BL,edmutilde)*gammagp(gp_nr)! Lagrange strain-type tensor
  END DO
  RETURN
END subroutine c3dtl8_emu



! 8-node brick element  based on total Lagrangian formulation	
subroutine c3dtl8_e2(ke,coord,D,ed,es)
  implicit none
  double precision                :: ke(:,:), coord(:,:), D(:,:,:)
  double precision                :: ed(:), es(:,:)

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8)
  double precision                :: BL0(6,24), BL(6,24)
  double precision                :: Dgp(6,6)
  DOUBLE PRECISION                :: STRESS(3,3)
  DOUBLE PRECISION                :: HE(9,24)
  DOUBLE PRECISION                :: A(6,9)
  DOUBLE PRECISION                :: dg(9)
  DOUBLE PRECISION                :: R(9,9)


  JT=MATMUL(DNR,transpose(COORD))

  KE=0D0
  DO GP_NR=1,NGP

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:24:3)=dNx(1,:)
    BL0(2,2:24:3)=dNx(2,:)
    BL0(3,3:24:3)=dNx(3,:)
    BL0(4,1:24:3)=dNx(2,:)
    BL0(4,2:24:3)=dNx(1,:)
    BL0(5,1:24:3)=dNx(3,:)
    BL0(5,3:24:3)=dNx(1,:)
    BL0(6,2:24:3)=dNx(3,:)
    BL0(6,3:24:3)=dNx(2,:)
  
    He=0d0
    He(1:3,1:24:3)=dNx(1:3,:)
    He(4:6,2:24:3)=dNx(1:3,:)
    He(7:9,3:24:3)=dNx(1:3,:)

! displacement gradient
    dg=matmul(He,ed)
    A=0d0
    A([1,4,5],1)=dg(1:3)
    A([2,4,6],2)=(/dg(2),dg(1),dg(3)/)
    A([3,5,6],3)=(/dg(3),dg(1),dg(2)/)
    A([1,4,5],4)=dg(4:6)
    A([2,4,6],5)=(/dg(5),dg(4),dg(6)/)
    A([3,5,6],6)=(/dg(6),dg(4),dg(5)/)
    A([1,4,5],7)=dg(7:9)
    A([2,4,6],8)=(/dg(8),dg(7),dg(9)/)
    A([3,5,6],9)=(/dg(9),dg(7),dg(8)/)

    BL=BL0+MATMUL(A,He)

    stress(1,:)=(/es(1,gp_nr), es(4,gp_nr), es(5,gp_nr)/)
    stress(2,:)=(/es(4,gp_nr), es(2,gp_nr), es(6,gp_nr)/)
    stress(3,:)=(/es(5,gp_nr), es(6,gp_nr), es(3,gp_nr)/)
    R=0D0
    R(1:3,1:3)=STRESS(:,:)
    R(4:6,4:6)=STRESS(:,:)
    R(7:9,7:9)=STRESS(:,:)

    Dgp=D(:,:,gp_nr)
  
    KE=KE+(MATMUL(TRANSPOSE(BL),MATMUL(Dgp,BL)) &
           +MATMUL(MATMUL(TRANSPOSE(He),R),He))*DETJ

  END DO
  RETURN
END subroutine c3dtl8_e2


subroutine c3dtl8_b(efb,coord,b)
  implicit none

  double precision                :: efb(:), coord(:,:), b(:,:)

  integer, parameter              :: NGP=8
  integer                         :: gp_nr

  double precision                :: JT(24,3), DetJ, JTinv(3,3)
  double precision                :: N(3,24), btmp(3)
 
  JT=MATMUL(DNR,transpose(COORD))

  efb=0D0
  DO GP_NR=1,NGP

    DETJ=det3(JT(index(:,gp_nr),:)) 

    N = 0d0
    N(1,1:24:3) = NR(GP_NR,:)
    N(2,2:24:3) = NR(GP_NR,:)
    N(3,3:24:3) = NR(GP_NR,:)
    
    bTMP = b(:,gp_nr)

    EFB=EFB+MATMUL(TRANSPOSE(N),bTMP)*detJ

  END DO
  RETURN

end subroutine c3dtl8_b


subroutine c3dtl8_bsens(efb,coord,b,arg)
  implicit none

  double precision                :: efb(:), coord(:,:), b(:,:), arg(:)

  integer, parameter              :: NGP=8
  integer                         :: gp_nr

  double precision                :: JT(24,3), DetJ, JTinv(3,3)
  double precision                :: N(3,24), btmp(3), efbgp, V2(8)
 
  JT=MATMUL(DNR,transpose(COORD))

  efb=0D0
  DO GP_NR=1,NGP

    DETJ=det3(JT(index(:,gp_nr),:)) 

    N = 0d0
    N(1,1:24:3) = NR(GP_NR,:)
    N(2,2:24:3) = NR(GP_NR,:)
    N(3,3:24:3) = NR(GP_NR,:)
    
    V2(:)=NR(GP_NR,:)
    bTMP = b(:,gp_nr)

    EFBgp=dot_product(MATMUL(TRANSPOSE(N),bTMP),arg)*detJ
    efb=efb+V2*efbgp

  END DO
  RETURN

end subroutine c3dtl8_bsens




end module elem_large_cont_3d
