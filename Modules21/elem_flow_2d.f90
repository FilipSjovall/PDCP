module elem_flow_2d

! Last rev.:
!  A. Dalklint 2020-04-29
!   - initial version 4-node brick element derived from 
!     elem_large_cont_2d 

! 4 node brick element, bilinear displacement interpolation
! 2x2 gauss integration points (in total 4), the weight for every
! gauss point is 1 for this scheme, and it therefore omitted in the code

! The total lagrangian formulation is implemented
use matrix_util, only: inv3, det3, inv2, det2

implicit none

double precision  xsi(4), eta(4), G1, Ifm(4)

double precision, allocatable                :: BL0_G(:,:,:,:)
DOUBLE PRECISION, allocatable                :: HE_G(:,:,:,:)
integer           index(2,4)
parameter        (G1=0.577350269189626D0)
parameter        (xsi=(/-1D0,  1D0, 1D0, -1D0/)*G1 )
parameter        (eta=(/-1D0, -1D0, 1D0,  1D0/)*G1 )
parameter        (Ifm=(/1d0,0d0,0d0,1d0/))
parameter        (index=[(/1,2/),(/3,4/),(/5,6/),(/7,8/)])
private xsi, eta, G1, Ifm, index

double precision  DNR(8,4) ! Maybe one should use DNR(4,8) instead for speed
! derivate of shape functions with respect to xsi
data             (DNR(1:8:2,1)= (ETA-1d0)/4D0)
data             (DNR(1:8:2,2)=-(ETA-1d0)/4D0)
data             (DNR(1:8:2,3)= (ETA+1d0)/4D0)
data             (DNR(1:8:2,4)=-(ETA+1d0)/4D0)
! derivate of shape functions with respect to eta
data             (DNR(2:8:2,1)= (XSI-1d0)/4D0)
data             (DNR(2:8:2,2)=-(XSI+1d0)/4D0)
data             (DNR(2:8:2,3)= (XSI+1d0)/4D0)
data             (DNR(2:8:2,4)=-(XSI-1d0)/4D0)
private DNR

double precision  NR(4,4) ! 4 Gauss points and 4 nodes
data             (NR(:,1)= (XSI-1d0)*(ETA-1d0)/4D0)
data             (NR(:,2)= -(XSI+1d0)*(ETA-1d0)/4D0)
data             (NR(:,3)= (XSI+1d0)*(ETA+1d0)/4D0)
data             (NR(:,4)= -(XSI-1d0)*(ETA+1d0)/4D0)
private NR

!------------------------------------------------------------------------------


interface fl2d4_e
  module procedure fl2d4_e2
end interface
private fl2d4_e2

interface fl2d4_m
  module procedure fl2d4_m2
end interface
private fl2d4_m2

interface fl2d4_vi
  module procedure fl2d4_vol1
  module procedure fl2d4_vol2
  module procedure fl2d4_vol3
  module procedure fl2d4_vol4
end interface
private fl2d4_vol1, fl2d4_vol2, fl2d4_vol3, fl2d4_vol4

interface fl2d4_gp
  module procedure fl2d4_gpval
end interface

interface fl2d4_bf
  module procedure fl2d4_bf1
  module procedure fl2d4_bf1x
  module procedure fl2d4_bf1a
  module procedure fl2d4_bf2
  module procedure fl2d4_bf3
end interface
private fl2d4_bf1, fl2d4_bf1a, fl2d4_bf2, fl2d4_bf3


!------------------------------------------------------------------------------
contains



subroutine fl2d4_e2(ke,coord,t,val)
  implicit none
  double precision                :: ke(:,:), coord(:,:)
  double precision                :: val, t

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4)

  JT=MATMUL(DNR,transpose(COORD))

  KE=0D0
  DO GP_NR=1,NGP

    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 

    KE=KE+MATMUL(TRANSPOSE(dNx),dNx)*DETJ*val*t

  END DO
  RETURN
END subroutine fl2d4_e2



subroutine fl2d4_m2(Me,coord,t,val)
  implicit none
  double precision                :: Me(:,:), coord(:,:)
  double precision                :: val, t

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,2), DetJ, JTinv(2,2)
  double precision                :: TMP(2,4), DNX(2,4)
  DOUBLE PRECISION                :: V(4,1)

  JT=MATMUL(DNR,transpose(COORD))

  ME=0D0
  DO GP_NR=1,NGP

    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 
    
    V(:,1)=NR(GP_NR,:)

    ME=ME+MATMUL(V,TRANSPOSE(V) )*DETJ*val*t
   
  END DO
  RETURN
END subroutine fl2d4_m2



! body force vector
subroutine fl2d4_bf1(fe,coord,t,gpval1)
  implicit none
  double precision                :: coord(:,:), fe(:), t
  double precision                :: gpval1(:)

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4)

  JT=MATMUL(DNR,transpose(COORD))

  FE=0D0
  DO GP_NR=1,NGP
    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 
  
    FE=FE+NR(GP_NR,:)*DETJ*gpval1(gp_nr)*t

  END DO
  RETURN
END subroutine fl2d4_bf1


! body force vector
subroutine fl2d4_bf1x(fe,coord,t,gpval1,MrX)
  implicit none
  double precision                :: coord(:,:), fe(:), t
  double precision                :: gpval1(:)

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR,ie,MrX

  double precision                :: JT(R2,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4)



  JT=MATMUL(DNR,transpose(COORD))

  FE=0D0
  DO GP_NR=1,NGP
    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:))   
    FE=FE+DETJ*gpval1(gp_nr)*t

  END DO
  RETURN
END subroutine fl2d4_bf1x



subroutine fl2d4_bf1a(fe,coord,t,val)
  implicit none
  double precision                :: coord(:,:), fe(:), t
  double precision                :: val

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4)

  JT=MATMUL(DNR,transpose(COORD))

  FE=0D0
  DO GP_NR=1,NGP
    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 
  
    FE=FE+NR(GP_NR,:)*DETJ*val*t

  END DO
  RETURN
END subroutine fl2d4_bf1a



subroutine fl2d4_bf2(fe,coord,t,gpval1,gpval2)
  implicit none
  double precision                :: coord(:,:), fe(:), t
  double precision                :: gpval1(:), gpval2(:)

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4)

  JT=MATMUL(DNR,transpose(COORD))

  FE=0D0
  DO GP_NR=1,NGP
    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 
  
    FE=FE+NR(GP_NR,:)*DETJ*gpval1(gp_nr)*gpval2(gp_nr)*t

  END DO
  RETURN
END subroutine fl2d4_bf2


subroutine fl2d4_bf3(fe,coord,t,ef,val)
  implicit none
  double precision                :: coord(:,:), fe(:), t
  double precision                :: ef(:,:), val

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4)
  DOUBLE PRECISION                :: F(3,3), Jvol

  JT=MATMUL(DNR,transpose(COORD))

  FE=0D0
  DO GP_NR=1,NGP

    F(1,:)=(/ef(1,GP_NR), ef(2,GP_NR), 0d0/)
    F(2,:)=(/ef(3,GP_NR), ef(4,GP_NR), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)

    Jvol=det2(F)      
    DETJ=det2(JT(index(:,gp_nr),:)) 

    FE=FE+Jvol*NR(GP_NR,:)*DETJ*val*t
   
  END DO

 
  RETURN
END subroutine fl2d4_bf3


! Extract gauss point values from nodal quantities
subroutine fl2d4_gpval(rhogp,rho)
  implicit none

  double precision                ::  rho(:)
  double precision                ::  rhogp(4)
 
  rhogp=MATMUL(NR,rho)

  RETURN
END subroutine fl2d4_gpval




! Calculate volume intergral
subroutine fl2d4_vol1(evol,coord,t)
  implicit none
  double precision                :: evol, coord(:,:), t

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR

  double precision                :: JT(R2,2), DetJ

  JT=MATMUL(DNR,transpose(COORD))

  evol=0D0
  DO GP_NR=1,NGP
    DETJ=det2(JT(index(:,gp_nr),:)) 
    evol=evol+DETJ*t
  END DO
  RETURN
END subroutine fl2d4_vol1


subroutine fl2d4_vol2(evol,coord,t,efin)
  implicit none
  double precision                :: evol, coord(:,:), efin(:,:), ef(4), t

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR

  double precision                :: JT(R2,2), DetJ, F(3,3), J

  JT=MATMUL(DNR,transpose(COORD))

  evol=0D0
  DO GP_NR=1,NGP

    ef=efin(:,gp_nr)

    F(1,:)=(/ef(1), ef(2), 0d0/)
    F(2,:)=(/ef(3), ef(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)
    J=det2(F)

    DETJ=det2(JT(index(:,gp_nr),:)) 
    evol=evol+DETJ*J*t
  END DO
  RETURN
END subroutine fl2d4_vol2


subroutine fl2d4_vol3(evol,coord,t,efin,gpval)
  implicit none
  double precision                :: evol, coord(:,:), efin(:,:), gpval(:), t

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR

  double precision                :: JT(R2,2), DetJ, F(3,3), J, ef(4)

  JT=MATMUL(DNR,transpose(COORD))

  evol=0D0
  DO GP_NR=1,NGP

    ef=efin(:,gp_nr)
    F(1,:)=(/ef(1), ef(2), 0d0/)
    F(2,:)=(/ef(3), ef(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)
    J=det2(F)

    DETJ=det2(JT(index(:,gp_nr),:)) 
    evol=evol+DETJ*J*gpval(GP_NR)*t
  END DO
  RETURN
END subroutine fl2d4_vol3


subroutine fl2d4_vol4(evol,coord,t,gpval)
  implicit none
  double precision                :: evol, coord(:,:), gpval(:), t

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR

  double precision                :: JT(R2,2), DetJ

  JT=MATMUL(DNR,transpose(COORD))

  evol=0D0
  DO GP_NR=1,NGP

    DETJ=det2(JT(index(:,gp_nr),:)) 
    evol=evol+DETJ*gpval(GP_NR)*t
  END DO
  RETURN
END subroutine fl2d4_vol4


end module elem_flow_2d
