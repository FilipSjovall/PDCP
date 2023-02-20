!This file contain: PLANI8E_PE, PLANI8F_PE, PLANI8M_PE and DEFORMATION8I_PE




!----------------------------------------------------------------------------------!
! NAME    : PLANI8E_PE                                                             !
!                                                                                   !
! PURPOSE : Calculate the stiffness matrix for a 8 node isoparametric              !
!           element in plane strain, with the theory for large deformations        !
!                                                                                   !
! INPUT   : coord       Node coordinates for the element [4 x 2]                   !
!                                                                                   !
!           D_el        The constitutive matrix for each integration               !
!                       point  [gp x 3 x 3]                                        !
!                                                                                   !
!           ed          Nodal displacment  [u1 u2 u3 u4 ...]                       !
!                                                                                   !
!           es          Second Piola-Kirchoff stress [gp x 3 x 3]                  ! 
!                                                                                   !
!           ngp         Number of integration points                               !
!                                                                                   !
!           t           Thickness of the element                                   !
!                                                                                   !
!                                                                                   !
! OUTPUT  : Ke          Element stiffness matrix  [16 x 16]                          !        
!                                                                                   !
!                                                                                   !
!----------------------------------------------------------------------------------!
!                                                                                   !
SUBROUTINE PLANI8E_AXI(Ke,coord,D_el,ed,es,ngp)

IMPLICIT NONE
INTEGER                                       :: r2, ngp, i 
DOUBLE PRECISION                              :: g1, g2, w1, w2,l11,l12,l22,l21,l33,detJ,xbar1
DOUBLE PRECISION, PARAMETER                   :: pi=3.14159265358979D0
DOUBLE PRECISION, DIMENSION(ngp,8)            :: N
DOUBLE PRECISION, DIMENSION(ngp)              :: eta, xsi,wp, w_1, w_2
DOUBLE PRECISION, DIMENSION(2*ngp,8),TARGET   :: dNr
DOUBLE PRECISION, DIMENSION(16,16)            :: Ke
DOUBLE PRECISION, DIMENSION(4,16)             :: BL0, BL1, BL
DOUBLE PRECISION, DIMENSION(5,16)             :: BNL
DOUBLE PRECISION, DIMENSION(2*ngp,2)          :: JT
DOUBLE PRECISION, DIMENSION(2,8)              :: dNx
DOUBLE PRECISION, DIMENSION(8,2)              :: coord
DOUBLE PRECISION, DIMENSION(5,5)              :: S
DOUBLE PRECISION, POINTER                     :: dNr_p(:,:)
DOUBLE PRECISION, DIMENSION(4,4)              :: D
DOUBLE PRECISION, DIMENSION(2,2)              :: JTinv
DOUBLE PRECISION, POINTER                     :: kappa_tan
DOUBLE PRECISION, DIMENSION(ngp,4,4)          :: D_el
DOUBLE PRECISION, DIMENSION(16)               :: ed
DOUBLE PRECISION, DIMENSION(ngp,3,3)          :: es




IF (ngp==1) THEN
    g1=0.0D0; w1=2.0D0
    xsi=(/g1/)
    eta=(/g1/)
    w_1=(/w1 /)
    w_2=(/w1 /)
ELSE IF (ngp==4) THEN
    g1=0.577350269189626D0; w1=1D0
    xsi=(/-g1, g1,-g1, g1/)
    eta=(/-g1, -g1,g1, g1/)
    w_1=(/w1, w1, w1, w1/)
    w_2=w_1
ELSE IF (ngp==9) THEN ! ????????????
    g1=0.774596669241483D0; g2=0.D0;
    w1=0.555555555555555D0; w2=0.888888888888888D0;
    xsi=(/-g1,-g2, g1,-g1, g2, g1,-g1, g2, g1/)
    eta=(/-g1,-g1,-g1, g2, g2, g2, g1, g1, g1/)
    w_1=(/w1, w2, w1, w1, w2, w1, w1, w2, w1/)
    w_2=(/w1, w1, w1, w2, w2, w2, w1, w1, w1/)
ELSE
    WRITE(*,*) 'Used number of integration points not implemented'
    
END IF
  wp=w_1*w_2


!--------- shape functions -----------------------------------

  N(:,1)=-(1D0-xsi)*(1D0-eta)*(1D0+xsi+eta)/4D0; N(:,5)=(1D0-xsi*xsi)*(1D0-eta)/2D0
  N(:,2)=-(1D0+xsi)*(1D0-eta)*(1D0-xsi+eta)/4D0; N(:,6)=(1D0+xsi)*(1D0-eta*eta)/2D0
  N(:,3)=-(1D0+xsi)*(1D0+eta)*(1D0-xsi-eta)/4D0; N(:,7)=(1D0-xsi*xsi)*(1D0+eta)/2D0
  N(:,4)=-(1D0-xsi)*(1D0+eta)*(1D0+xsi-eta)/4D0; N(:,8)=(1D0-xsi)*(1D0-eta*eta)/2D0

  dNr(1:2*ngp-1:2,1)=-(-(1D0-eta)*(1D0+xsi+eta)+(1D0-xsi)*(1D0-eta))/4D0
  dNr(1:2*ngp-1:2,2)=-( (1D0-eta)*(1D0-xsi+eta)-(1D0+xsi)*(1D0-eta))/4D0
  dNr(1:2*ngp-1:2,3)=-( (1D0+eta)*(1D0-xsi-eta)-(1D0+xsi)*(1D0+eta))/4D0
  dNr(1:2*ngp-1:2,4)=-(-(1D0+eta)*(1D0+xsi-eta)+(1D0-xsi)*(1D0+eta))/4D0
  dNr(1:2*ngp-1:2,5)=-xsi*(1D0-eta)
  dNr(1:2*ngp-1:2,6)=(1D0-eta*eta)/2D0
  dNr(1:2*ngp-1:2,7)=-xsi*(1D0+eta)
  dNr(1:2*ngp-1:2,8)=-(1D0-eta*eta)/2D0
  dNr(2:2*ngp:2,1)=-(-(1D0-xsi)*(1D0+xsi+eta)+(1D0-xsi)*(1D0-eta))/4D0
  dNr(2:2*ngp:2,2)=-(-(1D0+xsi)*(1D0-xsi+eta)+(1D0+xsi)*(1D0-eta))/4D0
  dNr(2:2*ngp:2,3)=-( (1D0+xsi)*(1D0-xsi-eta)-(1D0+xsi)*(1D0+eta))/4D0
  dNr(2:2*ngp:2,4)=-( (1D0-xsi)*(1D0+xsi-eta)-(1D0-xsi)*(1D0+eta))/4D0
  dNr(2:2*ngp:2,5)=-(1D0-xsi*xsi)/2D0
  dNr(2:2*ngp:2,6)=-eta*(1D0+xsi)
  dNr(2:2*ngp:2,7)=(1D0-xsi*xsi)/2D0
  dNr(2:2*ngp:2,8)=-eta*(1D0-xsi)
  

  
  JT=MATMUL(dNr,coord)


!--------- plane strain conditions -----------------------------
Ke=0D0
DO i=1,ngp
 
  detJ=JT(2*i-1,1)*JT(2*i,2)-JT(2*i-1,2)*JT(2*i,1)
  IF (detJ<0) THEN
    WRITE(*,*)'Jacobideterminant equal or less than zero!'
  END IF
  JTinv=1/detJ*RESHAPE((/JT(2*i,2), -JT(2*i,1),-JT(2*i-1,2),JT(2*i-1,1)/),(/2,2/))

  dNr_p=>dNr(2*i-1:2*i,:)

  dNx=MATMUL(JTinv,dNr_P)

    
  xbar1=SUM(N(i,:)*coord(:,1))
 
  BL0(1,1:16:2)=dNx(1,:)
  BL0(2,2:16:2)=dNx(2,:)
  BL0(3,1:16:2)=N(i,:)/xbar1 
  BL0(4,1:16:2)=dNx(2,:)  
  BL0(4,2:16:2)=dNx(1,:)   
   
  
  l11=SUM(dNx(1,:)*ed(1:16:2))
  l22=SUM(dNx(2,:)*ed(2:16:2))
  l12=SUM(dNx(2,:)*ed(1:16:2))
  l21=SUM(dNx(1,:)*ed(2:16:2))
  l33=SUM(N(i,:)*ed(1:16:2)/xbar1)
      
  BL1(1,:)=(/dNx(1,1)*l11, dNx(1,1)*l21, dNx(1,2)*l11, dNx(1,2)*l21,&
             dNx(1,3)*l11, dNx(1,3)*l21, dNx(1,4)*l11, dNx(1,4)*l21,&
             dNx(1,5)*l11, dNx(1,5)*l21, dNx(1,6)*l11, dNx(1,6)*l21,&
             dNx(1,7)*l11, dNx(1,7)*l21, dNx(1,8)*l11, dNx(1,8)*l21/)
  
  BL1(2,:)=(/dNx(2,1)*l12, dNx(2,1)*l22, dNx(2,2)*l12, dNx(2,2)*l22,&
             dNx(2,3)*l12, dNx(2,3)*l22, dNx(2,4)*l12, dNx(2,4)*l22,&
             dNx(2,5)*l12, dNx(2,5)*l22, dNx(2,6)*l12, dNx(2,6)*l22,&
             dNx(2,7)*l12, dNx(2,7)*l22, dNx(2,8)*l12, dNx(2,8)*l22/)

  BL1(3,:)=l33*BL0(3,:); 
  
  BL1(4,:)=(/(l11*dNx(2,1)+l12*dNx(1,1)),  (l21*dNx(2,1)+l22*dNx(1,1)),&
             (l11*dNx(2,2)+l12*dNx(1,2)),  (l21*dNx(2,2)+l22*dNx(1,2)),& 
             (l11*dNx(2,3)+l12*dNx(1,3)),  (l21*dNx(2,3)+l22*dNx(1,3)),& 
             (l11*dNx(2,4)+l12*dNx(1,4)),  (l21*dNx(2,4)+l22*dNx(1,4)),&
             (l11*dNx(2,5)+l12*dNx(1,5)),  (l21*dNx(2,5)+l22*dNx(1,5)),&
             (l11*dNx(2,6)+l12*dNx(1,6)),  (l21*dNx(2,6)+l22*dNx(1,6)),& 
             (l11*dNx(2,7)+l12*dNx(1,7)),  (l21*dNx(2,7)+l22*dNx(1,7)),& 
             (l11*dNx(2,8)+l12*dNx(1,8)),  (l21*dNx(2,8)+l22*dNx(1,8))/)                

 
 
  BNL(1,1:16:2)=dNx(1,:)
  BNL(2,1:16:2)=dNx(2,:)
  BNL(3,2:16:2)=dNx(1,:)
  BNL(4,2:16:2)=dNx(2,:)
  BNL(5,:)=BL0(3,:)  
 
      
      
  S(1:2,1:2)=es(i,1:2,1:2)
  S(3:5,3:5)=es(i,:,:)
  BL=BL0+BL1      
 
  D=D_el(i,1:4,1:4)

  Ke=Ke+(MATMUL(TRANSPOSE(BL),MATMUL(D,BL))+MATMUL(MATMUL(TRANSPOSE(BNL),S),BNL))*detJ*wp(i) *xbar1*2.D0*pi

END DO
RETURN
END




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!----------------------------------------------------------------------------------!
! NAME    : PLANI4F_PE                                                             !
!                                                                                   !
! PURPOSE : Calculate the force vector for a 4 node isoparametric                  !
!           element in plane strain, with the theory for large deformations        !
!                                                                                   !
! INPUT   : coord       Node coordinates for the element [4 x 2]                   !
!                                                                                   !
!           ed          Nodal displacment  [u1 u2 u3 u4 ...]                       !
!                                                                                   !
!           es          Second Piola-Kirchoff stress [gp x 3 x 3]                  ! 
!                                                                                   !
!           ngp         Number of integration points                               !
!                                                                                   !
!           t           Thickness of the element                                   !  
!                                                                                   !
!                                                                                   !
! OUTPUT  : fe          Element force vector  [1 x 8]                              !        
!                                                                                   !
!                                                                                   !
!----------------------------------------------------------------------------------!
!                                                                                   !
SUBROUTINE PLANI8F_AXI(fe,coord,ed,es,ngp)


IMPLICIT NONE
INTEGER                                       :: r2, ngp, i 
DOUBLE PRECISION                              :: g1, g2, w1, w2,l11,l12,l22,l21,l33,detJ,xbar1
DOUBLE PRECISION, PARAMETER                   :: pi=3.14159265358979D0
DOUBLE PRECISION, DIMENSION(ngp,8)            :: N
DOUBLE PRECISION, DIMENSION(ngp)              :: eta, xsi,wp, w_1, w_2, Phi
DOUBLE PRECISION, DIMENSION(2*ngp,8),TARGET   :: dNr
DOUBLE PRECISION, DIMENSION(4,16)             :: BL0, BL1, BL
DOUBLE PRECISION, DIMENSION(2*ngp,2)          :: JT
DOUBLE PRECISION, DIMENSION(2,8)              :: dNx
DOUBLE PRECISION, DIMENSION(8,2)              :: coord
DOUBLE PRECISION, POINTER                     :: dNr_p(:,:)
DOUBLE PRECISION, DIMENSION(2,2)              :: JTinv
DOUBLE PRECISION, DIMENSION(16)               :: ed
DOUBLE PRECISION, DIMENSION(ngp,3,3)          :: es
DOUBLE PRECISION, DIMENSION(16,1)             :: ef
DOUBLE PRECISION, DIMENSION(1,16)             :: fe


IF (ngp==1) THEN
    g1=0.0D0; w1=2.0D0
    xsi=(/g1/)
    eta=(/g1/)
    w_1=(/w1 /)
    w_2=(/w1 /)
ELSE IF (ngp==4) THEN
    g1=0.577350269189626D0; w1=1D0
    xsi=(/-g1, g1,-g1, g1/)
    eta=(/-g1, -g1,g1, g1/)
    w_1=(/w1, w1, w1, w1/)
    w_2=w_1
ELSE IF (ngp==9) THEN ! ????????????
    g1=0.774596669241483D0; g2=0.D0;
    w1=0.555555555555555D0; w2=0.888888888888888D0;
    xsi=(/-g1,-g2, g1,-g1, g2, g1,-g1, g2, g1/)
    eta=(/-g1,-g1,-g1, g2, g2, g2, g1, g1, g1/)
    w_1=(/w1, w2, w1, w1, w2, w1, w1, w2, w1/)
    w_2=(/w1, w1, w1, w2, w2, w2, w1, w1, w1/)
ELSE
    WRITE(*,*) 'Used number of integration points not implemented'
    
END IF
  wp=w_1*w_2


!--------- shape functions -----------------------------------

  N(:,1)=-(1D0-xsi)*(1D0-eta)*(1D0+xsi+eta)/4D0; N(:,5)=(1D0-xsi*xsi)*(1D0-eta)/2D0
  N(:,2)=-(1D0+xsi)*(1D0-eta)*(1D0-xsi+eta)/4D0; N(:,6)=(1D0+xsi)*(1D0-eta*eta)/2D0
  N(:,3)=-(1D0+xsi)*(1D0+eta)*(1D0-xsi-eta)/4D0; N(:,7)=(1D0-xsi*xsi)*(1D0+eta)/2D0
  N(:,4)=-(1D0-xsi)*(1D0+eta)*(1D0+xsi-eta)/4D0; N(:,8)=(1D0-xsi)*(1D0-eta*eta)/2D0

  dNr(1:2*ngp-1:2,1)=-(-(1D0-eta)*(1D0+xsi+eta)+(1D0-xsi)*(1D0-eta))/4D0
  dNr(1:2*ngp-1:2,2)=-( (1D0-eta)*(1D0-xsi+eta)-(1D0+xsi)*(1D0-eta))/4D0
  dNr(1:2*ngp-1:2,3)=-( (1D0+eta)*(1D0-xsi-eta)-(1D0+xsi)*(1D0+eta))/4D0
  dNr(1:2*ngp-1:2,4)=-(-(1D0+eta)*(1D0+xsi-eta)+(1D0-xsi)*(1D0+eta))/4D0
  dNr(1:2*ngp-1:2,5)=-xsi*(1D0-eta)
  dNr(1:2*ngp-1:2,6)=(1D0-eta*eta)/2D0
  dNr(1:2*ngp-1:2,7)=-xsi*(1D0+eta)
  dNr(1:2*ngp-1:2,8)=-(1D0-eta*eta)/2D0
  dNr(2:2*ngp:2,1)=-(-(1D0-xsi)*(1D0+xsi+eta)+(1D0-xsi)*(1D0-eta))/4D0
  dNr(2:2*ngp:2,2)=-(-(1D0+xsi)*(1D0-xsi+eta)+(1D0+xsi)*(1D0-eta))/4D0
  dNr(2:2*ngp:2,3)=-( (1D0+xsi)*(1D0-xsi-eta)-(1D0+xsi)*(1D0+eta))/4D0
  dNr(2:2*ngp:2,4)=-( (1D0-xsi)*(1D0+xsi-eta)-(1D0-xsi)*(1D0+eta))/4D0
  dNr(2:2*ngp:2,5)=-(1D0-xsi*xsi)/2D0
  dNr(2:2*ngp:2,6)=-eta*(1D0+xsi)
  dNr(2:2*ngp:2,7)=(1D0-xsi*xsi)/2D0
  dNr(2:2*ngp:2,8)=-eta*(1D0-xsi)
  
  
  JT=MATMUL(dNr,coord)


!--------- plane strain conditions -----------------------------
ef=0D0
DO i=1,ngp
 
  detJ=JT(2*i-1,1)*JT(2*i,2)-JT(2*i-1,2)*JT(2*i,1)
  IF (detJ<0) THEN
    WRITE(*,*)'Jacobideterminant equal or less than zero!'
  END IF
  JTinv=1/detJ*RESHAPE((/JT(2*i,2), -JT(2*i,1),-JT(2*i-1,2),JT(2*i-1,1)/),(/2,2/))

  dNr_p=>dNr(2*i-1:2*i,:)

  dNx=MATMUL(JTinv,dNr_P)

  xbar1=SUM(N(i,:)*coord(:,1))
 
  BL0(1,1:16:2)=dNx(1,:)
  BL0(2,2:16:2)=dNx(2,:)
  BL0(3,1:16:2)=N(i,:)/xbar1 
  BL0(4,1:16:2)=dNx(2,:)  
  BL0(4,2:16:2)=dNx(1,:)   
   
  
  l11=SUM(dNx(1,:)*ed(1:16:2))
  l22=SUM(dNx(2,:)*ed(2:16:2))
  l12=SUM(dNx(2,:)*ed(1:16:2))
  l21=SUM(dNx(1,:)*ed(2:16:2))
  l33=SUM(N(i,:)*ed(1:16:2)/xbar1)
      
  BL1(1,:)=(/dNx(1,1)*l11, dNx(1,1)*l21, dNx(1,2)*l11, dNx(1,2)*l21,&
             dNx(1,3)*l11, dNx(1,3)*l21, dNx(1,4)*l11, dNx(1,4)*l21,&
             dNx(1,5)*l11, dNx(1,5)*l21, dNx(1,6)*l11, dNx(1,6)*l21,&
             dNx(1,7)*l11, dNx(1,7)*l21, dNx(1,8)*l11, dNx(1,8)*l21/)
  
  BL1(2,:)=(/dNx(2,1)*l12, dNx(2,1)*l22, dNx(2,2)*l12, dNx(2,2)*l22,&
             dNx(2,3)*l12, dNx(2,3)*l22, dNx(2,4)*l12, dNx(2,4)*l22,&
             dNx(2,5)*l12, dNx(2,5)*l22, dNx(2,6)*l12, dNx(2,6)*l22,&
             dNx(2,7)*l12, dNx(2,7)*l22, dNx(2,8)*l12, dNx(2,8)*l22/)

  BL1(3,:)=l33*BL0(3,:); 
  
  BL1(4,:)=(/(l11*dNx(2,1)+l12*dNx(1,1)),  (l21*dNx(2,1)+l22*dNx(1,1)),&
             (l11*dNx(2,2)+l12*dNx(1,2)),  (l21*dNx(2,2)+l22*dNx(1,2)),& 
             (l11*dNx(2,3)+l12*dNx(1,3)),  (l21*dNx(2,3)+l22*dNx(1,3)),& 
             (l11*dNx(2,4)+l12*dNx(1,4)),  (l21*dNx(2,4)+l22*dNx(1,4)),&
             (l11*dNx(2,5)+l12*dNx(1,5)),  (l21*dNx(2,5)+l22*dNx(1,5)),&
             (l11*dNx(2,6)+l12*dNx(1,6)),  (l21*dNx(2,6)+l22*dNx(1,6)),& 
             (l11*dNx(2,7)+l12*dNx(1,7)),  (l21*dNx(2,7)+l22*dNx(1,7)),& 
             (l11*dNx(2,8)+l12*dNx(1,8)),  (l21*dNx(2,8)+l22*dNx(1,8))/)                

 
 
 
      
      
   BL=BL0+BL1      
   ef=ef+MATMUL(TRANSPOSE(BL),RESHAPE((/es(i,1,1),es(i,2,2),es(i,3,3),es(i,1,2)/),(/4,1/)))*detJ*wp(i)*xbar1*2.D0*pi

END DO
fe=transpose(ef)
RETURN
END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------------------!
! NAME    : PLANI8M_PE                                                             !
!                                                                                   !
! PURPOSE : Calculate the mass matrix for a 4 node isoparametric                   !
!           element in plane strain, with the theory for large deformations        !
!                                                                                   !
! INPUT   : coord       Node coordinates for the element [4 x 2]                   !
!                                                                                   !
!           rho         Density in the reference configuration                     ! 
!                                                                                   !
!           ngp         Number of integration points                               !
!                                                                                   !
!           t           Thickness of the element                                   !
!                                                                                   !
!                                                                                   !
! OUTPUT  : Me          Element mass matrix  [8 x 8]                               !        
!                                                                                   !
!                                                                                   !
!----------------------------------------------------------------------------------!
!                                                                                   !

SUBROUTINE PLANI8M_AXI(Me,coord,rho,ngp)

IMPLICIT NONE
INTEGER                                       :: ngp, i ,m
DOUBLE PRECISION                              :: g1, g2, w1, w2,detJ,rho,xbar1
DOUBLE PRECISION, PARAMETER                   :: pi=3.14159265358979D0
DOUBLE PRECISION, DIMENSION(ngp,8)            :: N
DOUBLE PRECISION, DIMENSION(16,16)            :: Me
DOUBLE PRECISION, DIMENSION(ngp)              :: eta, xsi,wp, w_1, w_2,vf
DOUBLE PRECISION, DIMENSION(2*ngp,8)          :: dNr
DOUBLE PRECISION, DIMENSION(2*ngp,2)          :: JT
DOUBLE PRECISION, DIMENSION(16,2)             :: Ntrans
DOUBLE PRECISION, DIMENSION(8,2)              :: coord
DOUBLE PRECISION, POINTER                     :: dNr_p(:,:)
DOUBLE PRECISION, DIMENSION(2,2)              :: JTinv



IF (ngp==1) THEN
    g1=0.0D0; w1=2.0D0
    xsi=(/g1/)
    eta=(/g1/)
    w_1=(/w1 /)
    w_2=(/w1 /)
ELSE IF (ngp==4) THEN
    g1=0.577350269189626D0; w1=1D0
    xsi=(/-g1, g1,-g1, g1/)
    eta=(/-g1, -g1,g1, g1/)
    w_1=(/w1, w1, w1, w1/)
    w_2=w_1
ELSE IF (ngp==9) THEN ! ????????????
    g1=0.774596669241483D0; g2=0.D0;
    w1=0.555555555555555D0; w2=0.888888888888888D0;
    xsi=(/-g1,-g2, g1,-g1, g2, g1,-g1, g2, g1/)
    eta=(/-g1,-g1,-g1, g2, g2, g2, g1, g1, g1/)
    w_1=(/w1, w2, w1, w1, w2, w1, w1, w2, w1/)
    w_2=(/w1, w1, w1, w2, w2, w2, w1, w1, w1/)
ELSE
    WRITE(*,*) 'Used number of integration points not implemented'
    
END IF
  wp=w_1*w_2


!--------- shape functions -----------------------------------
  N(:,1)=-(1D0-xsi)*(1D0-eta)*(1D0+xsi+eta)/4D0; N(:,5)=(1D0-xsi*xsi)*(1D0-eta)/2D0
  N(:,2)=-(1D0+xsi)*(1D0-eta)*(1D0-xsi+eta)/4D0; N(:,6)=(1D0+xsi)*(1D0-eta*eta)/2D0
  N(:,3)=-(1D0+xsi)*(1D0+eta)*(1D0-xsi-eta)/4D0; N(:,7)=(1D0-xsi*xsi)*(1D0+eta)/2D0
  N(:,4)=-(1D0-xsi)*(1D0+eta)*(1D0+xsi-eta)/4D0; N(:,8)=(1D0-xsi)*(1D0-eta*eta)/2D0

  dNr(1:2*ngp-1:2,1)=-(-(1D0-eta)*(1D0+xsi+eta)+(1D0-xsi)*(1D0-eta))/4D0
  dNr(1:2*ngp-1:2,2)=-( (1D0-eta)*(1D0-xsi+eta)-(1D0+xsi)*(1D0-eta))/4D0
  dNr(1:2*ngp-1:2,3)=-( (1D0+eta)*(1D0-xsi-eta)-(1D0+xsi)*(1D0+eta))/4D0
  dNr(1:2*ngp-1:2,4)=-(-(1D0+eta)*(1D0+xsi-eta)+(1D0-xsi)*(1D0+eta))/4D0
  dNr(1:2*ngp-1:2,5)=-xsi*(1D0-eta)
  dNr(1:2*ngp-1:2,6)=(1D0-eta*eta)/2D0
  dNr(1:2*ngp-1:2,7)=-xsi*(1D0+eta)
  dNr(1:2*ngp-1:2,8)=-(1D0-eta*eta)/2D0
  dNr(2:2*ngp:2,1)=-(-(1D0-xsi)*(1D0+xsi+eta)+(1D0-xsi)*(1D0-eta))/4D0
  dNr(2:2*ngp:2,2)=-(-(1D0+xsi)*(1D0-xsi+eta)+(1D0+xsi)*(1D0-eta))/4D0
  dNr(2:2*ngp:2,3)=-( (1D0+xsi)*(1D0-xsi-eta)-(1D0+xsi)*(1D0+eta))/4D0
  dNr(2:2*ngp:2,4)=-( (1D0-xsi)*(1D0+xsi-eta)-(1D0-xsi)*(1D0+eta))/4D0
  dNr(2:2*ngp:2,5)=-(1D0-xsi*xsi)/2D0
  dNr(2:2*ngp:2,6)=-eta*(1D0+xsi)
  dNr(2:2*ngp:2,7)=(1D0-xsi*xsi)/2D0
  dNr(2:2*ngp:2,8)=-eta*(1D0-xsi)

  
  
  JT=MATMUL(dNr,coord)


!--------- PLANE conditions -----------------------------
Me=0D0
DO i=1,ngp
 
  detJ=JT(2*i-1,1)*JT(2*i,2)-JT(2*i-1,2)*JT(2*i,1)
  IF (detJ<0) THEN
    WRITE(*,*)'Jacobideterminant equal or less than zero!'
  END IF
  xbar1=SUM(N(i,:)*coord(:,1))
  DO m=1,8
    Ntrans(m*2-1,1)=N(i,m)
    Ntrans(m*2,2)=N(i,m)
  END DO
  Me=Me+rho*matmul(Ntrans,transpose(Ntrans))*detJ*wp(i)*xbar1*2.D0*pi
END DO
RETURN
END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------------------!
! NAME    : DEFORMATION8I_PE                                                       !
!                                                                                   !
! PURPOSE : Calculate the deformation gradient  and Green-Lagrange strain tensor   !
!           for a 8 node isoparametric element with a mixed formulation            !
!           in plane strain                                                        ! 
!                                                                                   !
! INPUT   : Nodes      Node coordinates for all element [num nodes x 2]            !
!                                                                                   !
!           element    The element matrix   [elgrp n1 n2 n3 n4], [nelm x 5]        !
!                                                                                   !
!           ed          Nodal displacment  [u1 u2 u3 u4 ...]                       !
!                                                                                   !
!           ngp         Number of integration points                               !
!                                                                                   !
!           nelm        Number of elements                                         !
!                                                                                   !
!           ndof        Number of degrees of freedom                               !
!                                                                                   !
!                                                                                   !
!                                                                                   !
! OUTPUT  : F_out       Deformation gradient, [nelm x ngp x 3 x 3]                 !        
!                                                                                   !
!                                                                                   !
!----------------------------------------------------------------------------------!
!                                                                                   !
SUBROUTINE DEFORMATION8I_AXI(F_out,Jac,nodes,element,ed,ngp,nelm,ndof)


IMPLICIT NONE
INTEGER                                       :: r2, ngp, i ,j, nelm,ndof
DOUBLE PRECISION                              :: g1, g2, w1, w2,l11,l12,l22,l21,l33,detJ,xbar1
DOUBLE PRECISION, DIMENSION(ngp,8)            :: N
DOUBLE PRECISION, DIMENSION(ngp)              :: eta, xsi,wp, w_1, w_2, Phi
DOUBLE PRECISION, DIMENSION(2*ngp,8),TARGET   :: dNr
DOUBLE PRECISION, DIMENSION(2*ngp,2)          :: JT
DOUBLE PRECISION, DIMENSION(2,8)              :: dNx
DOUBLE PRECISION, DIMENSION(8,2)              :: coord
DOUBLE PRECISION, POINTER                     :: dNr_p(:,:)
DOUBLE PRECISION, DIMENSION(2,2)              :: JTinv
DOUBLE PRECISION, DIMENSION(nelm,16)          :: ed
DOUBLE PRECISION, DIMENSION(nelm,ngp,3,3)     :: F_out
DOUBLE PRECISION, DIMENSION(nelm,ngp)         :: Jac
DOUBLE PRECISION, DIMENSION(3,3)              :: F
INTEGER, DIMENSION(nelm,9)                    :: element
DOUBLE PRECISION, DIMENSION(ndof/2,2)         :: nodes

IF (ngp==1) THEN
    g1=0.0D0; w1=2.0D0
    xsi=(/g1/)
    eta=(/g1/)
    w_1=(/w1 /)
    w_2=(/w1 /)
ELSE IF (ngp==4) THEN
    g1=0.577350269189626D0; w1=1D0
    xsi=(/-g1, g1,-g1, g1/)
    eta=(/-g1, -g1,g1, g1/)
    w_1=(/w1, w1, w1, w1/)
    w_2=w_1
ELSE IF (ngp==9) THEN ! ????????????
    g1=0.774596669241483D0; g2=0.D0;
    w1=0.555555555555555D0; w2=0.888888888888888D0;
    xsi=(/-g1,-g2, g1,-g1, g2, g1,-g1, g2, g1/)
    eta=(/-g1,-g1,-g1, g2, g2, g2, g1, g1, g1/)
    w_1=(/w1, w2, w1, w1, w2, w1, w1, w2, w1/)
    w_2=(/w1, w1, w1, w2, w2, w2, w1, w1, w1/)
ELSE
    WRITE(*,*) 'Used number of integration points not implemented'
    
END IF
  wp=w_1*w_2


!--------- shape functions -----------------------------------
  N(:,1)=-(1D0-xsi)*(1D0-eta)*(1D0+xsi+eta)/4D0; N(:,5)=(1D0-xsi*xsi)*(1D0-eta)/2D0
  N(:,2)=-(1D0+xsi)*(1D0-eta)*(1D0-xsi+eta)/4D0; N(:,6)=(1D0+xsi)*(1D0-eta*eta)/2D0
  N(:,3)=-(1D0+xsi)*(1D0+eta)*(1D0-xsi-eta)/4D0; N(:,7)=(1D0-xsi*xsi)*(1D0+eta)/2D0
  N(:,4)=-(1D0-xsi)*(1D0+eta)*(1D0+xsi-eta)/4D0; N(:,8)=(1D0-xsi)*(1D0-eta*eta)/2D0

  dNr(1:2*ngp-1:2,1)=-(-(1D0-eta)*(1D0+xsi+eta)+(1D0-xsi)*(1D0-eta))/4D0
  dNr(1:2*ngp-1:2,2)=-( (1D0-eta)*(1D0-xsi+eta)-(1D0+xsi)*(1D0-eta))/4D0
  dNr(1:2*ngp-1:2,3)=-( (1D0+eta)*(1D0-xsi-eta)-(1D0+xsi)*(1D0+eta))/4D0
  dNr(1:2*ngp-1:2,4)=-(-(1D0+eta)*(1D0+xsi-eta)+(1D0-xsi)*(1D0+eta))/4D0
  dNr(1:2*ngp-1:2,5)=-xsi*(1D0-eta)
  dNr(1:2*ngp-1:2,6)=(1D0-eta*eta)/2D0
  dNr(1:2*ngp-1:2,7)=-xsi*(1D0+eta)
  dNr(1:2*ngp-1:2,8)=-(1D0-eta*eta)/2D0
  dNr(2:2*ngp:2,1)=-(-(1D0-xsi)*(1D0+xsi+eta)+(1D0-xsi)*(1D0-eta))/4D0
  dNr(2:2*ngp:2,2)=-(-(1D0+xsi)*(1D0-xsi+eta)+(1D0+xsi)*(1D0-eta))/4D0
  dNr(2:2*ngp:2,3)=-( (1D0+xsi)*(1D0-xsi-eta)-(1D0+xsi)*(1D0+eta))/4D0
  dNr(2:2*ngp:2,4)=-( (1D0-xsi)*(1D0+xsi-eta)-(1D0-xsi)*(1D0+eta))/4D0
  dNr(2:2*ngp:2,5)=-(1D0-xsi*xsi)/2D0
  dNr(2:2*ngp:2,6)=-eta*(1D0+xsi)
  dNr(2:2*ngp:2,7)=(1D0-xsi*xsi)/2D0
  dNr(2:2*ngp:2,8)=-eta*(1D0-xsi)

 

!--------- plane strain conditions -----------------------------
DO j=1,nelm
   coord=Nodes(element(j,2:9),:)
   JT=MATMUL(dNr,coord)

   DO i=1,ngp
 
      detJ=JT(2*i-1,1)*JT(2*i,2)-JT(2*i-1,2)*JT(2*i,1)
      IF (detJ<0) THEN
           WRITE(*,*)'Jacobideterminant equal or less than zero!'
      END IF
      JTinv=1D0/detJ*RESHAPE((/JT(2*i,2), -JT(2*i,1),-JT(2*i-1,2),JT(2*i-1,1)/),(/2,2/))

      dNr_p=>dNr(2*i-1:2*i,:)

      dNx=MATMUL(JTinv,dNr_P)
      xbar1=SUM(N(i,:)*Nodes(element(j,2:9),1))
   
  
      l11=SUM(dNx(1,:)*ed(j,1:16:2))
      l22=SUM(dNx(2,:)*ed(j,2:16:2))
      l12=SUM(dNx(2,:)*ed(j,1:16:2))
      l21=SUM(dNx(1,:)*ed(j,2:16:2))
      l33=SUM(N(i,:)*ed(j,1:16:2)/xbar1)
      
      F=RESHAPE((/1D0+l11, l12, 0D0,l21, 1D0+l22, 0D0,0D0, 0D0,  1D0+l33/),(/3,3/),ORDER=(/2,1/))
      F_out(j,i,:,:)=F
      Jac(j,i)=F(1,1)*F(2,2)*F(3,3) - F(1,1)*F(2,3)*F(3,2) &
          - F(2,1)*F(1,2)*F(3,3) + F(2,1)*F(1,3)*F(3,2) &
          + F(3,1)*F(1,2)*F(2,3) - F(3,1)*F(1,3)*F(2,2)

   END DO
END DO
RETURN
END

