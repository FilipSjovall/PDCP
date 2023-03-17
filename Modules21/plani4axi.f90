module plani4axi
!This file contain: PLANI4E_AXI, PLANI4F_AXI, elemvol, and DEFORMATION4I_AXI

CONTAINS 


!----------------------------------------------------------------------------------!
! NAME    : PLANI4E_AXI                                                            !
!										                                                     !
! PURPOSE : Calculate the stiffness matrix for a 4 node isoparametric              !
!           element in axisymmetri, with the theory for large deformations         !
!										                                                     !
! INPUT   : coord       Node coordinates for the element [4 x 2]                   !
!										                                                     !
!           D_el        The constitutive matrix for each integration               !
!                       point  [gp x 4 x 4]                                        !
!										                                                     !
!           ed          Nodal displacment  [u1 u2 u3 u4 ...]                       !
!										                                                     !
!           es          Second Piola-Kirchoff stress [gp x 3 x 3]                  ! 
!										                                                     !
!           ngp         Number of integration points                               !
!										                                                     !
!										                                                     !
! OUTPUT  : Ke          Element stiffness matrix  [8 x 8]                          !        
!										                                                     !
!										                                                     !
!----------------------------------------------------------------------------------!
!										                                                     !


SUBROUTINE PLANI4E_AXI(Ke,coord,D_el,ed,es,ngp)

IMPLICIT NONE
INTEGER                                       :: r2, ngp, i , gp
DOUBLE PRECISION                              :: g1, g2, w1, w2,l11,l12,l22,l21,l33,detJ,xbar1
DOUBLE PRECISION, PARAMETER                   :: pi=3.14159265358979D0 
DOUBLE PRECISION, DIMENSION(ngp,4)            :: N
DOUBLE PRECISION, DIMENSION(ngp)              :: eta, xsi,wp, w_1, w_2
DOUBLE PRECISION, DIMENSION(2*ngp,4),TARGET   :: dNr 
DOUBLE PRECISION, DIMENSION(8,8), intent(inout)              :: Ke
DOUBLE PRECISION, DIMENSION(4,8)              :: BL0, BL1, BL
DOUBLE PRECISION, DIMENSION(5,8)              :: BNL
DOUBLE PRECISION, DIMENSION(2*ngp,2)          :: JT
DOUBLE PRECISION, DIMENSION(2,4)              :: dNx
DOUBLE PRECISION, DIMENSION(4,2)              :: coord
DOUBLE PRECISION, DIMENSION(5,5)              :: S
DOUBLE PRECISION, POINTER                     :: D(:,:)
DOUBLE PRECISION, POINTER                     :: dNr_p(:,:)
DOUBLE PRECISION, DIMENSION(2,2)              :: JTinv
DOUBLE PRECISION, DIMENSION(ngp,4,4), TARGET  :: D_el
DOUBLE PRECISION, DIMENSION(4)                :: ex,ey
DOUBLE PRECISION, DIMENSION(8)                :: ed
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
    g1=0.774596699241483D0; g2=0.D0;
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
  N(:,1)=(1D0-xsi)*(1D0-eta)/4.D0;  N(:,2)=(1D0+xsi)*(1D0-eta)/4.D0
  N(:,3)=(1D0+xsi)*(1D0+eta)/4.D0;  N(:,4)=(1D0-xsi)*(1D0+eta)/4.D0

  dNr(1:2*ngp-1:2,1)=-(1D0-eta)/4.D0;     dNr(1:2*ngp-1:2,2)= (1D0-eta)/4.D0
  dNr(1:2*ngp-1:2,3)= (1D0+eta)/4.D0;     dNr(1:2*ngp-1:2,4)=-(1D0+eta)/4.D0
  dNr(2:2*ngp:2,1)=-(1D0-xsi)/4.D0;       dNr(2:2*ngp:2,2)=-(1D0+xsi)/4.D0
  dNr(2:2*ngp:2,3)= (1D0+xsi)/4.D0;       dNr(2:2*ngp:2,4)= (1D0-xsi)/4.D0

  
  
  JT=MATMUL(dNr,coord)


!--------- axisymmetric conditions -----------------------------
Bl0 = 0d0
Bl1 = 0d0
BL  = 0d0
BNL = 0d0
S   = 0d0
Ke=0D0
DO gp=1,ngp
 
  detJ=JT(2*gp-1,1)*JT(2*gp,2)-JT(2*gp-1,2)*JT(2*gp,1)
  IF (detJ<0) THEN
    WRITE(*,*)'Jacobideterminant equal or less than zero! (K-matrix)'
  END IF
  JTinv=1D0/detJ*RESHAPE((/JT(2*gp,2), -JT(2*gp,1),-JT(2*gp-1,2),JT(2*gp-1,1)/),(/2,2/))

  dNr_p=>dNr(2*gp-1:2*gp,:)

  dNx=MATMUL(JTinv,dNr_P)

  xbar1=SUM(N(gp,:)*coord(:,1))
  BL0(1,1:8:2)=dNx(1,:)
  BL0(2,2:8:2)=dNx(2,:)
  BL0(4,1:8:2)=N(gp,:)/xbar1 
  BL0(3,1:8:2)=dNx(2,:)  
  BL0(3,2:8:2)=dNx(1,:)	
	
  
  l11=SUM(dNx(1,:)*ed(1:8:2))
  l22=SUM(dNx(2,:)*ed(2:8:2))
  l12=SUM(dNx(2,:)*ed(1:8:2))
  l21=SUM(dNx(1,:)*ed(2:8:2))
  l33=SUM(N(gp,:)*ed(1:8:2)/xbar1)
		
  BL1(1,:)=(/dNx(1,1)*l11, dNx(1,1)*l21, dNx(1,2)*l11, dNx(1,2)*l21,&
            dNx(1,3)*l11, dNx(1,3)*l21, dNx(1,4)*l11, dNx(1,4)*l21/)
		
  BL1(2,:)=(/dNx(2,1)*l12, dNx(2,1)*l22, dNx(2,2)*l12, dNx(2,2)*l22,&
            dNx(2,3)*l12, dNx(2,3)*l22, dNx(2,4)*l12, dNx(2,4)*l22/)
  BL1(4,:)=l33*BL0(4,:); 

  BL1(3,:)=(/(l11*dNx(2,1)+l12*dNx(1,1)),  (l21*dNx(2,1)+l22*dNx(1,1)),&
            (l11*dNx(2,2)+l12*dNx(1,2)),  (l21*dNx(2,2)+l22*dNx(1,2)),& 
            (l11*dNx(2,3)+l12*dNx(1,3)),  (l21*dNx(2,3)+l22*dNx(1,3)),& 
            (l11*dNx(2,4)+l12*dNx(1,4)),  (l21*dNx(2,4)+l22*dNx(1,4))/)                

 
 
  BNL(1,1:8:2)=dNx(1,:)
  BNL(2,1:8:2)=dNx(2,:)
  BNL(3,2:8:2)=dNx(1,:)
  BNL(4,2:8:2)=dNx(2,:)
  BNL(5,:)=BL0(4,:)  

		
		
  S(1:2,1:2)=es(gp,1:2,1:2)
  S(3:5,3:5)=es(gp,:,:)
  BL=BL0+BL1
 
  !D=>D_el(i,1:4,1:4)
 	D=>D_el(1:4,1:4,gp)

  Ke=Ke+(MATMUL(TRANSPOSE(BL),MATMUL(D,BL)) + MATMUL(MATMUL(TRANSPOSE(BNL),S),BNL))*detJ*wp(gp)*xbar1*2.D0*pi 
  
  !print * , "Ke i plani4axi - loop" , Ke
  !if(gp.eq.1) then
  !  print * , "D" , D
  !  print * , "S" , S
  !  print * , "BL", BL
  !  print * , "BNL", BNL
  !  print * , "Ke i plani4axi" , Ke
  !endif

END DO


END subroutine PLANI4E_AXI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!----------------------------------------------------------------------------------!
! NAME    : PLANI4F_AXI                                                            !
!										                                                     !
! PURPOSE : Calculate the force vector for a 4 node isoparametric                  !
!           element in axisymmetri, with the theory for large deformations         !
!										                                                     !
! INPUT   : coord       Node coordinates for the element [4 x 2]                   !
!										                                                     !
!           ed          Nodal displacment  [u1 u2 u3 u4 ...]                       !
!										                                                     !
!           es          Second Piola-Kirchoff stress [gp x 3 x 3]                  ! 
!										                                                     !
!           ngp         Number of integration points                               !
!										                                                     !
!										                                                     !
! OUTPUT  : fe          Element foce vector  [1 x 8]                               !        
!										                                                     !
!										                                                     !
!----------------------------------------------------------------------------------!
!										                                                     !


SUBROUTINE PLANI4F_AXI(fe,coord,ed,es,ngp)


IMPLICIT NONE
INTEGER                                       :: r2, ngp, i , gp
DOUBLE PRECISION                              :: g1, g2, w1, w2,l11,l12,l22,l21,l33,detJ,xbar1
DOUBLE PRECISION, PARAMETER                   :: pi=3.14159265358979D0
DOUBLE PRECISION, DIMENSION(ngp,4)            :: N
DOUBLE PRECISION, DIMENSION(ngp)              :: eta, xsi,wp, w_1, w_2, Phi
DOUBLE PRECISION, DIMENSION(2*ngp,4),TARGET   :: dNr
DOUBLE PRECISION, DIMENSION(4,8)              :: BL0, BL1, BL
DOUBLE PRECISION, DIMENSION(5,8)              :: BNL
DOUBLE PRECISION, DIMENSION(2*ngp,2)          :: JT
DOUBLE PRECISION, DIMENSION(2,4)              :: dNx
DOUBLE PRECISION, DIMENSION(4,2)              :: coord
DOUBLE PRECISION, DIMENSION(5,5)              :: S
DOUBLE PRECISION, POINTER                     :: D(:,:),Cp(:,:)
DOUBLE PRECISION, POINTER                     :: dNr_p(:,:)
DOUBLE PRECISION, DIMENSION(2,2)              :: JTinv
DOUBLE PRECISION, DIMENSION(8)                :: ed
DOUBLE PRECISION, DIMENSION(ngp,3,3)          :: es
DOUBLE PRECISION, DIMENSION(8,1)              :: ef
DOUBLE PRECISION, DIMENSION(1,8)              :: fe

JT = 0D0

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
    g1=0.774596699241483D0; g2=0.D0;
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
  N(:,1)=(1D0-xsi)*(1D0-eta)/4.D0;  N(:,2)=(1D0+xsi)*(1D0-eta)/4.D0
  N(:,3)=(1D0+xsi)*(1D0+eta)/4.D0;  N(:,4)=(1D0-xsi)*(1D0+eta)/4.D0

  dNr(1:2*ngp-1:2,1)=-(1D0-eta)/4.D0;     dNr(1:2*ngp-1:2,2)= (1D0-eta)/4.D0
  dNr(1:2*ngp-1:2,3)= (1D0+eta)/4.D0;     dNr(1:2*ngp-1:2,4)=-(1D0+eta)/4.D0
  dNr(2:2*ngp:2,1)=-(1D0-xsi)/4.D0;       dNr(2:2*ngp:2,2)=-(1D0+xsi)/4.D0
  dNr(2:2*ngp:2,3)= (1D0+xsi)/4.D0;       dNr(2:2*ngp:2,4)= (1D0-xsi)/4.D0

  
  
  JT=MATMUL(dNr,coord)


!--------- axisymmetric conditions -----------------------------
ef=0D0
Bl0 = 0d0
Bl1 = 0d0
BL  = 0d0
BNL = 0d0
S   = 0d0
DO gp=1,ngp
 
  detJ=JT(2*gp-1,1)*JT(2*gp,2)-JT(2*gp-1,2)*JT(2*gp,1)
  IF (detJ<0) THEN
    WRITE(*,*)'Jacobideterminant equal or less than zero! (fe routine)'
  END IF
  JTinv=1/detJ*RESHAPE((/JT(2*gp,2), -JT(2*gp,1),-JT(2*gp-1,2),JT(2*gp-1,1)/),(/2,2/))

  dNr_p=>dNr(2*gp-1:2*gp,:)

  dNx=MATMUL(JTinv,dNr_P)

	xbar1=SUM(N(gp,:)*coord(:,1))

  BL0(1,1:8:2)=dNx(1,:)
  BL0(2,2:8:2)=dNx(2,:)
  BL0(4,1:8:2)=N(gp,:)/xbar1 
  BL0(3,1:8:2)=dNx(2,:)  
  BL0(3,2:8:2)=dNx(1,:)	
	
  
  l11=SUM(dNx(1,:)*ed(1:8:2))
  l22=SUM(dNx(2,:)*ed(2:8:2))
  l12=SUM(dNx(2,:)*ed(1:8:2))
  l21=SUM(dNx(1,:)*ed(2:8:2))
  l33=SUM(N(gp,:)*ed(1:8:2)/xbar1)
		
  BL1(1,:)=(/dNx(1,1)*l11, dNx(1,1)*l21, dNx(1,2)*l11, dNx(1,2)*l21,&
             dNx(1,3)*l11, dNx(1,3)*l21, dNx(1,4)*l11, dNx(1,4)*l21/)
		
  BL1(2,:)=(/dNx(2,1)*l12, dNx(2,1)*l22, dNx(2,2)*l12, dNx(2,2)*l22,&
             dNx(2,3)*l12, dNx(2,3)*l22, dNx(2,4)*l12, dNx(2,4)*l22/)
  BL1(4,:)=l33*BL0(4,:); 

  BL1(3,:)=(/(l11*dNx(2,1)+l12*dNx(1,1)),  (l21*dNx(2,1)+l22*dNx(1,1)),&
             (l11*dNx(2,2)+l12*dNx(1,2)),  (l21*dNx(2,2)+l22*dNx(1,2)),& 
             (l11*dNx(2,3)+l12*dNx(1,3)),  (l21*dNx(2,3)+l22*dNx(1,3)),& 
             (l11*dNx(2,4)+l12*dNx(1,4)),  (l21*dNx(2,4)+l22*dNx(1,4))/)                
 
  BNL(1,1:8:2)=dNx(1,:)
  BNL(2,1:8:2)=dNx(2,:)
  BNL(3,2:8:2)=dNx(1,:)
  BNL(4,2:8:2)=dNx(2,:)
  BNL(5,:)=BL0(4,:)  

		
		
  S(1:2,1:2)=es(gp,1:2,1:2)
  S(3:5,3:5)=es(gp,:,:)
		
  BL=BL0+BL1

  ! Changed to fit another format
  !ef=ef+MATMUL(TRANSPOSE(BL),RESHAPE((/es(gp,1,1),es(gp,2,2),es(gp,3,3),es(gp,1,2)/),(/4,1/)))*detJ*wp(gp)*xbar1*2.D0*pi

	ef=ef+MATMUL(TRANSPOSE(BL),RESHAPE((/es(gp,1,1),es(gp,2,2),es(gp,1,2),es(gp,3,3)/),(/4,1/)))*detJ*wp(gp)*xbar1*2.D0*pi

END DO
fe=transpose(ef)
RETURN
END



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------------------!
! NAME    : PLANI4M_PE                                                             !
!										                                                     !
! PURPOSE : Calculate the mass matrix for a 4 node isoparametric                   !
!           element in plane strain, with the theory for large deformations        !
!										                                                     !
! INPUT   : coord       Node coordinates for the element [4 x 2]                   !
!										                                                     !
!           rho         Density in the reference configuration                     ! 
!										                                                     !
!           ngp         Number of integration points                               !
!										                                                     !
!           t           Thickness of the element                                   !
!										                                                     !
!										                                                     !
! OUTPUT  : Me          Element mass matrix  [8 x 8]                               !        
!										                                                     !
!										                                                     !
!----------------------------------------------------------------------------------!
!										                                                     !

SUBROUTINE PLANI4M_AXI(Me,coord,rho,ngp)

IMPLICIT NONE
INTEGER                                       :: ngp, i ,m
DOUBLE PRECISION                              :: g1, g2, w1, w2,detJ,rho,xbar1
DOUBLE PRECISION, PARAMETER                   :: pi=3.14159265358979D0
DOUBLE PRECISION, DIMENSION(ngp,4)            :: N
DOUBLE PRECISION, DIMENSION(8,8)              :: Me
DOUBLE PRECISION, DIMENSION(ngp)              :: eta, xsi,wp, w_1, w_2,vf
DOUBLE PRECISION, DIMENSION(2*ngp,4)          :: dNr
DOUBLE PRECISION, DIMENSION(2*ngp,2)          :: JT
DOUBLE PRECISION, DIMENSION(8,2)              :: Ntrans
DOUBLE PRECISION, DIMENSION(4,2)              :: coord
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
    g1=0.774596699241483D0; g2=0.D0;
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
  N(:,1)=(1D0-xsi)*(1D0-eta)/4.D0;  N(:,2)=(1D0+xsi)*(1D0-eta)/4.D0
  N(:,3)=(1D0+xsi)*(1D0+eta)/4.D0;  N(:,4)=(1D0-xsi)*(1D0+eta)/4.D0

  dNr(1:2*ngp-1:2,1)=-(1D0-eta)/4.D0;     dNr(1:2*ngp-1:2,2)= (1D0-eta)/4.D0
  dNr(1:2*ngp-1:2,3)= (1D0+eta)/4.D0;     dNr(1:2*ngp-1:2,4)=-(1D0+eta)/4.D0
  dNr(2:2*ngp:2,1)=-(1D0-xsi)/4.D0;       dNr(2:2*ngp:2,2)=-(1D0+xsi)/4.D0
  dNr(2:2*ngp:2,3)= (1D0+xsi)/4.D0;       dNr(2:2*ngp:2,4)= (1D0-xsi)/4.D0

  
  
  JT=MATMUL(dNr,coord)


!--------- PLANE conditions -----------------------------
Me=0D0
DO i=1,ngp
 
  detJ=JT(2*i-1,1)*JT(2*i,2)-JT(2*i-1,2)*JT(2*i,1)
  IF (detJ<0) THEN
    WRITE(*,*)'Jacobideterminant equal or less than zero!'
  END IF
  xbar1=SUM(N(i,:)*coord(:,1))
  DO m=1,4
    Ntrans(m*2-1,1)=N(i,m)
    Ntrans(m*2,2)=N(i,m)
  END DO
  Me=Me+rho*matmul(Ntrans,transpose(Ntrans))*detJ*wp(i)*xbar1*2.D0*pi
END DO
RETURN
END






!----------------------------------------------------------------------------------!
! NAME    : DEFORMATION4I_AXI                                                      !
!										                                                     !
! PURPOSE : Calculate the deformation gradient  and Green-Lagrange strain tensor   !
!           for a 4 node isoparametric element with a mixed formulation            !
!           in axisymmetri,                                                        ! 
!										                                                     !
! INPUT   : Nodes      Node coordinates for all element [num nodes x 2]            !
!										                                                     !
!           element    The element matrix   [elgrp n1 n2 n3 n4], [nelm x 5]        !
!										                                                     !
!           ed          Nodal displacment  [u1 u2 u3 u4 ...]                       !
!										                                                     !
!           ngp         Number of integration points                               !
!										                                                     !
!           nelm        Number of elements                                         !
!										                                                     !
!           ndof        Number of degrees of freedom                               !
!										                                                     !
!										                                                     !
!										                                                     !
! OUTPUT  : F_out       Deformation gradient, [nelm x ngp x 3 x 3]                 !        
!										                                                     !
!           E_out       Green-Lagrange strain tensor,[nelm x ngp x 3 x 3]	        !
!										                                                     !
!	          J_out       Jacobian   [nelm x ngp]                                    !
!                                                                                  !
!										                                                     !
!----------------------------------------------------------------------------------!
!										                                                     !

SUBROUTINE DEFORMATION4I_AXI(F_out,J_out,nodes,element,ed,ngp,nelm,ndof)


IMPLICIT NONE
INTEGER                                       :: r2, ngp, i ,j, nelm,ndof
DOUBLE PRECISION                              :: g1, g2, w1, w2,l11,l12,l22,l21,l33,detJ,xbar1
DOUBLE PRECISION, DIMENSION(ngp,4)            :: N
DOUBLE PRECISION, DIMENSION(ngp)              :: eta, xsi,wp, w_1, w_2, Phi
DOUBLE PRECISION, DIMENSION(2*ngp,4),TARGET   :: dNr
DOUBLE PRECISION, DIMENSION(2*ngp,2)          :: JT
DOUBLE PRECISION, DIMENSION(2,4)              :: dNx
DOUBLE PRECISION, DIMENSION(4,2)              :: coord
DOUBLE PRECISION, POINTER                     :: dNr_p(:,:)
DOUBLE PRECISION, DIMENSION(2,2)              :: JTinv
DOUBLE PRECISION, DIMENSION(nelm,8)           :: ed
DOUBLE PRECISION, DIMENSION(nelm,ngp,3,3)     :: F_out,E_out
DOUBLE PRECISION, DIMENSION(nelm,ngp)         :: J_out
DOUBLE PRECISION, DIMENSION(3,3)              :: F
INTEGER, DIMENSION(nelm,5)                    :: element
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
    g1=0.774596699241483D0; g2=0.D0;
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
  N(:,1)=(1D0-xsi)*(1D0-eta)/4.D0;  N(:,2)=(1D0+xsi)*(1D0-eta)/4.D0
  N(:,3)=(1D0+xsi)*(1D0+eta)/4.D0;  N(:,4)=(1D0-xsi)*(1D0+eta)/4.D0

  dNr(1:2*ngp-1:2,1)=-(1D0-eta)/4.D0;     dNr(1:2*ngp-1:2,2)= (1D0-eta)/4.D0
  dNr(1:2*ngp-1:2,3)= (1D0+eta)/4.D0;     dNr(1:2*ngp-1:2,4)=-(1D0+eta)/4.D0
  dNr(2:2*ngp:2,1)=-(1D0-xsi)/4.D0;       dNr(2:2*ngp:2,2)=-(1D0+xsi)/4.D0
  dNr(2:2*ngp:2,3)= (1D0+xsi)/4.D0;       dNr(2:2*ngp:2,4)= (1D0-xsi)/4.D0

 

!--------- axisymmetric conditions -----------------------------
DO j=1,nelm
	 
   coord=Nodes(element(j,2:5),:)
	 JT=MATMUL(dNr,coord)

   DO i=1,ngp
 
      detJ=JT(2*i-1,1)*JT(2*i,2)-JT(2*i-1,2)*JT(2*i,1)
      IF (detJ<0) THEN
      	  WRITE(*,*)'Jacobideterminant equal or less than zero! (deformations)'
					!write(*,*) 'i,j,JT:' ,i,j,JT
      END IF
      JTinv=1D0/detJ*RESHAPE((/JT(2*i,2), -JT(2*i,1),-JT(2*i-1,2),JT(2*i-1,1)/),(/2,2/))

      dNr_p=>dNr(2*i-1:2*i,:)

      dNx=MATMUL(JTinv,dNr_P)
      xbar1=SUM(N(i,:)*Nodes(element(j,2:5),1))
	
  
      l11=SUM(dNx(1,:)*ed(j,1:8:2))
      l22=SUM(dNx(2,:)*ed(j,2:8:2))
      l12=SUM(dNx(2,:)*ed(j,1:8:2))
      l21=SUM(dNx(1,:)*ed(j,2:8:2))
      l33=SUM(N(i,:)*ed(j,1:8:2)/xbar1)
	
      F=RESHAPE((/1D0+l11, l12, 0D0, l21, 1D0+l22, 0D0, 0D0, 0D0,  1D0+l33/),(/3,3/),ORDER=(/2,1/))  

			F_out(j,i,:,:)=F
      
      ! Lagrangian deformation tensor
      !E_out(j,i,:,:)=1D0/2D0*(MATMUL(TRANSPOSE(F),F)-RESHAPE((/1D0,0D0,0D0,0D0,1D0,0D0,0D0,0D0,1D0/),(/3,3/)))
   	
      J_out(j,i)=F(1,1)*F(2,2)*F(3,3) - F(1,1)*F(2,3)*F(3,2)&
	       - F(2,1)*F(1,2)*F(3,3) + F(2,1)*F(1,3)*F(3,2) &
			 + F(3,1)*F(1,2)*F(2,3) - F(3,1)*F(1,3)*F(2,2)

   
   
   
   
   END DO
END DO
RETURN
END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine elemvol(Ve,coord,ir)

IMPLICIT NONE

INTEGER                                       :: ir, i, m,ngp
!INTEGER, PARAMTER  														:: ngp=ir*ir
DOUBLE PRECISION                              :: g1, g2, w1, w2,detJ,rho,xbar1,Ve,Ai
DOUBLE PRECISION, PARAMETER                   :: pi=3.141592653589793D0
DOUBLE PRECISION, DIMENSION(ir*ir,4)          :: N
DOUBLE PRECISION, DIMENSION(ir*ir)            :: eta, xsi,wp, w_1, w_2,vf
DOUBLE PRECISION, DIMENSION(2*ir**2,4)        :: dNr
DOUBLE PRECISION, DIMENSION(2*ir**2,2)        :: JT
DOUBLE PRECISION, DIMENSION(8,2)              :: Ntrans
DOUBLE PRECISION, DIMENSION(4,2)              :: coord
DOUBLE PRECISION, POINTER                     :: dNr_p(:,:)
DOUBLE PRECISION, DIMENSION(2,2)              :: JTinv

Ve  = 0d0
JT  = 0d0
ngp = ir*ir

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
    g1=0.774596699241483D0; g2=0.D0;
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
  N(:,1)=(1D0-xsi)*(1D0-eta)/4.D0;  N(:,2)=(1D0+xsi)*(1D0-eta)/4.D0
  N(:,3)=(1D0+xsi)*(1D0+eta)/4.D0;  N(:,4)=(1D0-xsi)*(1D0+eta)/4.D0

  dNr(1:2*ngp-1:2,1)=-(1D0-eta)/4.D0;     dNr(1:2*ngp-1:2,2)= (1D0-eta)/4.D0
  dNr(1:2*ngp-1:2,3)= (1D0+eta)/4.D0;     dNr(1:2*ngp-1:2,4)=-(1D0+eta)/4.D0
  dNr(2:2*ngp:2,1)=-(1D0-xsi)/4.D0;       dNr(2:2*ngp:2,2)=-(1D0+xsi)/4.D0
  dNr(2:2*ngp:2,3)= (1D0+xsi)/4.D0;       dNr(2:2*ngp:2,4)= (1D0-xsi)/4.D0

  
  
  JT=MATMUL(dNr,coord)
	DO i=1,ngp
  	detJ=JT(2*i-1,1)*JT(2*i,2)-JT(2*i-1,2)*JT(2*i,1)
  	IF (detJ<0) THEN
    	WRITE(*,*)'Jacobideterminant equal or less than zero!'
  	END IF
  	xbar1=SUM(N(i,:)*coord(:,1))
		!indx = (/ 2*i-1, 2*i /) 
		Ve   = Ve + detJ*wp(i)*xbar1*2d0*pi
	END DO

end subroutine elemvol


end module plani4axi
