module elem_large_cont_2d

! Last rev.:
!  M. Ristinmaa 2011-04-30
!   - Code from matlab implemented c2dxx3	
!  M. Ristinmaa 2011-05-03
!   - Coordinate array transposed, column is now node coordinates	
!------------------------------------------------------------------------------
! A. Dalklint 2020-04-28
!   - Implemented plane strain 4 node total lagragian elements
!
!   Module elem_large_cont_2d contains element subroutines 
!   large deformations.	



use matrix_util, only: inv2, det2, inv3

implicit none


double precision  xsi(4), eta(4), G1, Ifm(4)

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
contains


! 3-node element based on updated Lagrangian formulation
subroutine c2dul3_e(ke,coord,ep,D,ed,es)
  implicit none
  double precision                :: ke(:,:), coord(:,:), ep
  double precision                :: D(:,:), ed(:), es(:)
  
  double precision                :: ex(3), ey(3), Area, idetA, detA, t
  double precision                :: B(3,6), H(4,6), S(4,4) 
  double precision                :: dN1dx, dN2dx, dN3dx, dN1dy, dN2dy, dN3dy
  
  ex=coord(1,:)+ed(1:5:2)
  ey=coord(2,:)+ed(2:6:2)
  
  detA=ex(2)*ey(3)+ex(1)*ey(2)+ey(1)*ex(3)-ex(2)*ey(1)-ex(3)*ey(2)-ey(3)*ex(1)
  Area=detA*0.5D0
  idetA=1d0/detA

  t=ep
  
  dN1dx=idetA*(ey(2)-ey(3))
  dN2dx=idetA*(ey(3)-ey(1))
  dN3dx=idetA*(ey(1)-ey(2))

  dN1dy=idetA*(ex(3)-ex(2))
  dN2dy=idetA*(ex(1)-ex(3))
  dN3dy=idetA*(ex(2)-ex(1))

  B(1,:)=(/dN1dx,   0d0, dN2dx,   0d0, dN3dx,   0d0/)
  B(2,:)=(/  0d0, dN1dy,   0d0, dN2dy,   0d0, dN3dy/)  
  B(3,:)=(/dN1dy, dN1dx, dN2dy, dN2dx, dN3dy, dN3dx/)

  H(1,:)=(/dN1dx,   0d0, dN2dx,   0d0, dN3dx,   0d0/)
  H(2,:)=(/dN1dy,   0d0, dN2dy,   0d0, dN3dy,   0d0/)
  H(3,:)=(/  0d0, dN1dx,   0d0, dN2dx,   0d0, dN3dx/)
  H(4,:)=(/  0d0, dN1dy,   0d0, dN2dy,   0d0, dN3dy/)

  S(1,:)=(/es(1), es(3),   0d0,   0d0/)
  S(2,:)=(/es(3), es(2),   0d0,   0d0/)
  S(3,:)=(/  0d0,   0d0, es(1), es(3)/)
  S(4,:)=(/  0d0,   0d0, es(3), es(2)/)

  Ke=matmul(transpose(B),matmul(D,B))+matmul(transpose(H),matmul(S,H))
  Ke=Ke*Area*t

end subroutine c2dul3_e


! calculates the deformation gradient
subroutine c2dul3_d(F,coord,ed)
  implicit none
  double precision                :: f(:), coord(:,:)
  double precision                :: ed(:)
  
  double precision                :: ex(3), ey(3)
  double precision                :: H(4,6), dudx(4), detA, Area, idetA
  double precision                :: dN1dx, dN2dx, dN3dx, dN1dy, dN2dy, dN3dy
  
  ex=coord(1,:)
  ey=coord(2,:)
  
  detA=ex(2)*ey(3)+ex(1)*ey(2)+ey(1)*ex(3)-ex(2)*ey(1)-ex(3)*ey(2)-ey(3)*ex(1)
  Area=detA*0.5D0
  idetA=1d0/detA
  
  dN1dx=idetA*(ey(2)-ey(3))
  dN2dx=idetA*(ey(3)-ey(1))
  dN3dx=idetA*(ey(1)-ey(2))

  dN1dy=idetA*(ex(3)-ex(2))
  dN2dy=idetA*(ex(1)-ex(3))
  dN3dy=idetA*(ex(2)-ex(1))

  H(1,:)=(/dN1dx,   0d0, dN2dx,   0d0, dN3dx,   0d0/)
  H(2,:)=(/dN1dy,   0d0, dN2dy,   0d0, dN3dy,   0d0/)
  H(3,:)=(/  0d0, dN1dx,   0d0, dN2dx,   0d0, dN3dx/)
  H(4,:)=(/  0d0, dN1dy,   0d0, dN2dy,   0d0, dN3dy/)

  dudx=matmul(H,ed)
! du_1/dx_1, du_1/dx_2, du_2/dx_1, du_2/dx_2 

  
   f(1)=dudx(1)+1d0
   f(2)=dudx(2)
   f(3)=dudx(3)
   f(4)=dudx(4)+1d0
	
end subroutine c2dul3_d


subroutine c2dul3_f(fe,coord,ep,ed,es)
  implicit none
  double precision                :: fe(:), coord(:,:), ep
  double precision                :: ed(:), es(:)
  
  double precision                :: ex(3), ey(3), Area, idetA, detA, t
  double precision                :: B(3,6)
  double precision                :: dN1dx, dN2dx, dN3dx, dN1dy, dN2dy, dN3dy
  
  ex=coord(1,:)+ed(1:5:2)
  ey=coord(2,:)+ed(2:6:2)
  
  detA=ex(2)*ey(3)+ex(1)*ey(2)+ey(1)*ex(3)-ex(2)*ey(1)-ex(3)*ey(2)-ey(3)*ex(1)
  Area=detA*0.5D0
  idetA=1d0/detA

  t=ep
  
  dN1dx=idetA*(ey(2)-ey(3))
  dN2dx=idetA*(ey(3)-ey(1))
  dN3dx=idetA*(ey(1)-ey(2))

  dN1dy=idetA*(ex(3)-ex(2))
  dN2dy=idetA*(ex(1)-ex(3))
  dN3dy=idetA*(ex(2)-ex(1))


  B(1,:)=(/dN1dx,   0d0, dN2dx,   0d0, dN3dx,   0d0/)
  B(2,:)=(/  0d0, dN1dy,   0d0, dN2dy,   0d0, dN3dy/)  
  B(3,:)=(/dN1dy, dN1dx, dN2dy, dN2dx, dN3dy, dN3dx/)

  fe=matmul(transpose(B),es)*Area*t

end subroutine c2dul3_f


! 3-node element based on total Lagrangian formulation	
subroutine c2dtl3_e(ke,coord,ep,D,ed,es)
  implicit none
  double precision                :: ke(:,:), coord(:,:), ep
  double precision                :: D(3,3), ed(:), es(:)
  
  double precision                :: ex(3), ey(3), Area, idetA, detA, t
  double precision                :: B(3,6), H(4,6), S(4,4), A(3,4), A_tmp(4)
  double precision                :: dN1dx, dN2dx, dN3dx, dN1dy, dN2dy, dN3dy
  
  ex=coord(1,:)
  ey=coord(2,:)
  
  detA=ex(2)*ey(3)+ex(1)*ey(2)+ey(1)*ex(3)-ex(2)*ey(1)-ex(3)*ey(2)-ey(3)*ex(1)
  Area=detA*0.5D0
  idetA=1d0/detA

  t=ep
  
  dN1dx=idetA*(ey(2)-ey(3))
  dN2dx=idetA*(ey(3)-ey(1))
  dN3dx=idetA*(ey(1)-ey(2))

  dN1dy=idetA*(ex(3)-ex(2))
  dN2dy=idetA*(ex(1)-ex(3))
  dN3dy=idetA*(ex(2)-ex(1))


  B(1,:)=(/dN1dx,   0d0, dN2dx,   0d0, dN3dx,   0d0/)
  B(2,:)=(/  0d0, dN1dy,   0d0, dN2dy,   0d0, dN3dy/)  
  B(3,:)=(/dN1dy, dN1dx, dN2dy, dN2dx, dN3dy, dN3dx/)

  H(1,:)=(/dN1dx,   0d0, dN2dx,   0d0, dN3dx,   0d0/)
  H(2,:)=(/dN1dy,   0d0, dN2dy,   0d0, dN3dy,   0d0/)
  H(3,:)=(/  0d0, dN1dx,   0d0, dN2dx,   0d0, dN3dx/)
  H(4,:)=(/  0d0, dN1dy,   0d0, dN2dy,   0d0, dN3dy/)

  A_tmp=matmul(H,ed)

  A(1,:)=(/A_tmp(1),      0d0, A_tmp(3),      0d0/) 
  A(2,:)=(/     0d0, A_tmp(2),      0d0, A_tmp(4)/)
  A(3,:)=(/A_tmp(2), A_tmp(1), A_tmp(4), A_tmp(3)/)

  B=B+matmul(A,H)

  S(1,:)=(/es(1), es(3),   0d0,   0d0/)
  S(2,:)=(/es(3), es(2),   0d0,   0d0/)
  S(3,:)=(/  0d0,   0d0, es(1), es(3)/)
  S(4,:)=(/  0d0,   0d0, es(3), es(2)/)

  Ke=matmul(transpose(B),matmul(D,B))+ matmul(transpose(H),matmul(S,H))
  Ke=Ke*Area*t

end subroutine c2dtl3_e


! calculates the total deformation gradient
subroutine c2dtl3_d(F,coord,ed)
  implicit none
  double precision                :: f(:), coord(:,:)
  double precision                :: ed(:)
  
  double precision                :: ex(3), ey(3)
  double precision                :: H(4,6), dudx(4), detA, Area, idetA
  double precision                :: dN1dx, dN2dx, dN3dx, dN1dy, dN2dy, dN3dy
  
  ex=coord(1,:)
  ey=coord(2,:)
  
  detA=ex(2)*ey(3)+ex(1)*ey(2)+ey(1)*ex(3)-ex(2)*ey(1)-ex(3)*ey(2)-ey(3)*ex(1)
  Area=detA*0.5D0
  idetA=1d0/detA
  
  dN1dx=idetA*(ey(2)-ey(3))
  dN2dx=idetA*(ey(3)-ey(1))
  dN3dx=idetA*(ey(1)-ey(2))

  dN1dy=idetA*(ex(3)-ex(2))
  dN2dy=idetA*(ex(1)-ex(3))
  dN3dy=idetA*(ex(2)-ex(1))

  H(1,:)=(/dN1dx,   0d0, dN2dx,   0d0, dN3dx,   0d0/)
  H(2,:)=(/dN1dy,   0d0, dN2dy,   0d0, dN3dy,   0d0/)
  H(3,:)=(/  0d0, dN1dx,   0d0, dN2dx,   0d0, dN3dx/)
  H(4,:)=(/  0d0, dN1dy,   0d0, dN2dy,   0d0, dN3dy/)

  dudx=matmul(H,ed)
! du_1/dx_1, du_1/dx_2, du_2/dx_1, du_2/dx_2 

  
   f(1)=dudx(1)+1d0
   f(2)=dudx(2)
   f(3)=dudx(3)
   f(4)=dudx(4)+1d0
	
end subroutine c2dtl3_d

subroutine c2dtl3_f(fe,coord,ep,ed,es)
  implicit none
  double precision                :: fe(:), coord(:,:), ep
  double precision                :: ed(:), es(:)
  
  double precision                :: ex(3), ey(3), Area, idetA, detA, t
  double precision                :: B(3,6), H(4,6), S(4,4), A(3,4), A_tmp(4)
  double precision                :: dN1dx, dN2dx, dN3dx, dN1dy, dN2dy, dN3dy
  
  ex=coord(1,:)
  ey=coord(2,:)
  
  detA=ex(2)*ey(3)+ex(1)*ey(2)+ey(1)*ex(3)-ex(2)*ey(1)-ex(3)*ey(2)-ey(3)*ex(1)
  Area=detA*0.5D0
  idetA=1d0/detA

  t=ep
  
  dN1dx=idetA*(ey(2)-ey(3))
  dN2dx=idetA*(ey(3)-ey(1))
  dN3dx=idetA*(ey(1)-ey(2))

  dN1dy=idetA*(ex(3)-ex(2))
  dN2dy=idetA*(ex(1)-ex(3))
  dN3dy=idetA*(ex(2)-ex(1))


  B(1,:)=(/dN1dx,   0d0, dN2dx,   0d0, dN3dx,   0d0/)
  B(2,:)=(/  0d0, dN1dy,   0d0, dN2dy,   0d0, dN3dy/)  
  B(3,:)=(/dN1dy, dN1dx, dN2dy, dN2dx, dN3dy, dN3dx/)

  H(1,:)=(/dN1dx,   0d0, dN2dx,   0d0, dN3dx,   0d0/)
  H(2,:)=(/dN1dy,   0d0, dN2dy,   0d0, dN3dy,   0d0/)
  H(3,:)=(/  0d0, dN1dx,   0d0, dN2dx,   0d0, dN3dx/)
  H(4,:)=(/  0d0, dN1dy,   0d0, dN2dy,   0d0, dN3dy/)

  A_tmp=matmul(H,ed)

  A(1,:)=(/A_tmp(1),      0d0, A_tmp(3),      0d0/) 
  A(2,:)=(/     0d0, A_tmp(2),      0d0, A_tmp(4)/)
  A(3,:)=(/A_tmp(2), A_tmp(1), A_tmp(4), A_tmp(3)/)

  B=B+matmul(A,H)

  fe=matmul(transpose(B),es)*Area*t

end subroutine c2dtl3_f











! 4-node element based on total Lagrangian formulation	
! Only plane strain implemented
subroutine c2dtl4_e(ke,coord,t,D,ed,es)
  implicit none
  double precision                :: ke(:,:), coord(:,:), D(:,:,:)
  double precision                :: ed(:), es(:,:), t

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR,ie

  double precision                :: JT(8,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4)
  double precision                :: BL0(3,8), BL(3,8)
  double precision                :: Dgp(3,3)
  DOUBLE PRECISION                :: STRESS(2,2)
  DOUBLE PRECISION                :: HE(4,8)
  DOUBLE PRECISION                :: A(3,4)
  DOUBLE PRECISION                :: dg(4), Dinv(6,6)
  DOUBLE PRECISION                :: R(4,4), Dtmp(3,3), D1inv(3,3), D2inv(3,3)


  JT=MATMUL(DNR,transpose(COORD))

  KE=0D0
  DO GP_NR=1,NGP

    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:8:2)=dNx(1,:)
    BL0(2,2:8:2)=dNx(2,:)
    BL0(3,1:8:2)=dNx(2,:)
    BL0(3,2:8:2)=dNx(1,:)
  
    He=0d0
    He(1:2,1:8:2)=dNx(1:2,:)
    He(3:4,2:8:2)=dNx(1:2,:)

! displacement gradient
    dg=matmul(He,ed)
    A=0d0
    A([1,3],1)=dg(1:2)
    A([1,3],3)=dg(3:4)
    A([2,3],2)=(/dg(2),dg(1)/)
    A([2,3],4)=(/dg(4),dg(3)/)

    BL=BL0+MATMUL(A,He)

    ! es = [S11 S22 S33 S12]
    stress(1,:)=(/es(1,gp_nr), es(4,gp_nr)/)
    stress(2,:)=(/es(4,gp_nr), es(2,gp_nr)/)
    R = 0d0
    R(1:2,1:2)=STRESS(:,:)
    R(3:4,3:4)=STRESS(:,:)
    
    Dgp=D(:,:,gp_nr)
  
    KE=KE+(MATMUL(TRANSPOSE(BL),MATMUL(Dgp,BL)) &
           +MATMUL(MATMUL(TRANSPOSE(He),R),He))*DETJ*t

  END DO
  RETURN
end subroutine c2dtl4_e



subroutine c2dtl4_egamma(ke,coord,t,D,ed,es,gammagp)
  implicit none
  double precision                :: ke(:,:), coord(:,:), D(:,:,:), gammagp(:)
  double precision                :: ed(:), es(:,:), t

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR,ie

  double precision                :: JT(8,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4)
  double precision                :: BL0(3,8), BL(3,8)
  double precision                :: Dgp(3,3)
  DOUBLE PRECISION                :: STRESS(2,2)
  DOUBLE PRECISION                :: HE(4,8)
  DOUBLE PRECISION                :: A(3,4)
  DOUBLE PRECISION                :: dg(4), Dinv(6,6)
  DOUBLE PRECISION                :: R(4,4), Dtmp(3,3), D1inv(3,3), D2inv(3,3)


  JT=MATMUL(DNR,transpose(COORD))

  KE=0D0
  DO GP_NR=1,NGP

    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:8:2)=dNx(1,:)
    BL0(2,2:8:2)=dNx(2,:)
    BL0(3,1:8:2)=dNx(2,:)
    BL0(3,2:8:2)=dNx(1,:)
  
    He=0d0
    He(1:2,1:8:2)=dNx(1:2,:)
    He(3:4,2:8:2)=dNx(1:2,:)

! displacement gradient
    dg=matmul(He,ed)*gammagp(GP_NR)
    A=0d0
    A([1,3],1)=dg(1:2)
    A([1,3],3)=dg(3:4)
    A([2,3],2)=(/dg(2),dg(1)/)
    A([2,3],4)=(/dg(4),dg(3)/)

    BL=BL0+MATMUL(A,He)

    ! es = [S11 S22 S33 S12]
    stress(1,:)=(/es(1,gp_nr), es(4,gp_nr)/)
    stress(2,:)=(/es(4,gp_nr), es(2,gp_nr)/)
    R = 0d0
    R(1:2,1:2)=STRESS(:,:)
    R(3:4,3:4)=STRESS(:,:)
    
    Dgp=D(:,:,gp_nr)
  
    KE=KE+(MATMUL(TRANSPOSE(BL),MATMUL(Dgp,BL)) &
           +MATMUL(MATMUL(TRANSPOSE(He),R),He))*DETJ*t

  END DO
  RETURN
end subroutine c2dtl4_egamma



subroutine c2dtl4_egammaSens(ke,coord,t,D,ed,larg,rarg,es,gammagp)
  implicit none
  double precision                :: ke(:), coord(:,:), D(:,:,:), gammagp(:)
  double precision                :: ed(:), es(:,:), t, larg(:), rarg(:)

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR,ie

  double precision                :: JT(8,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4)
  double precision                :: BL0(3,8), BL(3,8)
  double precision                :: Dgp(3,3)
  DOUBLE PRECISION                :: STRESS(2,2)
  DOUBLE PRECISION                :: HE(4,8)
  DOUBLE PRECISION                :: A(3,4)
  DOUBLE PRECISION                :: dg(4), Dinv(6,6)
  DOUBLE PRECISION                :: R(4,4), Dtmp(3,3), V2(4), kegp


  JT=MATMUL(DNR,transpose(COORD))

  KE=0D0
  DO GP_NR=1,NGP

    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:8:2)=dNx(1,:)
    BL0(2,2:8:2)=dNx(2,:)
    BL0(3,1:8:2)=dNx(2,:)
    BL0(3,2:8:2)=dNx(1,:)
  
    He=0d0
    He(1:2,1:8:2)=dNx(1:2,:)
    He(3:4,2:8:2)=dNx(1:2,:)

! displacement gradient
    dg=matmul(He,ed)*gammagp(GP_NR)
    A=0d0
    A([1,3],1)=dg(1:2)
    A([1,3],3)=dg(3:4)
    A([2,3],2)=(/dg(2),dg(1)/)
    A([2,3],4)=(/dg(4),dg(3)/)

    BL=BL0+MATMUL(A,He)

    ! es = [S11 S22 S33 S12]
    stress(1,:)=(/es(1,gp_nr), es(4,gp_nr)/)
    stress(2,:)=(/es(4,gp_nr), es(2,gp_nr)/)
    R = 0d0
    R(1:2,1:2)=STRESS(:,:)
    R(3:4,3:4)=STRESS(:,:)
    
    Dgp=D(:,:,gp_nr)
    V2(:)=NR(GP_NR,:)
  
    KEGP=dot_product(larg,(matmul((MATMUL(TRANSPOSE(BL),MATMUL(Dgp,BL)) &
           +MATMUL(MATMUL(TRANSPOSE(He),R),He)),rarg)))*DETJ*t

    ke=ke+V2*KEGP

  END DO
  RETURN
end subroutine c2dtl4_egammaSens


! calculates the total deformation gradient
subroutine c2dtl4_d(dg,coord,ed)
  implicit none

  double precision                :: dg(:,:), coord(:,:),  ed(:)

  integer                         :: gp_nr
  integer, parameter              :: NGP=4

  double precision                :: JT(8,2), JTinv(2,2)
  double precision                :: dNx(2,4)
  double precision                :: He(4,8)

  JT=MATMUL(DNR,transpose(COORD))
   
  do GP_NR=1,ngp

   CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
   dNx=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
 
   He=0d0
   He(1:2,1:8:2)=dNx(1:2,:)
   He(3:4,2:8:2)=dNx(1:2,:)
     
   dg(:,gp_nr)=Ifm+MATMUL(HE,ED)

  end do

  return

end subroutine c2dtl4_d




subroutine c2dtl4_f(ef,coord,t,ed,es)
  implicit none

  double precision                :: ef(:), coord(:,:),  ed(:), es(:,:), t

  integer, parameter              :: NGP=4
  integer                         :: gp_nr

  double precision                :: JT(8,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4), esTMP(3)
  double precision                :: BL0(3,8), He(4,8)
  double precision                :: A(3,4), dg(4)
 
  JT=MATMUL(DNR,transpose(COORD))

  ef=0D0
  DO GP_NR=1,NGP

    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:8:2)=dNx(1,:)
    BL0(2,2:8:2)=dNx(2,:)
    BL0(3,1:8:2)=dNx(2,:)
    BL0(3,2:8:2)=dNx(1,:)
  
    He=0d0
    He(1:2,1:8:2)=dNx(1:2,:)
    He(3:4,2:8:2)=dNx(1:2,:)

! displacement gradient
    dg=matmul(He,ed)
    A=0d0
    A([1,3],1)=dg(1:2)
    A([1,3],3)=dg(3:4)
    A([2,3],2)=(/dg(2),dg(1)/)
    A([2,3],4)=(/dg(4),dg(3)/)
    
    ! es = [S11 S22 S33 S12]
    esTMP = es([1, 2, 4],gp_nr)

    EF=EF+MATMUL(TRANSPOSE(BL0+MATMUL(A,He)),esTMP)*detJ*t

  END DO
  RETURN

end subroutine c2dtl4_f



subroutine c2dtl4_fgamma(ef,coord,t,ed,es,gammagp)
  implicit none

  double precision                :: ef(:), coord(:,:),  ed(:), es(:,:), t, gammagp(:)

  integer, parameter              :: NGP=4
  integer                         :: gp_nr

  double precision                :: JT(8,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4), esTMP(3)
  double precision                :: BL0(3,8), He(4,8)
  double precision                :: A(3,4), dg(4)
 
  JT=MATMUL(DNR,transpose(COORD))

  ef=0D0
  DO GP_NR=1,NGP

    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:8:2)=dNx(1,:)
    BL0(2,2:8:2)=dNx(2,:)
    BL0(3,1:8:2)=dNx(2,:)
    BL0(3,2:8:2)=dNx(1,:)
  
    He=0d0
    He(1:2,1:8:2)=dNx(1:2,:)
    He(3:4,2:8:2)=dNx(1:2,:)

! displacement gradient
    dg=matmul(He,ed)*gammagp(GP_NR)
    A=0d0
    A([1,3],1)=dg(1:2)
    A([1,3],3)=dg(3:4)
    A([2,3],2)=(/dg(2),dg(1)/)
    A([2,3],4)=(/dg(4),dg(3)/)
   
    ! es = [S11 S22 S33 S12]
    esTMP = es([1, 2, 4],gp_nr)

    EF=EF+MATMUL(TRANSPOSE(BL0+MATMUL(A,He)),esTMP)*detJ*t

  END DO
  RETURN

end subroutine c2dtl4_fgamma





subroutine c2dtl4_fgammaSens(ef,coord,t,ed,edvarphi,es,gammagp)
  implicit none

  double precision                :: ef(:), coord(:,:),  ed(:), es(:,:), t, gammagp(:)
  double precision                :: edvarphi(:)

  integer, parameter              :: NGP=4
  integer                         :: gp_nr

  double precision                :: JT(8,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4), esTMP(3)
  double precision                :: BL0(3,8), He(4,8)
  double precision                :: A(3,4), dg(4), efgp, v2(4)
 
  JT=MATMUL(DNR,transpose(COORD))

  ef=0D0
  DO GP_NR=1,NGP

    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:8:2)=dNx(1,:)
    BL0(2,2:8:2)=dNx(2,:)
    BL0(3,1:8:2)=dNx(2,:)
    BL0(3,2:8:2)=dNx(1,:)
  
    He=0d0
    He(1:2,1:8:2)=dNx(1:2,:)
    He(3:4,2:8:2)=dNx(1:2,:)

! displacement gradient
    dg=matmul(He,ed)*gammagp(GP_NR)
    A=0d0
    A([1,3],1)=dg(1:2)
    A([1,3],3)=dg(3:4)
    A([2,3],2)=(/dg(2),dg(1)/)
    A([2,3],4)=(/dg(4),dg(3)/)
    
    V2=NR(GP_NR,:)

    ! es = [S11 S22 S33 S12]
    esTMP = es([1, 2, 4],gp_nr)

    efgp=dot_product(edvarphi,MATMUL(TRANSPOSE(BL0+MATMUL(A,He)),esTMP))*detJ*t
    ef=ef+V2*efgp

  END DO
  RETURN

end subroutine c2dtl4_fgammaSens


subroutine elmvol(evol,coord,t)
  implicit none
  double precision                :: evol, coord(:,:), t

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR

  double precision                :: JT(R2,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4)

  JT=MATMUL(DNR,transpose(COORD))

  evol=0D0
  DO GP_NR=1,NGP
    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 
    evol=evol+DETJ*t
  END DO
  RETURN
END subroutine elmvol


subroutine c2dtl4_emu2(coord,ed,edmutilde,Emutilde,dgmutilde,gammagp,ones)
  implicit none
  double precision                :: coord(:,:)
  double precision                :: ed(:), edmutilde(:), Emutilde(:,:), gammagp(:), ones(:)

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR

  double precision                :: JT(R2,2), DetJ, JTinv(2,2)
  double precision                :: TMP(2,4), DNX(2,4)
  double precision                :: BL0(3,8), BL(3,8)
  DOUBLE PRECISION                :: HE(4,8)
  DOUBLE PRECISION                :: A(3,4)
  DOUBLE PRECISION                :: dg(4), dgmutilde(:,:)

  JT=MATMUL(DNR,transpose(COORD))
  
  DO GP_NR=1,NGP

    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:8:2)=dNx(1,:)
    BL0(2,2:8:2)=dNx(2,:)
    BL0(3,1:8:2)=dNx(2,:)
    BL0(3,2:8:2)=dNx(1,:)
  
    He=0d0
    He(1:2,1:8:2)=dNx(1:2,:)
    He(3:4,2:8:2)=dNx(1:2,:)

! displacement gradient
    dg=matmul(He,ed)*gammagp(GP_NR)
    A=0d0
    A([1,3],1)=dg(1:2)
    A([1,3],3)=dg(3:4)
    A([2,3],2)=(/dg(2),dg(1)/)
    A([2,3],4)=(/dg(4),dg(3)/)

    BL=BL0+MATMUL(A,He)

    dgmutilde(:,GP_Nr)=matmul(He,edmutilde)
    Emutilde(:,GP_Nr) =matmul(BL,edmutilde)
  END DO
  RETURN
END subroutine c2dtl4_emu2


subroutine c2dtl4_emu(coord,ed,edmutilde,Emutilde,dgmutilde,gammagp)
  implicit none
  double precision                :: coord(:,:)
  double precision                :: ed(:), edmutilde(:), Emutilde(:,:), gammagp(:)

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR

  double precision                :: JT(R2,2), DetJ, JTinv(2,2)
  double precision                :: TMP(2,4), DNX(2,4)
  double precision                :: BL0(3,8), BL(3,8)
  DOUBLE PRECISION                :: HE(4,8)
  DOUBLE PRECISION                :: A(3,4)
  DOUBLE PRECISION                :: dg(4), dgmutilde(:,:)

  JT=MATMUL(DNR,transpose(COORD))
  
  DO GP_NR=1,NGP

    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:8:2)=dNx(1,:)
    BL0(2,2:8:2)=dNx(2,:)
    BL0(3,1:8:2)=dNx(2,:)
    BL0(3,2:8:2)=dNx(1,:)
  
    He=0d0
    He(1:2,1:8:2)=dNx(1:2,:)
    He(3:4,2:8:2)=dNx(1:2,:)

! displacement gradient
    dg=matmul(He,ed)*gammagp(GP_NR)
    A=0d0
    A([1,3],1)=dg(1:2)
    A([1,3],3)=dg(3:4)
    A([2,3],2)=(/dg(2),dg(1)/)
    A([2,3],4)=(/dg(4),dg(3)/)

    BL=BL0+MATMUL(A,He)


    dgmutilde(:,GP_Nr)=matmul(He,edmutilde)*gammagp(gp_nr)
    Emutilde(:,GP_Nr) =matmul(BL,edmutilde)*gammagp(gp_nr)
  END DO
  RETURN
END subroutine c2dtl4_emu




subroutine c2dtl4_b(efb,coord,t,b)
  implicit none

  double precision                :: efb(:), coord(:,:), b(:,:), t

  integer, parameter              :: NGP=4
  integer                         :: gp_nr

  double precision                :: JT(8,2), DetJ, JTinv(2,2)
  double precision                :: N(2,8), btmp(2)
 
  JT=MATMUL(DNR,transpose(COORD))

  efb=0D0
  DO GP_NR=1,NGP

    DETJ=det2(JT(index(:,gp_nr),:)) 

    N = 0d0
    N(1,1:8:2) = NR(GP_NR,:)
    N(2,2:8:2) = NR(GP_NR,:)
    
    bTMP = b(:,gp_nr)

    EFB=EFB+MATMUL(TRANSPOSE(N),bTMP)*detJ*t

  END DO
  RETURN

end subroutine c2dtl4_b


subroutine c2dtl4_bsens(efb,coord,t,b,arg)
  implicit none

  double precision                :: efb(:), coord(:,:), b(:,:), t, arg(:)

  integer, parameter              :: NGP=4
  integer                         :: gp_nr

  double precision                :: JT(8,2), DetJ, JTinv(2,2)
  double precision                :: N(2,8), btmp(2), V2(4), efbgp
 
  JT=MATMUL(DNR,transpose(COORD))

  efb=0D0
  DO GP_NR=1,NGP

    DETJ=det2(JT(index(:,gp_nr),:)) 

    N = 0d0
    N(1,1:8:2) = NR(GP_NR,:)
    N(2,2:8:2) = NR(GP_NR,:)
    
    V2=NR(GP_NR,:)

    bTMP = b(:,gp_nr)

    EFBgp=dot_product(MATMUL(TRANSPOSE(N),bTMP),arg)*detJ*t
    efb=efb+V2*efbgp
  
  END DO
  RETURN

end subroutine c2dtl4_bsens




subroutine c2dtl4_m(Me,coord,t,val,rhogp)
  implicit none
  double precision                :: Me(:,:), coord(:,:)
  double precision                :: val,t

  integer, parameter              :: NGP=4
  integer                         :: GP_NR,ie

  double precision                :: rhogp(:), rhogpCurr
  double precision                :: JT(8,2), DetJ
  DOUBLE PRECISION                :: N(8,2)

  JT=MATMUL(DNR,transpose(COORD))

  Me=0D0
  DO GP_NR=1,NGP    

    DETJ=det2(JT(index(:,gp_nr),:)) 
   
    N = 0d0 
    N(1:8:2,1)=NR(GP_NR,:)
    N(2:8:2,2)=NR(GP_NR,:)
    rhogpCurr = rhogp(GP_NR)

    Me=Me+rhogpCurr*MATMUL(N,TRANSPOSE(N))*DETJ*val*t

  END DO
  RETURN
END subroutine c2dtl4_m


subroutine dc2dtl4_m(fme,coord,t,val,rhogp,larg,rarg)
  implicit none
  double precision                :: fme(:), coord(:,:), larg(:), rarg(:)
  double precision                :: val, t

  integer, parameter              :: NGP=4
  integer                         :: GP_NR,ie

  double precision                :: rhogp(:),rhogpCurr
  double precision                :: JT(8,2), DetJ 
  DOUBLE PRECISION                :: N(8,2), V2(4), tmp2(8,8), tmp3

  JT=MATMUL(DNR,transpose(COORD))
  fme=0D0
  DO GP_NR=1,NGP    

    DETJ=det2(JT(index(:,gp_nr),:)) 
    
    N = 0d0
    N(1:8:2,1)=NR(GP_NR,:)
    N(2:8:2,2)=NR(GP_NR,:)
    rhogpCurr = rhogp(GP_NR)

    V2(:)=NR(GP_NR,:)
 
    tmp2 = rhogpCurr*MATMUL(N,TRANSPOSE(N))*DETJ*val*t
    tmp3 = dot_product(larg,MATMUL(tmp2,rarg))

    fme=fme+V2*tmp3

  END DO
  RETURN

END subroutine dc2dtl4_m




! 4-node element based on updated Lagrangian formulation	
! Only plane strain implemented
subroutine c2dul4_e(ke,coord,t,D,ed,es)
  implicit none
  double precision                :: ke(:,:), coord(:,:), D(:,:,:)
  double precision                :: ed(:), es(:,:), t

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR,ie

  double precision                :: JT(8,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4)
  double precision                :: BL0(3,8), BL(3,8)
  double precision                :: Dgp(3,3)
  DOUBLE PRECISION                :: STRESS(2,2)
  DOUBLE PRECISION                :: HE(4,8)
  DOUBLE PRECISION                :: A(3,4)
  DOUBLE PRECISION                :: dg(4), Dinv(6,6)
  DOUBLE PRECISION                :: R(4,4), Dtmp(3,3), D1inv(3,3), D2inv(3,3)
  double precision                :: ucoord(2,4)

  ucoord(1,:)=coord(1,:)+ed([1,3,5,7])
  ucoord(2,:)=coord(2,:)+ed([2,4,6,8])

  JT=MATMUL(DNR,transpose(ucoord))

  KE=0D0
  DO GP_NR=1,NGP

    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:8:2)=dNx(1,:)
    BL0(2,2:8:2)=dNx(2,:)
    BL0(3,1:8:2)=dNx(2,:)
    BL0(3,2:8:2)=dNx(1,:)
  
    He=0d0
    He(1:2,1:8:2)=dNx(1:2,:)
    He(3:4,2:8:2)=dNx(1:2,:)

    ! es = [S11 S22 S33 S12]
    stress(1,:)=(/es(1,gp_nr), es(4,gp_nr)/)
    stress(2,:)=(/es(4,gp_nr), es(2,gp_nr)/)
    R = 0d0
    R(1:2,1:2)=STRESS(:,:)
    R(3:4,3:4)=STRESS(:,:)
    
    Dgp=D(:,:,gp_nr)
  
    KE=KE+(MATMUL(TRANSPOSE(BL0),MATMUL(Dgp,BL0)) &
           +MATMUL(MATMUL(TRANSPOSE(He),R),He))*DETJ*t

  END DO
  RETURN
end subroutine c2dul4_e




subroutine c2dul4_esens(ke,coord,t,D,ed,larg,rarg,es)
  implicit none
  double precision                :: ke(:), coord(:,:), D(:,:,:)
  double precision                :: ed(:), es(:,:), t, larg(:), rarg(:)

  integer, parameter              :: NGP=4, R2=NGP*2
  integer                         :: GP_NR,ie

  double precision                :: JT(8,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4)
  double precision                :: BL0(3,8), BL(3,8)
  double precision                :: Dgp(3,3)
  DOUBLE PRECISION                :: STRESS(2,2)
  DOUBLE PRECISION                :: HE(4,8)
  DOUBLE PRECISION                :: A(3,4)
  DOUBLE PRECISION                :: dg(4), Dinv(6,6)
  DOUBLE PRECISION                :: R(4,4), Dtmp(3,3), V2(4), kegp
  double precision                :: ucoord(2,4)

  ucoord(1,:)=coord(1,:)+ed([1,3,5,7])
  ucoord(2,:)=coord(2,:)+ed([2,4,6,8])

  JT=MATMUL(DNR,transpose(ucoord))

  KE=0D0
  DO GP_NR=1,NGP

    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:8:2)=dNx(1,:)
    BL0(2,2:8:2)=dNx(2,:)
    BL0(3,1:8:2)=dNx(2,:)
    BL0(3,2:8:2)=dNx(1,:)
  
    He=0d0
    He(1:2,1:8:2)=dNx(1:2,:)
    He(3:4,2:8:2)=dNx(1:2,:)

    ! es = [S11 S22 S33 S12]
    stress(1,:)=(/es(1,gp_nr), es(4,gp_nr)/)
    stress(2,:)=(/es(4,gp_nr), es(2,gp_nr)/)
    R = 0d0
    R(1:2,1:2)=STRESS(:,:)
    R(3:4,3:4)=STRESS(:,:)
    
    Dgp=D(:,:,gp_nr)
    V2(:)=NR(GP_NR,:)
  
    KEGP=dot_product(larg,(matmul((MATMUL(TRANSPOSE(BL0),MATMUL(Dgp,BL0)) &
           +MATMUL(MATMUL(TRANSPOSE(He),R),He)),rarg)))*DETJ*t

    ke=ke+V2*KEGP

  END DO
  RETURN
end subroutine c2dul4_eSens



subroutine c2dul4_f(ef,coord,t,ed,es)
  implicit none

  double precision                :: ef(:), coord(:,:),  ed(:), es(:,:), t

  integer, parameter              :: NGP=4
  integer                         :: gp_nr

  double precision                :: JT(8,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4), esTMP(3)
  double precision                :: BL0(3,8), He(4,8)
  double precision                :: A(3,4), dg(4)
  double precision                :: ucoord(2,4)

  ucoord(1,:)=coord(1,:)+ed([1,3,5,7])
  ucoord(2,:)=coord(2,:)+ed([2,4,6,8])

  JT=MATMUL(DNR,transpose(ucoord))

  ef=0D0
  DO GP_NR=1,NGP

    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:8:2)=dNx(1,:)
    BL0(2,2:8:2)=dNx(2,:)
    BL0(3,1:8:2)=dNx(2,:)
    BL0(3,2:8:2)=dNx(1,:)
  
    ! es = [S11 S22 S33 S12]
    esTMP = es([1, 2, 4],gp_nr)

    EF=EF+MATMUL(TRANSPOSE(BL0),esTMP)*detJ*t

  END DO
  RETURN

end subroutine c2dul4_f






subroutine c2dul4_fSens(ef,coord,t,ed,edvarphi,es)
  implicit none

  double precision                :: ef(:), coord(:,:),  ed(:), es(:,:), t
  double precision                :: edvarphi(:)

  integer, parameter              :: NGP=4
  integer                         :: gp_nr

  double precision                :: JT(8,2), DetJ, JTinv(2,2)
  double precision                :: DNX(2,4), esTMP(3)
  double precision                :: BL0(3,8), He(4,8)
  double precision                :: A(3,4), dg(4), efgp, v2(4)
  double precision                :: ucoord(2,4)

  ucoord(1,:)=coord(1,:)+ed([1,3,5,7])
  ucoord(2,:)=coord(2,:)+ed([2,4,6,8])

  JT=MATMUL(DNR,transpose(ucoord))

  ef=0D0
  DO GP_NR=1,NGP

    CALL INV2(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det2(JT(index(:,gp_nr),:)) 

    BL0=0D0
    BL0(1,1:8:2)=dNx(1,:)
    BL0(2,2:8:2)=dNx(2,:)
    BL0(3,1:8:2)=dNx(2,:)
    BL0(3,2:8:2)=dNx(1,:)
  
    V2=NR(GP_NR,:)

    ! es = [S11 S22 S33 S12]
    esTMP = es([1, 2, 4],gp_nr)

    efgp=dot_product(edvarphi,MATMUL(TRANSPOSE(BL0),esTMP))*detJ*t
    ef=ef+V2*efgp

  END DO
  RETURN

end subroutine c2dul4_fSens



end module elem_large_cont_2d
