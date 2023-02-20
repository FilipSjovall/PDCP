module elem_large_cont_3d_inv

! Last rev.:
!  M. Ristinmaa 2011-10-07
!   - initial version 
!  M. Ristinmaa 2011-10-13
!   - modified to use global parameters
!------------------------------------------------------------------------------

!
!   Module elem_large_cont_pot_3d_inv contains element subroutines 
!   potential field to large deformations coupling
!   for the inverse motion problem

! 8 node brick element, trilinear displacement interpolation
! 2x2x2 gauss integration points (in total 8), the weight for every
! gauss point is one for this scheme

use matrix_util, only: inv3, det3
use mater_large, only: getF

implicit none



! This is an attempt to predefine some vectors, to obtain speed
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

! local degrees of freedom
integer, parameter,dimension(24):: u_list=(/ 1, 2, 3, & 
                                             5, 6, 7, &
                                             9,10,11, &
                                            13,14,15, &
                                            17,18,19, &
                                            21,22,23, &
                                            25,26,27, &
                                            29,30,31/)
integer, parameter,dimension(8) :: p_list=(/4,8,12,16,20,24,28,32/)
private u_list, p_list

private det3, inv3

interface c3dinv8p0_f
  module procedure c3dinv8_f
end interface

	
!------------------------------------------------------------------------------
contains


! 8-node brick element for use in a inverse motion formulation	
subroutine c3dinv8_e(ke,coord,duu)
  implicit none
  double precision                :: ke(:,:), coord(:,:)
  double precision                :: duu(:,:,:)

  integer, parameter              :: NGP=8
  integer                         :: gp_nr

  double precision                :: JT(24,3), DetJ, JTinv(3,3)
  double precision                :: Bp(3,8)
  double precision                :: Bu(6,24)
  DOUBLE PRECISION                :: HE(9,24)

  JT=MATMUL(DNR,transpose(coord))

  Ke=0D0
  do GP_NR=1,NGP

    call INV3(JTinv,JT(INDEX(:,gp_nr),:))
    Bp=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:)) !dNx
    DETJ=det3(JT(index(:,gp_nr),:)) 

    Bu=0D0
    Bu(1,1:24:3)=Bp(1,:)
    Bu(2,2:24:3)=Bp(2,:)
    Bu(3,3:24:3)=Bp(3,:)
    Bu(4,1:24:3)=Bp(2,:)
    Bu(4,2:24:3)=Bp(1,:)
    Bu(5,1:24:3)=Bp(3,:)
    Bu(5,3:24:3)=Bp(1,:)
    Bu(6,2:24:3)=Bp(3,:)
    Bu(6,3:24:3)=Bp(2,:)
  
    He=0d0
    He(1:3,1:24:3)=Bp(1:3,:)
    He(4:6,2:24:3)=Bp(1:3,:)
    He(7:9,3:24:3)=Bp(1:3,:)

    ke=ke+(MATMUL(TRANSPOSE(Bu),MATMUL(Duu(:,:,gp_nr),He)))*DETJ

  end do

  return   
end subroutine c3dinv8_e


subroutine c3dinv8_f(ef,coord,es)
  implicit none

  double precision                :: ef(:), coord(:,:), es(:,:)

  integer, parameter              :: NGP=8
  integer                         :: gp_nr

  double precision                :: JT(24,3), JTinv(3,3)
  double precision                :: detJ
  double precision                :: Bu(6,24), Bp(3,8)

  JT=MATMUL(DNR,transpose(coord))

  ef=0d0
  do GP_NR=1,NGP
    call INV3(JTinv,JT(INDEX(:,gp_nr),:))
    Bp=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:)) !dNx
    detJ=det3(JT(index(:,gp_nr),:)) 

    Bu=0D0
    Bu(1,1:24:3) =Bp(1,:)
    Bu(2,2:24:3) =Bp(2,:)
    Bu(3,3:24:3) =Bp(3,:)
    Bu(4,1:24:3) =Bp(2,:)
    Bu(4,2:24:3) =Bp(1,:)
    Bu(5,1:24:3) =Bp(3,:)
    Bu(5,3:24:3) =Bp(1,:)
    Bu(6,2:24:3) =Bp(3,:)
    Bu(6,3:24:3) =Bp(2,:)

    ef=ef+MATMUL(TRANSPOSE(Bu),es(:,gp_nr))*detJ

  end do

  return
end subroutine c3dinv8_f



! calculates the total inverse deformation gradient
subroutine c3dinv8_d(dg,coord,ed)
  implicit none

  double precision                :: dg(:,:), coord(:,:),  ed(:)

  integer, parameter              :: NGP=8
  integer                         :: gp_nr

  double precision                :: JT(24,3), JTinv(3,3)
  double precision                :: DNX(3,8), He(9,24)
  

  JT=MATMUL(DNR,transpose(COORD))
   
  do GP_NR=1,NGP
    call INV3(JTinv,JT(INDEX(:,gp_nr),:))
    dNx=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:)) 

! Written such that dg is obtained directly
    He=0D0
    He(1:3,1:24:3)=dNx(1:3,:)
    He(4:6,2:24:3)=dNx(1:3,:)
    He(7:9,3:24:3)=dNx(1:3,:)
     
! Deformation gradient      
    dg(:,gp_nr)=Ifm+MATMUL(HE,ED)

  end do
  
  return
end subroutine c3dinv8_d


! Use of u/p elements
subroutine cp3dinv8p0_e(ke,coord,duu,dup,dut,dpu,dpp,ed,th)
  implicit none
  double precision                :: ke(:,:), coord(:,:)
  double precision                :: duu(:,:,:), dup(:,:,:), dut(:)
  double precision                :: dpu(:,:,:), dpp(:,:,:)
  double precision                :: ed(:), th

  integer, parameter              :: NGP=8
  integer                         :: gp_nr

  double precision                :: JT(24,3), JT0(24,3), DetJ, DetJ0, JTinv(3,3)
  double precision                :: Bp(3,8), Bp0(3,8)
  double precision                :: Bu(6,24), Bu0(6,24)
!  DOUBLE PRECISION                :: STRESS(3,3)
  DOUBLE PRECISION                :: HE(9,24), PBu(1,24), JPBu(1,24)
  double precision                :: kuu(24,24), Kup(24,8), Kpu(8,24), Kpp(8,8)
  double precision                :: f(3,3), Volume0, Volume, j
  double precision                :: ucoord(3,8), dg(9)

  double precision, parameter     :: P(1,6)=(/1d0,1d0,1d0,0d0,0d0,0d0/)


  ucoord(1,:)=coord(1,:)+ed([1,5, 9,13,17,21,25,29])
  ucoord(2,:)=coord(2,:)+ed([2,6,10,14,18,22,26,30])
  ucoord(3,:)=coord(3,:)+ed([3,7,11,15,19,23,27,31])

  JT=MATMUL(DNR,transpose(coord))
  JT0=MATMUL(DNR,transpose(ucoord))

  Kuu=0D0
  Kup=0D0
  Kpu=0D0
  Kpp=0D0
  PBu=0d0
  JPBu=0d0
  Volume=0d0
  Volume0=0d0

  do GP_NR=1,NGP

    call INV3(JTinv,JT(INDEX(:,gp_nr),:))
    Bp=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:)) !dNx
    DETJ=det3(JT(index(:,gp_nr),:)) 
    DETJ0=det3(JT0(index(:,gp_nr),:)) 

    call INV3(JTinv,JT0(INDEX(:,gp_nr),:))
    Bp0=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:)) !dNx

    Bu=0D0
    Bu(1,1:24:3)=Bp(1,:)
    Bu(2,2:24:3)=Bp(2,:)
    Bu(3,3:24:3)=Bp(3,:)
    Bu(4,1:24:3)=Bp(2,:)
    Bu(4,2:24:3)=Bp(1,:)
    Bu(5,1:24:3)=Bp(3,:)
    Bu(5,3:24:3)=Bp(1,:)
    Bu(6,2:24:3)=Bp(3,:)
    Bu(6,3:24:3)=Bp(2,:)

    Bu0=0D0
    Bu0(1,1:24:3)=Bp0(1,:)
    Bu0(2,2:24:3)=Bp0(2,:)
    Bu0(3,3:24:3)=Bp0(3,:)
    Bu0(4,1:24:3)=Bp0(2,:)
    Bu0(4,2:24:3)=Bp0(1,:)
    Bu0(5,1:24:3)=Bp0(3,:)
    Bu0(5,3:24:3)=Bp0(1,:)
    Bu0(6,2:24:3)=Bp0(3,:)
    Bu0(6,3:24:3)=Bp0(2,:)

    He=0d0
    He(1:3,1:24:3)=Bp(1:3,:)
    He(4:6,2:24:3)=Bp(1:3,:)
    He(7:9,3:24:3)=Bp(1:3,:)

! Note that dg is the invers deformation gradient     
    dg=Ifm+MATMUL(HE,ED(u_list))
    f=getF(dg)
    j=det3(f)

! Note that detJ0=detJ*j
    PBu=PBu+matmul(P,Bu)*DetJ
    JPBu=JPBu+matmul(P,Bu0)*DetJ*j

    kuu=kuu+(MATMUL(TRANSPOSE(Bu),MATMUL(Duu(:,:,gp_nr),He)))*DETJ
    kup=kup-(MATMUL(TRANSPOSE(Bu),MATMUL(Dup(:,:,gp_nr),Bp)))*DETJ
    kpu=kpu-(MATMUL(TRANSPOSE(Bp),MATMUL(Dpu(:,:,gp_nr),He)))*DETJ
    kpp=kpp+(MATMUL(TRANSPOSE(Bp),MATMUL(Dpp(:,:,gp_nr),Bp)))*DETJ

    Volume=Volume+DetJ
    Volume0=Volume0+DetJ0

  end do
! write(*,*)'diff volume',Volume/Volume0, th
  Ke(u_list,u_list)=kuu-matmul(transpose(PBu),JPBu)*dut(1)/Volume0*th
  Ke(u_list,p_list)=kup
  Ke(p_list,u_list)=kpu
  Ke(p_list,p_list)=kpp

  return   
end subroutine cp3dinv8p0_e


! calculates the total inverse deformation gradient and the Electric field
! and the element averaged volume jacobian
subroutine cp3dinv8p0_d(dg,de,th,coord,ed)
  implicit none

  double precision                :: dg(:,:), de(:,:), coord(:,:),  ed(:), th

  integer, parameter              :: NGP=8
  integer                         :: gp_nr

  double precision                :: JT(24,3), JT0(24,3), JTinv(3,3)
  double precision                :: DNX(3,8), He(9,24)
  
!  double precision                :: invF(3,3)
  double precision                :: p_val(8), ucoord(3,8)
  double precision                :: DetJ, DetJ0
  double precision                :: Volume, Volume0

! current nodal values
  ucoord(1,:)=coord(1,:)+ed([1,5, 9,13,17,21,25,29])
  ucoord(2,:)=coord(2,:)+ed([2,6,10,14,18,22,26,30])
  ucoord(3,:)=coord(3,:)+ed([3,7,11,15,19,23,27,31])


  p_val=ed(p_list)

  JT=MATMUL(DNR,transpose(COORD))
  JT0=MATMUL(DNR,transpose(ucoord))

  Volume=0d0
  Volume0=0d0
   
  do GP_NR=1,NGP
    call INV3(JTinv,JT(INDEX(:,gp_nr),:))
    dNx=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:)) 
    DETJ=det3(JT(index(:,gp_nr),:)) 
    DETJ0=det3(JT0(index(:,gp_nr),:)) 

    He=0D0
    He(1:3,1:24:3)=dNx(1:3,:)
    He(4:6,2:24:3)=dNx(1:3,:)
    He(7:9,3:24:3)=dNx(1:3,:)
     
! Deformation gradient      
    dg(:,gp_nr)=Ifm+MATMUL(HE,ED(u_list))

! Spatial gradient of scalar field
    de(:,gp_nr)=-matmul(dNx,p_val)

    Volume=Volume+DetJ
    Volume0=Volume0+DetJ0

  end do

  th=Volume/Volume0
  
  return
end subroutine cp3dinv8p0_d

 

end module elem_large_cont_3d_inv
