module elem_large_mix_th_2d

! Last rev.:
!  M. Ristinmaa 2012-01-12
!   - Initial code _f subroutine
!  M. Ristinmaa 2012-01-15
!   - _d subroutine included
!  M. Ristinmaa 2012-01-17
!   - bugg, inv3 and det3 were used.
!   - bugg, p_list
!  M. Ristinmaa 2013-04-14
!  - _bn subroutine included
!------------------------------------------------------------------------------	

!
!   Module elem_large_cont_2d_th contains element subroutines 
!   for 2-D mixture formulation (displacement - pressure)
!   The balance equations than can be solved are defined by the equation set
!
!    int(Bu^T s dA) - int(Nu^T t dS) =0
!    int(Bp^T h dA) + int(Np^T g dA) - int(Np^T q dS)=0
!
!   where the possible constitutive assumptions are
!
!    h=h(grad(u),p,grad(p))
!    g=g(grad(u),p)
!    q=constant
!    s=s(grad(u),p)
!    t=constant
!
!   - Q2Q1 2D Taylor-Hood element (8 node for displacements and 4 for pressure)
!  
!   nodal numbering: corner node Counter Clock Wise, then mid node CCW
!   number of dof=(3,3,3,3,2,2,2,2) total 20 dofs
!   3 for corner nodes (displacements and pressure)
!   2 for mid nodes    (displacements)
!
!   For the coordinate mapping the 8 nodes are used for both the displacement
!   and the pressure.

use matrix_util, only: inv2, det2

implicit none
integer, parameter :: nr_gp=9, nr_p=4, nr_u=8, dim=2
integer, parameter :: nrgpd=nr_gp*dim

double precision   :: xsi(nr_gp), eta(nr_gp), w_1(nr_gp), w_2(nr_gp)
double precision   :: G1, G2, W1, W2
double precision   :: Ifm(4)
integer            :: index(dim,nr_gp)

parameter        (G1=0.774596669241483D0)
parameter        (G2=0.D0)
parameter        (W1=0.555555555555555D0)
parameter        (W2=0.888888888888888D0)

! same order as Abaqus
parameter        (xsi=(/-g1, g2, g1,-g1, g2, g1,-g1, g2, g1/))
parameter        (eta=(/-g1,-g1,-g1, g2, g2, g2, g1, g1, g1/))
parameter        (w_1=(/w1, w2, w1, w1, w2, w1, w1, w2, w1/))
parameter        (w_2=(/w1, w1, w1, w2, w2, w2, w1, w1, w1/))

parameter        (Ifm=(/1d0,0d0,0d0,1d0/))

parameter        (index=[(/1,2/),(/3,4/),(/5,6/),(/7,8/),(/9,10/),(/11,12/), &
                         (/13,14/),(/15,16/),(/17,18/)])

! list of dof for element (positions)
integer, parameter,dimension(16):: u_list=(/ 1, 2, & 
                                             4, 5, &
                                             7, 8, &
                                             10, 11, &
                                             13, 14, &
                                             15, 16, &
                                             17, 18, &
                                             19, 20/)
integer, parameter,dimension(4) :: p_list=(/3,6,9,12/)

! pressure field
! shape functions
double precision  Np(nr_gp,nr_p) 
data  (Np(:,1)=(1D0-xsi)*(1D0-eta)/4.D0)
data  (Np(:,2)=(1D0+xsi)*(1D0-eta)/4.D0)
data  (Np(:,3)=(1D0+xsi)*(1D0+eta)/4.D0)
data  (Np(:,4)=(1D0-xsi)*(1D0+eta)/4.D0)
double precision  DNpR(nrgpd,nr_p) 
! derivate of shape functions for p with respect to xsi
data  (dNpR(1:nrgpd-1:2,1)=-(1D0-eta)/4D0)
data  (dNpR(1:nrgpd-1:2,2)= (1D0-eta)/4D0)
data  (dNpR(1:nrgpd-1:2,3)= (1D0+eta)/4D0)
data  (dNpR(1:nrgpd-1:2,4)=-(1D0+eta)/4D0)
! derivate of shape functions for p with respect to eta
data  (dNpR(2:nrgpd:2,1)=-(1D0-xsi)/4D0)
data  (dNpR(2:nrgpd:2,2)=-(1D0+xsi)/4D0)
data  (dNpR(2:nrgpd:2,3)= (1D0+xsi)/4D0)
data  (dNpR(2:nrgpd:2,4)= (1D0-xsi)/4D0)
	
! diplacement field
! shape functions
 double precision  Nu(nr_gp,nr_u) 
data  (Nu(:,1)=-(1D0-xsi)*(1D0-eta)*(1D0+xsi+eta)/4D0)
data  (Nu(:,2)=-(1D0+xsi)*(1D0-eta)*(1D0-xsi+eta)/4D0)
data  (Nu(:,3)=-(1D0+xsi)*(1D0+eta)*(1D0-xsi-eta)/4D0)
data  (Nu(:,4)=-(1D0-xsi)*(1D0+eta)*(1D0+xsi-eta)/4D0)
data  (Nu(:,5)=(1D0-xsi*xsi)*(1D0-eta)/2D0)
data  (Nu(:,6)=(1D0+xsi)*(1D0-eta*eta)/2D0)
data  (Nu(:,7)=(1D0-xsi*xsi)*(1D0+eta)/2D0)
data  (Nu(:,8)=(1D0-xsi)*(1D0-eta*eta)/2D0)
double precision  DNuR(nrgpd,nr_u) 
! derivate of shape functions for u with respect to xsi
!data  (dNuR(1:nrgpd-1:2,1)=-(-(1D0-eta)*(1D0+xsi+eta)+(1D0-xsi)*(1D0-eta))/4D0)
!data  (dNuR(1:nrgpd-1:2,2)=-( (1D0-eta)*(1D0-xsi+eta)-(1D0+xsi)*(1D0-eta))/4D0)
!data  (dNuR(1:nrgpd-1:2,3)=-( (1D0+eta)*(1D0-xsi-eta)-(1D0+xsi)*(1D0+eta))/4D0)
!data  (dNuR(1:nrgpd-1:2,4)=-(-(1D0+eta)*(1D0+xsi-eta)+(1D0-xsi)*(1D0+eta))/4D0)
!data  (dNuR(1:nrgpd-1:2,5)=-xsi*(1D0-eta))
!data  (dNuR(1:nrgpd-1:2,6)= (1D0-eta*eta)/2D0)
!data  (dNuR(1:nrgpd-1:2,7)=-xsi*(1D0+eta))
!data  (dNuR(1:nrgpd-1:2,8)=-(1D0-eta*eta)/2D0)

data  (dNuR(1:nrgpd-1:2,1)=(1d0/4d0*(1d0-eta))*(1d0+xsi+eta)-(1d0/4d0*(1d0-xsi))*(1d0-eta))
data  (dNuR(1:nrgpd-1:2,2)=-(1d0/4d0*(1d0-eta))*(1d0-xsi+eta)+(1d0/4d0*(1d0+xsi))*(1d0-eta))
data  (dNuR(1:nrgpd-1:2,3)=-(1d0/4d0*(1d0+eta))*(1d0-xsi-eta)+(1d0/4d0*(1d0+xsi))*(1d0+eta))
data  (dNuR(1:nrgpd-1:2,4)=(1d0/4d0*(1d0+eta))*(1d0+xsi-eta)-(1d0/4d0*(1d0-xsi))*(1d0+eta))
data  (dNuR(1:nrgpd-1:2,5)=-xsi*(1d0-eta))
data  (dNuR(1:nrgpd-1:2,6)=1d0/2d0-(1d0/2d0)*eta*eta)
data  (dNuR(1:nrgpd-1:2,7)=-xsi*(1d0+eta))
data  (dNuR(1:nrgpd-1:2,8)=-1d0/2d0+(1d0/2d0)*eta*eta)


! derivate of shape functions for u with respect to xsi
!data  (dNuR(2:nrgpd:2,1)=-(-(1D0-xsi)*(1D0+xsi+eta)+(1D0-xsi)*(1D0-eta))/4D0)
!data  (dNuR(2:nrgpd:2,2)=-(-(1D0+xsi)*(1D0-xsi+eta)+(1D0+xsi)*(1D0-eta))/4D0)
!data  (dNuR(2:nrgpd:2,3)=-( (1D0+xsi)*(1D0-xsi-eta)-(1D0+xsi)*(1D0+eta))/4D0)
!data  (dNuR(2:nrgpd:2,4)=-( (1D0-xsi)*(1D0+xsi-eta)-(1D0-xsi)*(1D0+eta))/4D0)
!data  (dNuR(2:nrgpd:2,5)=-(1D0-xsi*xsi)/2D0)
!data  (dNuR(2:nrgpd:2,6)=-eta*(1D0+xsi))
!data  (dNuR(2:nrgpd:2,7)= (1D0-xsi*xsi)/2D0)
!data  (dNuR(2:nrgpd:2,8)=-eta*(1D0-xsi))


data  (dNuR(2:nrgpd:2,1)=(1d0/4d0*(1d0-xsi))*(1d0+xsi+eta)-(1d0/4d0*(1d0-xsi))*(1d0-eta))
data  (dNuR(2:nrgpd:2,2)=(1d0/4d0*(1d0+xsi))*(1d0-xsi+eta)-(1d0/4d0*(1d0+xsi))*(1d0-eta))
data  (dNuR(2:nrgpd:2,3)=-(1d0/4d0*(1d0+xsi))*(1d0-xsi-eta)+(1d0/4d0*(1d0+xsi))*(1d0+eta))
data  (dNuR(2:nrgpd:2,4)=-(1d0/4d0*(1d0-xsi))*(1d0+xsi-eta)+(1d0/4d0*(1d0-xsi))*(1d0+eta))
data  (dNuR(2:nrgpd:2,5)=-1d0/2d0+(1d0/2d0)*xsi*xsi)
data  (dNuR(2:nrgpd:2,6)=-(1d0+xsi)*eta)
data  (dNuR(2:nrgpd:2,7)=1d0/2d0-(1d0/2d0)*xsi*xsi)
data  (dNuR(2:nrgpd:2,8)=-(1d0-xsi)*eta)


!-------------------------------------------------------------------------------
contains
	
! Taylor-Hood element Q2Q1
! Update Lagrangian setting

! Element stiffness
subroutine m2dpulq2q1_e(ke,coord,D,ed,es,ep)
  implicit none
  double precision                :: ke(:,:), coord(:,:)
  double precision                :: D(:,:), ed(:), es(:), ep(:)
  
	
end subroutine m2dpulq2q1_e

! calculates the total deformation gradient and the pressure gradient
subroutine m2dpulq2q1_d(dg,pg,dp,coord,ed)
  implicit none

  double precision                :: dg(:,:), pg(:), dp(:,:)
  double precision                :: coord(:,:), ed(:)

  double precision                :: JT(nrgpd,2), JT0(nrgpd,2)
  double precision                :: JTinv(2,2), JT0inv(2,2)
  double precision                :: dNux(2,8), dNpx(2,4)

  double precision                :: ucoord(2,8), He(4,16)
  integer                         :: igp

! current nodal coordinates
  ucoord(1,:)=coord(1,:)+ed([1,4,7,10,13,15,17,19])
  ucoord(2,:)=coord(2,:)+ed([2,5,8,11,14,16,18,20])

  JT=MATMUL(dNuR,transpose(ucoord))
  JT0=MATMUL(dNuR,transpose(coord)) 
   
  do igp=1,nr_gp
  
    call inv2(JTinv,JT(index(:,igp),:))
    call inv2(JT0inv,JT0(index(:,igp),:))

! displacements
    dNuX=matmul(JT0inv,dNuR(index(:,igp),:))
! pressure
    dNpx=matmul(JTinv,dNpR(index(:,igp),:))

! deformation gradient
    He=0D0
    He(1:2,1:16:2)=dNuX 
    He(3:4,2:16:2)=dNuX 

    dg(:,igp)=Ifm+matmul(He,ed(u_list))

! could be done as, speed should be checked
!   dg(1,igp)=dot_product(Dnux(igp,1),ed([1,4,7,10,13,15,17,19])))
!   dg(2,igp)=dot_product(Dnux(igp,1),ed([1,4,7,10,13,15,17,19])))
!   dg(3,igp)=dot_product(Dnux(igp,2),ed([2,5,8,11,14,16,18,20])))
!   dg(4,igp)=dot_product(Dnux(igp,2),ed([2,5,8,11,14,16,18,20])))

! pressure gradient
    pg(igp)=dot_product(Np(igp,:),ed(p_list))
    dp(:,igp)=matmul(dNpx,ed(p_list))
    
  end do
  
  return
end subroutine m2dpulq2q1_d
 

! internal force vector
subroutine m2dpulq2q1_f(fe,coord,ed,es,h,g)
  implicit none

  double precision                :: fe(:), coord(:,:),  ed(:)
  double precision                :: es(:,:), h(:,:), g(:)

  double precision                :: JT(nrgpd,2), JTinv(2,2), detJ
  double precision                :: dNux(2,8), dNpx(2,4), Bu(3,16)

  double precision                :: efu(16), efp(4)
  double precision                :: ucoord(2,8)

  integer                         :: igp

! corrdinates in deformed configuration
  ucoord(1,:)=coord(1,:)+ed([1,4,7,10,13,15,17,19])
  ucoord(2,:)=coord(2,:)+ed([2,5,8,11,14,16,18,20])

  JT=matmul(dNuR,transpose(ucoord))

  efu=0D0
  efp=0D0
  do igp=1,nr_gp
  
    detJ=det2(JT(index(:,igp),:)) 
    call inv2(JTinv,JT(index(:,igp),:))

! displacements
    dNuX=matmul(JTinv,dNuR(index(:,igp),:))
! pressure
    dNpx=matmul(JTinv,dNpR(index(:,igp),:))

! strain displacement matrix
    Bu=0D0
    Bu(1,1:16:2) =dNuX(1,:)
    Bu(2,2:16:2) =dNuX(2,:)
    Bu(3,1:16:2) =dNuX(2,:)
    Bu(3,2:16:2) =dNuX(1,:)

    efu=efu+matmul(transpose(Bu),es([1,2,4],igp))*detJ*w_1(igp)*w_2(igp)
    efp=efp+(matmul(transpose(dNpx),h(:,igp)) &
             +Np(igp,:)*g(igp))*detJ*w_1(igp)*w_2(igp)
             
  end do

! element force vector
  fe(u_list)=efu
  fe(p_list)=efp

  return
end subroutine m2dpulq2q1_f

! internal force vector
subroutine m2dpulq2q1_bn(fe,coord,ed,ea,rho)
  implicit none

  double precision                :: fe(:), coord(:,:),  ed(:)
  double precision                :: ea(:), rho

  double precision                :: detJ, JT(nrgpd,2)
  double precision                :: ag(16)

  double precision                :: efu(16), efp(4)
  double precision                :: ucoord(2,8), N2(2,16)

  integer                         :: igp

! coordinates in undeformed configuration
  ucoord(1,:)=coord(1,:)!+ed([1,4,7,10,13,15,17,19])
  ucoord(2,:)=coord(2,:)!+ed([2,5,8,11,14,16,18,20])

  JT=matmul(dNuR,transpose(ucoord))

  efu=0D0
  ag=ea(u_list)   ! mechanical dofs

  do igp=1,nr_gp
  
    detJ=det2(JT(index(:,igp),:)) 
   
    N2(1,:)=(/Nu(igp,1), 0d0, Nu(igp,2), 0d0, Nu(igp,3), 0d0, Nu(igp,4), &
              0d0, Nu(igp,5), 0d0, Nu(igp,6), 0d0, Nu(igp,7), 0d0, Nu(igp,8), 0d0 /)
    N2(2,:)=(/0d0, Nu(igp,1), 0d0, Nu(igp,2), 0d0, Nu(igp,3), 0d0, Nu(igp,4), &
              0d0, Nu(igp,5), 0d0, Nu(igp,6), 0d0, Nu(igp,7), 0d0, Nu(igp,8) /)
 
    efu=efu+matmul(transpose(N2),matmul(N2,ag))*detJ*w_1(igp)*w_2(igp)*rho

  end do

! element force vector
  fe(u_list)=efu
  fe(p_list)=0d0

  return
end subroutine m2dpulq2q1_bn


end module elem_large_mix_th_2d


