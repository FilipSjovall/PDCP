module elem_flow_3d

! Last rev.:
!  M. Ristinmaa 2014-04-22
!   - initial version 8-node brick element derived from 
!     elem_large_cont_3d 

!   Module elem_flow_3d contains element subroutines 
!   diffusion where also large deformation of continuum can be considered

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

!------------------------------------------------------------------------------

interface fl3d8_e
  module procedure fl3d8_e1
  module procedure fl3d8_e2
end interface
private fl3d8_e1, fl3d8_e2

interface fl3d8_m
  module procedure fl3d8_m1
  module procedure fl3d8_m2
end interface
private fl3d8_m1, fl3d8_m2

interface fl3d8_vi
  module procedure fl3d8_vol1
  module procedure fl3d8_vol2
  module procedure fl3d8_vol3
  module procedure fl3d8_vol4
end interface
private fl3d8_vol1, fl3d8_vol2, fl3d8_vol3, fl3d8_vol4

interface fl3d8_gp
  module procedure fl3d8_gpval
end interface

interface fl3d8_bf
  module procedure fl3d8_bf1
  module procedure fl3d8_bf1x
  module procedure fl3d8_bf1a
  module procedure fl3d8_bf2
  module procedure fl3d8_bf3
end interface
private fl3d8_bf1, fl3d8_bf1a, fl3d8_bf2, fl3d8_bf3

!!!!!!! ADDED BY ANNA 2019-01-04 !!!!!!!!!!!!
interface fl3d8_volPres
  module procedure fun
  module procedure etafun
end interface
private fun, etafun

!------------------------------------------------------------------------------
contains


! WEIGHTED BTB-MATRIX
subroutine fl3d8_e1(ke,coord,ef,scaleval)
  implicit none
  double precision                :: ke(:,:), coord(:,:)
  double precision                :: ef(:,:), scaleval

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: DNX(3,8)
  DOUBLE PRECISION                :: F(3,3), C(3,3), iC(3,3), Jvol

  JT=MATMUL(DNR,transpose(COORD))

  KE=0D0
  DO GP_NR=1,NGP

    F(1,:)=(/ef(1,GP_NR), ef(2,GP_NR), ef(3,GP_NR)/)
    F(2,:)=(/ef(4,GP_NR), ef(5,GP_NR), ef(6,GP_NR)/)
    F(3,:)=(/ef(7,GP_NR), ef(8,GP_NR), ef(9,GP_NR)/)

    Jvol=det3(F)   
    C=matmul(transpose(F),F)
    call inv3(iC,C)

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 

    KE=KE+Jvol*MATMUL(TRANSPOSE(dNx),MATMUL(IC,dNx))*DETJ*scaleval

  END DO
  RETURN
END subroutine fl3d8_e1


subroutine fl3d8_e2(ke,coord,scaleval)
  implicit none
  double precision                :: ke(:,:), coord(:,:)
  double precision                :: scaleval

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: DNX(3,8)

  JT=MATMUL(DNR,transpose(COORD))

  KE=0D0
  DO GP_NR=1,NGP

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 

    KE=KE+MATMUL(TRANSPOSE(dNx),dNx)*DETJ*scaleval

  END DO
  RETURN
END subroutine fl3d8_e2


! MASS MATRIX
subroutine fl3d8_m1(Me,coord,ef,val)
  implicit none
  double precision                :: Me(:,:), coord(:,:)
  double precision                :: ef(:,:), val

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8)
  DOUBLE PRECISION                :: F(3,3), Jvol, V(8,1)

  JT=MATMUL(DNR,transpose(COORD))

  ME=0D0
  DO GP_NR=1,NGP
    F(1,:)=(/ef(1,GP_NR), ef(2,GP_NR), ef(3,GP_NR)/)
    F(2,:)=(/ef(4,GP_NR), ef(5,GP_NR), ef(6,GP_NR)/)
    F(3,:)=(/ef(7,GP_NR), ef(8,GP_NR), ef(9,GP_NR)/)

    Jvol=det3(F)      

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 
    
    V(:,1)=NR(GP_NR,:)

    ME=ME+Jvol* MATMUL(V,TRANSPOSE(V) )*DETJ*val
   
  END DO
  RETURN
END subroutine fl3d8_m1


subroutine fl3d8_m2(Me,coord,val)
  implicit none
  double precision                :: Me(:,:), coord(:,:)
  double precision                :: val

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8)
  DOUBLE PRECISION                :: V(8,1)

  JT=MATMUL(DNR,transpose(COORD))

  ME=0D0
  DO GP_NR=1,NGP

    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 
    
    V(:,1)=NR(GP_NR,:)

    ME=ME+MATMUL(V,TRANSPOSE(V) )*DETJ*val
   
  END DO
  RETURN
END subroutine fl3d8_m2


subroutine ShapeFun(ShapeF)
  implicit none
! Ship the shape functions to the main program 
  integer, parameter              :: NGP=8 
  integer                         :: GP_NR
  double precision                :: ShapeF(8,8)
  
  DO GP_NR=1,NGP
    ShapeF(GP_NR,:)=NR(GP_NR,:)
  END DO
  RETURN
END subroutine ShapeFun



! body force vector
subroutine fl3d8_bf1(fe,coord,gpval1)
  implicit none
  double precision                :: coord(:,:), fe(:)
  double precision                :: gpval1(:)

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8)

  JT=MATMUL(DNR,transpose(COORD))

  FE=0D0
  DO GP_NR=1,NGP
    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 
  
    FE=FE+NR(GP_NR,:)*DETJ*gpval1(gp_nr)

  END DO
  RETURN
END subroutine fl3d8_bf1


! body force vector
subroutine fl3d8_bf1x(fe,coord,gpval1,MrX)
  implicit none
  double precision                :: coord(:,:), fe(:)
  double precision                :: gpval1(:)

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR,ie,MrX

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8)



  JT=MATMUL(DNR,transpose(COORD))

  FE=0D0
  DO GP_NR=1,NGP
    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:))   
    FE=FE+DETJ*gpval1(gp_nr)

  END DO
  RETURN
END subroutine fl3d8_bf1x



subroutine fl3d8_bf1a(fe,coord,val)
  implicit none
  double precision                :: coord(:,:), fe(:)
  double precision                :: val

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8)

  JT=MATMUL(DNR,transpose(COORD))

  FE=0D0
  DO GP_NR=1,NGP
    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 
  
    FE=FE+NR(GP_NR,:)*DETJ*val

  END DO
  RETURN
END subroutine fl3d8_bf1a



subroutine fl3d8_bf2(fe,coord,gpval1,gpval2)
  implicit none
  double precision                :: coord(:,:), fe(:)
  double precision                :: gpval1(:), gpval2(:)

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8)

  JT=MATMUL(DNR,transpose(COORD))

  FE=0D0
  DO GP_NR=1,NGP
    CALL INV3(JTinv,JT(INDEX(:,gp_nr),:))
    DNX=MATMUL(JTinv,DNR(INDEX(:,gp_nr),:))
    DETJ=det3(JT(index(:,gp_nr),:)) 
  
    FE=FE+NR(GP_NR,:)*DETJ*gpval1(gp_nr)*gpval2(gp_nr)

  END DO
  RETURN
END subroutine fl3d8_bf2


subroutine fl3d8_bf3(fe,coord,ef,val)
  implicit none
  double precision                :: coord(:,:), fe(:)
  double precision                :: ef(:,:), val

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR,ie

  double precision                :: JT(R2,3), DetJ, JTinv(3,3)
  double precision                :: TMP(3,8), DNX(3,8)
  DOUBLE PRECISION                :: F(3,3), Jvol, V(8,1)

  JT=MATMUL(DNR,transpose(COORD))

  FE=0D0
  DO GP_NR=1,NGP
    F(1,:)=(/ef(1,GP_NR), ef(2,GP_NR), ef(3,GP_NR)/)
    F(2,:)=(/ef(4,GP_NR), ef(5,GP_NR), ef(6,GP_NR)/)
    F(3,:)=(/ef(7,GP_NR), ef(8,GP_NR), ef(9,GP_NR)/)

    Jvol=det3(F)      
    DETJ=det3(JT(index(:,gp_nr),:)) 

    FE=FE+Jvol*NR(GP_NR,:)*DETJ*val
   
  END DO

 
  RETURN
END subroutine fl3d8_bf3


! Extract gauss point values from nodal quantities
subroutine fl3d8_gpval(rhogp,rho)
  implicit none

  double precision                ::  rho(:)
  double precision                ::  rhogp(8)
 
  rhogp=MATMUL(NR,rho)

  RETURN
END subroutine fl3d8_gpval


! Calculate volume intergral
subroutine fl3d8_vol1(evol,coord)
  implicit none
  double precision                :: evol, coord(:,:)

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR

  double precision                :: JT(R2,3), DetJ

  JT=MATMUL(DNR,transpose(COORD))

  evol=0D0
  DO GP_NR=1,NGP
    DETJ=det3(JT(index(:,gp_nr),:)) 
    evol=evol+DETJ
  END DO
  RETURN
END subroutine fl3d8_vol1


subroutine fl3d8_vol2(evol,coord,efin)
  implicit none
  double precision                :: evol, coord(:,:), efin(:,:), ef(9)

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR

  double precision                :: JT(R2,3), DetJ, F(3,3), J

  JT=MATMUL(DNR,transpose(COORD))

  evol=0D0
  DO GP_NR=1,NGP

    ef=efin(:,gp_nr)
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
    J=det3(F)

    DETJ=det3(JT(index(:,gp_nr),:)) 
    evol=evol+DETJ*J
  END DO
  RETURN
END subroutine fl3d8_vol2


subroutine fl3d8_vol3(evol,coord,efin,gpval)
  implicit none
  double precision                :: evol, coord(:,:), efin(:,:), gpval(:)

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR

  double precision                :: JT(R2,3), DetJ, F(3,3), J, ef(9)

  JT=MATMUL(DNR,transpose(COORD))

  evol=0D0
  DO GP_NR=1,NGP

    ef=efin(:,gp_nr)
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
    J=det3(F)

    DETJ=det3(JT(index(:,gp_nr),:)) 
    evol=evol+DETJ*J*gpval(GP_NR)
  END DO
  RETURN
END subroutine fl3d8_vol3


subroutine fl3d8_vol4(evol,coord,gpval)
  implicit none
  double precision                :: evol, coord(:,:), gpval(:)

  integer, parameter              :: NGP=8, R2=NGP*3
  integer                         :: GP_NR

  double precision                :: JT(R2,3), DetJ

  JT=MATMUL(DNR,transpose(COORD))

  evol=0D0
  DO GP_NR=1,NGP

    DETJ=det3(JT(index(:,gp_nr),:)) 
    evol=evol+DETJ*gpval(GP_NR)
  END DO
  RETURN
END subroutine fl3d8_vol4


subroutine etafun(eta,x1,x2,edrho,coord,enod,B,pdeFilter,rhotilde,rho)
  implicit none
  double precision                 ::  x1, x2, pp, error, edrho(:,:), B, coord(:,:), output1, output2, eta, prod, rhotilde(:), rho(:)
  integer 								  ::  enod(:,:)
  logical								  ::  pdeFilter

  pp = (x1 + x2)/2d0
  call fun(output1,B,edrho,coord,enod,pp,pdeFilter,rhotilde,rho)
  error = abs(output1)
  
  do while (error.gt.1D-7)
  		call fun(output1,B,edrho,coord,enod,x1,pdeFilter,rhotilde,rho)
  		call fun(output2,B,edrho,coord,enod,pp,pdeFilter,rhotilde,rho)
  		prod = output1*output2
    	if (prod.lt.0d0) then
    		x2 = pp
    	else
       	x1 = pp       
    	endif
      pp = (x1 + x2)/2d0 
      call fun(output1,B,edrho,coord,enod,pp,pdeFilter,rhotilde,rho)
      error = abs(output1)
      !write(*,*) 'Error: ', error
  enddo
  eta = pp

  return
end subroutine etafun

subroutine fun(output,B,edrho,coord,enod,eta,pdeFilter,rhotilde,rho)
  implicit none
  double precision                 ::  B, eta, rhogp(8), Hfun, output, Scalar(8), edrho(:,:), tmp, coord(:,:), rhotilde(:), rho(:)
  integer								  ::  ngp, igp, nelm, ie, enod(:,:)
  logical								  ::  pdeFilter
  
  ngp = 8
  Scalar = 0d0
  nelm = size(edrho,2)
  output = 0d0

  do ie=1,nelm   
  	  if (pdeFilter) then
  	  		call fl3d8_gp(rhogp,edrho(:,ie))
  	  else
  	  		rhogp = rhotilde(ie)
  	  endif
     !!!!!!! IS IT CORRECT???!!!!!!!!!
     do igp=1,ngp
         Scalar(igp)= Hfun(rhogp(igp),B,eta) - rho(ie)!rhogp(igp)   
     enddo
	  call fl3d8_vi(tmp,coord(:,enod(:,ie)),Scalar)
	  output = output + tmp
  enddo

  return
end subroutine fun

end module elem_flow_3d
