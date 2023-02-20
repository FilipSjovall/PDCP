module mater_J2iso

! last modified
! M. Wallin 2011-05-16
!  - Finite strain isotropic hardening implemented
! M. Wallin 2011-05-17
!  - init_J2_iso added
!  - accept_J2_iso added
! M. Wallin 2011-05-17
!  - material parameters added to J2iso_init
! M. Wallin 2011-05-17
!  - functionality J2iso_getVal added
! M. Ristinmaa 2011-05-19
!  - Changed arguments J2iso_getVal J2iso_init
!  - Private declaration of global variables
!  - Bugg J2iso_accept, New and Old was declared in routine
!  - Added routine for elements with one gauss point	
! M. Ristinmaa 2011-09-14
!  - changed in J2iso_accept, update wrong
! M. Ristinmaa 2011-10-24
!  - calls used for mixed fields u/p implemented.
! ------------------------------------------------------------------------------
!   Ref. Simo and Miehe, Associative coupled thermoplasticity, 
!   CMAME (1992), 98, pp. 41-104

use mater_large

implicit none

double precision, allocatable       :: Isv(:,:,:,:)
double precision, allocatable       :: dgn(:,:,:)
integer                             :: New, Old
double precision                    :: ep(6)
private Isv, dgn, New, Old, ep
	

interface J2iso_getVal
  module procedure J2iso_getValScalar
  module procedure J2iso_getValVector
end interface
private J2iso_getValScalar, J2iso_getValVector


interface J2iso_init
  module procedure init_J2_iso1
end interface
private init_J2_iso1


interface J2iso_accept
  module procedure accept_J2_iso1
  module procedure accept_J2_iso2
  module procedure accept_J2_iso3
end interface
private accept_J2_iso1


interface dJ2iso
  module procedure dJ2_iso0
  module procedure dJ2_iso1
  module procedure dJ2_iso2
  module procedure dJ2_iso1_up
end interface
private dJ2_iso0, dJ2_iso1, dJ2_iso2


interface J2iso
  module procedure J2_iso0
  module procedure J2_iso1
  module procedure J2_iso2
  module procedure J2_iso1_up
end interface
private J2_iso0, J2_iso1, J2_iso2


!-------------------------------------------------------------------------------
contains


subroutine J2iso_getValScalar(stype,out)
  implicit none
  integer                         :: ierr
  character(len=*)                :: stype
  double precision, allocatable   :: out(:,:)
 
  if (stype.eq.'alpha') then
    out=Isv(7,:,:,new)
  elseif (stype.eq.'plastic') then
    out=0d0
    where(Isv(8,:,:,new).gt.0d0) out=1d0   !dl
  else
    stop "Variable does not exist"
  end if
  
  return
end subroutine J2iso_getValScalar


subroutine J2iso_getValVector(stype,out)
  implicit none
  integer                         :: ierr
  character(len=*)                :: stype
  double precision, allocatable   :: out(:,:,:)

  
  if (stype.eq.'be') then
    out=Isv(1:6,:,:,new)
  else
    stop "Variable does not exist"
  end if
  
  return
end subroutine J2iso_getValVector

subroutine init_J2_iso1(mp,Nelm,Ngp)
  implicit none
  integer                         :: Nisv, Ngp, Nelm, ierr
  double precision                :: mp(6)

  Nisv=9
  allocate(Isv(Nisv,Ngp,nelm,2), stat=ierr)   
  New=1
  Old=2
  ep=mp

  Isv=0D0
  Isv(1,:,:,:)=1D0
  Isv(2,:,:,:)=1D0
  Isv(3,:,:,:)=1D0

  allocate(dgn(9,Ngp,nelm), stat=ierr)
  dgn=0d0
  dgn((/1,5,9/),:,:)=1d0

  return
end subroutine init_J2_iso1


subroutine accept_J2_iso1
  implicit none

! This must be done otherwise the tangent stiffness
! in the first iteration is calculated at the n-1 state
  Isv(:,:,:,Old)=Isv(:,:,:,New)

  if (New.eq.1) then
    New=2
    Old=1
  else
    New=1
    Old=2
  endif

  return
end subroutine accept_J2_iso1

! one gauss point per element
subroutine accept_J2_iso2(dg)
  implicit none
  double precision                :: dg(:,:)

  double precision                :: ddg(9)
  integer                         :: ie, ig

! This must be done otherwise the tangent stiffness
! in the first iteration is calculated at the n-1 state
  Isv(:,:,:,Old)=Isv(:,:,:,New)

  if (New.eq.1) then
    New=2
    Old=1
  else
    New=1
    Old=2
  endif
 
  if (size(dg,1).eq.9) then
    dgn(:,1,:)=dg(:,:)
  else
    dgn(1,1,:)=dg(1,:)
    dgn(2,1,:)=dg(2,:)
    dgn(4,1,:)=dg(3,:)
    dgn(5,1,:)=dg(4,:)
  end if

  return
end subroutine accept_J2_iso2

! more than one gauss point per element
subroutine accept_J2_iso3(dg)
  implicit none
  double precision                :: dg(:,:,:)

  double precision                :: ddg(9)
  integer                         :: ie, ig

! This must be done otherwise the tangent stiffness
! in the first iteration is calculated at the n-1 state
  Isv(:,:,:,Old)=Isv(:,:,:,New)

  if (New.eq.1) then
    New=2
    Old=1
  else
    New=1
    Old=2
  endif

  if (size(dg,1).eq.9) then
    dgn=dg
  else
    dgn(1,:,:)=dg(1,:,:)
    dgn(2,:,:)=dg(2,:,:)
    dgn(4,:,:)=dg(3,:,:)
    dgn(5,:,:)=dg(4,:,:)
  end if

  return
end subroutine accept_J2_iso3


subroutine J2_iso1_up(stype,stress,ef_New,th,ie)
  implicit none
  double precision                :: stress(:,:), ef_New(:,:), th
  character(len=*)                :: stype

  double precision                :: p
  integer                         :: gp, nr_gp,ie

  nr_gp=size(ef_New,2)

  do gp=1,nr_gp
    call J2_iso2('Cauchy', stress(:,gp), Isv(:,gp,ie,old), Isv(:,gp,ie,new), & 
                 ep, ef_New(:,gp) ,ie,gp)
                 
    p=(stress(1,gp)+stress(2,gp)+stress(3,gp))/3d0
    stress(:,gp)=stress(:,gp)  &
               +(ep(1)/2d0*(th-1d0/th)-p)*(/1d0,1d0,1d0,0d0,0d0,0d0/)
  enddo

  return
end subroutine J2_iso1_up




subroutine J2_iso0(stype,stress,ef_New,ie)
  implicit none
  double precision                :: stress(:), ef_New(:)
  character(len=*)                :: stype
  integer                         :: gp, nr_gp,ie


  call J2_iso2(stype, stress(:), Isv(:,1,ie,old), Isv(:,1,ie,new), & 
                 ep, ef_New(:),ie,1)

  return
end subroutine J2_iso0


subroutine J2_iso1(stype,stress,ef_New,ie)
  implicit none
  double precision                :: stress(:,:), ef_New(:,:)
  character(len=*)                :: stype
  integer                         :: gp, nr_gp,ie

  nr_gp=size(ef_New,2)

  do gp=1,nr_gp
    call J2_iso2(stype, stress(:,gp), Isv(:,gp,ie,old), Isv(:,gp,ie,new), & 
                 ep, ef_New(:,gp),ie,gp)
  enddo

  return
end subroutine J2_iso1


subroutine J2_iso2(stype,stress,Isv_Old,Isv_New,ep,ef_New,ie,gp)
  implicit none
  double precision                :: stress(:), Isv_Old(:), Isv_New(:)
  double precision                :: ep(:), ef_New(:)
  character(len=*)                :: stype  
  integer                         :: ie, gp
  
  double precision                :: E, v, kappa, mu, sigy0, H, delta, y_oo
  double precision                :: Sy, tr_stst, f_trial
  double precision                :: dalpha_o, dalpha, f1, dSy, df
  double precision                :: J_new, Uptheta
  double precision                :: alpha_old, alpha_new, detf_rel
  double precision                :: tr_be_bar_trial, Ie_new, mu_bar, dl
  double precision                :: id(3,3), S(3,3), tmp33(3,3)
  double precision                :: F_Old(3,3), F_New(3,3), iF_Old(3,3)
  double precision                :: iF_New(3,3), f_rel(3,3), be_bar_old(3,3)
  double precision                :: f_rel_bar(3,3), temp(3,3)
  double precision                :: be_bar_trial(3,3), s_new(3,3)
  double precision                :: be_bar_new(3,3), Kirch_new(3,3)
  double precision                :: s_trial(3,3), n(3,3) 


  kappa=ep(1)
  mu=ep(2)
  sigy0=ep(3)
  H=ep(4)
  delta=ep(5)
  y_oo=ep(6)

  F_New=getF(ef_New)
  F_Old=getF(dgn(:,gp,ie))
  Id   =getI()
  call inv3(iF_Old,F_Old)
  call inv3(iF_New,F_New)
  f_rel=MATMUL(F_new,iF_old)

  detf_rel=det3(f_rel)
  f_rel_bar=detf_rel**(-1D0/3D0)*f_rel

  be_bar_old=RESHAPE((/Isv_Old(1),    Isv_Old(4),    Isv_Old(5),&
                       Isv_Old(4),    Isv_Old(2),    Isv_Old(6),&
                       Isv_Old(5),    Isv_Old(6),    Isv_Old(3)/), &
                        (/3,3/),ORDER=(/2,1/))
  
  alpha_old=Isv_Old(7)
 
  temp=MATMUL(f_rel_bar,be_bar_old)
  be_bar_trial=MATMUL(temp,TRANSPOSE(f_rel_bar))
		
  tr_be_bar_trial=tr3(be_bar_trial)

  s_trial=mu*dev3(be_bar_trial)
	 
  Sy=sigy0+h*alpha_old+y_oo*(1-exp(-delta*alpha_old))

  tmp33=matmul(s_trial,s_trial)
  tr_stst=tr3(tmp33)

  f_trial=SQRT(tr_stst)-SQRT(2D0/3D0)*Sy

  IF (f_trial<0D0) THEN
     s_new=s_trial
     alpha_new=alpha_old
     be_bar_new=be_bar_trial
     Ie_new=1D0/3D0*tr_be_bar_trial
     mu_bar=Ie_new*mu
     dl=0D0
  ELSE
     Ie_new=1.D0/3.D0*tr_be_bar_trial
     mu_bar=Ie_new*mu
	
     dalpha=0.5D-3	
     dalpha_o=1D-2	
		
     DO WHILE (ABS(dalpha-dalpha_o)>1D-15)
        dalpha_o=dalpha
        Sy=sigy0+h*(alpha_old+dalpha)+y_oo*(1D0-exp(-delta*(alpha_old+dalpha)))
        f1=sqrt(2.D0/3.D0)*Sy-SQRT(tr_stst) &
            +SQRT(2.D0/3.D0)*mu*dalpha*tr_be_bar_trial
        dSy=h+y_oo*delta*exp(-delta*(alpha_old+dalpha))
        df=SQRT(2.D0/3.D0)*dSy+SQRT(2.D0/3.D0)*mu*tr_be_bar_trial
        dalpha=dalpha_o-f1/df
     END DO
	
     dl=sqrt(3D0/2D0)*dalpha  			
     n=s_trial/SQRT(tr_stst)	
     s_new=s_trial-2D0*mu_bar*dl*n	
     alpha_new=alpha_old+sqrt(2D0/3D0)*dl	
		
  END IF	
	
  J_new=det3(F_new)	 
  Uptheta=kappa/2D0*(J_new**2D0-1D0)/J_new
	
  Kirch_new=J_new*Uptheta*id+s_new
  be_bar_new=s_new/mu+Ie_new*id	
		
  IsV_New(1:6)=(/be_bar_new(1,1),be_bar_new(2,2),be_bar_new(3,3), &
                 be_bar_new(1,2),be_bar_new(1,3),be_bar_new(2,3)/)
  Isv_New(7)=alpha_new
  Isv_New(8)=dl

  if (stype.eq.'Cauchy') then
    S=Kirch_new/J_New
  elseif (stype.eq.'2ndPiola') then
    S=MATMUL(MATMUL(iF_new,Kirch_new),TRANSPOSE(iF_new))
  else
    stop 'stype not implemented'
  endif

  if (size(ef_New).eq.4) then
    stress=(/S(1,1),S(2,2),S(3,3),S(1,2)/)
  else
    stress=(/S(1,1),S(2,2),S(3,3),S(1,2),S(1,3),S(2,3)/)
  endif
   
  return
end subroutine J2_iso2


subroutine dJ2_iso1_up(stype,D,ddwdjj,ef_New,th,ie)
  implicit none
  character(len=*)                :: stype
  double precision                :: D(:,:,:), ef_New(:,:), th, ddwdjj(:)

  double precision                :: F_m(3,3), Ja, Id(3,3)
  double precision                :: ddudjj, dudj, tmp1, tmp2
  integer                         :: gp,nr_gp, ie, ii, jj, i, j, k, l
  integer, parameter              :: d_list(2,6)=[(/1,1/),(/2,2/),(/3,3/),&
                                                  (/1,2/),(/1,3/),(/2,3/)]

  double precision :: hD(6,6), hDth

  nr_gp=size(D,3)
  Id   =getI()

  do gp=1,nr_gp
    F_m=getF(ef_New(:,gp))
    Ja=det3(F_m)
    call dJ2_iso2('ul',D(:,:,gp),ep,ef_New(:,gp), &
                  Isv(:,gp,ie,Old), Isv(:,gp,ie,New),ie,gp)

    ddudjj=ep(1)/2d0*(1d0+1d0/Ja**2d0)
    dudj=ep(1)/2D0*(Ja-1D0/Ja)
    tmp1=dudj
    tmp2=ddudjj*Ja+dudj
    do ii=1,6
      i=d_list(1,ii)
      j=d_list(2,ii)
      do jj=1,6
        k=d_list(1,jj)
        l=d_list(2,jj)
        D(ii,jj,gp)=D(ii,jj,gp)-tmp2*Id(i,j)*Id(k,l) &
                   +tmp1*(Id(i,k)*Id(j,l)+Id(i,l)*Id(j,k))
      end do
    end do  

    ddwdjj(gp)=ep(1)/2d0*(1d0+1d0/th**2d0)

  enddo


  return
end subroutine dJ2_iso1_up


subroutine dJ2_iso0(stype,D,ef_New,ie)
  implicit none
  character(len=*)                :: stype
  double precision                :: D(:,:), ef_New(:)
  integer                         :: gp,nr_gp, ie

  call dJ2_iso2(stype,D(:,:),ep,ef_New(:), &
                  Isv(:,1,ie,Old), Isv(:,1,ie,New),ie,1)

  return
end subroutine dJ2_iso0


subroutine dJ2_iso1(stype,D,ef_New,ie)
  implicit none
  character(len=*)                :: stype
  double precision                :: D(:,:,:), ef_New(:,:)
  integer                         :: gp,nr_gp, ie

  nr_gp=size(D,3)

  do gp=1,nr_gp
    call dJ2_iso2(stype,D(:,:,gp),ep,ef_New(:,gp), &
                  Isv(:,gp,ie,Old), Isv(:,gp,ie,New),ie,gp)
  enddo

  return
end subroutine dJ2_iso1


subroutine dJ2_iso2(stype,Dout,ep,ef_New,Isv_Old,Isv_New,ie,gp)
  implicit none
  character(len=*)                :: stype
  double precision                :: ep(:), ef_New(:)
  double precision                :: Isv_New(:), Isv_Old(:), Dout(:,:) 
  integer                         :: ie, gp
  
  double precision                :: E, v, kappa, mu, sigy0, H, delta, y_oo, J 
  double precision                :: F(3,3), C(3,3), Id(3,3), s(3,3)
  double precision                :: invC(3,3), invF(3,3)
  double precision                :: alpha, dl, tr_be, mu_bar, tr_ss
  double precision                :: tr_be_bar_trial

  double precision                :: n1(3,3), n_times_n(3,3), n1_hat(3,3)
  double precision                :: temp(3,3), F_Old(3,3)
  double precision                :: N_matrix(1,6),  Unit(6,6), D_ref(6,6)
  double precision                :: Kron(1,6), be(3,3), tmp33(3,3)
  double precision                :: N_dyad_N(6,6)
  double precision                :: hprime, beta0, beta1, beta2, beta3
  double precision                :: beta4, tr_stst
  double precision                :: f_rel(3,3), detf_rel, invF_old(3,3)
  double precision                :: f_rel_bar(3,3), be_bar_old(3,3)
  double precision                :: be_bar_trial(3,3), strial(3,3)
  double precision                :: Final(6,6), D(6,6)
  double precision                :: one_dyad_one(6,6), N_dyad_Kron(6,6)
  double precision                :: C_bar(6,6), C1(6,6), N_hat_matrix(1,6)

  double precision                :: dg(size(ef_New))
  

  kappa=ep(1)
  mu=ep(2)
  sigy0=ep(3)
  H=ep(4)
  delta=ep(5)
  y_oo=ep(6)


  dg=ef_New
  F    =getF(ef_New)
  F_Old=getF(dgn(:,gp,ie))
  Id   =getI()

  J  =det3(F)
  C=MATMUL(TRANSPOSE(F),F)   
  call inv3s(invC,C)
  call inv3(invF,F)
  call inv3(invF_old,F_Old)

  f_rel=MATMUL(F,invF_old)

  detf_rel=det3(f_rel)
  f_rel_bar=detf_rel**(-1D0/3D0)*f_rel

  be_bar_old=RESHAPE((/Isv_Old(1),    Isv_Old(4),    Isv_Old(5),&
                       Isv_Old(4),    Isv_Old(2),    Isv_Old(6),&
                       Isv_Old(5),    Isv_Old(6),    Isv_Old(3)/), &
                        (/3,3/),ORDER=(/2,1/))
  temp=MATMUL(f_rel_bar,be_bar_old)
  be_bar_trial=MATMUL(temp,TRANSPOSE(f_rel_bar))
		
  tr_be_bar_trial=tr3(be_bar_trial)

  strial=mu*dev3(be_bar_trial)

  be=RESHAPE((/Isv_New(1),    Isv_New(4),    Isv_New(5),&
               Isv_New(4),    Isv_New(2),    Isv_New(6),&
               Isv_New(5),    Isv_New(6),    Isv_New(3)/), &
                (/3,3/),ORDER=(/2,1/))
 
  alpha=Isv_New(7)
  dl   =Isv_New(8)

  tr_be=tr3(be)
  mu_bar=mu*tr_be/3D0

  s     =mu*dev3(be)
  
  tmp33=matmul(s,s)
  tr_ss=tr3(tmp33) 
   
  tmp33=matmul(strial,strial)
  tr_stst=tr3(tmp33)

  IF (sqrt(tr_ss)==0) THEN 
  	n1=id
  ELSE
  	n1=s/sqrt(tr_ss)
  END IF

  n_times_n=MATMUL(n1,n1)
  n1_hat=dev3(n_times_n)

  Unit=0D0
  Unit(1,1)=1D0
  Unit(2,2)=1D0
  Unit(3,3)=1D0
  Unit(4,4)=.5D0
  Unit(5,5)=.5D0
  Unit(6,6)=.5D0
 
  Kron=RESHAPE((/1D0,  1D0,  1D0,  0D0,   0D0,  0D0 /),(/1,6/))
  one_dyad_one=0D0
  one_dyad_one(1:3,1:3)=1D0

  N_matrix=RESHAPE((/n1(1,1), n1(2,2),  n1(3,3),    n1(1,2), &
                    n1(1,3),     n1(2,3) /),(/1,6/))
  N_hat_matrix=RESHAPE((/n1_hat(1,1), n1_hat(2,2), n1_hat(3,3), &
                      n1_hat(1,2),  n1_hat(1,3),  n1_hat(2,3) /),(/1,6/))
   
  N_dyad_N=MATMUL(TRANSPOSE(N_matrix),N_matrix)

  N_dyad_Kron=MATMUL(TRANSPOSE(N_matrix),Kron)+MATMUL(TRANSPOSE(Kron),N_matrix)

  Final=MATMUL(TRANSPOSE(N_matrix),N_hat_matrix)   

  C_bar= 2D0*mu_bar*(Unit-1D0/3D0*one_dyad_one) &
        -2D0/3D0*sqrt(tr_stst)*N_dyad_Kron
  C1=kappa*J**2D0*one_dyad_one-kappa*(J**2D0-1D0)*Unit+C_bar 

  IF (sqrt(tr_stst)>0) THEN
     hprime=h+y_oo*delta*exp(-delta*(alpha))
     beta0=1D0+hprime/(3D0*mu_bar)
     beta1=2D0*mu_bar*dl/sqrt(tr_stst)
     beta2=(1D0-1D0/beta0)*2D0/3D0*sqrt(tr_stst)/mu_bar*dl
     beta3=1D0/beta0-beta1+beta2
     beta4=(1D0/beta0-beta1)*sqrt(tr_stst)/mu_bar
  END IF

  IF (dl>0D0) THEN
     D=C1-beta1*C_bar-2D0*mu_bar*beta3*N_dyad_N-2D0*mu_bar*beta4*Final
  ELSE
     D=C1
  END IF

  if (size(ef_New).eq.4) then
    Dout(1,:)=(/D(1,1), D(2,2), D(1,4)/)
    Dout(2,:)=(/D(2,1), D(2,2), D(2,4)/)
    Dout(3,:)=(/D(4,1), D(4,2), D(4,4)/)
  else
    Dout=D
  endif 

  if (stype.eq.'ul') then
    Dout=Dout/J
  elseif (stype.eq.'tl') then
    call pullback(Dout,Dout,ef_new)
  else
    stop 'stype not implemented'
  endif


end subroutine dJ2_iso2

function det3(F)
  implicit none
  double precision                 :: F(:,:), det3

   det3= F(1,1)*F(2,2)*F(3,3)-F(1,1)*F(2,3)*F(3,2)  &
        -F(2,1)*F(1,2)*F(3,3)+F(2,1)*F(1,3)*F(3,2)  &
        +F(3,1)*F(1,2)*F(2,3)-F(3,1)*F(1,3)*F(2,2)
  return
end function det3
  




function getI
  implicit none
  double precision                :: getI(3,3)
 
  getI(1,:)=(/1D0, 0D0, 0D0/)
  getI(2,:)=(/0D0, 1D0, 0D0/)
  getI(3,:)=(/0D0, 0D0, 1d0/)

  return
end function getI




end module mater_J2iso
