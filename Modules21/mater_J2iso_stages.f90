module mater_J2iso_stages

! last modified
! 
! M. Ristinmaa 2012-02-08
!  - copied from mater_J2_iso
!  - introducing stages, requires that variables are stored for the different
!    stages
!  - deformation gradient is stored in dgn when stress calculation is carried out
! ------------------------------------------------------------------------------

use mater_large
use some_constants

implicit none

double precision, allocatable       :: Isv(:,:,:,:)
double precision, allocatable       :: Isv_init(:)
double precision, allocatable       :: dgn(:,:,:,:)
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


interface dJ2iso
  module procedure dJ2_iso0
  module procedure dJ2_iso1
  module procedure dJ2_iso2
end interface
private dJ2_iso0, dJ2_iso1, dJ2_iso2

interface diJ2iso
  module procedure diJ2_iso1
  module procedure diJ2_iso2
end interface
private diJ2_iso1, diJ2_iso2


interface J2iso
  module procedure J2_iso0
  module procedure J2_iso1
  module procedure J2_iso2
end interface
private J2_iso0, J2_iso1, J2_iso2


!-------------------------------------------------------------------------------
contains


subroutine state_save(filename)
  implicit none
  character(len=*)                :: filename
! don't declear nvis global
    integer                         :: nr_gp, nelem, nstate, nisv

! Can be used to make binary file system-independent
!  INQUIRE(IOLENGTH=length) specimen, list, of, items. Require DIRECT access
  open (1, file = filename, access = 'sequential', status = 'REPLACE')

  nisv=size(Isv,1)
  nr_gp =size(Isv,2)
  nelem =size(Isv,3)
  nstate=size(Isv,4)

  write(1,*) nr_gp, nisv, nelem, nstate
  write(1,*) ep
  write(1,*) Isv
  write(1,*) dgn
  write(1,*) Isv_init

  close(1)

  return
end subroutine state_save


subroutine state_load(filename)
  implicit none
  character(len=*)                :: filename
! don't declear nvis global
  integer                         :: nr_gp,  nelem, nstate, ierr, nisv

! Can be used to make binary file system-independent
!  INQUIRE(IOLENGTH=length) specimen, list, of, items 
  open (1, file = filename, access = 'sequential', status = 'old')

  read(1,*) nr_gp, nisv, nelem, nstate

! Allocate structures
  allocate(dgn(9,nr_gp,nelem,nstate), stat=ierr)
  allocate(Isv(nisv,nr_gp,nelem,nstate), stat=ierr)
  allocate(Isv_init(nisv), stat=ierr)

  read(1,*) ep
  read(1,*) Isv
  read(1,*) dgn
  read(1,*) Isv_init

  close(1)

  return
end subroutine state_load

! Get data from material routine

function state_numbers()
  implicit none
  integer                         :: state_numbers

  state_numbers=size(Isv,4)

  return
end function state_numbers




subroutine state_copy(ist,istn)
  implicit none
  integer          :: ist,istn

! Copy the state
  Isv(:,:,:,ist)=Isv(:,:,:,istn)
  dgn(:,:,:,ist)=dgn(:,:,:,istn)

  return
end subroutine state_copy


subroutine J2iso_getValScalar(stype,out,stage)
  implicit none
  character(len=*)                :: stype
  double precision, allocatable   :: out(:,:)
  integer                         :: stage
 
  if (stype.eq.'alpha') then
    out=Isv(7,:,:,stage)
  elseif (stype.eq.'plastic') then
    out=0d0
    where(Isv(8,:,:,stage).gt.0d0) out=1d0   !dl
  else
    stop "Variable does not exist"
  end if
  
  return
end subroutine J2iso_getValScalar


subroutine J2iso_getValVector(stype,out,stage)
  implicit none
  character(len=*)                :: stype
  double precision, allocatable   :: out(:,:,:)
  integer                         :: stage

  
  if (stype.eq.'be') then
    out=Isv(1:6,:,:,stage)
  else
    stop "Variable does not exist"
  end if
  
  return
end subroutine J2iso_getValVector


subroutine init_J2_iso1(mp,Nsta,Nelm,Ngp)
  implicit none
  integer                         :: Nisv, Ngp, Nelm, Nsta, ierr
  double precision                :: mp(6)

  Nisv=9
  allocate(Isv(Nisv,Ngp,nelm,Nsta), stat=ierr)   
  allocate(Isv_init(Nisv), stat=ierr)   
  New=1
  Old=2
  ep=mp

  Isv=0D0
  Isv(1,:,:,:)=1D0
  Isv(2,:,:,:)=1D0
  Isv(3,:,:,:)=1D0

  Isv_init=0D0
  Isv_init(1)=1D0
  Isv_init(2)=1D0
  Isv_init(3)=1D0

  allocate(dgn(9,Ngp,nelm,Nsta), stat=ierr)
  dgn=0d0
  dgn((/1,5,9/),:,:,:)=1d0

  return
end subroutine init_J2_iso1


subroutine J2_iso0(stype,stress,ef_New,ie,stage)
  implicit none
  double precision                :: stress(:), ef_New(:)
  character(len=*)                :: stype
  integer                         :: ie, stage

  if (stage.eq.1) then
    call J2_iso2(stype, stress(:), Isv_init, Isv(:,1,ie,stage), & 
                 ep, ef_New(:),stage,ie,1)
  else
    call J2_iso2(stype, stress(:), Isv(:,1,ie,stage-1), Isv(:,1,ie,stage), & 
                 ep, ef_New(:),stage,ie,1)
  end if

  return
end subroutine J2_iso0


subroutine J2_iso1(stype,stress,ef_New,ie,stage)
  implicit none
  double precision                :: stress(:,:), ef_New(:,:)
  character(len=*)                :: stype
  integer                         :: gp, nr_gp,ie, stage

  nr_gp=size(ef_New,2)

  if (stage.eq.1) then
    do gp=1,nr_gp
      call J2_iso2(stype, stress(:,gp), Isv_init, Isv(:,gp,ie,stage), & 
                 ep, ef_New(:,gp),stage,ie,gp)
    enddo
  else
    do gp=1,nr_gp
      call J2_iso2(stype, stress(:,gp), Isv(:,gp,ie,stage-1), Isv(:,gp,ie,stage), & 
                 ep, ef_New(:,gp),stage,ie,gp)
    enddo
  end if

  return
end subroutine J2_iso1


subroutine J2_iso2(stype,stress,Isv_Old,Isv_New,ep,ef_New,stage,ie,gp)
  implicit none
  double precision                :: stress(:), Isv_Old(:), Isv_New(:)
  double precision                :: ep(:), ef_New(:)
  character(len=*)                :: stype  
  integer                         :: ie, gp, stage
  
  double precision                :: kappa, mu, sigy0, H, delta, y_oo
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
  double precision                :: s_trial(3,3), n(3,3), invF(3,3)

! First store the deformation gradient for this stage
  dgn(:,gp,ie,stage)=ef_New

! Recall material parameters
  kappa=ep(1)
  mu=ep(2)
  sigy0=ep(3)
  H=ep(4)
  delta=ep(5)
  y_oo=ep(6)

!
  F_New=getF(ef_New)
! Check if a inverse deformation gradient is provided as input
  if (stype.eq.'Cauchy-inv') then
    call inv3(invF,F_New)
    F_New=invF
    dgn(:,gp,ie,stage)=(/F_New(1,1),F_New(1,2),F_New(1,3), &
                         F_New(2,1),F_New(2,2),F_New(2,3), &
                         F_New(3,1),F_New(3,2),F_New(3,3)/)
  end if

  if (stage.ne.1) then
    F_Old=getF(dgn(:,gp,ie,stage-1))
  else
    F_Old=getI()
  end if

  Id   =getI()
  call inv3(iF_Old,F_Old)
  call inv3(iF_New,F_New)
  f_rel=MATMUL(F_new,iF_old)

!
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

  if ((stype.eq.'Cauchy').or.(stype.eq.'Cauchy-inv'))  then
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


subroutine dJ2_iso0(stype,D,ef_New,ie,stage)
  implicit none
  character(len=*)                :: stype
  double precision                :: D(:,:), ef_New(:)
  integer                         :: ie,stage

  if (stage.eq.1) then
    call dJ2_iso2(stype,D(:,:),ep,ef_New(:), &
                  Isv_init, Isv(:,1,ie,stage),stage,ie,1)
  else
    call dJ2_iso2(stype,D(:,:),ep,ef_New(:), &
                  Isv(:,1,ie,stage-1), Isv(:,1,ie,stage),stage,ie,1)
  end if

  return
end subroutine dJ2_iso0


subroutine dJ2_iso1(stype,D,ef_New,ie,stage)
  implicit none
  character(len=*)                :: stype
  double precision                :: D(:,:,:), ef_New(:,:)
  integer                         :: gp,nr_gp, ie, stage

  nr_gp=size(D,3)

  if (stage.eq.1) then
    do gp=1,nr_gp
      call dJ2_iso2(stype,D(:,:,gp),ep,ef_New(:,gp), &
                  Isv_init, Isv(:,gp,ie,stage),stage,ie,gp)
    enddo
  else
    do gp=1,nr_gp
      call dJ2_iso2(stype,D(:,:,gp),ep,ef_New(:,gp), &
                  Isv(:,gp,ie,stage-1), Isv(:,gp,ie,stage),stage,ie,gp)
    enddo
  end if

  return
end subroutine dJ2_iso1


subroutine dJ2_iso2(stype,Dout,ep,ef_New,Isv_Old,Isv_New,stage,ie,gp)
  implicit none
  character(len=*)                :: stype
  double precision                :: ep(:), ef_New(:)
  double precision                :: Isv_New(:), Isv_Old(:), Dout(:,:) 
  integer                         :: ie, gp, stage
  
  double precision                :: kappa, mu, sigy0, H, delta, y_oo, J 
  double precision                :: F(3,3), C(3,3), Id(3,3), s(3,3)
  double precision                :: invC(3,3), invF(3,3)
  double precision                :: alpha, dl, tr_be, mu_bar, tr_ss
  double precision                :: tr_be_bar_trial

  double precision                :: n1(3,3), n_times_n(3,3), n1_hat(3,3)
  double precision                :: temp(3,3), F_Old(3,3)
  double precision                :: N_matrix(1,6),  Unit(6,6)
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
  if (stage.ne.1) then
    F_Old=getF(dgn(:,gp,ie,stage-1))
  else
    F_old=getI()
  end if

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

 
subroutine diJ2_iso1(stype,Duu,es,ef,ie,ist)
  implicit none
  double precision                :: Duu(:,:,:)
  double precision                :: es(:,:)
  double precision                :: ef(:,:)
  character(len=*)                :: stype
  integer                         :: gp,nr_gp, ie, ist

!  double precision :: dgx(9,8), ddx(3,8), pert, esx(6,8)
 
  nr_gp=size(Duu,3)

  do gp=1,nr_gp
    call diJ2_iso2(stype,Duu(:,:,gp), &
                    es(:,gp),ef(:,gp),ie,gp,ist)
  enddo

!  pert=1d-8
!  do ii=1,9
!    dgx=ef(:,:)
!    dgx(ii,:)=dgx(ii,:)+pert
!    esx=es
!    call vis_update('Cauchy-inv',dgx,dt,ie,ist)
!    call neohooke_vis_elect('Cauchy-inv',esx,ddx,dgx,e(:,:),ie,ist)
!    do gp=1,8
!       duu(:,ii,gp)=(esx(:,gp)-es(:,gp))/pert
!    enddo
!  enddo

  return
end subroutine diJ2_iso1



subroutine diJ2_iso2(stype,Duu,es,ef,ie,gp,ist)
  implicit none
  double precision                :: Duu(:,:)
  double precision                :: es(:), ef(:)
  character(len=*)                :: stype
  integer                         :: ie, gp, ist
  
  double precision                :: f(3,3), fi(3,3)
  integer                         :: i, j, k, l, ii, jj
  
  double precision                :: Duux(6,6), duun(6,9)
  double precision                :: efp(9), s(3,3)

  integer, parameter              :: d_list(2,6)=[(/1,1/),(/2,2/),(/3,3/),&
                                                  (/1,2/),(/1,3/),(/2,3/)]
  integer, parameter              :: f_list(2,9)=[(/1,1/),(/1,2/),(/1,3/),&
                                                  (/2,1/),(/2,2/),(/2,3/),&
                                                  (/3,1/),(/3,2/),(/3,3/)]
  double precision :: pert, esx1(6),esx2(6),  efx(9)

! Note it is the inverse deformation gradient that is in the 
! call arguments
  if (stype.ne.'inv') stop

  if (size(ef).eq.4) then
    f(1,:)=(/ef(1), ef(2), zero/)
    f(2,:)=(/ef(3), ef(4), zero/)
    f(3,:)=(/ zero,  zero,  one/)
  else
    f(1,:)=(/ef(1), ef(2), ef(3)/)
    f(2,:)=(/ef(4), ef(5), ef(6)/)
    f(3,:)=(/ef(7), ef(8), ef(9)/)
  endif

  s(1,:)=(/es(1),es(4),es(5)/)
  s(2,:)=(/es(4),es(2),es(6)/)
  s(3,:)=(/es(5),es(6),es(3)/)

  call inv3(fi,f)

  do ii=1,6
    i=d_list(1,ii)
    j=d_list(2,ii)
    do jj=1,9
      k=f_list(1,jj)
      l=f_list(2,jj)
      duun(ii,jj)=s(i,j)*fi(l,k)-fi(i,k)*s(l,j)-s(i,l)*fi(j,k)
    end do
  end do  


  efp=(/fi(1,1),fi(1,2),fi(1,3),fi(2,1),fi(2,2),fi(2,3),fi(3,1),fi(3,2),fi(3,3)/)
  if (ist.eq.1) then
    call dJ2_iso2('ul',Duux,ep,efp,Isv_init,Isv(:,gp,ie,ist),ist,ie,gp)
  else
    call dJ2_iso2('ul',Duux,ep,efp,Isv(:,gp,ie,ist-1),Isv(:,gp,ie,ist),ist,ie,gp)
  end if

!  call dneohooke2('ul',Duux,Dup,Dpux,Dpp, &
!                    efp,e,intVar(gp,:,:,ie,ist),dt)

! Note that the matrices dup and dpp found from the updated Lagrange
! scheme are the same as for the inverse motion problem, only duu and dpu need
! to be evaluated

! Duu matrix for inverse motion
! order 11 12 13 21 22 23 31 32 33
  duun(:,1)=duun(:,1)-duux(:,1)*fi(1,1)-duux(:,4)*fi(2,1)-duux(:,5)*fi(3,1)
  duun(:,2)=duun(:,2)-duux(:,4)*fi(1,1)-duux(:,2)*fi(2,1)-duux(:,6)*fi(3,1)
  duun(:,3)=duun(:,3)-duux(:,5)*fi(1,1)-duux(:,6)*fi(2,1)-duux(:,3)*fi(3,1)

  duun(:,4)=duun(:,4)-duux(:,1)*fi(1,2)-duux(:,4)*fi(2,2)-duux(:,5)*fi(3,2)
  duun(:,5)=duun(:,5)-duux(:,4)*fi(1,2)-duux(:,2)*fi(2,2)-duux(:,6)*fi(3,2)
  duun(:,6)=duun(:,6)-duux(:,5)*fi(1,2)-duux(:,6)*fi(2,2)-duux(:,3)*fi(3,2)

  duun(:,7)=duun(:,7)-duux(:,1)*fi(1,3)-duux(:,4)*fi(2,3)-duux(:,5)*fi(3,3)
  duun(:,8)=duun(:,8)-duux(:,4)*fi(1,3)-duux(:,2)*fi(2,3)-duux(:,6)*fi(3,3)
  duun(:,9)=duun(:,9)-duux(:,5)*fi(1,3)-duux(:,6)*fi(2,3)-duux(:,3)*fi(3,3)


  duu=duun


! pertubation method
! This pertubation gives better convergence than tha
! above analytical tangent
! DEBUGGING REQUIRED
!
if (1.eq.0) then
! if (ie.eq.1.and.gp.eq.1)  write(*,*)'duu_ana ',duu
  pert=1e-10
  esx1=es
  do ii=1,9
    efx=ef
    efx(ii)=efx(ii)+pert    
!    if (ist.eq.1) then
!      call update2('Cauchy-inv',efx, dt, Cv_init, Cv_gp)
!    else
!      call update2('Cauchy-inv',efx, dt, intVar(gp,:,:,ie,ist-1), Cv_gp)
!    end if
!    call neohooke2('Cauchy-inv',esx2,Cv_gp,ddx,efx,e)
    duu(:,ii)=(esx2-esx1)/pert
  enddo
!if (ie.eq.1.and.gp.eq.1)   write(*,*)'duu_per ',duu
!if (ie.eq.1.and.gp.eq.1) pause
endif

  return
end subroutine diJ2_iso2









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




end module mater_J2iso_stages
