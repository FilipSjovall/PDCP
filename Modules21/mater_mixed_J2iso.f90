module mater_mixed_J2iso

! last modified
! M. Ristinmaa 2012-01-15
!  - initial implementation
! -----------------------------------------------------------------------------

!use mater_large
use matrix_util, only: det3, inv3s, inv3
use mater_large, only: pushforward
use mater_hyperel !, only: neohooke_init, neohooke
use mater_J2iso, only: J2iso_init, J2iso, J2iso_accept, J2iso_getVal

implicit none
private det3

double precision                  :: K_mp, G_mp, RT_mp, rhofo_mp
double precision                  :: Nso_mp, Ko_mp, kappa_mp, pf_init
private K_mp, G_mp, RT_mp, rhofo_mp, Nso_mp, Ko_mp, kappa_mp, pf_init

integer                           :: nelm_save, ngp_save
private nelm_save, ngp_save

double precision, allocatable     :: pgn(:,:), dgn(:,:,:)
private pgn, dgn

interface matermix_getval
  module procedure matermix_getval1
  module procedure matermix_getval2
end interface

private matermix_getval1, matermix_getval2

interface dmatermix
  module procedure dneohooke1_mix
  module procedure dneohooke2_mix
end interface

private dneohooke1_mix, dneohooke2_mix

interface matermix_hyper
  module procedure neohooke1_mix
  module procedure neohooke2_mix
end interface

private neohooke1_mix, neohooke2_mix

!------------------------------------------------------------------------------

contains
 


subroutine matermix_init(mp,nelm,ngp)
  implicit none
  double precision                :: mp(:)
  integer                         :: nelm, ngp

  integer                         :: ie, ig, ierr

  nelm_save=nelm
  ngp_save=ngp

  K_mp=mp(1)
  G_mp=mp(2)
!  call neohooke_init(mp([1,2]))

  RT_mp=mp(3)
  rhofo_mp=mp(4)
  Nso_mp=mp(5)
  Ko_mp=mp(6)
  kappa_mp=mp(7)

  allocate(pgn(ngp,nelm),stat=ierr)
  allocate(dgn(9,ngp,nelm),stat=ierr)
  
  pf_init=RT_mp*rhofo_mp

  pgn=0d0 !RT_mp*rhofo_mp

! Note that the full deformation gradient is stored even for the
! 2D case.
  dgn=0d0
  dgn(1,:,:)=1d0
  dgn(5,:,:)=1d0
  dgn(9,:,:)=1d0

! Initiate plasticity routine
  call J2iso_init(mp([1,2,8,9,10,11]),nelm,9)

  return
end subroutine matermix_init


subroutine accept(dg,pg)
  implicit none
  double precision                :: dg(:,:,:), pg(:,:)

  dgn([1,2,4,5],:,:)=dg
  pgn=pg

  call J2iso_accept(dg)

end subroutine accept



subroutine matermix_getval1(stype,val)
  implicit none
  character(len=*)                :: stype
  double precision, allocatable   :: val(:,:)

  double precision                :: F(3,3), J
  integer                         :: ie, ig

  val=0d0
  if (stype.eq.'ns') then

    do ie=1,nelm_save
      do ig=1,ngp_save
! deformation gradient
        F(1,:)=(/dgn(1,ig,ie), dgn(2,ig,ie), dgn(3,ig,ie)/)
        F(2,:)=(/dgn(4,ig,ie), dgn(5,ig,ie), dgn(6,ig,ie)/)
        F(3,:)=(/dgn(7,ig,ie), dgn(8,ig,ie), dgn(9,ig,ie)/)
        J=det3(F)
! volume fractions for solid
        val(ig,ie)=Nso_mp/J
      end do
    end do
  elseif (stype.eq.'pore_pressure') then
    val=pgn+pf_init   
  elseif (stype.eq.'plastic') then
    call J2iso_getVal('plastic',val)
  else
    stop 'Not implemented - get val'
  end if

  return
end subroutine matermix_getval1

subroutine matermix_getval2(stype,val)
  implicit none
  character(len=*)                :: stype
  double precision                :: val(:,:,:)

  double precision                :: F(3,3), J

  double precision, allocatable   :: dp(:,:), h(:,:)
  double precision                :: g, dt
  integer                         :: ie, ig, ierr

  val=0d0
  if (stype.eq.'Terzaghi stress') then

    if (size(val,1).eq.4) then !2D
      allocate(h(2,ngp_save),stat=ierr)
      allocate(dp(2,ngp_save),stat=ierr)
      h=0d0
      dp=0d0
    else
      stop 'please implement'
    end if
    dt=1d0;
! found by calling stress routine with zero pore pressure
    do ie=1,nelm_save
      do ig=1,ngp_save
        call neohooke2_mix('Cauchy',val(:,ig,ie),h(:,ig),g,dgn(:,ig,ie), &
                         0d0,dp(:,ig),dt,ig,ie)
      end do
    end do
  else
    stop 'Not implemented - get val'
  end if

  return
end subroutine matermix_getval2

! Call for several gauss points
subroutine neohooke1_mix(stype,stress,h,g,dg,pg,dp,dt,ie)
  implicit none
  double precision                :: stress(:,:), h(:,:), g(:)
  double precision                :: dg(:,:), pg(:), dp(:,:), dt
  integer                         :: ie
  character(len=*)                :: stype
  integer                         :: gp, nr_gp

  nr_gp=size(dg,2)

  call J2iso('Cauchy',stress,dg,ie)
  do gp=1,nr_gp
    call neohooke2_mix(stype,stress(:,gp),h(:,gp),g(gp),dg(:,gp),pg(gp),dp(:,gp), &
                   dt,gp,ie)
  enddo

  return
end subroutine neohooke1_mix

! Call for one gauss point
subroutine neohooke2_mix(stype,stress,h,g,dg,pg,dp,dt,gp,ie)
  implicit none
  double precision                :: stress(:), h(:), g, dg(:), pg, dp(:), dt
  double precision, allocatable   :: es(:,:)
  integer                         :: gp, ie, ierr
  character(len=*)                :: stype
  
  double precision                :: F(3,3), C(3,3), iC(3,3), id(3,3)
  double precision                :: J, trC, S(3,3)
  double precision                :: invF(3,3), Fn(3,3), dundx(3,3), dudx(3,3)
  double precision                :: ns, nf, K, pf, Jn, Kf, rhof
  
! NOTE WE NEED A DET3 ROUTINE WHICH TAKES dg DIRECTLY
! CHANGE MATRIX_UTIL

  if (size(dg).eq.4) then
    F(1,:)=(/dg(1), dg(2), 0d0/)
    F(2,:)=(/dg(3), dg(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)
    allocate(es(4,8), stat=ierr) 
  else
    F(1,:)=(/dg(1), dg(2), dg(3)/)
    F(2,:)=(/dg(4), dg(5), dg(6)/)
    F(3,:)=(/dg(7), dg(8), dg(9)/)
    allocate(es(6,8), stat=ierr) 
  endif
  J=det3(F)

! volume fractions
  ns=Nso_mp/J
  nf=1-ns

! total pore pressure
  pf=pg+pf_init

! fluid density
  rhof=pf/RT_mp

! fluid stiffness
  Kf=rhof*RT_mp
!write(*,*)'Kf ',Kf,pf,rhof

! Calculate effective stress
!  call J2iso('Cauchy',es,dg,ie)
!  stress=es(:,gp)    
!  call neohooke('Cauchy',stress,dg)

! Terzagi stress = Effective stress minus pore pressure
  stress=ns*stress
  stress(1)=stress(1)-pf
  stress(2)=stress(2)-pf
  stress(3)=stress(3)-pf

! interaction in linear momentum
  K=Ko_mp*((J-Nso_mp)/(1d0-Nso_mp))**kappa_mp
  h=nf*rhof*K*dp

! other part
  call inv3(invF,F)
  Fn(1,:)=(/dgn(1,gp,ie), dgn(2,gp,ie), dgn(3,gp,ie)/)
  Fn(2,:)=(/dgn(4,gp,ie), dgn(5,gp,ie), dgn(6,gp,ie)/)
  Fn(3,:)=(/dgn(7,gp,ie), dgn(8,gp,ie), dgn(9,gp,ie)/)
  Jn=det3(Fn)

! Darcys law, note sign due to element definition
  g=rhof*nf/Kf*(pg-pgn(gp,ie))/dt &
   +rhof*(J-Jn)/dt/J

  return
end subroutine neohooke2_mix


! Material tangent stiffness for a pure displacement formulation

subroutine dneohooke1_mix(stype,D,ef)
  implicit none
  double precision                :: D(:,:,:), ef(:,:)
  character(len=*)                :: stype
  integer                         :: gp,nr_gp

  nr_gp=size(D,3)

  do gp=1,nr_gp
    call dneohooke2_mix(stype,D(:,:,gp),ef(:,gp))
  enddo

  return
end subroutine dneohooke1_mix


subroutine dneohooke2_mix(stype,D,ef)
  implicit none
  double precision                :: D(:,:), ef(:)
  character(len=*)                :: stype
  
  double precision                :: E, v, K, G
  double precision                :: F(3,3), C(3,3), iC(3,3), Id(3,3)
  double precision                :: J, trC, a1, a2, a3, Ds(21)
  integer                         :: i, jj, kk, l, in(21,4), el
  
  K=K_mp
  G=G_mp

 
  if (size(ef).eq.4) then
    F(1,:)=(/ef(1), ef(2), 0d0/)
    F(2,:)=(/ef(3), ef(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)
  else
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
  endif
  
  C=matmul(transpose(F),F)
  call inv3s(iC,C)
  J=det3(F)
  Id=0d0
  Id(1,1)=1d0
  Id(2,2)=1d0
  Id(3,3)=1d0
  trC=C(1,1)+C(2,2)+C(3,3)

  a1=K*J**2d0+2d0/9d0*G*J**(-2d0/3d0)*trC
  a2=2d0*G/3d0*J**(-2d0/3d0)
  a3=G/3d0*J**(-2d0/3d0)*trC-K/2d0*(J**2d0-1d0)

  if (size(ef).eq.4) then
    in(1,:)=(/1, 1, 1, 1/)
    in(2,:)=(/1, 1, 2, 2/)
    in(3,:)=(/1, 1, 1, 2/)
    in(4,:)=(/2, 2, 2, 2/)
    in(5,:)=(/2, 2, 1, 2/)
    in(6,:)=(/1, 2, 1, 2/)    

    do el=1,6
      i=in(el,1)
      jj=in(el,2)
      kk=in(el,3)
      l=in(el,4)

      Ds(el)=a1*iC(i,jj)*iC(kk,l)-a2*(Id(i,jj)*iC(kk,l)+iC(i,jj)*Id(kk,l)) &
       +a3*(iC(i,kk)*iC(jj,l)+iC(i,l)*iC(jj,kk))
    enddo

    D(1,:)=(/Ds(1), Ds(2), Ds(3)/)
    D(2,:)=(/Ds(2), Ds(4), Ds(5)/)
    D(3,:)=(/Ds(3), Ds(5), Ds(6)/)

  else
    in(1,:) =(/1, 1, 1, 1/)
    in(2,:) =(/1, 1, 2, 2/)
    in(3,:) =(/1, 1, 3, 3/)
    in(4,:) =(/1, 1, 1, 2/)
    in(5,:) =(/1, 1, 1, 3/)
    in(6,:) =(/1, 1, 2, 3/)    
    in(7,:) =(/2, 2, 2, 2/)    
    in(8,:) =(/2, 2, 3, 3/)    
    in(9,:) =(/2, 2, 1, 2/)    
    in(10,:)=(/2, 2, 1, 3/)    
    in(11,:)=(/2, 2, 2, 3/)    
    in(12,:)=(/3, 3, 3, 3/)    
    in(13,:)=(/3, 3, 1, 2/)    
    in(14,:)=(/3, 3, 1, 3/)    
    in(15,:)=(/3, 3, 2, 3/)    
    in(16,:)=(/1, 2, 1, 2/)    
    in(17,:)=(/1, 2, 1, 3/)    
    in(18,:)=(/1, 2, 2, 3/)    
    in(19,:)=(/1, 3, 1, 3/)    
    in(20,:)=(/1, 3, 2, 3/)    
    in(21,:)=(/2, 3, 2, 3/)    

    do el=1,21
      i=in(el,1)
      jj=in(el,2)
      kk=in(el,3)
      l=in(el,4)

      Ds(el)=a1*iC(i,jj)*iC(kk,l)-a2*(Id(i,jj)*iC(kk,l)+iC(i,jj)*Id(kk,l)) &
       +a3*(iC(i,kk)*iC(jj,l)+iC(i,l)*iC(jj,kk))
    enddo

    D(1,:)=(/Ds(1),  Ds(2),  Ds(3),  Ds(4),  Ds(5),  Ds(6)/)
    D(2,:)=(/Ds(2),  Ds(7),  Ds(8),  Ds(9), Ds(10), Ds(11)/)
    D(3,:)=(/Ds(3),  Ds(8), Ds(12), Ds(13), Ds(14), Ds(15)/)
    D(4,:)=(/Ds(4),  Ds(9), Ds(13), Ds(16), Ds(17), Ds(18)/)
    D(5,:)=(/Ds(5), Ds(10), Ds(14), Ds(17), Ds(19), Ds(20)/)
    D(6,:)=(/Ds(6), Ds(11), Ds(15), Ds(18), Ds(20), Ds(21)/)


  endif 

  if (stype.eq.'ul') then
    call pushforward(D,D,ef,'j')
  elseif (stype.eq.'tl') then
  else
    stop 'stype not implemented'
  endif


end subroutine dneohooke2_mix


end module mater_mixed_J2iso
