module mater_hyperel

! last modified
! M. Ristinmaa 2011-05-06
!  - hyper elastic material implemented
! M. Ristinmaa 2011-05-06
!  - pushforward and pullback implemented
!  - upddefgr implemented
! M. Ristinmaa 2011-05-06
!  - pushforward, pullback and updefgr debugged for 2-d
!  - 3-d parts missing in pushforward and pullback for matrices
! M. Ristinmaa 2011-05-06
!  - stype introduced
! M. Ristinmaa 2011-05-16
!  - introduced new module hyper elastic material derived from mater_large
! M. Ristinmaa 2011-05-23
!  - Initiation of material parameters
! M. Ristinmaa 2011-09-07
!  - Introduced isochoric volumetric split for use with mixed elements
! M. Ristinmaa 2011-09-12
!  - Bugg dneohooke2 - missed pushforward in 3D
!  - Debugged neohooeXup and dneohookeXup
! M. Ristinmaa 2011-09-12
!  - Changed order stress and strain components
!    [11 22 33 12 23 13]
! ------------------------------------------------------------------------------


use mater_large
use matrix_util

implicit none

double precision                  :: Kmp, Gmp
private Kmp, Gmp

interface dneohooke
  module procedure dneohooke1
  module procedure dneohooke1gamma
  module procedure dneohooke2
  module procedure dneohooke2gamma
  module procedure dneohooke1up
  module procedure dneohooke2up
  module procedure dneohooke3up
end interface

private dneohooke1, dneohooke2


interface dneohookesens
  module procedure dneohooke1sens
  module procedure dneohooke2sens
end interface

private dneohooke1sens, dneohooke2sens


interface dfat
  module procedure dfat1
  module procedure dfat2
end interface

private dfat1, dfat2


interface neohooke
  module procedure neohooke1
  module procedure neohooke1gamma
  module procedure neohooke2
  module procedure neohooke2gamma
  module procedure neohooke1up
  module procedure neohooke2up
  module procedure neohooke3up
  module procedure neohooke1energy
  module procedure neohooke2energy
end interface

private neohooke1, neohooke2
private neohooke1energy, neohooke2energy
!private inv3, inv3s, det3

!-------------------------------------------------------------------------------

contains

subroutine neohooke_init(mp)
  implicit none
  double precision                :: mp(:)
  Kmp=mp(1)
  Gmp=mp(2)

  return
end subroutine neohooke_init

! Material routines for a pure displacement formulation

! Call for several gauss points
subroutine neohooke1(stype,stress,ef)
  implicit none
  double precision                :: stress(:,:), ef(:,:)
  character(len=*)                :: stype
  integer                         :: gp, nr_gp

  nr_gp=size(ef,2)

  do gp=1,nr_gp
    call neohooke2(stype,stress(:,gp), ef(:,gp))
  enddo

  return
end subroutine neohooke1


subroutine neohooke1gamma(stype,stress,ef,efgamma)
  implicit none
  double precision                :: stress(:,:), ef(:,:), efgamma(:,:)
  character(len=*)                :: stype
  integer                         :: gp, nr_gp

  nr_gp=size(ef,2)

  do gp=1,nr_gp
    call neohooke2gamma(stype,stress(:,gp), ef(:,gp),efgamma(:,gp))
  enddo

  return
end subroutine neohooke1gamma


! Call for one gauss point
subroutine neohooke2(stype,stress,ef)
  implicit none
  double precision                :: stress(:), ef(:)
  character(len=*)                :: stype
  
  double precision                :: E, v, D1
  double precision                :: F(3,3), C(3,3), iC(3,3), id(3,3)
  double precision                :: J, trC, S(3,3)

  D1=2D0/Kmp

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

  S=Kmp/2d0*(J**2d0-1d0)*iC+Gmp*J**(-2d0/3d0)*(Id-trC/3d0*iC)
  if (1.eq.2) then
     S=2d0/D1*(J**2d0-J)*iC+Gmp*J**(-2d0/3d0)*(Id-trC/3d0*iC)
  endif

 

  if (stype.eq.'Cauchy') then
    S=matmul(F,matmul(S,transpose(F)))/J
  elseif (stype.eq.'2ndPiola') then
  else
    stop 'stype not implemented'
  endif
  if (size(ef).eq.4) then
    stress=(/S(1,1),S(2,2),S(3,3),S(1,2)/)
  else
    stress=(/S(1,1),S(2,2),S(3,3),S(1,2),S(1,3),S(2,3)/)
  endif

  return
end subroutine neohooke2


subroutine neohooke2gamma(stype,stress,ef,efgamma)
  implicit none
  double precision                :: stress(:), ef(:), efgamma(:)
  character(len=*)                :: stype
  
  double precision                :: E, v, D1
  double precision                :: F(3,3), C(3,3), iC(3,3), id(3,3)
  double precision                :: J, trC, S(3,3)

  D1=2D0/Kmp

  if (size(ef).eq.4) then
    F(1,:)=(/efgamma(1), efgamma(2), 0d0/)
    F(2,:)=(/efgamma(3), efgamma(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)
  else
    F(1,:)=(/efgamma(1), efgamma(2), efgamma(3)/)
    F(2,:)=(/efgamma(4), efgamma(5), efgamma(6)/)
    F(3,:)=(/efgamma(7), efgamma(8), efgamma(9)/)
  endif
  
  C=matmul(transpose(F),F)
  call inv3s(iC,C)
  J=det3(F)
  Id=0d0
  Id(1,1)=1d0
  Id(2,2)=1d0
  Id(3,3)=1d0
  trC=C(1,1)+C(2,2)+C(3,3)

  S=Kmp/2d0*(J**2d0-1d0)*iC+Gmp*J**(-2d0/3d0)*(Id-trC/3d0*iC)
  if (1.eq.2) then
     S=2d0/D1*(J**2d0-J)*iC+Gmp*J**(-2d0/3d0)*(Id-trC/3d0*iC)
  endif

 

  if (stype.eq.'Cauchy') then
     if (size(ef).eq.4) then
       F(1,:)=(/ef(1), ef(2), 0d0/)
       F(2,:)=(/ef(3), ef(4), 0d0/)
       F(3,:)=(/  0d0,   0d0, 1d0/)
     else
       F(1,:)=(/ef(1), ef(2), ef(3)/)
       F(2,:)=(/ef(4), ef(5), ef(6)/)
       F(3,:)=(/ef(7), ef(8), ef(9)/)
     endif
     S=matmul(F,matmul(S,transpose(F)))/J
  elseif (stype.eq.'2ndPiola') then
  else
    stop 'stype not implemented'
  endif
  if (size(ef).eq.4) then
    stress=(/S(1,1),S(2,2),S(3,3),S(1,2)/)
  else
    stress=(/S(1,1),S(2,2),S(3,3),S(1,2),S(1,3),S(2,3)/)
  endif

  return
end subroutine neohooke2gamma



! Calculate strain energy!
subroutine neohooke1energy(energy,ef,remove)
  implicit none
  double precision                :: energy(:), ef(:,:)
!  character(len=*)                :: stype
  integer                         :: gp, nr_gp
  logical								 :: remove

  nr_gp=size(ef,2)
  remove = .false.
  
  do gp=1,nr_gp
    call neohooke2energy(energy(gp), ef(:,gp), remove)
  enddo

  return
end subroutine neohooke1energy


! Call for one gauss point
subroutine neohooke2energy(energy,ef,remove)
  implicit none
  double precision                ::  ef(:), energy
!  character(len=*)                :: stype
 
  double precision                :: E, v, D1
  double precision                :: F(3,3), C(3,3), iC(3,3), id(3,3)
  double precision                :: J, trC, S(3,3)
  logical								 :: remove
  
  D1=2D0/Kmp

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
!  Id=0d0
!  Id(1,1)=1d0
!  Id(2,2)=1d0
!  Id(3,3)=1d0
  trC=C(1,1)+C(2,2)+C(3,3)

!  S=Kmp/2d0*(J**2d0-1d0)*iC+Gmp*J**(-2d0/3d0)*(Id-trC/3d0*iC)

  if (J.lt.1d-1) then
  	  write(*,*) 'det(F) < 1d-1, J = ', J
  	  remove = .true.
  endif

  energy=Kmp/2d0*(0.5D0*(J**2d0-1D0)-dlog(J))+Gmp/2D0*(J**(-2d0/3d0)*trC-3D0)

  return
end subroutine neohooke2energy




! Routines for a mixed u-p-J formulation

! Call for several gauss points
subroutine neohooke1up(stype,stress,ef,th)
  implicit none
  double precision                :: stress(:,:), ef(:,:), th(:)
  character(len=*)                :: stype
  integer                         :: gp, nr_gp

  nr_gp=size(ef,2)

  do gp=1,nr_gp
    call neohooke2up(stype,stress(:,gp), ef(:,gp),th(gp))
  enddo

  return
end subroutine neohooke1up

! Call for one gauss point
subroutine neohooke2up(stype,stress,ef,th)
  implicit none
  double precision                :: stress(:), ef(:), th
  character(len=*)                :: stype
  
  double precision                :: E, v
  double precision                :: F(3,3), C(3,3), iC(3,3), id(3,3)
  double precision                :: J, trC, S(3,3), b(3,3), trb, p
  
  if (size(ef).eq.4) then
    F(1,:)=(/ef(1), ef(2), 0d0/)
    F(2,:)=(/ef(3), ef(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)
  else
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
  endif
  
  J=det3(F)
  Id=0d0
  Id(1,1)=1d0
  Id(2,2)=1d0
  Id(3,3)=1d0

  p=Kmp/2d0*(th-1d0/th)

  if (stype.eq.'Cauchy') then
    b=J**(-2d0/3d0)*matmul(F,transpose(F))
    trb=b(1,1)+b(2,2)+b(3,3)
    b=b-trb/3d0*Id
    S=Gmp*b/J+p*id
  elseif (stype.eq.'2ndPiola') then
    C=matmul(transpose(F),F)
    call inv3s(iC,C)
    trC=C(1,1)+C(2,2)+C(3,3)
    S=p*J*iC+Gmp*J**(-2d0/3d0)*(Id-trC/3d0*iC)
  else
    stop 'stype not implemented'
  endif
  if (size(ef).eq.4) then
    stress=(/S(1,1),S(2,2),S(3,3),S(1,2)/)
  else
    stress=(/S(1,1),S(2,2),S(3,3),S(1,2),S(1,3),S(2,3)/)
  endif

  return
end subroutine neohooke2up

! Call for constant th
subroutine neohooke3up(stype,stress,ef,th)
  implicit none
  double precision                :: stress(:,:), ef(:,:), th
  character(len=*)                :: stype
  integer                         :: gp, nr_gp

  nr_gp=size(ef,2)

  do gp=1,nr_gp
    call neohooke2up(stype,stress(:,gp), ef(:,gp),th)
  enddo

  return
end subroutine neohooke3up



! Material tangent stiffness for a pure displacement formulation

subroutine dneohooke1(stype,D,ef)
  implicit none
  double precision                :: D(:,:,:), ef(:,:)
  character(len=*)                :: stype
  integer                         :: gp,nr_gp

  nr_gp=size(D,3)

  do gp=1,nr_gp
    call dneohooke2(stype,D(:,:,gp),ef(:,gp))
  enddo

  return
end subroutine dneohooke1

subroutine dneohooke1gamma(stype,D,ef,efgamma)
  implicit none
  double precision                :: D(:,:,:), ef(:,:), efgamma(:,:)
  character(len=*)                :: stype
  integer                         :: gp,nr_gp

  nr_gp=size(D,3)

  do gp=1,nr_gp
    call dneohooke2gamma(stype,D(:,:,gp),ef(:,gp),efgamma(:,gp))
  enddo

  return
end subroutine dneohooke1gamma


subroutine dneohooke2(stype,D,ef)
  implicit none
  double precision                :: D(:,:), ef(:)
  character(len=*)                :: stype
  
  double precision                :: E, v, K, G
  double precision                :: F(3,3), C(3,3), iC(3,3), Id(3,3)
  double precision                :: J, trC, a1, a2, a3, Ds(21), D1
  integer                         :: i, jj, kk, l, in(21,4), el
  
  K=Kmp
  G=Gmp
  D1=2D0/K
 
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

  if (1.eq.2) then
    a1=2D0/D1*(2D0*J**2D0-J)+2d0/9d0*G*J**(-2d0/3d0)*trC
    a3=G/3d0*J**(-2d0/3d0)*trC-2D0/D1*(J**2d0-J)
  endif


  if (size(ef).eq.4) then
  ! if (1.eq.0) then
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


end subroutine dneohooke2



subroutine dneohooke2gamma(stype,D,ef,efgamma)
  implicit none
  double precision                :: D(:,:), ef(:), efgamma(:)
  character(len=*)                :: stype
  
  double precision                :: E, v, K, G
  double precision                :: F(3,3), C(3,3), iC(3,3), Id(3,3)
  double precision                :: J, trC, a1, a2, a3, Ds(21), D1
  integer                         :: i, jj, kk, l, in(21,4), el
  
  K=Kmp
  G=Gmp
  D1=2D0/K
 
  if (size(ef).eq.4) then
    F(1,:)=(/efgamma(1), efgamma(2), 0d0/)
    F(2,:)=(/efgamma(3), efgamma(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)
  else
    F(1,:)=(/efgamma(1), efgamma(2), efgamma(3)/)
    F(2,:)=(/efgamma(4), efgamma(5), efgamma(6)/)
    F(3,:)=(/efgamma(7), efgamma(8), efgamma(9)/)
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

  if (1.eq.2) then
    a1=2D0/D1*(2D0*J**2D0-J)+2d0/9d0*G*J**(-2d0/3d0)*trC
    a3=G/3d0*J**(-2d0/3d0)*trC-2D0/D1*(J**2d0-J)
  endif


  if (size(ef).eq.4) then
  ! if (1.eq.0) then
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


end subroutine dneohooke2gamma



! Material tangent stiffness for a mixed displacement formulation

subroutine dneohooke1up(stype,D,Dth,ef,th)
  implicit none
  double precision                :: D(:,:,:), Dth(:), ef(:,:), th(:)
  character(len=*)                :: stype
  integer                         :: gp,nr_gp

  nr_gp=size(D,3)

  do gp=1,nr_gp
    call dneohooke2up(stype,D(:,:,gp), Dth(gp),ef(:,gp),th(gp))
  enddo

  return
end subroutine dneohooke1up


! Call when constant th in element
subroutine dneohooke3up(stype,D,Dth,ef,th)
  implicit none
  double precision                :: D(:,:,:), Dth(:), ef(:,:), th
  character(len=*)                :: stype
  integer                         :: gp,nr_gp

  nr_gp=size(D,3)

  do gp=1,nr_gp
    call dneohooke2up(stype,D(:,:,gp), Dth(gp),ef(:,gp),th)
  enddo

  return
end subroutine dneohooke3up


subroutine dneohooke2up(stype,D,Dth,ef,th)
  implicit none
  double precision                :: D(:,:), Dth, ef(:), th
  character(len=*)                :: stype
  
  double precision                :: F(3,3), b(3,3), bt(3,3), Id(3,3)
  double precision                :: J, trb, a1, a2, a3, Ds(21)
  double precision                :: tmp, tmp1, tmp2, tmp3, taud(3,3)
  integer                         :: i, jj, kk, l, in(21,4), el
  
  if (size(ef).eq.4) then
    F(1,:)=(/ef(1), ef(2), 0d0/)
    F(2,:)=(/ef(3), ef(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)
  else
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
  endif
  
  b=matmul(F,transpose(F))
  J=det3(F)
  Id=0d0
  Id(1,1)=1d0
  Id(2,2)=1d0
  Id(3,3)=1d0
  trb=b(1,1)+b(2,2)+b(3,3)
  
  taud=Gmp*J**(-2d0/3d0)*(b-1d0/3d0*trb*Id)

if (1.eq.1) then
  D(1,:)=(/2d0*taud(1,1),taud(1,1)+taud(2,2),taud(1,1)+taud(3,3), &
         taud(1,2),taud(1,3),taud(2,3)/)
  D(2,:)=(/taud(2,2)+taud(1,1),2d0*taud(2,2),taud(2,2)+taud(3,3), &
         taud(1,2),taud(1,3),taud(2,3)/)
  D(3,:)=(/taud(3,3)+taud(1,1),taud(3,3)+taud(2,2),2d0*taud(3,3), &
         taud(1,2),taud(1,3),taud(2,3)/)
  D(4,:)=(/taud(1,2),taud(1,2),taud(1,2),0d0,0d0,0d0/)
  D(5,:)=(/taud(1,3),taud(1,3),taud(1,3),0d0,0d0,0d0/)
  D(6,:)=(/taud(2,3),taud(2,3),taud(2,3),0d0,0d0,0d0/)
  D=D*(-1d0/J*2d0/3d0)

  tmp=-Gmp*trb*J**(-2d0/3d0)/J*2d0/3d0
  tmp1=-4d0/6d0*tmp
  tmp2=1d0/3d0*tmp
  tmp3=-0.5d0*tmp
  D(1,:)=D(1,:)+(/tmp1, tmp2, tmp2,  0d0,  0d0,  0d0/)
  D(2,:)=D(2,:)+(/tmp2, tmp1, tmp2,  0d0,  0d0,  0d0/)
  D(3,:)=D(3,:)+(/tmp2, tmp2, tmp1,  0d0,  0d0,  0d0/)
  D(4,:)=D(4,:)+(/ 0d0,  0d0,  0d0, tmp3,  0d0,  0d0/)
  D(5,:)=D(5,:)+(/ 0d0,  0d0,  0d0,  0d0, tmp3,  0d0/)
  D(6,:)=D(6,:)+(/ 0d0,  0d0,  0d0,  0d0,  0d0, tmp3/)
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


    tmp1=-1d0/J*2d0/3d0
    tmp2=tmp1*Gmp*trb*J**(-2d0/3d0)
    do el=1,21
      i=in(el,1)
      jj=in(el,2)
      kk=in(el,3)
      l=in(el,4)

      Ds(el)=tmp1*(taud(i,jj)*Id(kk,l)+Id(i,jj)*taud(kk,l)) &
            +tmp2*(1d0/3d0*Id(i,jj)*Id(kk,l) &
       -1d0/2d0*(Id(i,kk)*Id(jj,l)+Id(i,l)*Id(jj,kk)))
    enddo

    D(1,:)=(/Ds(1),  Ds(2),  Ds(3),  Ds(4),  Ds(5),  Ds(6)/)
    D(2,:)=(/Ds(2),  Ds(7),  Ds(8),  Ds(9), Ds(10), Ds(11)/)
    D(3,:)=(/Ds(3),  Ds(8), Ds(12), Ds(13), Ds(14), Ds(15)/)
    D(4,:)=(/Ds(4),  Ds(9), Ds(13), Ds(16), Ds(17), Ds(18)/)
    D(5,:)=(/Ds(5), Ds(10), Ds(14), Ds(17), Ds(19), Ds(20)/)
    D(6,:)=(/Ds(6), Ds(11), Ds(15), Ds(18), Ds(20), Ds(21)/)

!    write(*,*)'D(1,:) ',D(1,:)
!    write(*,*)'D(2,:) ',D(2,:)
!    write(*,*)'D(3,:) ',D(3,:)
!    write(*,*)'D(4,:) ',D(4,:)
!    write(*,*)'D(5,:) ',D(5,:)
!    write(*,*)'D(6,:) ',D(6,:)
!stop
endif

  if (stype.eq.'ul') then
  elseif (stype.eq.'tl') then
    call pullback(D,D,ef,'j')
  else
    stop 'stype not implemented'
  endif

  Dth=1d0/2d0*Kmp*(1d0+1d0/th/th)

  return
end subroutine dneohooke2up


! The xxxxxsens routines is used for topology optimzation



!subroutine dneohooke1sens(stype,D,ed,Q,ef,Emu,efmutilde,sigtilde,Test,gfun)
subroutine dneohooke1sens(stype,D,Q,ef,Emu,efmutilde,sigtilde,SD)
  implicit none
  double precision                :: D(:,:,:), Q(:,:),ef(:,:) ,Emu(:,:),efmutilde(:,:),sigtilde(:,:), SD(:,:)
!  double precision                :: ed(:), rho, Test(:), gfun
  double precision                ::  rho, gfun
  character(len=*)                :: stype
  integer                         :: gp,nr_gp

  nr_gp=size(D,3)

!  Test=0D0
  do gp=1,nr_gp
!    call dneohooke2sens(stype,D(:,:,gp),ed,Q(:,gp),ef(:,gp),Emu(:,gp),efmutilde(:,gp),sigtilde(:,gp),Test(gp),gfun)
    call dneohooke2sens(stype,D(:,:,gp),Q(:,gp),ef(:,gp),Emu(:,gp),efmutilde(:,gp),sigtilde(:,gp),SD(:,gp))
  enddo

  return
end subroutine dneohooke1sens

!subroutine dneohooke2sens(stype,D,ed,Q,ef,Emu,efmutilde,sigtilde,Test,gfun)
subroutine dneohooke2sens(stype,D,Q,ef,Emu,efmutilde,sigtilde,SD)
  implicit none
  double precision                :: D(:,:), Q(:),ef(:), rho, Emu(:), efmutilde(:), sigtilde(:)
  double precision                :: tmp33(3,3), Q2(6),Test, gfun, SD(6), SS(6)
  character(len=*)                :: stype

  double precision                :: E, v, K, G
  double precision                :: F(3,3), C(3,3), iC(3,3), Id(3,3), Emusq(3,3), Qsq(3,3)
  double precision                :: a1, a2, a3, BigA1, BigA2, BigA3, Emutr
 
  double precision                :: da1(3,3), da2(3,3), da3(3,3), dBigA1(3,3), dBigA2(3,3), dBigA3(3,3)
  double precision                :: iC_Emu(3,3), iC_Emutr

  double precision                :: J, trC, D1
  integer                         :: i, jj, kk, l,  el

  Id     =0d0
  Id(1,1)=1d0
  Id(2,2)=1d0
  Id(3,3)=1d0


! Scaling is performed outside material routine
  K=Kmp
  G=Gmp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 
  D1=2D0/K
   
  Emusq(1,:)=(/Emu(1),     Emu(4)/2D0, Emu(5)/2D0 /)
  Emusq(2,:)=(/Emu(4)/2D0, Emu(2),     Emu(6)/2D0 /)
  Emusq(3,:)=(/Emu(5)/2D0, Emu(6)/2D0, Emu(3) /)
  
  if (size(ef).eq.4) then
    F(1,:)=(/ef(1), ef(2), 0d0/)
    F(2,:)=(/ef(3), ef(4), 0d0/)
    F(3,:)=(/ 0d0,   0d0,  1d0/)
  else
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
  endif
  
  C=matmul(transpose(F),F)
  call inv3s(iC,C)
  J        = det3(F)
  trC      = C(1,1)+C(2,2)+C(3,3)

  iC_Emu   = matmul(iC,Emusq)
  iC_Emutr = iC_Emu(1,1)+iC_Emu(2,2)+iC_Emu(3,3)
  Emutr    = Emu(1)+Emu(2)+Emu(3)

  BigA1    = iC_Emutr**2d0
  BigA2    = 2D0*Emutr*iC_Emutr
  tmp33    = 2D0*matmul(iC_Emu,iC_Emu) 
  BigA3    = tmp33(1,1)+tmp33(2,2)+tmp33(3,3)

  a1 = K*J**2d0+2d0/9d0*G*J**(-2d0/3d0)*trC
  a2 = -2d0*G/3d0*J**(-2d0/3d0)
  a3 = -K/2d0*(J**2d0-1d0)+G/3d0*J**(-2d0/3d0)*trC

  tmp33  = matmul(iC_Emu,iC)
  dBigA1 = -4D0*iC_Emutr*tmp33
  dBigA2 = -4D0*Emutr*tmp33
  dBigA3 = -8D0*matmul(tmp33,transpose(iC_Emu))

  da1 = (2D0*K*J**2D0-4D0*G*J**(-2D0/3D0)/27D0*trC)*iC+4D0*G*J**(-2d0/3D0)/9D0*Id
  da2 = 4D0*G/9D0*J**(-2D0/3D0)*iC
  da3 =-(K*J**2D0+2D0*J**(-2D0/3D0)*G/9D0*trC)*iC+2D0*J**(-2D0/3D0)*G/3D0*Id

  if (1.eq.1) then
    a1=2D0/D1*(2D0*J**2D0-J)+2d0/9d0*G*J**(-2d0/3d0)*trC
    a3=G/3d0*J**(-2d0/3d0)*trC-2D0/D1*(J**2d0-J)

    da1 = (2D0/D1*(4D0*J**2D0-J)-4D0*G*J**(-2D0/3D0)/27D0*trC)*iC+4D0*G*J**(-2d0/3D0)/9D0*Id
    da3 =-(2D0/D1*(2D0*J**2D0-J)+2D0*J**(-2D0/3D0)*G/9D0*trC)*iC+2D0*J**(-2D0/3D0)*G/3D0*Id
  endif



!  Test=a1*BigA1+a2*BigA2+a3*BigA3

  Qsq = BigA1*da1+BigA2*da2+BigA3*da3+dBigA1*a1+dBigA2*a2+dBigA3*a3

  SD=(/Qsq(1,1), Qsq(2,2), Qsq(3,3), Qsq(1,2), Qsq(1,3), Qsq(2,3)/)

  Q2(1)=efmutilde(1)**2D0+efmutilde(4)**2D0+efmutilde(7)**2D0
  Q2(2)=efmutilde(2)**2D0+efmutilde(5)**2D0+efmutilde(8)**2D0
  Q2(3)=efmutilde(3)**2D0+efmutilde(6)**2D0+efmutilde(9)**2D0

  Q2(4)=efmutilde(1)*efmutilde(2)+efmutilde(4)*efmutilde(5)+efmutilde(7)*efmutilde(8)
  Q2(5)=efmutilde(1)*efmutilde(3)+efmutilde(4)*efmutilde(6)+efmutilde(7)*efmutilde(9)
  Q2(6)=efmutilde(2)*efmutilde(3)+efmutilde(5)*efmutilde(6)+efmutilde(8)*efmutilde(9)

  Q2(4:6)=2D0*Q2(4:6)
!  Q2=matmul(D,Q2)
  SS=matmul(D,Q2)

  Q=SD+SS

  sigtilde=matmul(D,Emu)

end subroutine dneohooke2sens


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dfat1(GradFAT,ef,Targ,Aarg)
  implicit none
  double precision                :: GradFAT(:,:),ef(:,:) ,Targ(:,:), Aarg(:,:)
  double precision                :: rho, gfun
!  character(len=*)               :: stype
  integer                         :: gp,nr_gp

  nr_gp=size(ef,2)

!  Test=0D0
  do gp=1,nr_gp
!    call dneohooke2sens(stype,D(:,:,gp),ed,Q(:,gp),ef(:,gp),Targ(:,gp),efmutilde(:,gp),sigtilde(:,gp),Test(gp),gfun)
!    call dneohooke2sens(stype,D(:,:,gp),Q(:,gp),ef(:,gp),Targ(:,gp),efmutilde(:,gp),Aarg(:,gp),efnu(:,gp),sigtilde(:,gp),GradFAT(:,gp))
!    call dneohooke2sens(stype,D(:,:,gp),Q(:,gp),ef(:,gp),Targ(:,gp),efmutilde(:,gp),Aarg(:,gp),efnu(:,gp),GradFAT(:,gp))
!    call dfat2(D(:,:,gp),Q(:,gp),ef(:,gp),Targ(:,gp),efmutilde(:,gp),Aarg(:,gp),efnu(:,gp))
     call dfat2(GradFAT(:,gp),ef(:,gp),Targ(:,gp),Aarg(:,gp))
  enddo

  return
end subroutine dfat1

!subroutine dneohooke2sens(stype,D,ed,Q,ef,Targtilde,efmutilde,sigtilde,Test,gfun)
!subroutine dneohooke2sens(stype,D,Q,ef,Targtilde,efmutilde, Aarg,efnu,GradFAT)
subroutine dfat2(GradFAT, ef,Targ, Aarg)

  implicit none
  double precision                :: ef(:), rho, Targ(:), Aarg(:)
  double precision                :: tmp33(3,3),Test, gfun, GradFAT(:)
!  character(len=*)               :: stype

  double precision                :: E, v, K, G
  double precision                :: F(3,3), C(3,3), iC(3,3), Id(3,3), Targsq(3,3), GradFAT_sq(3,3)
  double precision                :: a1, a2, a3, BigA1, BigA2, BigA3
 
  double precision                :: da1(3,3), da2(3,3), da3(3,3), dBigA1(3,3), dBigA2(3,3), dBigA3(3,3)
  double precision                :: iC_Targ(3,3), iC_Targtr, iC_Aarg(3,3), iC_Aargtr

  double precision                :: J, trC, D1
  integer                         :: i, jj, kk, l,  el

  double precision                :: Aargsq(3,3), tmp33_A(3,3), tmp33_T(3,3), Aargtr, Targtr
 
  Id     =0d0
  Id(1,1)=1d0
  Id(2,2)=1d0
  Id(3,3)=1d0

! Scaling is performed outside material routine
  K=Kmp
  G=Gmp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  D1=2D0/K
   
  if (size(Targ).eq.3) then
     Targsq(1,:)=(/Targ(1),     Targ(3)/2D0, 0D0 /)
     Targsq(2,:)=(/Targ(3)/2D0, Targ(2),     0D0 /)
     Targsq(3,:)=(/0D0,         0D0,         0d0 /) 
     
     Aargsq(1,:)=(/Aarg(1),     Aarg(3)/2D0, 0D0 /)
     Aargsq(2,:)=(/Aarg(3)/2D0, Aarg(2),     0D0 /)
     Aargsq(3,:)=(/0D0,         0D0,         0d0 /)
  else
     Targsq(1,:)=(/Targ(1),     Targ(4)/2D0, Targ(5)/2D0 /)
     Targsq(2,:)=(/Targ(4)/2D0, Targ(2),     Targ(6)/2D0 /)
     Targsq(3,:)=(/Targ(5)/2D0, Targ(6)/2D0, Targ(3) /)
     
     Aargsq(1,:)=(/Aarg(1),     Aarg(4)/2D0, Aarg(5)/2D0 /)
     Aargsq(2,:)=(/Aarg(4)/2D0, Aarg(2),     Aarg(6)/2D0 /)
     Aargsq(3,:)=(/Aarg(5)/2D0, Aarg(6)/2D0, Aarg(3) /)
  endif


  if (size(ef).eq.4) then
    F(1,:)=(/ef(1), ef(2), 0d0/)
    F(2,:)=(/ef(3), ef(4), 0d0/)
    F(3,:)=(/ 0d0,   0d0,  1d0/)
  else
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
  endif
  
  C=matmul(transpose(F),F)
  call inv3s(iC,C)
!  J        = det3(F)
  J        = dabs(det3(F))
  trC      = C(1,1)+C(2,2)+C(3,3)

  iC_Targ= matmul(iC,Targsq)
  iC_Aarg= matmul(iC,Aargsq)

  iC_Targtr = iC_Targ(1,1)+iC_Targ(2,2)+iC_Targ(3,3)
  iC_Aargtr = iC_Aarg(1,1)+iC_Aarg(2,2)+iC_Aarg(3,3)

  if (size(Targ).eq.3) then
      Targtr = Targ(1)+Targ(2)
      Aargtr = Aarg(1)+Aarg(2)
  else
      Targtr = Targ(1)+Targ(2)+Targ(3)
      Aargtr = Aarg(1)+Aarg(2)+Aarg(3)
  endif   


!  BigA1 = iC_Targtr**2d0
  BigA1  = iC_Targtr*iC_Aargtr

!  BigA2    = 2D0*Targtr*iC_Targtr
  BigA2    = Aargtr*iC_Targtr+Targtr*iC_Aargtr

!  tmp33    = 2D0*matmul(iC_Targ,iC_Targ)
  tmp33    = 2D0*matmul(iC_Aarg,iC_Targ)  
  BigA3    = tmp33(1,1)+tmp33(2,2)+tmp33(3,3)

  a1 = K*J**2d0+2d0/9d0*G*J**(-2d0/3d0)*trC
  a2 = -2d0*G/3d0*J**(-2d0/3d0)
  a3 = -K/2d0*(J**2d0-1d0)+G/3d0*J**(-2d0/3d0)*trC

  tmp33_T = matmul(iC_Targ,iC)  ! =iC* Targ*iC
  tmp33_A  = matmul(iC_Aarg,iC)  ! =iC* Aarg*iC

!  dBigA1 = -4D0*iC_Targtr*tmp33
!  dBigA2 = -4D0*Targtr*tmp33
!  dBigA3 = -8D0*matmul(tmp33,transpose(iC_Targ))

  dBigA1 = -2D0*(iC_Targtr*tmp33_A+iC_Aargtr*tmp33_T)
  dBigA2 = -2D0*(Targtr*tmp33_A+Aargtr*tmp33_T)
  dBigA3 = -4D0*(matmul(tmp33_T,transpose(iC_Aarg))+matmul(tmp33_A,transpose(iC_Targ)))

  da1 = (2D0*K*J**2D0-4D0*G*J**(-2D0/3D0)/27D0*trC)*iC+4D0*G*J**(-2d0/3D0)/9D0*Id
  da2 = 4D0*G/9D0*J**(-2D0/3D0)*iC
  da3 =-(K*J**2D0+2D0*J**(-2D0/3D0)*G/9D0*trC)*iC+2D0*J**(-2D0/3D0)*G/3D0*Id

  if (1.eq.2) then
    a1=2D0/D1*(2D0*J**2D0-J)+2d0/9d0*G*J**(-2d0/3d0)*trC
    a3=G/3d0*J**(-2d0/3d0)*trC-2D0/D1*(J**2d0-J)

    da1 = (2D0/D1*(4D0*J**2D0-J)-4D0*G*J**(-2D0/3D0)/27D0*trC)*iC+4D0*G*J**(-2d0/3D0)/9D0*Id
    da3 =-(2D0/D1*(2D0*J**2D0-J)+2D0*J**(-2D0/3D0)*G/9D0*trC)*iC+2D0*J**(-2D0/3D0)*G/3D0*Id
  endif

  GradFAT_sq = BigA1*da1+BigA2*da2+BigA3*da3+dBigA1*a1+dBigA2*a2+dBigA3*a3
  
  if (size(gradFAT,1).eq.6) then
     GradFAT=(/GradFAT_sq(1,1), GradFAT_sq(2,2), GradFAT_sq(3,3), GradFAT_sq(1,2), GradFAT_sq(1,3), GradFAT_sq(2,3)/)
  else
     GradFAT=(/GradFAT_sq(1,1), GradFAT_sq(2,2), GradFAT_sq(1,2)/)
  endif
end subroutine dfat2




end module mater_hyperel
