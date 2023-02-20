module mater_hyperel_elect

! last modified
! M. Ristinmaa 2011-09-20
!  - initial code hyper elector-elastic material
! M. Ristinmaa 2011-09-22
!  - Elastic part debugged, Cm part debugged, good convergence
! M. Ristinmaa 2011-10-07
!  - Tangent stiffness for inverse motion problem stype='inv' implemented
!  - Stress calculation using inverse deformation gradient stype='Cauchy-inv'
!    implemented for use in inverse motion problem
! ------------------------------------------------------------------------------

! Eo part is not checked. DEGUGGING required

use mater_large
use some_constants
use matrix_util

implicit none

double precision                  :: mpK, mpMy, mpCm, mpCe, mpEo
private mpK, mpMy, mpCm, mpCe, mpEo

interface dneohooke_elect
  module procedure dneohooke1
  module procedure dneohooke2a
  module procedure dneohooke1up
  module procedure dneohooke2up
end interface

private dneohooke1, dneohooke2a

interface dineohooke_elect
  module procedure dineohooke1
  module procedure dineohooke2
  module procedure dineohooke1up
  module procedure dineohooke2up
end interface

private dineohooke1, dineohooke2, dineohooke2b

interface neohooke_elect
  module procedure neohooke1
  module procedure neohooke2
  module procedure neohooke_up1
  module procedure neohooke_up2
end interface

private neohooke1, neohooke2

!-------------------------------------------------------------------------------

contains

subroutine neohooke_elect_init(mp)
  implicit none
  double precision                :: mp(:)
  mpK =mp(1)
  mpMy=mp(2)
  mpCm=mp(3)
  mpCe=mp(4)

! This part needs a debugging
  mpEo=zero
!  mpEo=mp(5)

  return
end subroutine neohooke_elect_init

! Call for several gauss points
subroutine neohooke1(stype,stress,edisp,ef,de)
  implicit none
  double precision                :: stress(:,:), edisp(:,:), ef(:,:), de(:,:)
  character(len=*)                :: stype
  integer                         :: gp, nr_gp

  nr_gp=size(ef,2)

  do gp=1,nr_gp
    call neohooke2(stype,stress(:,gp), edisp(:,gp),ef(:,gp), de(:,gp))
  enddo

  return
end subroutine neohooke1

! Call for one gauss point
subroutine neohooke2(stype,stress,edisp,ef,e)
  implicit none
  double precision                :: stress(:), edisp(:), ef(:), e(:)
  character(len=*)                :: stype
  
  double precision                :: F(3,3), b(3,3), id(3,3), invF(3,3)
  double precision                :: J, Jm53, tr_b, S(3,3)
  double precision                :: hat_e(3), tilde_e(3), tr_e, tr_hat_e
  double precision                :: dy_e(3,3), dy_hat_e(3,3)
  
  if (size(ef).eq.4) then
    F(1,:)=(/ef(1), ef(2), zero/)
    F(2,:)=(/ef(3), ef(4), zero/)
    F(3,:)=(/ zero,  zero,  one/)
  else
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
  endif

! Check if a inverse deformation gradient is provided as input
  if (stype.eq.'Cauchy-inv') then
    call inv3(invF,F)
    F=invF
  end if

  b=matmul(F,transpose(F))
  J=det3(F)

  Jm53=J**(-five*onethird)
  Id=zero
  Id(1,1)=one
  Id(2,2)=one
  Id(3,3)=one
  tr_b=b(1,1)+b(2,2)+b(3,3)

! Assumes 3 components

  hat_e=matmul(b,e)
  tilde_e=matmul(b,hat_e)
  tr_e=dot_product(e,e)
  tr_hat_e=dot_product(hat_e,hat_e)

! Dyadic products
  dy_hat_e=spread(hat_e,2,3)*spread(hat_e,1,3)
  dy_e=spread(e,2,3)*spread(e,1,3)

! Cauchy stress tensor
  S=mpK*(J-one)*Id+mpMy*Jm53*(b-onethird*tr_b*Id) &
    +two*mpCm*Jm53*(dy_hat_e-onethird*tr_hat_e*Id) &
    +mpEo*(dy_e-onehalf*tr_e*Id)

! Electric displacement (spatial)
  edisp=-two/J*mpCe*hat_e-two*mpCm*Jm53*tilde_e+mpEo*e

  if (stype.eq.'Cauchy') then
  elseif (stype.eq.'Cauchy-inv') then
  elseif (stype.eq.'2ndPiola') then
    stop 'stype not implemented'
  else
    stop 'stype not implemented'
  endif
  if (size(f).eq.4) then
    stress=(/S(1,1),S(2,2),S(3,3),S(1,2)/)
  else
    stress=(/S(1,1),S(2,2),S(3,3),S(1,2),S(1,3),S(2,3)/)
  endif

  return
end subroutine neohooke2


subroutine dneohooke1(stype,Duu,Dup,Dpu,Dpp,ef,e)
  implicit none
  double precision                :: Duu(:,:,:), Dup(:,:,:)
  double precision                :: Dpu(:,:,:), Dpp(:,:,:)
  double precision                :: ef(:,:), e(:,:)
  character(len=*)                :: stype
  integer                         :: gp,nr_gp

  nr_gp=size(Duu,3)

  do gp=1,nr_gp
    call dneohooke2a(stype,Duu(:,:,gp),Dup(:,:,gp),Dpu(:,:,gp),Dpp(:,:,gp), &
                    ef(:,gp),e(:,gp))
  enddo

  return
end subroutine dneohooke1


!subroutine dneohooke2(stype,Duu,Dup,Dpu,Dpp,ef,e)
!  implicit none
!  double precision                :: Duu(:,:), Dup(:,:), Dpu(:,:), Dpp(:,:)
!  double precision                :: ef(:), e(:)
!  character(len=*)                :: stype
!  
! Check if inverse motion problem
!  if (stype.eq.'inv') then
!    call dneohooke2b(stype,Duu,Dup,Dpu,Dpp,ef,e)
!    call dneohooke2c(stype,Duu,Dup,Dpu,Dpp,ef,e)
!  else
!    call dneohooke2a(stype,Duu,Dup,Dpu,Dpp,ef,e)
!  end if
!
!  return
!end subroutine dneohooke2


subroutine dneohooke2a(stype,Duu,Dup,Dpu,Dpp,ef,e)
  implicit none
  double precision                :: Duu(:,:), Dup(:,:), Dpu(:,:), Dpp(:,:)
  double precision                :: ef(:), e(:)
  character(len=*)                :: stype
  
  double precision                :: F(3,3), b(3,3), Id(3,3)
  double precision                :: Ja, Jm23, tr_b
  integer                         :: i, j, k, l, ii, jj
  
  double precision                :: hat_e(3), tilde_e(3), tr_e, tr_hat_e
  double precision                :: c1, c2, c3, c4, c5 ,c6, c7

  integer, parameter              :: d_list(2,6)=[(/1,1/),(/2,2/),(/3,3/),&
                                                  (/1,2/),(/1,3/),(/2,3/)]

 
  if (size(ef).eq.4) then
    F(1,:)=(/ef(1), ef(2), zero/)
    F(2,:)=(/ef(3), ef(4), zero/)
    F(3,:)=(/ zero,  zero,  one/)
  else
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
  endif
  
  b=matmul(F,transpose(F))
  tr_b=b(1,1)+b(2,2)+b(3,3)
  Ja=det3(F)
  Jm23=Ja**(-twothird)

  Id=zero
  Id(1,1)=one
  Id(2,2)=one
  Id(3,3)=one

  hat_e=matmul(b,e)
  tilde_e=matmul(b,hat_e)
  tr_e=dot_product(e,e)
  tr_hat_e=dot_product(hat_e,hat_e)

  c1=two*mpK*Ja*(Ja-onehalf)+two*oneninth*mpMy*Jm23*tr_b &
    +four*oneninth*mpCm*Jm23*tr_hat_e-Ja*onehalf*mpEo*tr_e
  c2=-twothird*mpMy*Jm23
  c3=onethird*mpMy*Jm23*tr_b+twothird*mpCm*Jm23*tr_hat_e &
    +onehalf*mpEo*Ja*tr_e-mpK*Ja*(Ja-one)
  c4=-four*onethird*mpCm*Jm23
  c5=mpEo*Ja
  c6=-c5
  c7=two*mpCm*Jm23

  if (size(ef).eq.9) then
  do ii=1,6
    i=d_list(1,ii)
    j=d_list(2,ii)

    do jj=1,6
      k=d_list(1,jj)
      l=d_list(2,jj)

      Duu(ii,jj)=c1*Id(i,j)*Id(k,l)+c2*(b(i,j)*Id(k,l)+Id(i,j)*b(k,l)) &
                +c3*(Id(i,k)*Id(j,l)+Id(i,l)*Id(j,k)) &
                +c4*(hat_e(i)*hat_e(j)*Id(k,l)+Id(i,j)*hat_e(k)*hat_e(l)) &
                +c5*(e(i)*e(j)*Id(k,l)+Id(i,j)*e(k)*e(l)) &
                +c6*(Id(i,k)*e(j)*e(l)+Id(i,l)*e(j)*e(k) &
                    +Id(j,k)*e(i)*e(l)+Id(j,l)*e(i)*e(k))
    end do

    do k=1,3
      Dup(ii,k)=c7*(b(i,k)*hat_e(j)+hat_e(i)*b(j,k)  &
                    -twothird*Id(i,j)*tilde_e(k))  &
               +c5*(Id(i,k)*e(j)+Id(j,k)*e(i)-Id(i,j)*e(k))
    end do
  end do
  else
    stop 'stype not implemented'
  end if

  Duu=Duu/Ja
  Dup=Dup/Ja
  Dpu=-transpose(Dup)
  Dpp=(-two*mpCe*b-2*mpCm*Jm23*matmul(b,b)+c5*Id)/Ja


  if (stype.eq.'ul') then
!    call pushforward(D,D,ef,'j')
  elseif (stype.eq.'tl') then
    stop 'stype not implemented'
  else
    stop 'stype not implemented'
  endif

  return
end subroutine dneohooke2a


! Material routines for use in u/p formulation


! Call for several gauss points
subroutine neohooke_up1(stype,stress,edisp,ef,de,th)
  implicit none
  double precision                :: stress(:,:), edisp(:,:), ef(:,:), de(:,:), th
  character(len=*)                :: stype
  integer                         :: gp, nr_gp

  double precision                :: p

  nr_gp=size(ef,2)

! This only works for Cauchy and Cauchy-inv
  if ((stype.ne.'Cauchy').and.(stype.ne.'Cauchy-inv')) then
     stop "Only Cauchy and Cauchy-inv implemented"
  endif

  do gp=1,nr_gp
     call neohooke_up2(stype,stress(:,gp),edisp(:,gp),ef(:,gp),de(:,gp),th)
  enddo

  return
end subroutine neohooke_up1

subroutine neohooke_up2(stype,stress,edisp,ef,de,th)
  implicit none
  double precision                :: stress(:), edisp(:), ef(:), de(:), th
  character(len=*)                :: stype
  integer                         :: gp, nr_gp

  double precision                :: p

! This only works for Cauchy and Cauchy-inv
  if ((stype.ne.'Cauchy').and.(stype.ne.'Cauchy-inv')) then
     stop "Only Cauchy and Cauchy-inv implemented"
  endif

    call neohooke2(stype,stress(:),edisp(:),ef(:),de(:))
                 
    p=(stress(1)+stress(2)+stress(3))/3d0
    stress(:)=stress(:)  &
               +(mpK*(th-1d0)-p)*(/1d0,1d0,1d0,0d0,0d0,0d0/)

  return
end subroutine neohooke_up2


subroutine dneohooke1up(stype,Duu,Dup,Dtt,Dpu,Dpp,ef,e,th)
  implicit none
  double precision                :: Duu(:,:,:), Dup(:,:,:)
  double precision                :: Dpu(:,:,:), Dpp(:,:,:)
  double precision                :: Dtt(:)
  double precision                :: ef(:,:), e(:,:), th
  character(len=*)                :: stype
  double precision                :: F_m(3,3), Ja, Id(3,3)
  double precision                :: ddudjj, dudj, tmp1, tmp2
  integer                         :: gp,nr_gp, ie, ii, jj, i, j, k, l
  integer, parameter              :: d_list(2,6)=[(/1,1/),(/2,2/),(/3,3/),&
                                                  (/1,2/),(/1,3/),(/2,3/)]

  double precision                :: p

  nr_gp=size(Duu,3)
  Id   =getI()

  do gp=1,nr_gp
 
    call dneohooke2up('ul',Duu(:,:,gp),Dup(:,:,gp),Dtt(gp),Dpu(:,:,gp), &
                           Dpp(:,:,gp),ef(:,gp),e(:,gp),th)
  enddo

  return
end subroutine dneohooke1up


subroutine dneohooke2up(stype,Duu,Dup,Dtt,Dpu,Dpp,ef,e,th)
  implicit none
  double precision                :: Duu(:,:), Dup(:,:)
  double precision                :: Dpu(:,:), Dpp(:,:)
  double precision                :: Dtt
  double precision                :: ef(:), e(:), th
  character(len=*)                :: stype
  double precision                :: F_m(3,3), Ja, Id(3,3)
  double precision                :: ddudjj, dudj, tmp1, tmp2
  integer                         :: gp,nr_gp, ie, ii, jj, i, j, k, l
  integer, parameter              :: d_list(2,6)=[(/1,1/),(/2,2/),(/3,3/),&
                                                  (/1,2/),(/1,3/),(/2,3/)]

  double precision                :: p

  Id   =getI()

  F_m=getF(ef(:))
  Ja=det3(F_m)
  call dneohooke2a('ul',Duu(:,:),Dup(:,:),Dpu(:,:), &
                          Dpp(:,:),ef(:),e(:))

  ddudjj=mpK
  dudj=mpK*(Ja-1D0)
  tmp1=dudj
  tmp2=ddudjj*Ja+dudj
  do ii=1,6
    i=d_list(1,ii)
    j=d_list(2,ii)
    do jj=1,6
      k=d_list(1,jj)
      l=d_list(2,jj)
      Duu(ii,jj)=Duu(ii,jj)-tmp2*Id(i,j)*Id(k,l) &
                   +tmp1*(Id(i,k)*Id(j,l)+Id(i,l)*Id(j,k))
    end do
  end do  
 
  dtt=mpK

  return
end subroutine dneohooke2up


! Material tangent stiffnesses for the inverse motion problem
subroutine dineohooke2b(stype,Duu,Dup,Dpu,Dpp,ef,e)
  implicit none
  double precision                :: Duu(:,:), Dup(:,:), Dpu(:,:), Dpp(:,:)
  double precision                :: ef(:), e(:)
  character(len=*)                :: stype
  
  double precision                :: f(3,3), c(3,3), ic(3,3), Id(3,3), icc(3,3)
  double precision                :: Ja, tr_ic
  integer                         :: i, j, k, l, ii, jj
  
  double precision                :: hat_e(3), tilde_e(3), tr_e, tr_hat_e
  double precision                :: q1, q2, q3, q4
  double precision                :: A(6,9), Duux(6,6), Dpux(3,6)

  integer, parameter              :: d_list(2,6)=[(/1,1/),(/2,2/),(/3,3/),&
                                                  (/1,2/),(/1,3/),(/2,3/)]
  integer, parameter              :: f_list(2,9)=[(/1,1/),(/1,2/),(/1,3/),&
                                                  (/2,1/),(/2,2/),(/2,3/),&
                                                  (/3,1/),(/3,2/),(/3,3/)]

  if (size(ef).eq.4) then
    f(1,:)=(/ef(1), ef(2), zero/)
    f(2,:)=(/ef(3), ef(4), zero/)
    f(3,:)=(/ zero,  zero,  one/)
  else
    f(1,:)=(/ef(1), ef(2), ef(3)/)
    f(2,:)=(/ef(4), ef(5), ef(6)/)
    f(3,:)=(/ef(7), ef(8), ef(9)/)
  endif

  c=matmul(transpose(f),f)
  call inv3s(ic,c)
  icc=matmul(ic,ic)

  tr_ic=ic(1,1)+ic(2,2)+ic(3,3)
  ja=det3(f)

  Id=zero
  Id(1,1)=one
  Id(2,2)=one
  Id(3,3)=one

  hat_e=matmul(ic,e)
  tilde_e=matmul(ic,hat_e)
  tr_e=dot_product(e,e)
  tr_hat_e=dot_product(hat_e,hat_e)

  q1=ja*mpCe
  q2=ja**(five*onethird)*mpCm
  q3=ja**(five*onethird)*mpMy
  q4=-onehalf*mpK/ja-5d0/18d0*q3*tr_ic-5d0/9d0*q2*tr_hat_e

  if (size(ef).eq.9) then
  do ii=1,6
    i=d_list(1,ii)
    j=d_list(2,ii)

    do jj=1,6
      k=d_list(1,jj)
      l=d_list(2,jj)

      Duux(ii,jj)=q4*Id(i,j)*ic(k,l)+5d0/6d0*q3*ic(i,j)*ic(k,l) &
                -onehalf*q3*(ic(i,k)*ic(l,j)+ic(i,l)*ic(k,j)) &
                +1d0/3d0*q3*Id(i,j)*icc(k,l) & 
                +5d0/3d0*q2*hat_e(i)*hat_e(j)*ic(k,l) &
                +2d0/3d0*q2*(tilde_e(k)*hat_e(l)+hat_e(k)*tilde_e(l))*Id(i,j) &
                -q2*(ic(i,k)*hat_e(l)*hat_e(j)+ic(i,l)*hat_e(k)*hat_e(j) &
                    +ic(j,k)*hat_e(l)*hat_e(i)+ic(j,l)*hat_e(k)*hat_e(i))

    end do

    do k=1,3

      Dpux(k,ii)=-q1*hat_e(k)*ic(i,j)-5d0/3d0*q2*tilde_e(k)*ic(i,j) &
               +q1*(ic(k,i)*hat_e(j)+ic(k,j)*hat_e(i)) &
               +q2*(ic(k,i)*tilde_e(j)+ic(k,j)*tilde_e(i) &
                   +icc(k,i)*hat_e(j)+icc(k,j)*hat_e(i))
      Dup(ii,k)=two*q2*(ic(i,k)*hat_e(j)+hat_e(i)*ic(j,k) &
                   -twothird*tilde_e(k)*Id(i,j)) &
                +mpEo*(Id(i,k)*e(j)+e(i)*Id(j,k)-e(k)*Id(i,j))

    end do
  end do
  else
    stop 'stype not implemented'
  end if

  Dpp=-two*q1*ic-two*q2*icc+mpEo*Id

  A=0d0
  A([1,4,5],1)=ef(1:3)
  A([2,4,6],2)=(/ef(2),ef(1),ef(3)/)
  A([3,5,6],3)=(/ef(3),ef(1),ef(2)/)
  A([1,4,5],4)=ef(4:6)
  A([2,4,6],5)=(/ef(5),ef(4),ef(6)/)
  A([3,5,6],6)=(/ef(6),ef(4),ef(5)/)
  A([1,4,5],7)=ef(7:9)
  A([2,4,6],8)=(/ef(8),ef(7),ef(9)/)
  A([3,5,6],9)=(/ef(9),ef(7),ef(8)/)
  A=2d0*A
   
  Duu=matmul(Duux,A)
  Dpu=matmul(Dpux,A)

  return
end subroutine dineohooke2b


! Material tangent stiffnesses for the inverse motion problem
! based on the usual material tangent matrices in an updated
! Lagrange scheme

subroutine dineohooke1(stype,Duu,Dup,Dpu,Dpp,es,d,ef,e)
  implicit none
  double precision                :: Duu(:,:,:), Dup(:,:,:)
  double precision                :: Dpu(:,:,:), Dpp(:,:,:)
  double precision                :: es(:,:), d(:,:)
  double precision                :: ef(:,:), e(:,:)
  character(len=*)                :: stype
  integer                         :: gp,nr_gp

  nr_gp=size(Duu,3)

  do gp=1,nr_gp
    call dineohooke2(stype,Duu(:,:,gp),Dup(:,:,gp),Dpu(:,:,gp),Dpp(:,:,gp), &
                    es(:,gp),d(:,gp),ef(:,gp),e(:,gp))
  enddo

  return
end subroutine dineohooke1


subroutine dineohooke2(stype,Duu,Dup,Dpu,Dpp,es,dd,ef,e)
  implicit none
  double precision                :: Duu(:,:), Dup(:,:), Dpu(:,:), Dpp(:,:)
  double precision                :: es(:), dd(:), ef(:), e(:)
  character(len=*)                :: stype
  
  double precision                :: f(3,3), fi(3,3)
  integer                         :: i, j, k, l, ii, jj
  
  double precision                :: Duux(6,6), Dpux(3,6), duun(6,9), dpun(3,9)
  double precision                :: efn(3), efp(9), s(3,3)

  integer, parameter              :: d_list(2,6)=[(/1,1/),(/2,2/),(/3,3/),&
                                                  (/1,2/),(/1,3/),(/2,3/)]
  integer, parameter              :: f_list(2,9)=[(/1,1/),(/1,2/),(/1,3/),&
                                                  (/2,1/),(/2,2/),(/2,3/),&
                                                  (/3,1/),(/3,2/),(/3,3/)]

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

  efn(1)=e(1)*fi(1,1)+e(2)*fi(2,1)+e(3)*fi(3,1)
  efn(2)=e(1)*fi(1,2)+e(2)*fi(2,2)+e(3)*fi(3,2)
  efn(3)=e(1)*fi(1,3)+e(2)*fi(2,3)+e(3)*fi(3,3)

  efp=(/fi(1,1),fi(1,2),fi(1,3),fi(2,1),fi(2,2),fi(2,3),fi(3,1),fi(3,2),fi(3,3)/)
  call dneohooke_elect('ul',duux,dup,dpux,dpp,efp,e)
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


! due to the coupling part
  duun(:,1)=duun(:,1)-dup(:,1)*efn(1)
  duun(:,2)=duun(:,2)-dup(:,2)*efn(1)
  duun(:,3)=duun(:,3)-dup(:,3)*efn(1)

  duun(:,4)=duun(:,4)-dup(:,1)*efn(2)
  duun(:,5)=duun(:,5)-dup(:,2)*efn(2)
  duun(:,6)=duun(:,6)-dup(:,3)*efn(2)

  duun(:,7)=duun(:,7)-dup(:,1)*efn(3)
  duun(:,8)=duun(:,8)-dup(:,2)*efn(3)
  duun(:,9)=duun(:,9)-dup(:,3)*efn(3)

! Dpu matrix for inverse motion
  do i=1,3
    do jj=1,9
      k=f_list(1,jj)
      l=f_list(2,jj)
      dpun(i,jj)=dd(i)*fi(l,k)-fi(i,k)*dd(l)-dpp(i,l)*efn(k)
    end do
  end do

! order 11 12 13 21 22 23 31 32 33
  dpun(:,1)=dpun(:,1)-dpux(:,1)*fi(1,1)-dpux(:,4)*fi(2,1)-dpux(:,5)*fi(3,1)
  dpun(:,2)=dpun(:,2)-dpux(:,4)*fi(1,1)-dpux(:,2)*fi(2,1)-dpux(:,6)*fi(3,1)
  dpun(:,3)=dpun(:,3)-dpux(:,5)*fi(1,1)-dpux(:,6)*fi(2,1)-dpux(:,3)*fi(3,1)

  dpun(:,4)=dpun(:,4)-dpux(:,1)*fi(1,2)-dpux(:,4)*fi(2,2)-dpux(:,5)*fi(3,2)
  dpun(:,5)=dpun(:,5)-dpux(:,4)*fi(1,2)-dpux(:,2)*fi(2,2)-dpux(:,6)*fi(3,2)
  dpun(:,6)=dpun(:,6)-dpux(:,5)*fi(1,2)-dpux(:,6)*fi(2,2)-dpux(:,3)*fi(3,2)

  dpun(:,7)=dpun(:,7)-dpux(:,1)*fi(1,3)-dpux(:,4)*fi(2,3)-dpux(:,5)*fi(3,3)
  dpun(:,8)=dpun(:,8)-dpux(:,4)*fi(1,3)-dpux(:,2)*fi(2,3)-dpux(:,6)*fi(3,3)
  dpun(:,9)=dpun(:,9)-dpux(:,5)*fi(1,3)-dpux(:,6)*fi(2,3)-dpux(:,3)*fi(3,3)

  duu=duun
  dpu=dpun

  return
end subroutine dineohooke2


subroutine dineohooke1up(stype,Duu,Dup,Dut,Dpu,Dpp,es,d,ef,e,th)
  implicit none
  double precision                :: Duu(:,:,:), Dup(:,:,:), Dut(:)
  double precision                :: Dpu(:,:,:), Dpp(:,:,:)
  double precision                :: es(:,:), d(:,:)
  double precision                :: ef(:,:), e(:,:), th
  character(len=*)                :: stype
  integer                         :: gp,nr_gp

  nr_gp=size(Duu,3)

  do gp=1,nr_gp
    call dineohooke2up(stype,Duu(:,:,gp),Dup(:,:,gp),Dut(gp),Dpu(:,:,gp),Dpp(:,:,gp), &
                    es(:,gp),d(:,gp),ef(:,gp),e(:,gp),th)
  enddo

  return
end subroutine dineohooke1up


subroutine dineohooke2up(stype,Duu,Dup,Dut,Dpu,Dpp,es,dd,ef,e,th)
  implicit none
  double precision                :: Duu(:,:), Dup(:,:), Dut, Dpu(:,:), Dpp(:,:)
  double precision                :: es(:), dd(:), ef(:), e(:), th
  character(len=*)                :: stype
  
  double precision                :: f(3,3), fi(3,3)
  integer                         :: i, j, k, l, ii, jj
  
  double precision                :: Duux(6,6), Dpux(3,6), duun(6,9), dpun(3,9)
  double precision                :: efn(3), efp(9), s(3,3), trs

  integer, parameter              :: d_list(2,6)=[(/1,1/),(/2,2/),(/3,3/),&
                                                  (/1,2/),(/1,3/),(/2,3/)]
  integer, parameter              :: f_list(2,9)=[(/1,1/),(/1,2/),(/1,3/),&
                                                  (/2,1/),(/2,2/),(/2,3/),&
                                                  (/3,1/),(/3,2/),(/3,3/)]
  double precision :: esx1(6), esx2(6), ddx(3), efx(9), pert, ex(3)

  if (size(ef).eq.4) then
    f(1,:)=(/ef(1), ef(2), zero/)
    f(2,:)=(/ef(3), ef(4), zero/)
    f(3,:)=(/ zero,  zero,  one/)
  else
    f(1,:)=(/ef(1), ef(2), ef(3)/)
    f(2,:)=(/ef(4), ef(5), ef(6)/)
    f(3,:)=(/ef(7), ef(8), ef(9)/)
  endif

! Deviatoric stress
  trs=(es(1)+es(2)+es(3))/3d0
  s(1,:)=(/es(1)-trs,es(4)    ,es(5)/)
  s(2,:)=(/es(4)    ,es(2)-trs,es(6)/)
  s(3,:)=(/es(5)    ,es(6)    ,es(3)-trs/)

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

  efn(1)=e(1)*fi(1,1)+e(2)*fi(2,1)+e(3)*fi(3,1)
  efn(2)=e(1)*fi(1,2)+e(2)*fi(2,2)+e(3)*fi(3,2)
  efn(3)=e(1)*fi(1,3)+e(2)*fi(2,3)+e(3)*fi(3,3)

  efp=(/fi(1,1),fi(1,2),fi(1,3),fi(2,1),fi(2,2),fi(2,3),fi(3,1),fi(3,2),fi(3,3)/)
  call dneohooke_elect('ul',duux,dup,dut,dpux,dpp,efp,e,th)
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


! due to the coupling part
  duun(:,1)=duun(:,1)-dup(:,1)*efn(1)
  duun(:,2)=duun(:,2)-dup(:,2)*efn(1)
  duun(:,3)=duun(:,3)-dup(:,3)*efn(1)

  duun(:,4)=duun(:,4)-dup(:,1)*efn(2)
  duun(:,5)=duun(:,5)-dup(:,2)*efn(2)
  duun(:,6)=duun(:,6)-dup(:,3)*efn(2)

  duun(:,7)=duun(:,7)-dup(:,1)*efn(3)
  duun(:,8)=duun(:,8)-dup(:,2)*efn(3)
  duun(:,9)=duun(:,9)-dup(:,3)*efn(3)

! Dpu matrix for inverse motion
  do i=1,3
    do jj=1,9
      k=f_list(1,jj)
      l=f_list(2,jj)
      dpun(i,jj)=dd(i)*fi(l,k)-fi(i,k)*dd(l)-dpp(i,l)*efn(k)
    end do
  end do

! order 11 12 13 21 22 23 31 32 33
  dpun(:,1)=dpun(:,1)-dpux(:,1)*fi(1,1)-dpux(:,4)*fi(2,1)-dpux(:,5)*fi(3,1)
  dpun(:,2)=dpun(:,2)-dpux(:,4)*fi(1,1)-dpux(:,2)*fi(2,1)-dpux(:,6)*fi(3,1)
  dpun(:,3)=dpun(:,3)-dpux(:,5)*fi(1,1)-dpux(:,6)*fi(2,1)-dpux(:,3)*fi(3,1)

  dpun(:,4)=dpun(:,4)-dpux(:,1)*fi(1,2)-dpux(:,4)*fi(2,2)-dpux(:,5)*fi(3,2)
  dpun(:,5)=dpun(:,5)-dpux(:,4)*fi(1,2)-dpux(:,2)*fi(2,2)-dpux(:,6)*fi(3,2)
  dpun(:,6)=dpun(:,6)-dpux(:,5)*fi(1,2)-dpux(:,6)*fi(2,2)-dpux(:,3)*fi(3,2)

  dpun(:,7)=dpun(:,7)-dpux(:,1)*fi(1,3)-dpux(:,4)*fi(2,3)-dpux(:,5)*fi(3,3)
  dpun(:,8)=dpun(:,8)-dpux(:,4)*fi(1,3)-dpux(:,2)*fi(2,3)-dpux(:,6)*fi(3,3)
  dpun(:,9)=dpun(:,9)-dpux(:,5)*fi(1,3)-dpux(:,6)*fi(2,3)-dpux(:,3)*fi(3,3)

  duu=duun
  dpu=dpun

! pertubation method
! Tangents have been check they are correct
if (1.eq.0) then
!  write(*,*)'duu_ana ',duu
  pert=1e-8  
  trs=(es(1)+es(2)+es(3))/3d0
  esx1=es
  esx1(1:3)=esx1(1:3)-trs
  do ii=1,9
    efx=ef
    efx(ii)=efx(ii)+pert    
    call neohooke_elect('Cauchy-inv',esx2,ddx,efx,e,th) 
    trs=(es(1)+es(2)+es(3))/3d0
    esx2(1:3)=esx2(1:3)-trs
    duu(:,ii)=(esx2-esx1)/pert
  enddo
!  write(*,*)'duu_per ',duu

!  write(*,*)'dup_ana ',dup
  pert=1e-8  
  trs=(es(1)+es(2)+es(3))/3d0
  esx1=es
  esx1(1:3)=esx1(1:3)-trs
  do ii=1,3
    ex=e
    ex(ii)=ex(ii)+pert    
    call neohooke_elect('Cauchy-inv',esx2,ddx,ef,ex,th) 
    trs=(es(1)+es(2)+es(3))/3d0
    esx2(1:3)=esx2(1:3)-trs
    dup(:,ii)=(esx2-esx1)/pert
  enddo
!  write(*,*)'dup_per ',dup

!  write(*,*)'dpu_ana ',dpu
  pert=1e-8  
  do ii=1,9
    efx=ef
    efx(ii)=efx(ii)+pert    
    call neohooke_elect('Cauchy-inv',esx2,ddx,efx,e,th) 
    trs=(es(1)+es(2)+es(3))/3d0
    dpu(:,ii)=(ddx-dd)/pert
  enddo
!  write(*,*)'dpu_per ',dpu

!  write(*,*)'dpp_ana ',dpp
  pert=1e-8  
  do ii=1,3
    ex=e
    ex(ii)=ex(ii)+pert    
    call neohooke_elect('Cauchy-inv',esx2,ddx,ef,ex,th) 
    dpp(:,ii)=(ddx-dd)/pert
  enddo
!  write(*,*)'dpp_per ',dpp
endif

  return
end subroutine dineohooke2up

end module mater_hyperel_elect
