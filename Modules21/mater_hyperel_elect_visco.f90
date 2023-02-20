module mater_hyperel_elect_visco

! last modified
! A. Ask 2011-10-14
!  - initial code hyper electro-visco-elastic material. 
!    Used electro-elastic module as starting point
! M. Ristinmaa 2011-10-14
!  - copy of intVar in accept routine
! ------------------------------------------------------------------------------

! ToDo
! move vis_update into stress call ??

use mater_large
use some_constants
use matrix_util

implicit none

double precision                  :: mpK, mpMy, mpCm, mpCe, mpEo
integer                           :: nvis, old, new
double precision, allocatable     :: mpBvis(:), mpGvis(:), intVar(:,:,:,:,:)
private mpK, mpMy, mpCm, mpCe, mpEo, nvis, mpBvis, mpGvis, intVar, old, new

interface dneohooke_vis_elect
  module procedure dneohooke1
  module procedure dneohooke2
end interface

private dneohooke1, dneohooke2

interface neohooke_vis_elect
  module procedure neohooke1
  module procedure neohooke2
end interface

private neohooke1, neohooke2

interface vis_update
  module procedure update1
  module procedure update2
end interface

private update1, update2

interface getIntVar
  module procedure getIntVar1
  module procedure getIntVar2
end interface

private getIntVar1, getIntVar2

interface dineohooke_vis_elect
  module procedure dineohooke1
  module procedure dineohooke2
end interface

private inv3

!-------------------------------------------------------------------------------
contains

subroutine neohooke_vis_elect_init(mp, nelem, nr_gp)
  implicit none
  double precision                :: mp(:)
  integer                         :: nelem, nr_gp, ierr, ie, iv, igp
  mpK =mp(1)
  mpMy=mp(2)
  mpCm=mp(3)
  mpCe=mp(4)
  ! This part needs debugging
  mpEo=zero
!  mpEo=mp(5)
  nvis = mp(6)
  allocate(mpBvis(nvis), stat=ierr)
  allocate(mpGvis(nvis), stat=ierr)
  
  mpBvis = mp(7:size(mp):2)
  mpGvis = mp(8:size(mp):2)
  
  ! allcoate internal variables (hidden from user!)
  allocate(intVar(nr_gp,nvis,6,nelem,2), stat=ierr)
  ! initiate internal variables
  do ie = 1,nelem
      do iv = 1,nvis
          do igp = 1,nr_gp
              intVar(igp, iv, :, ie, 1) = (/1d0, 1d0, 1d0, 0d0, 0d0, 0d0/)
              intVar(igp, iv, :, ie, 2) = (/1d0, 1d0, 1d0, 0d0, 0d0, 0d0/)
          end do
      end do
  end do
  old = 1
  new = 2

  return
end subroutine neohooke_vis_elect_init

! Call for several gauss points
subroutine neohooke1(stype,stress,edisp,ef,de,ie)
  implicit none
  double precision                :: stress(:,:), edisp(:,:), ef(:,:), de(:,:)
  character(len=*)                :: stype
  integer                         :: ie, gp, nr_gp

  nr_gp=size(ef,2)

  do gp=1,nr_gp
    call neohooke2(stype,stress(:,gp), intVar(gp,:,:,ie,new), edisp(:,gp),ef(:,gp), de(:,gp))
  enddo

  return
end subroutine neohooke1

! Call for one gauss point
subroutine neohooke2(stype,stress,intVar_gp,edisp,ef,e)
  implicit none
  double precision                :: stress(:), intVar_gp(:,:), edisp(:), ef(:), e(:)
  character(len=*)                :: stype
  
  double precision                :: F(3,3), b(3,3), id(3,3), invF(3,3)
  double precision                :: J, Jm53, tr_b, S(3,3), Cv(3,3), invCv(3,3), CC(3,3), tr_CC
  double precision                :: hat_e(3), tilde_e(3), tr_e, tr_hat_e
  double precision                :: dy_e(3,3), dy_hat_e(3,3)
  integer                         :: iv
  
if (stype.ne.'debugg') then
  if (size(ef).eq.4) then
    F(1,:)=(/ef(1), ef(2), zero/)
    F(2,:)=(/ef(3), ef(4), zero/)
    F(3,:)=(/ zero,  zero,  one/)
  else
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
  endif

! No inverse motion with viscosity (as of yet)

! Check if a inverse deformation gradient is provided as input
  if (stype.eq.'Cauchy-inv') then
    call inv3(invF,F)
    F=invF
  end if


  b=matmul(F,transpose(F))
  J=det3(F)

else
    F(1,:)=(/ef(1), ef(4), ef(6)/)
    F(2,:)=(/ef(4), ef(2), ef(5)/)
    F(3,:)=(/ef(6), ef(5), ef(3)/)
   call inv3(b,F)
   J=dsqrt(det3(b))
!write(*,*)'b',b
!write(*,*)'J',J

end if

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
! Viscous contribution
  do iv = 1,nvis
      Cv(1,:) = (/intVar_gp(iv,1), intVar_gp(iv,4), intVar_gp(iv,5)/)
      Cv(2,:) = (/intVar_gp(iv,4), intVar_gp(iv,2), intVar_gp(iv,6)/)
      Cv(3,:) = (/intVar_gp(iv,5), intVar_gp(iv,6), intVar_gp(iv,3)/)
      call inv3(invCv,Cv)
      CC = matmul(matmul(transpose(F),F),invCv)
      tr_CC = CC(1,1)+CC(2,2)+CC(3,3)
      S = S+mpBvis(iv)*mpMy*Jm53*(matmul(matmul(F,invCv),transpose(F))-onethird*tr_CC*Id)
  end do

! Electric displacement (spatial)
  edisp=-two/J*mpCe*hat_e-two*mpCm*Jm53*tilde_e+mpEo*e

  if (stype.eq.'Cauchy') then
  elseif (stype.eq.'Cauchy-inv') then
  elseif (stype.eq.'debugg') then
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


subroutine dneohooke1(stype,Duu,Dup,Dpu,Dpp,ef,e,ie,dt)
  implicit none
  double precision                :: Duu(:,:,:), Dup(:,:,:)
  double precision                :: Dpu(:,:,:), Dpp(:,:,:)
  double precision                :: ef(:,:), e(:,:),dt
  character(len=*)                :: stype
  integer                         :: ie, gp,nr_gp

  nr_gp=size(Duu,3)

  do gp=1,nr_gp
    call dneohooke2(stype,Duu(:,:,gp),Dup(:,:,gp),Dpu(:,:,gp),Dpp(:,:,gp), &
                    ef(:,gp),e(:,gp),intVar(gp,:,:,ie,new),dt)
  enddo

  return
end subroutine dneohooke1

subroutine dneohooke2(stype,Duu,Dup,Dpu,Dpp,ef,e,intVar_gp,dt)
  implicit none
  double precision                :: Duu(:,:), Dup(:,:), Dpu(:,:), Dpp(:,:)
  double precision                :: ef(:), e(:), intVar_gp(:,:),dt
  character(len=*)                :: stype

! No inverse motion with viscosity

! Check if inverse motion problem
!  if (stype.eq.'inv') then
!    call dneohooke2b(stype,Duu,Dup,Dpu,Dpp,ef,e)
!  else
!    call dneohooke2a(stype,Duu,Dup,Dpu,Dpp,ef,e)
!  end if
  call dneohooke2a(stype,Duu,Dup,Dpu,Dpp,ef,e,intVar_gp,dt)
  return
end subroutine dneohooke2


subroutine dneohooke2a(stype,Duu,Dup,Dpu,Dpp,ef,e,intVar_gp,dt)
  implicit none
  double precision                :: Duu(:,:), Dup(:,:), Dpu(:,:), Dpp(:,:)
  double precision                :: ef(:), e(:), intVar_gp(:,:), dt
  character(len=*)                :: stype
  
  double precision                :: F(3,3), b(3,3), Id(3,3), Cbar(3,3), CCC(3,3) 
  double precision                :: Cv(3,3), invCv(3,3), CC(3,3), FCF(3,3), FCv(3,3)
  double precision                :: Ja, Jm23, tr_b, Ds(21), tr_CC, det_Cv
  integer                         :: i, j, k, l, ii, jj, iv
  
  double precision                :: hat_e(3), tilde_e(3), tr_e, tr_hat_e
  double precision                :: dy_e(3,3), dy_hat_e(3,3)
  double precision                :: c1, c2, c3, c4, c5 ,c6, c7
  double precision                :: Duv(6,6), Km(6,6), Am(6,6), Bm(6,6)

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
  Cbar = Jm23*transpose(b)

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
  Duv = 0d0
  do iv=1,nvis
    Am = 0d0
    Bm = 0d0
    Cv(1,:) = (/intVar_gp(iv,1), intVar_gp(iv,4), intVar_gp(iv,5)/)
    Cv(2,:) = (/intVar_gp(iv,4), intVar_gp(iv,2), intVar_gp(iv,6)/)
    Cv(3,:) = (/intVar_gp(iv,5), intVar_gp(iv,6), intVar_gp(iv,3)/)
    call inv3(invCv,Cv) 
    det_Cv = det3(Cv)
    FCF = matmul(matmul(F,invCv),transpose(F)) 
    CC = matmul(matmul(transpose(F),F),invCv)
    CCC = matmul(invCv,CC)
    FCv = matmul(F, invCv)
    tr_CC = CC(1,1)+CC(2,2)+CC(3,3)
    Km = 0d0
    call kmat2(Cv, invCv, Cbar, det_Cv, mpBvis(iv), mpGvis(iv), dt, tr_CC, Km)
    do ii=1,6
      i=d_list(1,ii)
      j=d_list(2,ii)
      do jj=1,6
        k=d_list(1,jj)
        l=d_list(2,jj)
        
        Duv(ii,jj) = Duv(ii,jj)+c2*mpBvis(iv)*(Id(i,j)*FCF(k,l)+Id(k,l)*FCF(i,j) &
                   -onethird*tr_CC*Id(i,j)*Id(k,l)-onehalf*tr_CC*(Id(i,k)*Id(j,l)+Id(i,l)*Id(j,k)))
        Am(ii,jj) = -two*mpBvis(iv)*mpMy*Jm23*(onehalf*(FCv(i,k)*FCv(j,l)+FCv(i,l)*FCv(j,k))&
                  -onethird*Id(i,j)*CCC(k,l))
        Bm(ii,jj) = Jm23*onehalf*(F(k,i)*F(l,j)+F(k,j)*F(l,i))-onethird*Cbar(i,j)*Id(k,l)
      end do
    end do
    Am(:,4:6) = two*Am(:,4:6)
    Bm(4:6,:) = two*Bm(4:6,:)
    Duv = Duv + matmul(matmul(Am,Km),Bm)
  end do
  else
    stop 'stype not implemented'
  end if
  Duu = Duu + Duv
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


! Material tangent stiffnesses for the inverse motion problem
subroutine dneohooke2b(stype,Duu,Dup,Dpu,Dpp,ef,e)
  implicit none
  double precision                :: Duu(:,:), Dup(:,:), Dpu(:,:), Dpp(:,:)
  double precision                :: ef(:), e(:)
  character(len=*)                :: stype
  
  double precision                :: f(3,3), c(3,3), ic(3,3), Id(3,3), icc(3,3)
  double precision                :: Ja, Jm23, tr_ic, Ds(21)
  integer                         :: i, j, k, l, ii, jj
  
  double precision                :: hat_e(3), tilde_e(3), tr_e, tr_hat_e
  double precision                :: dy_e(3,3), dy_hat_e(3,3)
  double precision                :: q1, q2, q3, q4

  double precision                :: stress1(6), edisp1(3)
  double precision                :: stress2(6), edisp2(3), cm(6)

  integer, parameter              :: d_list(2,6)=[(/1,1/),(/2,2/),(/3,3/),&
                                                  (/1,2/),(/1,3/),(/2,3/)]
 
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
  call inv3(ic,c)
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

      Duu(ii,jj)=q4*Id(i,j)*ic(k,l)+5d0/6d0*q3*ic(i,j)*ic(k,l) &
                -onehalf*q3*(ic(i,k)*ic(l,j)+ic(i,l)*ic(k,j)) &
                +1d0/3d0*q3*Id(i,j)*icc(k,l) & 
                +5d0/3d0*q2*hat_e(i)*hat_e(j)*ic(k,l) &
                +2d0/3d0*q2*(tilde_e(k)*hat_e(l)+hat_e(k)*tilde_e(l))*Id(i,j) &
                -q2*(ic(i,k)*hat_e(l)*hat_e(j)+ic(i,l)*hat_e(k)*hat_e(j) &
                    +ic(j,k)*hat_e(l)*hat_e(i)+ic(j,l)*hat_e(k)*hat_e(i))

    end do

    do k=1,3

      Dpu(k,ii)=-q1*hat_e(k)*ic(i,j)-5d0/3d0*q2*tilde_e(k)*ic(i,j) &
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
 
  return
end subroutine dneohooke2b

! Update - element level (several gauss points)
subroutine update1(stype,ef, dt, ie)
  implicit none
  double precision      :: ef(:,:), dt, Cv_gp(nvis,6)
  character(len=*)                :: stype
  integer               :: ie, nr_gp, gp
  
  nr_gp = size(ef,2)
  if(nvis>0) then
    do gp=1,nr_gp
      call update2(stype,ef(:,gp), dt, intVar(gp,:,:,ie,old), Cv_gp)
      intVar(gp,:,:,ie,new) = Cv_gp
    end do 
  end if
  return
end subroutine update1

! Update - gauss point level
subroutine update2(stype,ef, dt, Cv0_gp, Cv_gp)
  implicit none
  double precision      :: ef(:), dt, Cv0_gp(:,:), Cv_gp(:,:)
  character(len=*)                :: stype
  double precision      :: gvec(6), fvec(6), Cv0(6), Cv(3,3) 
  double precision      :: F(3,3), K(6,6), Cbar(3,3), Jvol, Jcoeff
  double precision      :: CC(3,3), tr_CC, invCv(3,3),det_Cv, vpar
  double precision      :: Cbar_vec(6), Cvn(6), dCv(6), res ,invF(3,3)
  integer               :: iv, ipiv(6), info, itemp
  double precision, parameter   :: Id(6) = (/1d0, 1d0, 1d0, 0d0, 0d0, 0d0/)
  
  double precision, external :: DNRM2
  
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


  Jvol = det3(F)
  Jcoeff = Jvol**(-2d0*onethird)
  Cbar = Jcoeff*matmul(transpose(F),F)
  Cbar_vec = (/Cbar(1,1), Cbar(2,2), Cbar(3,3), Cbar(1,2), Cbar(1,3), Cbar(2,3)/)
  
  do iv=1,nvis
    vpar = onehalf*mpBvis(iv)*mpMy*mpGvis(iv)
    Cv0 = Cv0_gp(iv,:)  
    Cv = Cbar
    det_Cv = det3(Cv)
    call inv3(invCv, Cv)
    CC = matmul(Cbar, invCv)
    tr_CC = CC(1,1) + CC(2,2) + CC(3,3)
    Cvn = Cbar_vec
    fvec = vpar*(Cbar_vec-onethird*tr_CC*Cvn)
    gvec = Cvn-Cv0-dt*fvec+1d2*(det_Cv-1d0)*Id
    res = DNRM2(6, gvec, 1)
    do while(res>1d-6)
      call kmat1(Cv, invCv, Cbar, det_Cv, mpBvis(iv), mpGvis(iv), dt, tr_CC, K)
      dCv = -gvec
      call dgesv(6,1,K,6,ipiv,dCv,6,info)
      Cvn = Cvn+dCv
      Cv(1,:) = (/Cvn(1), Cvn(4), Cvn(5)/)
      Cv(2,:) = (/Cvn(4), Cvn(2), Cvn(6)/) 
      Cv(3,:) = (/Cvn(5), Cvn(6), Cvn(3)/)
      det_Cv = det3(Cv)
      call inv3(invCv, Cv)
      CC = matmul(Cbar, invCv)
      tr_CC = CC(1,1) + CC(2,2) + CC(3,3)
      fvec = vpar*(Cbar_vec-onethird*tr_CC*Cvn)
      gvec = Cvn-Cv0-dt*fvec+1d2*(det_Cv-1d0)*Id
      res = DNRM2(6, gvec, 1)
    end do
    Cv_gp(iv,:) = Cvn
  end do
  return
end subroutine update2

subroutine kmat1(Cv, invCv, Cbar, det_Cv, bpar, gpar, dt, c0, Km)
! Use in update routine
  implicit none
  double precision              :: Cv(:,:), invCv(:,:), Cbar(:,:), det_Cv
  double precision              :: bpar, gpar, dt, c0, Km(6,6)
  double precision              :: c1, Id(3,3), CC(3,3)
  integer                       :: ii, jj, i, j, k, l
  integer, parameter            :: d_list(2,6)=[(/1,1/),(/2,2/),(/3,3/),&
                                                (/1,2/),(/1,3/),(/2,3/)]
  Id(1,:) = (/1d0, 0d0, 0d0/)
  Id(2,:) = (/0d0, 1d0, 0d0/)
  Id(3,:) = (/0d0, 0d0, 1d0/)                                      
  c1 = onehalf*onethird*dt*bpar*gpar*mpMy 
  CC = matmul(matmul(invCv,Cbar), invCv) 
  Km = 0d0                                 
  do ii = 1,6
    i = d_list(1,ii)
    j = d_list(2,ii)
    do jj = 1,6
      k = d_list(1,jj)
      l = d_list(2,jj)
      Km(ii,jj) = onehalf*(1d0+c1*c0)*(Id(i,k)*Id(j,l)+Id(i,l)*Id(j,k))&
                -c1*Cv(i,j)*CC(k,l)+1d2*det_Cv*Id(i,j)*invCv(k,l)
    end do
  end do
  Km(:,4:6) = Km(:,4:6)*2d0
  return
end subroutine kmat1

subroutine kmat2(Cv, invCv, Cbar, det_Cv, bpar, gpar, dt, c0, Km)
! Use to calculate Dats
  implicit none
  double precision              :: Cv(:,:), invCv(:,:), Cbar(:,:), det_Cv
  double precision              :: bpar, gpar, dt, c0, Km(6,6), Kiso(6,6)
  double precision              :: c1, Id(3,3), CC(3,3)
  integer                       :: ii, jj, i, j, k, l, ipiv(6), info
  integer, parameter            :: d_list(2,6)=[(/1,1/),(/2,2/),(/3,3/),&
                                                (/1,2/),(/1,3/),(/2,3/)]
  double precision, external :: DNRM2
  
  Id(1,:) = (/1d0, 0d0, 0d0/)
  Id(2,:) = (/0d0, 1d0, 0d0/)
  Id(3,:) = (/0d0, 0d0, 1d0/)                                      
  c1 = onehalf*onethird*dt*bpar*gpar*mpMy 
  CC = matmul(matmul(invCv,Cbar), invCv) 
  Km = 0d0                 
  Kiso = 0d0                  
  do ii = 1,6
    i = d_list(1,ii)
    j = d_list(2,ii)
    do jj = 1,6
      k = d_list(1,jj)
      l = d_list(2,jj)
      Km(ii,jj) = onehalf*(1d0+c1*c0)*(Id(i,k)*Id(j,l)+Id(i,l)*Id(j,k))&
                -c1*Cv(i,j)*CC(k,l)+1d2*det_Cv*Id(i,j)*invCv(k,l)
      Kiso(ii,jj) = dt*onehalf*bpar*gpar*mpMy*(-onethird*invCv(k,l)*Cv(i,j)&
                  +onehalf*(id(i,k)*id(j,l)+id(i,l)*id(j,k)))
    end do
  end do
  Km(:,4:6) = Km(:,4:6)*2d0
  CALL DGESV(6, 6, Km, 6, ipiv, Kiso, 6, info)
  Km = Kiso
  return
end subroutine kmat2

subroutine state_accept()
  implicit none

! This copy must be done otherwise the tangent stiffness
! in the first iteration is calculated at the n-1 state
  intVar(:,:,:,:,Old)=intVar(:,:,:,:,New)

  if(new==1) then
    if(old==1) then
      stop 'Error new = old in state_accept'
    end if
    new = 2
    old = 1
  else if(new==2) then
    if(old==2) then
      stop 'Error new = old in state_accept'
    end if
    new = 1
    old = 2
  end if
  return
end subroutine state_accept

subroutine getIntVar1(ie, Cv)
! Get internal variables in one element
  implicit none
  integer               :: ie
  double precision      :: Cv(:,:,:)
  
  Cv = intVar(:,:,:,ie,new)
  return
end subroutine getIntVar1

subroutine getIntVar2(ie, Cv)
! Get internal variables for many (possibly all) elements
  implicit none
  integer               :: ie(:)
  double precision      :: Cv(:,:,:,:)

  Cv = intVar(:,:,:,ie,new)
  return
end subroutine getIntVar2

!('inv',Duu,Dup,Dpu,Dpp,es(:,:,ie),dd(:,:,ie),dg(:,:,ie),ef(:,:,ie),ie,dt)
subroutine dineohooke1(stype,Duu,Dup,Dpu,Dpp,es,d,ef,e,ie,dt)
  implicit none
  double precision                :: Duu(:,:,:), Dup(:,:,:)
  double precision                :: Dpu(:,:,:), Dpp(:,:,:)
  double precision                :: es(:,:), d(:,:)
  double precision                :: ef(:,:), e(:,:), dt
  character(len=*)                :: stype
  integer                         :: gp,nr_gp, ie, ii

  double precision :: dgx(9,8), ddx(3,8), pert, esx(6,8)

  nr_gp=size(Duu,3)

  do gp=1,nr_gp
    call dineohooke2(stype,Duu(:,:,gp),Dup(:,:,gp),Dpu(:,:,gp),Dpp(:,:,gp), &
                    es(:,gp),d(:,gp),ef(:,gp),e(:,gp),ie,dt,gp)
  enddo

  pert=1d-8
  do ii=1,9
    dgx=ef(:,:)
    dgx(ii,:)=dgx(ii,:)+pert
    esx=es
    call vis_update('Cauchy-inv',dgx,dt,ie)
    call neohooke_vis_elect('Cauchy-inv',esx,ddx,dgx,e(:,:),ie)
    do gp=1,8
       duu(:,ii,gp)=(esx(:,gp)-es(:,gp))/pert
    enddo
  enddo

  return
end subroutine dineohooke1



subroutine dineohooke2(stype,Duu,Dup,Dpu,Dpp,es,dd,ef,e,ie,dt,gp)
  implicit none
  double precision                :: Duu(:,:), Dup(:,:), Dpu(:,:), Dpp(:,:)
  double precision                :: es(:), dd(:), ef(:), e(:), dt
  character(len=*)                :: stype
  integer                         :: ie, gp
  
  double precision                :: f(3,3), fi(3,3)
  integer                         :: i, j, k, l, ii, jj
  
  double precision                :: Duux(6,6), Dpux(3,6), duun(6,9), dpun(3,9)
  double precision                :: efn(3), efp(9), s(3,3)

  integer, parameter              :: d_list(2,6)=[(/1,1/),(/2,2/),(/3,3/),&
                                                  (/1,2/),(/1,3/),(/2,3/)]
  integer, parameter              :: f_list(2,9)=[(/1,1/),(/1,2/),(/1,3/),&
                                                  (/2,1/),(/2,2/),(/2,3/),&
                                                  (/3,1/),(/3,2/),(/3,3/)]
  double precision :: pert, esx1(6),esx2(6),  efx(9), Cv_gp(nvis,6)
  double precision :: ex(3), ddx(3)


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
!  call dneohooke_elect('ul',duux,dup,dpux,dpp,efp,e)
!  call dneohooke_vis_elect('ul',Duux,Dup,Dpux,Dpp,efp,e,ie,dt)
  call dneohooke2('ul',Duux,Dup,Dpux,Dpp, &
                    efp,e,intVar(gp,:,:,ie,new),dt)
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
if (1.eq.0) then
!  write(*,*)'duu_ana ',duu
  pert=1e-8
  esx1=es
  do ii=1,9
    efx=ef
    efx(ii)=efx(ii)+pert    
    call update2('Cauchy-inv',efx, dt, intVar(gp,:,:,ie,old), Cv_gp)
    call neohooke2('Cauchy-inv',esx2,Cv_gp,ddx,efx,e)
!    call neohooke_elect('Cauchy-inv',esx2,ddx,efx,e,th) 
    duu(:,ii)=(esx2-esx1)/pert
  enddo
!  write(*,*)'duu_per ',duu
endif

if (1.eq.0) then
  write(*,*)'dup_ana ',dup
  pert=1e-8  
  esx1=es
  do ii=1,3
    ex=e
    ex(ii)=ex(ii)+pert    
!    call neohooke_elect('Cauchy-inv',esx2,ddx,ef,ex,th) 
    call neohooke2('Cauchy-inv',esx2,intVar(gp,:,:,ie,new),ddx,ef,ex)
    dup(:,ii)=(esx2-esx1)/pert
  enddo
  write(*,*)'dup_per ',dup

endif

if (1.eq.0) then
!  write(*,*)'dpu_ana ',dpu
  pert=1e-8  
  do ii=1,9
    efx=ef
    efx(ii)=efx(ii)+pert    
!    call neohooke_elect('Cauchy-inv',esx2,ddx,efx,e,th) 
    call update2('Cauchy-inv',efx, dt, intVar(gp,:,:,ie,old), Cv_gp)
    call neohooke2('Cauchy-inv',esx2,Cv_gp,ddx,efx,e)
    dpu(:,ii)=(ddx-dd)/pert
  enddo
!  write(*,*)'dpu_per ',dpu
endif

if (1.eq.0) then
!  write(*,*)'dpp_ana ',dpp
  pert=1e-8  
  do ii=1,3
    ex=e
    ex(ii)=ex(ii)+pert    
!    call neohooke_elect('Cauchy-inv',esx2,ddx,ef,ex,th) 
    call neohooke2('Cauchy-inv',esx2,intVar(gp,:,:,ie,new),ddx,ef,ex)
    dpp(:,ii)=(ddx-dd)/pert
  enddo
!  write(*,*)'dpp_per ',dpp
endif


  return
end subroutine dineohooke2




end module mater_hyperel_elect_visco
