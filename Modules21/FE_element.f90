module FE_element

! Last rev.:
! A. Dalklint 2020-04-28
!
! ------------------------------------------------------------------
! Module FE_element contains generic interfaces for 3D/2D elements
! Total lagrangian implementation
! ------------------------------------------------------------------

use elem_large_cont_3d
use elem_large_cont_2d
use elem_flow_3d
use elem_flow_2d

implicit none

interface elem_stiff
  module procedure stiff
end interface
private stiff


interface elem_stiff_gamma
  module procedure stiff_gamma
end interface
private stiff_gamma


interface elem_stiff_gamma_sens
  module procedure stiff_gamma_sens
end interface
private stiff_gamma_sens

interface elem_mass
  module procedure mass
end interface
private mass

interface elem_mass_sens
  module procedure mass_sens
end interface
private mass_sens

interface elem_force
  module procedure force
end interface
private force


interface elem_force_gamma
  module procedure force_gamma
end interface
private force_gamma


interface elem_force_gamma_sens
  module procedure force_gamma_sens
end interface
private force_gamma_sens


interface elem_flow_stiff
  module procedure stiff_flow
end interface
private stiff_flow


interface elem_flow_mass
  module procedure mass_flow
end interface
private mass_flow


interface elem_dg
  module procedure defgrad
end interface
private defgrad


interface elem_gp
  module procedure gp
end interface
private gp


interface elem_flow_bf
  module procedure bf_flow_val
  module procedure bf_flow_vec
end interface
private bf_flow_val, bf_flow_vec


interface elem_bf
  module procedure bf
end interface
private bf


interface elem_bf_sens
  module procedure bf_sens
end interface
private bf_Sens



interface elem_vol
  module procedure vol
  module procedure vol_vec 
end interface
private vol, vol_vec


interface elem_Etilde1
  module procedure Etilde1
end interface
private Etilde1


interface elem_Etilde2
  module procedure Etilde2 
end interface
private Etilde2


!------------------------------------------------------------------------------
contains


!------------------------------------------------------------------------------
!                         Stiffness matrices
!------------------------------------------------------------------------------

subroutine stiff(ctype,ke,coord,t,D,ed,es)
   implicit none
   double precision                :: ke(:,:), coord(:,:), D(:,:,:)
   double precision                :: ed(:), es(:,:), t
   character(2)                    :: ctype

   if (ctype.eq.'tl') then
      if (t.gt.0d0) then
         call c2dtl4_e(ke,coord,t,D,ed,es) 
      else
         call c3dtl8_e(ke,coord,D,ed,es) 
      endif
   elseif (ctype.eq.'ul') then
      if (t.gt.0d0) then
         call c2dul4_e(ke,coord,t,D,ed,es) 
      else
         call c3dul8_e(ke,coord,D,ed,es) 
      endif
   else
      write(*,*) 'Error in Newton solver: You must choose between tl or ul!'
      stop
   endif

end subroutine stiff


subroutine stiff_gamma(ctype,ke,coord,t,D,ed,es,gammagp)
   implicit none
   double precision                :: ke(:,:), coord(:,:), D(:,:,:)
   double precision                :: ed(:), es(:,:), t, gammagp(:)
   character(2)                    :: ctype

   if (ctype.eq.'tl') then
      if (t.gt.0d0) then
         call c2dtl4_egamma(ke,coord,t,D,ed,es,gammagp) 
      else
         call c3dtl8_egamma(ke,coord,D,ed,es,gammagp)
      endif
   elseif (ctype.eq.'ul') then
      if (t.gt.0d0) then
         call c2dul4_e(ke,coord,t,D,ed,es) 
      else
         call c3dul8_e(ke,coord,D,ed,es) 
      endif
   else
      write(*,*) 'Error in Newton solver: You must choose between tl or ul!'
      stop
   endif

end subroutine stiff_gamma


subroutine stiff_gamma_sens(ctype,ke,coord,t,D,ed,larg,rarg,es,gammagp)
   implicit none
   double precision                :: ke(:), coord(:,:), D(:,:,:)
   double precision                :: ed(:), es(:,:), t, gammagp(:)
   double precision                :: larg(:), rarg(:)
   character(2)                    :: ctype

   if (ctype.eq.'tl') then
      if (t.gt.0d0) then
         call c2dtl4_egammaSens(ke,coord,t,D,ed,larg,rarg,es,gammagp) 
      else
         call c3dtl8_egammaSens(ke,coord,D,ed,larg,rarg,es,gammagp)
      endif
   elseif (ctype.eq.'ul') then
      if (t.gt.0d0) then
         call c2dul4_eSens(ke,coord,t,D,ed,larg,rarg,es) 
      else
         call c3dul8_eSens(ke,coord,D,ed,larg,rarg,es)
      endif
   else
      write(*,*) 'Error in Newton solver: You must choose between tl or ul!'
      stop
   endif



end subroutine stiff_gamma_sens



subroutine stiff_flow(ke,coord,t,val)
   implicit none
   double precision                :: ke(:,:), coord(:,:)
   double precision                :: val, t

   if (t.gt.0d0) then
      call fl2d4_e(ke,coord,t,val)
   else
      call fl3d8_e(ke,coord,val)
   endif

end subroutine stiff_flow

!------------------------------------------------------------------------------
!                         Mass matrices
!------------------------------------------------------------------------------

subroutine mass(me,coord,t,val,rhogp)
   implicit none
   double precision                :: me(:,:), coord(:,:)
   double precision                :: rhogp(:), t, val

   if (t.gt.0d0) then
      call c2dtl4_m(me,coord,t,val,rhogp) 
   else
      call c3dtl8_m(me,coord,val,rhogp) 
   endif

end subroutine mass


subroutine mass_sens(fme,coord,t,val,rhogp,larg,rarg)
   implicit none
   double precision                :: fme(:), coord(:,:)
   double precision                :: rhogp(:), t, val
   double precision                :: larg(:), rarg(:)

   if (t.gt.0d0) then
      call dc2dtl4_m(fme,coord,t,val,rhogp,larg,rarg) 
   else
      call dc3dtl8_m(fme,coord,val,rhogp,larg,rarg) 
   endif
end subroutine mass_sens



!------------------------------------------------------------------------------
!                          Deformation gradient
!------------------------------------------------------------------------------

subroutine defgrad(dg,coord,ed)
   implicit none
   double precision                :: dg(:,:), coord(:,:), ed(:)

   ! If 3D
   if (size(coord,1).eq.3) then
      call c3dtl8_d(dg,coord,ed) 
      ! If 2D
   elseif (size(coord,1).eq.2) then
      call c2dtl4_d(dg,coord,ed) 
   endif

end subroutine defgrad


!------------------------------------------------------------------------------
!                            Mass matrices
!------------------------------------------------------------------------------


subroutine mass_flow(me,coord,t,val)
   implicit none
   double precision                :: me(:,:), coord(:,:)
   double precision                :: val, t

   if (t.gt.0d0) then
      call fl2d4_m(me,coord,t,val)
   else
      call fl3d8_m(me,coord,val)
   endif

end subroutine mass_flow


!------------------------------------------------------------------------------
!                             Force vectors 
!------------------------------------------------------------------------------



subroutine force(ctype,fe,coord,t,ed,es)
   implicit none
   double precision                :: fe(:), coord(:,:)
   double precision                :: ed(:), es(:,:), t
   character(2)                    :: ctype

   if (ctype.eq.'tl') then
      if (t.gt.0d0) then
         call c2dtl4_f(fe,coord,t,ed,es) 
      else
         call c3dtl8_f(fe,coord,ed,es) 
      endif
   elseif (ctype.eq.'ul') then
      if (t.gt.0d0) then
         call c2dul4_f(fe,coord,t,ed,es) 
      else
         call c3dul8_f(fe,coord,ed,es) 
      endif
   else
      write(*,*) 'Error in Newton solver: You must choose between tl or ul!'
      stop
   endif



end subroutine force



subroutine force_gamma(ctype,fe,coord,t,ed,es,gammagp)
   implicit none
   double precision                :: fe(:), coord(:,:)
   double precision                :: ed(:), es(:,:), t, gammagp(:)
   character(2)                    :: ctype

   if (ctype.eq.'tl') then
      if (t.gt.0d0) then
         call c2dtl4_fgamma(fe,coord,t,ed,es,gammagp) 
      else
         call c3dtl8_fgamma(fe,coord,ed,es,gammagp)
      endif
   elseif (ctype.eq.'ul') then
      if (t.gt.0d0) then
         call c2dul4_f(fe,coord,t,ed,es) 
      else
         call c3dul8_f(fe,coord,ed,es) 
      endif
   else
      write(*,*) 'Error in Newton solver: You must choose between tl or ul!'
      stop
   endif

end subroutine force_gamma


subroutine force_gamma_sens(ctype,fe,coord,t,ed,larg,es,gammagp)
   implicit none
   double precision                :: fe(:), coord(:,:)
   double precision                :: ed(:), es(:,:), t, gammagp(:)
   double precision                :: larg(:)
   character(2)                    :: ctype

   if (ctype.eq.'tl') then
      if (t.gt.0d0) then
         call c2dtl4_fgammaSens(fe,coord,t,ed,larg,es,gammagp) 
      else
         call c3dtl8_fgammaSens(fe,coord,ed,larg,es,gammagp)
      endif
   elseif (ctype.eq.'ul') then
      if (t.gt.0d0) then
         call c2dul4_fSens(fe,coord,t,ed,larg,es) 
      else
         call c3dul8_fSens(fe,coord,ed,larg,es)
      endif
   else
      write(*,*) 'Error in Newton solver: You must choose between tl or ul!'
      stop
   endif

end subroutine force_gamma_sens

!------------------------------------------------------------------------------
!                         Gauss point values 
!------------------------------------------------------------------------------

subroutine gp(rhogp,rho)
   implicit none
   double precision                ::  rho(:)
   double precision                ::  rhogp(:)

   if (size(rhogp).eq.8) then
      call fl3d8_gp(rhogp,rho) 
   elseif (size(rhogp).eq.4) then
      call fl2d4_gp(rhogp,rho) 
   else
      write(*,*) 'elem_gp is not implement for this nbr of gps'
      stop
   endif

end subroutine gp



!------------------------------------------------------------------------------
!                             Body forces - flow element
!------------------------------------------------------------------------------

subroutine bf_flow_val(fe,coord,t,val)
   implicit none
   double precision                :: coord(:,:), fe(:)
   double precision                :: val, t

   if (t.gt.0d0) then
      call fl2d4_bf(fe,coord,t,val) 
   else
      call fl3d8_bf(fe,coord,val) 
   endif

end subroutine bf_flow_val


subroutine bf_flow_vec(fe,coord,t,vec)
   implicit none
   double precision                :: coord(:,:), fe(:)
   double precision                :: vec(:), t

   if (t.gt.0d0) then
      call fl2d4_bf(fe,coord,t,vec) 
   else
      call fl3d8_bf(fe,coord,vec) 
   endif

end subroutine bf_flow_vec

!------------------------------------------------------------------------------
!                             Body forces
!------------------------------------------------------------------------------

subroutine bf(fe,coord,t,b)
   implicit none
   double precision                :: coord(:,:), fe(:), b(:,:)
   double precision                :: t

   if (t.gt.0d0) then
      call c2dtl4_b(fe,coord,t,b) 
   else 
      call c3dtl8_b(fe,coord,b) 
   endif

end subroutine bf


subroutine bf_sens(fe,coord,t,b,arg)
   implicit none
   double precision                :: coord(:,:), fe(:), b(:,:), arg(:)
   double precision                :: t

   if (t.gt.0d0) then
      call c2dtl4_bsens(fe,coord,t,b,arg) 
   else
      call c3dtl8_bsens(fe,coord,b,arg) 
   endif

end subroutine bf_sens


!------------------------------------------------------------------------------
!                             Element volumes 
!------------------------------------------------------------------------------


subroutine vol(evol,coord,t)
   implicit none
   double precision                :: evol, coord(:,:), t

   if (t.gt.0d0) then
      call fl2d4_vi(evol,coord,t) 
   else
      call fl3d8_vi(evol,coord) 
   endif

end subroutine vol 


subroutine vol_vec(evol,coord,t,vec)
   implicit none
   double precision                :: evol, coord(:,:), t, vec(:)

   if (t.gt.0d0) then
      call fl2d4_vi(evol,coord,t,vec) 
   else
      call fl3d8_vi(evol,coord,vec) 
   endif

end subroutine vol_vec 

!------------------------------------------------------------------------------
!                     Largrangian strain type quantity
!------------------------------------------------------------------------------

subroutine Etilde1(coord,ed,edmutilde,Emutilde,dgmutilde,gammagp)
   implicit none
   double precision                :: coord(:,:), dgmutilde(:,:)
   double precision                :: ed(:), edmutilde(:), Emutilde(:,:), gammagp(:)

   if (size(coord,1).eq.3) then
      call c3dtl8_emu(coord,ed,edmutilde,Emutilde,dgmutilde,gammagp)
   elseif (size(coord,1).eq.2) then
      call c2dtl4_emu(coord,ed,edmutilde,Emutilde,dgmutilde,gammagp)
   endif

end subroutine Etilde1 


subroutine Etilde2(coord,ed,edmutilde,Emutilde,dgmutilde,gammagp,ones)
   implicit none
   double precision                :: coord(:,:), dgmutilde(:,:), ones(:)
   double precision                :: ed(:), edmutilde(:), Emutilde(:,:), gammagp(:)

   if (size(coord,1).eq.3) then
      call c3dtl8_emu2(coord,ed,edmutilde,Emutilde,dgmutilde,gammagp,ones)
   elseif (size(coord,1).eq.2) then
      call c2dtl4_emu2(coord,ed,edmutilde,Emutilde,dgmutilde,gammagp,ones)
   endif

end subroutine Etilde2



end module FE_element





