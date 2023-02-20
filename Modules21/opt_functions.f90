
!############################################################
! The module contains various small functions used in the
! optimization.
!############################################################
! Created: 2020-12-10
! Anna Dalklint
! Lund University
!############################################################
module opt_functions

!$ use omp_lib

contains

!=======================================================================

!                 Choose aggregation function

!  arg(:): Vector containing the variables subject to the norm
!  num_arg: Number of variables 

!=======================================================================
subroutine aggfun(aggtype,val,arg,num_arg,pAgg)
   implicit none
   double precision                 :: val, arg(:), pAgg
   integer                          :: num_arg
   character(len=*)                 :: aggtype
   
   if (aggtype.eq.'Pnorm') then
      call pfun(val,arg,num_arg,pAgg)
   else 
      stop 'aggfun: This type of aggregation function has not been implemented!'
   endif

   return
end subroutine aggfun



!=======================================================================

!                      Calculates the p-norm

!=======================================================================
subroutine pfun(val,arg,num_arg,pAgg)
   implicit none
   double precision                :: val, arg(:), pAgg
   integer                         :: i, num_arg

   val = 0d0
  
   do i=1,num_arg
      ! Should I use absolute value?!?!
      val = val + arg(i)**pAgg
   enddo
  
   val = val**(1d0/pAgg)

   return
end subroutine pfun

!=======================================================================

!                 Derivative of aggregation function

!=======================================================================
subroutine Aggprim(aggtype,dAggCurr,arg,dfdarg,dargdrho,num_arg,pAgg)
   implicit none
   double precision                 :: dAggCurr(:), arg(:), dfdarg(:), dargdrho(:,:), pAgg
   integer                          :: num_arg
   character(len=*)                 :: aggtype
   
   if (aggtype.eq.'Pnorm') then
       call pAggprim(dAggCurr,arg,dfdarg,dargdrho,num_arg,pAgg)
   else 
     stop 'Aggprim: This type of aggregation function has not been implemented!'
   endif

   return
end subroutine Aggprim


!=======================================================================

!                 Derivative of p-norm function

!  dfdarg: If chain rule must be used to obtain dargdrho (if not applicable 
!          set to 1d0)

!=======================================================================
subroutine pAggprim(dAggCurr,arg,dfdarg,dargdrho,num_arg,pAgg)
   implicit none
   double precision     :: val, dAggCurr(:), arg(:), dfdarg(:), dargdrho(:,:), pAgg
   integer              :: i, num_arg

   dAggCurr = 0d0
   val = 0d0
   do i=1,num_arg
     dAggCurr = dAggCurr + arg(i)**(pAgg-1d0)*dfdarg(i)*dargdrho(i,:)
   enddo
  
   call pfun(val,arg,num_arg,pAgg)
   dAggCurr = val**(1d0-pAgg)*dAggCurr
  
   return
end subroutine pAggprim


!=======================================================================

!                 Calculates the gamma(=Heaviside) function

!=======================================================================
subroutine gammafun(val, rhogpCurr, pcurr, Hproj, B1, rho0)
   implicit none
   double precision         :: val2, rhogpCurr, val, pcurr, B1, rho0
   integer                   :: Hproj

   val = 0d0
   val2 = 0d0

   if (Hproj.eq.1) then
      call Hfun(val2,rhogpCurr)
   else
      val2 = rhogpCurr
   endif

   val = 0D0
   if (Gscaling.eq.1) then
      val = (dtanh(B1*rho0)+dtanh(B1*(val2**pcurr-rho0)))/(dtanh(B1*rho0)+dtanh(B1*(0.1D1-rho0)))   
   else
      val=1D0
   endif

   return
end subroutine gammafun



!=======================================================================

!      Calculates the derivative of the gamma(=Heaviside)-function

!=======================================================================
subroutine gammaprimfun(val,rhogpCurr, pcurr,Hproj,B1,rho0)
   implicit none
   double precision                 :: val2, val3, rhogpCurr, val, pcurr, B1, rho0
   integer                          :: Hproj

   val = 0d0
   val2 = 0d0
   val3 = 0d0

   if (Hproj.eq.1) then
      call Hfun(val2,rhogpCurr)
   else
      val2 = rhogpCurr
   endif

   val=0D0
   if (Gscaling.eq.1) then
      val = (0.1D1-dtanh(B1*(val2**pcurr-rho0))**2D0) *  & 
                       B1*pcurr*val2**(pcurr-1d0)/(dtanh(B1*rho0) + dtanh(B1 * (0.1D1 - rho0)))

      if (Hproj.eq.1) then
         call Hprimfun(val3, rhogpCurr)
         val = val*val3
      endif
   else
      val=0D0
   endif
   
   return
end subroutine gammaprimfun


!=======================================================================

!                 Choose penalization function

!=======================================================================
subroutine gfun(pentype,val,rhogpcurr,pcurr)
   implicit none
   double precision                 :: val, rhogpCurr, pCurr
   character(len=*)                 :: pentype
   
   if (pentype.eq.'SIMP') then
     call simpfun(val,rhogpcurr,pcurr)
  elseif (pentype.eq.'RAMP') then
     call rampfun(val,rhogpcurr,pcurr)
   else
     stop 'This type of aggregation function has not been implemented!'
   endif

   return
end subroutine gfun



!=======================================================================

!                    Calculates the Xi-function
!                            (SIMP)

!=======================================================================
subroutine simpfun(val,rhogpCurr, pcurr, Hproj, delta0)
   implicit none
   double precision                 :: rhogpCurr, val2, val, pcurr, delta0
   integer                          :: Hproj
   
   val2 = 0d0
   val = 0d0

   if (Hproj.eq.1) then
      call Hfun(val2,rhogpCurr)
   else
      val2 = rhogpCurr
   endif

   ! If Helmholtz filter is used xval can be negative
   if (val2.ge.0D0) then
      val = (1D0-delta0)*val2**pcurr+delta0
   else
      val = delta0
   endif

   return
end subroutine simpfun


!=======================================================================

!                    Calculates the Xi-function
!                            (RAMP)

!=======================================================================
subroutine rampfun(val,rhogpCurr, pcurr, Hproj, delta0)
   implicit none
   double precision                 :: rhogpCurr, val2, val, pcurr, delta0
   integer                          :: Hproj
   
   val2 = 0d0
   val = 0d0

   if (Hproj.eq.1) then
      call Hfun(val2,rhogpCurr)
   else
      val2 = rhogpCurr
   endif

   ! If Helmholtz filter is used xval can be negative
   if (val2.ge.0D0) then
      val = delta0 + val2*(1d0-delta0)/(1d0+pcurr*(1d0-val2))
   else
      val = delta0
   endif

   return
end subroutine rampfun


 


!=======================================================================

!                 Derivative of Xi-function

!=======================================================================
subroutine gfunprim(pentype,val,rhogpCurr, pcurr)
   implicit none
   double precision                 :: val, rhogpCurr, pcurr
   character(len=*)                 :: pentype
   
   if (pentype.eq.'SIMP') then
     call simpfunprim(val,rhogpcurr,pcurr)
   elseif (pentype.eq.'RAMP') then
     call rampfunprim(val,rhogpcurr,pcurr)
   else
     stop 'This type of aggregation function has not been implemented!'
   endif

   return
end subroutine gfunprim




!=======================================================================

!           Calculates the derivative of Xi-function
!                            (SIMP)

!=======================================================================
subroutine simpfunprim(val,rhogpCurr,pcurr,Hproj,delta0, B, eta)
   implicit none
   double precision          :: rhogpCurr, val2, val3, val, pcurr, delta0, B, eta
   integer                   :: Hproj

   val = 0d0
   val2 = 0d0
   val3 = 0d0

   if (Hproj.eq.1) then
      call Hfun(val2, rhogpCurr, B, eta)
      call Hprimfun(val3, rhogpCurr, B, eta)
      val=pcurr*(1D0-delta0)*val2**(pcurr-1D0)*val3
   else
      val=pcurr*(1D0-delta0)*rhogpCurr**(pcurr-1D0)
   endif

   if (rhogpCurr.lt.0D0)  val=0D0

   return
end subroutine simpfunprim




!=======================================================================

!           Calculates the derivative of Xi-function
!                            (RAMP)

!=======================================================================
subroutine rampfunprim(val,rhogpCurr,pcurr,Hproj,delta0,B,eta)
   implicit none
   double precision          :: rhogpCurr, val2, val3, val, pcurr, delta0, B, eta
   integer                   :: Hproj

   val = 0d0
   val2 = 0d0
   val3 = 0d0

   if (Hproj.eq.1) then
      call Hfun(val2, rhogpCurr, B, eta)
      call Hprimfun(val3, rhogpCurr, B, eta)
      val=((1D0-delta0)*(1d0+pcurr)/((1d0+pcurr*(1d0-val2))**2D0))*val3
   else
      val=(1D0-delta0)*(1d0+pcurr)/((1d0+pcurr*(1d0-rhogpCurr))**2D0)
   endif

   if (rhogpCurr.lt.0D0)  val=0D0

   return
end subroutine rampfunprim






!=======================================================================

!                 Calculates the Heaviside function

!=======================================================================
subroutine Hfun(val, rhogpCurr, B, eta)
   implicit none
   double precision                 ::  rhogpCurr, val, B, eta
      
   val = 0d0  
   val = (dtanh(B*eta)+dtanh(B*(rhogpCurr-eta)))/(dtanh(B*eta)+dtanh(B*(0.1D1-eta)))  
   
   return
end subroutine Hfun


    
!=======================================================================

!        Calculates the derivative of the Heaviside function

!=======================================================================   
subroutine Hprimfun(val, rhogpCurr, B, eta)
   implicit none
   double precision                 ::  rhogpCurr, val, B, eta
    
   val = 0d0
   val = B*(0.1D1-dtanh(B*(rhogpCurr-eta))**2D0)/(dtanh(B*eta)+dtanh(B*(0.1D1-eta)))
   
   return
end subroutine Hprimfun




!=======================================================================

!          Calculates the deformation gradient with gamma

!=======================================================================
subroutine calcdggamma(dggamma,dge,gammagpe)
   implicit none
   double precision                    :: dggamma(:,:), dge(:,:), gammagpe(:)
   integer                             :: igp
   double precision                    :: unitVec(eldef)
   
   UnitVec   =0D0
   if (FE_dim.eq.'3D') then
      UnitVec(1)=1d0
      UnitVec(5)=1d0
      UnitVec(9)=1d0 
   elseif (FE_dim.eq.'2D') then
      UnitVec(1)=1d0
      UnitVec(4)=1d0
   endif
   
   do igp=1,ngp   
     dggamma(:,igp)=UnitVec+(dge(:,igp)-UnitVec)*gammagpe(igp)
   enddo

end subroutine calcdggamma


end module opt_functions
