module mater_large

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
!  - removed hyper elastic material derived from mater_large
! M. Wallin 2011-05-16
!  - General 3D pushforward implemented
! M. Wallin 2011-05-16
!  - General 2D/3D pullback implemented
! M. Ristinmaa 2011-05-23
!  - changed order of arguments in upddefgr
! ------------------------------------------------------------------------------
! M. Wallin 2011-09-16
!  - padelog, padeexp, dpadelog, dpadeexp
! M. Wallin 2011-09-16
!  - BAR_DYAD implemented
! M. Wallin 2011-09-16
!  - inv3, inv3s, det3 moved to matrix_util
! M. Wallin 2011-09-16
!  - getF moved to mater_large


! Define pullback and pushforward

use matrix_util

implicit none

interface reldefgr
  module procedure reldefgr1
  module procedure reldefgr2
  module procedure reldefgr3
end interface

private reldefgr1, reldefgr2, reldefgr3

interface pushforward
  module procedure pushforward1
  module procedure pushforward2
  module procedure pushforward3
  module procedure pushforward4
  module procedure pushforward5
end interface

private pushforward1, pushforward2, pushforward3, pushforward4, pushforward5

interface pullback
  module procedure pullback1
  module procedure pullback2
  module procedure pullback3
  module procedure pullback4
  module procedure pullback5
end interface

private pullback1, pullback2, pullback3, pullback4, pullback5




private  det3, getI

!-------------------------------------------------------------------------------

contains

subroutine upddefgr(dgn,ddg,dg)
  implicit none
  double precision                :: dgn(:), dg(:), ddg(:)

  if (size(dg).eq.4) then
      dgn(1) = dg(1) * ddg(1) + ddg(2) * dg(3)
      dgn(2) = ddg(1) * dg(2) + dg(4) * ddg(2)
      dgn(3) = dg(1) * ddg(3) + ddg(4) * dg(3)
      dgn(4) = ddg(3) * dg(2) + dg(4) * ddg(4)
  else         
      dgn(1) = dg(1) * ddg(1) + dg(4) * ddg(2) + ddg(3) * dg(7)
      dgn(2) = ddg(1) * dg(2) + ddg(2) * dg(5) + ddg(3) * dg(8)
      dgn(3) = dg(3) * ddg(1) + ddg(2) * dg(6) + ddg(3) * dg(9)
      dgn(4) = ddg(4) * dg(1) + ddg(5) * dg(4) + ddg(6) * dg(7)
      dgn(5) = dg(2) * ddg(4) + ddg(5) * dg(5) + ddg(6) * dg(8)
      dgn(6) = ddg(4) * dg(3) + ddg(5) * dg(6) + ddg(6) * dg(9)
      dgn(7) = ddg(7) * dg(1) + ddg(8) * dg(4) + ddg(9) * dg(7)
      dgn(8) = ddg(7) * dg(2) + ddg(8) * dg(5) + ddg(9) * dg(8)
      dgn(9) = ddg(7) * dg(3) + ddg(8) * dg(6) + ddg(9) * dg(9)

  endif

  return
end subroutine upddefgr

subroutine reldefgr1(ddg,dg,dgn)
  implicit none
  double precision                :: dgn(:,:,:), dg(:,:,:), ddg(:,:,:)

  integer                         :: nelm, ie

  nelm=size(dgn,3)
  do ie=1,nelm
    call reldefgr(ddg(:,:,ie),dg(:,:,ie),dgn(:,:,ie))
  end do

end subroutine reldefgr1

subroutine reldefgr2(ddg,dg,dgn)
  implicit none
  double precision                :: dgn(:,:), dg(:,:), ddg(:,:)

  integer                         :: ngp, igp

  ngp=size(dgn,2)
  do igp=1,ngp
    call reldefgr(ddg(:,igp),dg(:,igp),dgn(:,igp))
  end do

end subroutine reldefgr2


subroutine reldefgr3(ddg,dg,dgn)
  implicit none
  double precision                :: dgn(:), dg(:), ddg(:)

  double precision   :: t1, t2, t3, t5, t6, t8, t10, t11, t13, t15, t16, t17
  double precision   :: t19, t20, t22, t23, t25, t28, t30, t33, t36, t39, t45
  double precision   :: t49, t53, t50, t59, t62, t65, t69, t72, t75, t93, t96
  double precision   :: t99

  if (size(dg).eq.4) then
     stop "not implemented"
  else         
      t1 = dg(1)
      t2 = dgn(5)
      t3 = dgn(9)
      t5 = dgn(6)
      t6 = dgn(8)
      t8 = t2*t3-t5*t6
      t10 = dgn(1)
      t11 = t10*t2
      t13 = t10*t5
      t15 = dgn(4)
      t16 = dgn(2)
      t17 = t15*t16
      t19 = dgn(3)
      t20 = t15*t19
      t22 = dgn(7)
      t23 = t22*t16
      t25 = t22*t19
      t28 = 1/(t11*t3-t13*t6-t17*t3+t20*t6+t23*t5-t25*t2)
      t30 = dg(2)
      t33 = t15*t3-t5*t22
      t36 = dg(3)
      t39 = -t15*t6+t2*t22
      t45 = t16*t3-t19*t6
      t49 = t10*t3-t25
      t53 = t10*t6-t23
      t59 = -t16*t5+t19*t2
      t62 = t13-t20
      t65 = t11-t17
      t69 = dg(4)
      t72 = dg(5)
      t75 = dg(6)
      t93 = dg(7)
      t96 = dg(8)
      t99 = dg(9)
      ddg(1) = t1*t8*t28-t30*t33*t28-t36*t39*t28
      ddg(2) = -t1*t45*t28+t30*t49*t28-t36*t53*t28
      ddg(3) = -t1*t59*t28-t30*t62*t28+t36*t65*t28
      ddg(4) = t69*t8*t28-t72*t33*t28-t75*t39*t28
      ddg(5) = -t69*t45*t28+t72*t49*t28-t75*t53*t28
      ddg(6) = -t69*t59*t28-t72*t62*t28+t75*t65*t28
      ddg(7) = t93*t8*t28-t96*t33*t28-t99*t39*t28
      ddg(8) = -t93*t45*t28+t96*t49*t28-t99*t53*t28
      ddg(9) = -t93*t59*t28-t96*t62*t28+t99*t65*t28

  endif

  return
end subroutine reldefgr3



subroutine pushforward5(res,var,ef,jtype)
  implicit none
  double precision               :: res(:,:,:), var(:,:,:), ef(:,:)
  integer                        :: ii, order
  character(len=*), optional     :: jtype
 
  do ii=1,size(var,3)
   if (present(jtype)) then
     call pushforward4(res(:,:,ii),var(:,:,ii),ef(:,ii),jtype)
   else
     call pushforward4(res(:,:,ii),var(:,:,ii),ef(:,ii))
   endif
  enddo
    
  return
end subroutine pushforward5




subroutine pushforward4(res,var,ef,jtype)
  implicit none
  double precision               :: res(:,:), var(:,:), ef(:)
  integer                        :: ii
  double precision               :: F(3,3), M3(3,3), j, M6(6,6)
  character(len=*), optional     :: jtype
 
!  M_pqij C_ijkl M_stkl
   if ((size(var,1).eq.3).and.(size(var,2).eq.3)) then
    F(1,:)=(/ef(1), ef(2), 0d0/)
    F(2,:)=(/ef(3), ef(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)
      
    M3(1,:)=(/F(1,1)*F(1,1), F(1,2)*F(1,2), 2d0*F(1,1)*F(1,2) /)
    M3(2,:)=(/F(2,1)*F(2,1), F(2,2)*F(2,2), 2d0*F(2,1)*F(2,2) /)
    M3(3,:)=(/F(1,1)*F(2,1), F(1,2)*F(2,2), F(1,1)*F(2,2)+F(1,2)*F(2,1) /)
    res=matmul(M3,matmul(var,transpose(M3)))
    if (present(jtype)) then 
       j=det3(F)
       res=res/j
    endif

   elseif ((size(var,1).eq.6).and.(size(var,2).eq.6)) then
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)

    M6(1,:)=(/F(1,1)*F(1,1), F(1,2)*F(1,2), F(1,3)*F(1,3), 2D0*F(1,1)*F(1,2), 2D0*F(1,1)*F(1,3), 2D0*F(1,2)*F(1,3) /)
    M6(2,:)=(/F(2,1)*F(2,1), F(2,2)*F(2,2), F(2,3)*F(2,3), 2D0*F(2,1)*F(2,2), 2D0*F(2,1)*F(2,3), 2D0*F(2,3)*F(2,2) /)
    M6(3,:)=(/F(3,1)*F(3,1), F(3,2)*F(3,2), F(3,3)*F(3,3), 2D0*F(3,2)*F(3,1), 2D0*F(3,3)*F(3,1), 2D0*F(3,3)*F(3,2) /)          
    M6(4,:)=(/F(1,1)*F(2,1), F(1,2)*F(2,2), F(1,3)*F(2,3), F(1,1)*F(2,2)+ F(1,2)*F(2,1), F(1,1)*F(2,3)+F(1,3)*F(2,1), F(1,2)*F(2,3)+F(1,3)*F(2,2) /)
    M6(5,:)=(/F(1,1)*F(3,1), F(1,2)*F(3,2), F(1,3)*F(3,3), F(1,1)*F(3,2)+ F(1,2)*F(3,1), F(1,1)*F(3,3)+F(1,3)*F(3,1), F(1,2)*F(3,3)+F(1,3)*F(3,2) /)
    M6(6,:)=(/F(2,1)*F(3,1), F(2,2)*F(3,2), F(2,3)*F(3,3), F(2,1)*F(3,2)+ F(2,2)*F(3,1), F(2,1)*F(3,3)+F(2,3)*F(3,1), F(2,2)*F(3,3)+F(2,3)*F(3,2) /)

    res=matmul(M6,matmul(var,transpose(M6)))
    if (present(jtype)) then 
       j=det3(F)
       res=res/j
    endif
   else
     stop 'pushforward dimensions does not fit'
   endif

  return
end subroutine pushforward4



subroutine pushforward3(res,var,ef,jtype)
  implicit none
  double precision               :: res(:,:,:), var(:,:,:), ef(:,:,:)
  integer                        :: ii, order
  character(len=*), optional     :: jtype
 
  do ii=1,size(var,3)
   if (present(jtype)) then
     call pushforward2(res(:,:,ii),var(:,:,ii),ef(:,:,ii),jtype)
   else
     call pushforward2(res(:,:,ii),var(:,:,ii),ef(:,:,ii))
   endif
  enddo
    
  return
end subroutine pushforward3


subroutine pushforward2(res,var,ef,jtype)
  implicit none
  double precision               :: res(:,:), var(:,:), ef(:,:)
  integer                        :: ii
  double precision               :: F(3,3), M3(3,3)
  character(len=*), optional     :: jtype
 
    do ii=1,size(var,2)
      if (present(jtype)) then
       call pushforward1(res(:,ii),var(:,ii),ef(:,ii),jtype)
     else
       call pushforward1(res(:,ii),var(:,ii),ef(:,ii))
     endif
    enddo

  return
end subroutine pushforward2

subroutine pushforward1(res,var,ef,jtype)
  implicit none
  double precision               :: res(:), var(:), ef(:)
  double precision               :: F(3,3), bv(3,3), br(3,3), j
  character(len=*), optional     :: jtype

    
  if (size(ef).eq.4) then
    F(1,:)=(/ef(1), ef(2), 0d0/)
    F(2,:)=(/ef(3), ef(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)
    if (present(jtype)) then 
       j=det3(F)
    else
       j=1d0
    endif
    bv(1,:)=(/var(1), var(4),    0d0/)
    bv(2,:)=(/var(4), var(2),    0d0/)
    bv(3,:)=(/  0d0,     0d0, var(3)/)
    br=matmul(F,matmul(bv,transpose(F)))/j
    res(1)=br(1,1)
    res(2)=br(2,2)
    res(3)=br(3,3)
    res(4)=br(1,2)
  else
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
    if (present(jtype)) then 
       j=det3(F)
    else
       j=1d0
    endif
    bv(1,:)=(/var(1), var(4), var(5)/)
    bv(2,:)=(/var(4), var(2), var(6)/)
    bv(3,:)=(/var(5), var(6), var(3)/)
    br=matmul(F,matmul(bv,transpose(F)))/j
    res(1)=br(1,1)
    res(2)=br(2,2)
    res(3)=br(3,3)
    res(4)=br(1,2)
    res(5)=br(1,3)
    res(6)=br(2,3)
  endif

  return
end subroutine pushforward1


subroutine pullback5(res,var,ef,jtype)
  implicit none
  double precision               :: res(:,:,:), var(:,:,:), ef(:,:)
  integer                        :: ii
  character(len=*), optional     :: jtype
 
  do ii=1,size(var,3)
   if (present(jtype)) then
     call pullback4(res(:,:,ii),var(:,:,ii),ef(:,ii),jtype)
   else
     call pullback4(res(:,:,ii),var(:,:,ii),ef(:,ii))
   endif
  enddo
    
  return
end subroutine pullback5

subroutine pullback4(res,var,ef,jtype)
  implicit none
  double precision               :: res(:,:), var(:,:), ef(:)
  integer                        :: ii
  character(len=*), optional     :: jtype
  double precision               :: F(3,3), M3(3,3), j, M6(6,6), Fi(3,3)

!  M_pqij C_ijkl M_stkl
   if ((size(var,1).eq.3).and.(size(var,2).eq.3)) then
    F(1,:)=(/ef(1), ef(2), 0d0/)
    F(2,:)=(/ef(3), ef(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)

    call inv3(Fi,F)          

    M3(1,:)=(/Fi(1,1)*Fi(1,1), Fi(1,2)*Fi(1,2), 2d0*Fi(1,1)*Fi(1,2) /)
    M3(2,:)=(/Fi(2,1)*Fi(2,1), Fi(2,2)*Fi(2,2), 2d0*Fi(2,1)*Fi(2,2) /)
    M3(3,:)=(/Fi(1,1)*Fi(2,1), Fi(1,2)*Fi(2,2), Fi(1,1)*Fi(2,2)+Fi(1,2)*Fi(2,1) /)
    res=matmul(M3,matmul(var,transpose(M3)))
    if (present(jtype)) then 
       j=det3(F)
       res=res*j
    endif

   elseif ((size(var,1).eq.6).and.(size(var,2).eq.6)) then
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
    
    call inv3(Fi,F)        
    
    M6(1,:)=(/Fi(1,1)*Fi(1,1), Fi(1,2)*Fi(1,2), Fi(1,3)*Fi(1,3), 2D0*Fi(1,1)*Fi(1,2), 2D0*Fi(1,1)*Fi(1,3), 2D0*Fi(1,2)*Fi(1,3) /)
    M6(2,:)=(/Fi(2,1)*Fi(2,1), Fi(2,2)*Fi(2,2), Fi(2,3)*Fi(2,3), 2D0*Fi(2,1)*Fi(2,2), 2D0*Fi(2,1)*Fi(2,3), 2D0*Fi(2,3)*Fi(2,2) /)
    M6(3,:)=(/Fi(3,1)*Fi(3,1), Fi(3,2)*Fi(3,2), Fi(3,3)*Fi(3,3), 2D0*Fi(3,2)*Fi(3,1), 2D0*Fi(3,3)*Fi(3,1), 2D0*Fi(3,3)*Fi(3,2) /)          
    M6(4,:)=(/Fi(1,1)*Fi(2,1), Fi(1,2)*Fi(2,2), Fi(1,3)*Fi(2,3), Fi(1,1)*Fi(2,2)+ Fi(1,2)*Fi(2,1), Fi(1,1)*Fi(2,3)+Fi(1,3)*Fi(2,1), Fi(1,2)*Fi(2,3)+Fi(1,3)*Fi(2,2) /)
    M6(5,:)=(/Fi(1,1)*Fi(3,1), Fi(1,2)*Fi(3,2), Fi(1,3)*Fi(3,3), Fi(1,1)*Fi(3,2)+ Fi(1,2)*Fi(3,1), Fi(1,1)*Fi(3,3)+Fi(1,3)*Fi(3,1), Fi(1,2)*Fi(3,3)+Fi(1,3)*Fi(3,2) /)
    M6(6,:)=(/Fi(2,1)*Fi(3,1), Fi(2,2)*Fi(3,2), Fi(2,3)*Fi(3,3), Fi(2,1)*Fi(3,2)+ Fi(2,2)*Fi(3,1), Fi(2,1)*Fi(3,3)+Fi(2,3)*Fi(3,1), Fi(2,2)*Fi(3,3)+Fi(2,3)*Fi(3,2) /)

    res=matmul(M6,matmul(var,transpose(M6)))
    if (present(jtype)) then 
       j=det3(F)
       res=res*j
    endif

   else
     stop 'pushforward dimensions does not fit'
   endif
    
  return
end subroutine pullback4


subroutine pullback3(res,var,ef,jtype)
  implicit none
  double precision               :: res(:,:,:), var(:,:,:), ef(:,:,:)
  integer                        :: ii
  character(len=*), optional     :: jtype
 
  do ii=1,size(var,3)
   if (present(jtype)) then
     call pullback2(res(:,:,ii),var(:,:,ii),ef(:,:,ii),jtype)
   else
     call pullback2(res(:,:,ii),var(:,:,ii),ef(:,:,ii))
   endif
  enddo
    
  return
end subroutine pullback3


subroutine pullback2(res,var,ef,jtype)
  implicit none
  double precision               :: res(:,:), var(:,:), ef(:,:)
  integer                        :: ii
  character(len=*), optional     :: jtype
 
    do ii=1,size(var,2)
     if (present(jtype)) then
       call pullback1(res(:,ii),var(:,ii),ef(:,ii),jtype)
     else
       call pullback1(res(:,ii),var(:,ii),ef(:,ii))
     endif
    enddo
    
  return
end subroutine pullback2



subroutine pullback1(res,var,ef,jtype)
  implicit none
  double precision               :: res(:), var(:), ef(:)
  double precision               :: F(3,3), bv(3,3), br(3,3), Finv(3,3), j
  character(len=*), optional     :: jtype


  if (size(ef).eq.4) then
    F(1,:)=(/ef(1), ef(2), 0d0/)
    F(2,:)=(/ef(3), ef(4), 0d0/)
    F(3,:)=(/  0d0,   0d0, 1d0/)
    call inv3(Finv,F)
    if (present(jtype)) then 
       j=det3(F)
    else
       j=1d0
    endif
    bv(1,:)=(/var(1), var(4),    0d0/)
    bv(2,:)=(/var(4), var(2),    0d0/)
    bv(3,:)=(/  0d0,     0d0, var(3)/)
    br=matmul(Finv,matmul(bv,transpose(Finv)))*j
    res(1)=br(1,1)
    res(2)=br(2,2)
    res(3)=br(3,3)
    res(4)=br(1,2)
  elseif (size(ef).eq.9) then
    F(1,:)=(/ef(1), ef(2), ef(3)/)
    F(2,:)=(/ef(4), ef(5), ef(6)/)
    F(3,:)=(/ef(7), ef(8), ef(9)/)
    if (present(jtype)) then 
       j=det3(F)
    else
       j=1d0
    endif
    call inv3(Finv,F)
    bv(1,:)=(/var(1), var(4), var(5)/)
    bv(2,:)=(/var(4), var(2), var(6)/)
    bv(3,:)=(/var(5), var(6), var(3)/)
    br=matmul(Finv,matmul(bv,transpose(Finv)))*j
    res(1)=br(1,1)
    res(2)=br(2,2)
    res(3)=br(3,3)
    res(4)=br(1,2)
    res(5)=br(1,3)
    res(6)=br(2,3)
  else
    stop 'size of deformation gradient wrong'
  endif

  return
end subroutine pullback1




 


SUBROUTINE PADELOG(M, LOGM)

DOUBLE PRECISION, DIMENSION(3,3) :: M, LOGM, A, T, INVT, Y, INVY, EYE
REAL		                 :: KM, KN
DOUBLE PRECISION,PARAMETER       :: ZERO=0., ONE=1., THREE=3., &
				    ONETHIRD=ONE/THREE 
				
   DO KM=1,3
      DO KN=1,3
         IF(KN.EQ.KM) THEN
            EYE(KM,KN)=ONE
         ELSE
   	    EYE(KM,KN)=ZERO
         END IF
      END DO
   END DO				

   T=M+EYE
!   CALL SYMINV(T,INVT)
   CALL inv3s(INVT,T)
   T=M-EYE
   A=MATMUL(T,INVT)
   Y=EYE-ONETHIRD*MATMUL(A,A)
!   CALL SYMINV(Y,INVY)
   CALL INV3s(INVY,Y)
   LOGM=2D0*MATMUL(A,INVY)

   END SUBROUTINE 

SUBROUTINE DPADELOG(C, S)

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3,3) :: C, EYE, T, A, Y, E, B, D,&
                                    F, M, N, INVT, INVY
DOUBLE PRECISION, DIMENSION(6,2) :: IJ				    
DOUBLE PRECISION, DIMENSION(6,6) :: S
DOUBLE PRECISION                 :: S1, S2, S3, S4, S5
REAL                             :: KM,KN, I, J, K, L
DOUBLE PRECISION, PARAMETER      :: ZERO=0., ONE=1., TWO=2., THREE=3.,&
 				    FOUR=4., TWELVE=12.,&
                                    HALF=ONE/TWO, THIRD=ONE/THREE, &
				    FOURTH=ONE/FOUR, TWELFTH=ONE/TWELVE

   DO KM=1,3
      DO KN=1,3
         IF(KN.EQ.KM) THEN
            EYE(KM,KN)=ONE
         ELSE
            EYE(KM,KN)=ZERO
         END IF
     END DO
   END DO

   IJ(1,1)=1 
   IJ(1,2)=1 
   IJ(2,1)=2 
   IJ(2,2)=2 
   IJ(3,1)=3 
   IJ(3,2)=3
   IJ(4,1)=1 
   IJ(4,2)=2 
   IJ(5,1)=1 
   IJ(5,2)=3 
   IJ(6,1)=2 
   IJ(6,2)=3

   T=C+EYE
!   CALL SYMINV(T,INVT)
   CALL INV3S(INVT,T)
   T=C-EYE
   A=MATMUL(T,INVT)
   Y=EYE-THIRD*MATMUL(A,A)
!   CALL SYMINV(Y,INVY)
   CALL INV3S(INVY,Y)
   T=INVT
   Y=INVY
   E=MATMUL(A,Y)
   B=MATMUL(E,A)
   D=MATMUL(T,Y)
   F=MATMUL(T,E)
   M=MATMUL(T,B)
   N=MATMUL(B,A)

   DO KM=1,6
       I=IJ(KM,1)
       J=IJ(KM,2)
       DO KN=KM,6
          K=IJ(KN,1)
          L=IJ(KN,2)
          S1=(FOURTH*(EYE(I,K)-A(I,K))+TWELFTH*(B(I,K)-N(I,K)))*D(L,J)+&
            (FOURTH*(EYE(I,L)-A(I,L))+TWELFTH*(B(I,L)-N(I,L)))*D(K,J)
          S2=(-FOURTH*T(I,L)+TWELFTH*(F(I,L)-M(I,L)))*E(K,J)+&
            (-FOURTH*T(I,K)+TWELFTH*(F(I,K)-M(I,K)))*E(L,J)
          S3=TWELFTH*(E(I,K)-B(I,K))*F(L,J)+TWELFTH*(E(I,L)-B(I,L))*F(K,J)
          S4=FOURTH*(Y(K,J)*T(I,L)+Y(L,J)*T(I,K))
          S5=TWELFTH*((M(I,L)*Y(J,K)+M(I,K)*Y(J,L))-(F(I,L)*B(J,K)+F(I,K)*B(J,L)))
          S(KM,KN)=S1+S2+S3+S4+S5
       END DO
   END DO

   DO KM=1,6 
     DO KN=1,KM-1
        S(KM,KN)=S(KN,KM)
     END DO
   END DO

   S=2D0*S   
   S(:,4:6)=2D0*S(:,4:6)


END SUBROUTINE DPADELOG



SUBROUTINE PADEEXP(IN, EXPM)

DOUBLE PRECISION, DIMENSION(3,3) :: IN, V, VINV, H, KR, EXPM

  KR     =getI()
  V      =KR-0.5D0*IN
!  call  SYMINV(V,VINV)
  call  INV3(VINV,V)
  H      =KR+0.5D0*IN
  EXPM   =MATMUL(VINV,H)
 
END SUBROUTINE 


SUBROUTINE DPADEEXP(IN, DEXPM)

DOUBLE PRECISION, DIMENSION(3,3) :: IN, EXPM, M, P, IM, TMP33, KR
DOUBLE PRECISION, DIMENSION(6,6) :: DEXPM, TMP66_1, TMP66_2

  KR     =getI()
  M      =KR-0.5D0*IN
  P      =KR+0.5D0*IN
!  call   SYMINV(M,IM)
  call   INV3(IM,M)
  TMP33  =MATMUL(P,IM)+KR
  
  CALL BAR_DYAD(TMP66_1,IM,TMP33)
  CALL BAR_DYAD(TMP66_2,TMP33,IM)

  DEXPM=0.5D0*(TMP66_1+TMP66_2)

  DEXPM=0.5D0*DEXPM 

END SUBROUTINE 


! extract deformation gradient
function getF(ef)
  implicit none
  double precision                :: getF(3,3), ef(:)
 
  if (size(ef).eq.4) then
    getF(1,:)=(/ef(1), ef(2), 0d0/)
    getF(2,:)=(/ef(3), ef(4), 0d0/)
    getF(3,:)=(/  0d0,   0d0, 1d0/)
  else    
    getF(1,:)=(/ef(1), ef(2), ef(3)/)
    getF(2,:)=(/ef(4), ef(5), ef(6)/)
    getF(3,:)=(/ef(7), ef(8), ef(9)/)
  endif
  return
end function getF



SUBROUTINE BAR_DYAD(OUT,A,B)

DOUBLE PRECISION, DIMENSION(3,3) :: A,B
DOUBLE PRECISION, DIMENSION(6,6) :: OUT
INTEGER                          :: M,N

        M=1
        N=1
        OUT(1,:)=(/A(M,1)*B(N,1), A(M,2)*B(N,2), A(M,3)*B(N,3), A(M,1)*B(N,2)+ A(M,2)*B(N,1), A(M,1)*B(N,3)+A(M,3)*B(N,1), A(M,2)*B(N,3)+A(M,3)*B(N,2) /)
        M=2
        N=2
        OUT(2,:)=(/A(M,1)*B(N,1), A(M,2)*B(N,2), A(M,3)*B(N,3), A(M,1)*B(N,2)+ A(M,2)*B(N,1), A(M,1)*B(N,3)+A(M,3)*B(N,1), A(M,2)*B(N,3)+A(M,3)*B(N,2) /)
        M=3
        N=3
        OUT(3,:)=(/A(M,1)*B(N,1), A(M,2)*B(N,2), A(M,3)*B(N,3), A(M,1)*B(N,2)+ A(M,2)*B(N,1), A(M,1)*B(N,3)+A(M,3)*B(N,1), A(M,2)*B(N,3)+A(M,3)*B(N,2) /)
        M=1
        N=2
        OUT(4,:)=(/A(M,1)*B(N,1), A(M,2)*B(N,2), A(M,3)*B(N,3), A(M,1)*B(N,2)+ A(M,2)*B(N,1), A(M,1)*B(N,3)+A(M,3)*B(N,1), A(M,2)*B(N,3)+A(M,3)*B(N,2) /)
        M=3
        N=1
        OUT(5,:)=(/A(M,1)*B(N,1), A(M,2)*B(N,2), A(M,3)*B(N,3), A(M,1)*B(N,2)+ A(M,2)*B(N,1), A(M,1)*B(N,3)+A(M,3)*B(N,1), A(M,2)*B(N,3)+A(M,3)*B(N,2) /)
        M=2
        N=3
        OUT(6,:)=(/A(M,1)*B(N,1), A(M,2)*B(N,2), A(M,3)*B(N,3), A(M,1)*B(N,2)+ A(M,2)*B(N,1), A(M,1)*B(N,3)+A(M,3)*B(N,1), A(M,2)*B(N,3)+A(M,3)*B(N,2) /)


END SUBROUTINE BAR_DYAD



end module mater_large
