module matrix_util

! last modified
! M. Wallin 2011-09-16
! det3 moved to matrix_util
! getI moved to matrix_util
! inv3 moved to matrix_util
! inv3s moved to matrix_util
! getI  moved to matrix_util
! dev3  moved to matrix_util
! tr 3  moved to matrix_util
! eye  moved to matrix_util
! inverse  moved to matrix_util
! norm  moved to matrix_util


implicit none


!-------------------------------------------------------------------------------

contains

subroutine inv2(res,Ca)
  implicit none
  double precision                :: Ca(:,:), res(:,:)

  double precision                :: t4

  t4 = 0.1D1 / (Ca(1,1) * Ca(2,2) - Ca(1,2) * Ca(2,1))
  res(1,1) = Ca(2,2) * t4
  res(1,2) = -Ca(1,2) * t4
  res(2,1) = -Ca(2,1) * t4
  res(2,2) = Ca(1,1) * t4

  return
end subroutine inv2

! inverse 3x3 matrix
subroutine inv3(res,Ca)
  implicit none
  double precision                :: Ca(:,:), res(:,:)

  double precision                :: t4, t6, t8, t10, t12, t14, t17
! Maple code
  t4 = Ca(1,1)*Ca(2,2)
  t6 = Ca(1,1)*Ca(2,3)
  t8 = Ca(1,2)*Ca(2,1)
  t10 = Ca(1,3)*Ca(2,1)
  t12 = Ca(1,2)*Ca(3,1)
  t14 = Ca(1,3)*Ca(3,1)
  t17 = 1d0/(t4*Ca(3,3)-t6*Ca(3,2)-t8*Ca(3,3)+t10*Ca(3,2) &
         +t12*Ca(2,3)-t14*Ca(2,2))

  res(1,1) = (Ca(2,2)*Ca(3,3)-Ca(2,3)*Ca(3,2))*t17
  res(1,2) = -(Ca(1,2)*Ca(3,3)-Ca(1,3)*Ca(3,2))*t17
  res(1,3) = (Ca(1,2)*Ca(2,3)-Ca(1,3)*Ca(2,2))*t17
  res(2,1) = -(Ca(2,1)*Ca(3,3)-Ca(2,3)*Ca(3,1))*t17
  res(2,2) = (Ca(1,1)*Ca(3,3)-t14)*t17
  res(2,3) = -(t6-t10)*t17
  res(3,1) = (Ca(2,1)*Ca(3,2)-Ca(2,2)*Ca(3,1))*t17
  res(3,2) = -(Ca(1,1)*Ca(3,2)-t12)*t17
  res(3,3) = (t4-t8)*t17

  return
end subroutine inv3


! inverse of 3x3 symmetric matrix
subroutine inv3s(res,cs)
  implicit none
  double precision                :: res(:,:), cs(:,:)

  double precision                :: t2, t4, t7, t9, t12, t15, t20, t24, t30

  t2 = Cs(2,3)**2d0
  t4 = Cs(1,1)*Cs(2,2)
  t7 = Cs(1,2)**2d0
  t9 = Cs(1,2)*Cs(1,3)
  t12 = Cs(1,3)**2d0
  t15 = 1d0/(t4*Cs(3,3)-Cs(1,1)*t2-t7*Cs(3,3)+2*t9*Cs(2,3)-t12*Cs(2,2))
  t20 = (Cs(1,2)*Cs(3,3)-Cs(1,3)*Cs(2,3))*t15
  t24 = (Cs(1,2)*Cs(2,3)-Cs(1,3)*Cs(2,2))*t15
  t30 = (Cs(1,1)*Cs(2,3)-t9)*t15
  res(1,1) = (Cs(2,2)*Cs(3,3)-t2)*t15
  res(1,2) = -t20
  res(1,3) = t24
  res(2,1) = -t20
  res(2,2) = (Cs(1,1)*Cs(3,3)-t12)*t15
  res(2,3) = -t30
  res(3,1) = t24
  res(3,2) = -t30
  res(3,3) = (t4-t7)*t15

  return
end subroutine inv3s
 
function det2(F)
  implicit none
  double precision                 :: F(:,:), det2

   det2 = F(1,1) * F(2,2) - F(1,2) * F(2,1)
  return
end function det2

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

function dev3(in)
  implicit none
  double precision                :: dev3(3,3), in(3,3), tr, id(3,3)
  id=getI()
  tr=in(1,1)+in(2,2)+in(3,3)
  dev3=in-tr*id/3D0 
  return
end function dev3


function tr3(in)
  implicit none
  double precision                :: tr3, in(3,3)
  tr3=in(1,1)+in(2,2)+in(3,3)
  return
end function tr3


SUBROUTINE eye(UNIT,N)
   
   IMPLICIT NONE
   INTEGER                                    :: N, I, J
   DOUBLE PRECISION,    DIMENSION(N,N)        :: UNIT

   DO I=1,N
      DO J=1,N
         IF (I==J) THEN
            UNIT(I,J)=1D0
	 ELSE
	    UNIT(I,J)=0D0
	 END IF
      END DO
   END DO   
   
  

RETURN
END SUBROUTINE EYE


SUBROUTINE INVERSE(A,AINVERSE,N)
   IMPLICIT NONE
   INTEGER                                    :: N, INFO1, INFO2, LWORK
   DOUBLE PRECISION,    DIMENSION(N,N)        :: A, AINVERSE
   INTEGER, DIMENSION(N)                      :: IPIV
   DOUBLE PRECISION,    DIMENSION(N)          :: WORK   

 !  stop 'not implemented'
   LWORK=N
   AINVERSE=A
!   CALL DGETRF(N,N,AINVERSE,N,IPIV,INFO1)
!   CALL DGETRI(N,AINVERSE,N,IPIV,WORK,LWORK,INFO2)
RETURN
END subroutine Inverse


SUBROUTINE norm(VEC,NORM_VEC,N)
   
       IMPLICIT NONE
       INTEGER                                   :: N, I
       DOUBLE PRECISION,    DIMENSION(:)         :: VEC
       DOUBLE PRECISION                          :: NORM_VEC

       NORM_VEC=0D0
       DO I=1,N
          NORM_VEC=NORM_VEC+VEC(I)*VEC(I)
       END DO
       NORM_VEC=DSQRT(NORM_VEC)
 
RETURN
END subroutine norm

! Cartesian basis
subroutine cross_product(a,b,vec_out)
   implicit none
   double precision           :: a(:), b(:), vec_out(:)

   vec_out = 0d0
   vec_out(1) = a(2)*b(3)-a(3)*b(2)
   vec_out(2) = -(a(1)*b(3)-a(3)*b(1))
   vec_out(3) = a(1)*b(2)-a(2)*b(1)

end subroutine cross_product


end module matrix_util
