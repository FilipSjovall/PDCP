module fem_element

! last modified
! M. Ristinmaa 2011-04-12
! -------------------------------------------------------------------------------------

!
! More elements should be included
!
! IMPORTANT SHOULD WE USE EX EY FORMAT OR COORD

implicit none

interface flw2te
   module procedure flw2te_KeFe
   module procedure flw2te_Ke
end interface


!----------------------------------------------------------------------

contains

subroutine flw2tm(me,fe,x,y,t)
  implicit none
  double precision, intent(in)  :: x(3), y(3)
  double precision, intent(out) :: me(3,3), fe(3)
  double precision              :: B(2,3), A, t

  double precision              :: N(3,3), gw, Ni(3)
  integer                       :: ip, i, j, status



  me = 0.0d0
  fe = 0.0d0

  A=0.5d0*(x(2)*y(3)+x(1)*y(2)+y(1)*x(3)-x(2)*y(1)-x(3)*y(2)-y(3)*x(1))

  N=reshape((/0.5d0,   0d0, 0.5d0, &
              0.5d0, 0.5d0,   0d0, &
                0d0, 0.5d0, 0.5d0/),[3,3],order=[2,1])
  gw=1d0/3d0*t*A;  

  do ip=1,3
    fe=fe+N(:,ip)*gw;
    me(1,:)=me(1,:)+N(1,ip)*N(:,ip)*gw
    me(2,:)=me(2,:)+N(2,ip)*N(:,ip)*gw
    me(3,:)=me(3,:)+N(3,ip)*N(:,ip)*gw
  end do

return
end subroutine flw2tm



subroutine flw2te_Ke(Ke,ex,ey,ep,D)
  implicit none
  double precision, intent(in)  :: ex(3), ey(3), ep, D(2,2)
  double precision, intent(out) :: Ke(3,3)
  double precision              :: eq, fe(3)

  Ke = 0.0d0
  fe = 0.0d0
  eq = 0.0d0

  call flw2te_KeFe(Ke,fe,ex,ey,ep,D,eq)

return
end subroutine flw2te_Ke

subroutine flw2te_KeFe(Ke,fe,ex,ey,ep,D,eq)
  implicit none
  double precision, intent(in)  :: ex(3), ey(3), ep, D(2,2), eq
  double precision, intent(out) :: Ke(3,3), fe(3)
  double precision              :: B(2,3), A, t, Q

  t = ep
  Q = eq
  Ke = 0.0d0
  fe = 0.0d0

  A = 0.5*((ex(2)*ey(3)-ex(3)*ey(2)) - (ex(1)*ey(3)-ex(3)*ey(1)) + (ex(1)*ey(2)-ex(2)*ey(1)))

  B(1,:) = (/ ey(2)-ey(3), ey(3)-ey(1), ey(1)-ey(2) /)
  B(2,:) = (/ ex(3)-ex(2), ex(1)-ex(3), ex(2)-ex(1) /)
  B = 0.5d0/A * B

  Ke = A*t*matmul(transpose(B),matmul(D,B))
  fe = Q*A*t/3.0d0*fe

return
end subroutine flw2te_KeFe



subroutine plante(Ke,x,y,ep,D)
  implicit none
  double precision        :: Ke(6,6)
  double precision        :: x(3), y(3), ep(2), D(3,3)

  double precision        :: detA, idetA, t, A
  double precision        :: B(3,6)

  detA=x(2)*y(3)+x(1)*y(2)+y(1)*x(3)-x(2)*y(1)-x(3)*y(2)-y(3)*x(1)
  A=detA*0.5D0
  idetA=1d0/detA

  t=ep(2);

  B=reshape((/y(2)-y(3),0d0      ,y(3)-y(1),0d0      ,y(1)-y(2),0d0, &
            0d0      ,x(3)-x(2),0d0      ,x(1)-x(3),0d0      ,x(2)-x(1), &
            x(3)-x(2),y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),y(1)-y(2)/), &
           [3,6],order=[2,1])*idetA

  Ke=matmul(matmul(transpose(B),D),B)*A*t;

return
end subroutine plante


end module fem_element
