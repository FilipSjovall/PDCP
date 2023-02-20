module elem_small_cont_2d

! Last rev.:
!  H. Hallberg, 2011-04-26
!  M. Ristinmaa 2011-04-27
!   - Change order of variables in call
!   - s-command only calculates the strains
!  M. Ristinmaa 2011-04-28
!   - Renamed module, renamed elements inspired by Abaqus
!  M. Ristinmaa 2011-04-29
!   - included c2d3 element
!  M. Ristinmaa 2011-04-30
!   - renamed module lin to small
!  M. Ristinmaa 2012-06-11
!   - Bugg, Hooke plane strain change to 3x3 matrix according to manual
!------------------------------------------------------------------------------


!
!   Module elem_lin_cont_2d contains element subroutines for
!   an isoparametric quadrilateral 4-node 2D
!   element. Formulation for small deformations.
!
 
implicit none
 
! Interfaces
! ----------

interface c2d3_e
   module procedure plante
end interface

interface c2d3_s 
   module procedure plants
end interface

interface c2d3_f
   module procedure plantf
end interface

interface c2d4_e
   module procedure plani4e_KeFe1
   module procedure plani4e_KeFe
   module procedure plani4e_Ke1
   module procedure plani4e_Ke2
end interface
private :: plani4e_KeFe, plani4e_Ke1, plani4e_KeFe1, plani4e_Ke2

interface c2d4_s
   module procedure plani4s_eset
   module procedure plani4s_eseteci
end interface
private :: plani4s_eset, plani4s_eseteci

private det2, inv2

!------------------------------------------------------------------------------
contains


subroutine plante(Ke,elcoord,ep,D)
  implicit none
  double precision        :: Ke(6,6), elcoord(:,:)
  double precision        :: x(3), y(3), ep, D(3,3)

  double precision        :: detA, idetA, t, A
  double precision        :: B(3,6)

  x=elcoord(1,:)
  y=elcoord(2,:)
  detA=x(2)*y(3)+x(1)*y(2)+y(1)*x(3)-x(2)*y(1)-x(3)*y(2)-y(3)*x(1)
  A=detA*0.5D0
  idetA=1d0/detA

  t=ep

  B=reshape((/y(2)-y(3),0d0      ,y(3)-y(1),0d0      ,y(1)-y(2),0d0, &
            0d0      ,x(3)-x(2),0d0      ,x(1)-x(3),0d0      ,x(2)-x(1), &
            x(3)-x(2),y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),y(1)-y(2)/), &
           [3,6],order=[2,1])*idetA

  Ke=matmul(matmul(transpose(B),D),B)*A*t

return
end subroutine plante

subroutine plants(ee,elcoord,ed)
  implicit none
  double precision        :: Ke(6,6), elcoord(2,3)
  double precision        :: x(3), y(3), ee(3), ed(6)

  double precision        :: detA, idetA, t, A
  double precision        :: B(3,6)

  x=elcoord(1,:)
  y=elcoord(2,:)
  detA=x(2)*y(3)+x(1)*y(2)+y(1)*x(3)-x(2)*y(1)-x(3)*y(2)-y(3)*x(1)
  A=detA*0.5D0
  idetA=1d0/detA


  B=reshape((/y(2)-y(3),0d0      ,y(3)-y(1),0d0      ,y(1)-y(2),0d0, &
            0d0      ,x(3)-x(2),0d0      ,x(1)-x(3),0d0      ,x(2)-x(1), &
            x(3)-x(2),y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),y(1)-y(2)/), &
           [3,6],order=[2,1])*idetA

  ee=matmul(B,ed)

return
end subroutine plants

subroutine plantf(fe,elcoord,ep,es)
  implicit none
  double precision        :: fe(6), elcoord(2,3), es(3)
  double precision        :: x(3), y(3), ep

  double precision        :: detA, idetA, t, A
  double precision        :: B(3,6)

  x=elcoord(1,:)
  y=elcoord(2,:)
  detA=x(2)*y(3)+x(1)*y(2)+y(1)*x(3)-x(2)*y(1)-x(3)*y(2)-y(3)*x(1)
  A=detA*0.5D0
  idetA=1d0/detA

  t=ep;

  B=reshape((/y(2)-y(3),0d0      ,y(3)-y(1),0d0      ,y(1)-y(2),0d0, &
            0d0      ,x(3)-x(2),0d0      ,x(1)-x(3),0d0      ,x(2)-x(1), &
            x(3)-x(2),y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),y(1)-y(2)/), &
           [3,6],order=[2,1])*idetA

  fe=matmul(transpose(B),es)*A*t

return
end subroutine plantf


subroutine plani4e_Ke1(Ke,elcoord,ep,D)
   implicit none

   double precision, intent(in)  :: D(:,:), elcoord(:,:), ep
   double precision, intent(out) :: Ke(8,8)
   double precision              :: eq(2), fe(8)
   double precision              :: Dg(3,3,8)

   integer                       :: ii
   eq = 0D0

   do ii=1,8
     Dg(:,:,ii)=D
   enddo
   call plani4e_KeFe(Ke,fe,elcoord,ep,Dg,eq)

   return

end subroutine plani4e_Ke1

subroutine plani4e_Ke2(Ke,elcoord,ep,D)
   implicit none

   double precision, intent(in)  :: D(:,:,:), elcoord(:,:), ep
   double precision, intent(out) :: Ke(8,8)
   double precision              :: eq(2), fe(8)
   double precision              :: Dg(3,3,8)

   eq = 0D0

   call plani4e_KeFe(Ke,fe,elcoord,ep,D,eq)

   return

end subroutine plani4e_Ke2

subroutine plani4e_KeFe1(Ke,fe,elcoord,ep,D,eq)
   implicit none
   double precision, intent(in)  :: D(:,:), eq(:), elcoord(:,:), ep
   double precision, intent(out) :: Ke(8,8), fe(8)

   double precision              :: Dg(3,3,8)

   integer                       :: ii

   do ii=1,8
     Dg(:,:,ii)=D
   enddo
   call plani4e_KeFe(Ke,fe,elcoord,ep,Dg,eq)

   return

end subroutine plani4e_KeFe1

subroutine plani4e_KeFe(Ke,fe,elcoord,ep,D,eq)
   implicit none

   double precision, intent(in)  :: D(:,:,:), eq(:), elcoord(:,:), ep
   double precision, intent(out) :: Ke(8,8), fe(8)
   double precision              :: Be(3,8), b(2), t, g1, g2, w1, w2, coord(2,4) 
   double precision              :: detJ, N2(2,8), JTinv(2,2), dNx(2,4), Dm(3,3)
   double precision, allocatable :: gp(:,:), w(:,:), xsi(:), eta(:)
   double precision, allocatable :: N(:,:), dNr(:,:), JT(:,:), wp(:)
   integer                       :: ptype, ir, ngp, r2, i, j, indx(2)

   t          = ep     ! Element thickness
   ir         = 2       ! Integration order (1, 2 and 3 implemented)
   ngp        = ir*ir   ! Number of Gauss points
   b          = eq      ! Uniform element load [N/m²]
   coord(1,:) = elcoord(1,:) ! x-coordinates of element nodes
   coord(2,:) = elcoord(2,:) ! y-coordinates of element nodes
   Ke         = 0D0
   fe         = 0D0

   ! Gauss points
   ! ------------
    allocate(gp(4,2),w(4,2),xsi(4),eta(4),N(4,4),wp(4))
    g1 = 0.577350269189626D0
    w1 = 1D0
    gp(:,1) = (/ -g1, g1, -g1, g1 /)
    gp(:,2) = (/ -g1, -g1, g1, g1 /)
    w(:,1)  = (/ w1, w1, w1, w1 /)
    w(:,2)  = (/ w1, w1, w1, w1 /)
    wp  = w(:,1)*w(:,2)      
    xsi = gp(:,1)
    eta = gp(:,2)

   r2  = ngp*2
   allocate(dNr(r2,4),JT(r2,2))
   dNr = 0D0
   JT  = 0D0
   N   = 0D0

   ! Shape functions
   ! ---------------
   do i=1,size(N,1)
      N(i,1) = (1-xsi(i))*(1-eta(i))/4
      N(i,2) = (1+xsi(i))*(1-eta(i))/4
      N(i,3) = (1+xsi(i))*(1+eta(i))/4
      N(i,4) = (1-xsi(i))*(1+eta(i))/4
   end do
   do i=1,size(gp,1)
      dNr(i*2-1,1) = -(1-eta(i))/4
      dNr(i*2-1,2) =  (1-eta(i))/4
      dNr(i*2-1,3) =  (1+eta(i))/4
      dNr(i*2-1,4) = -(1+eta(i))/4
      dNr(i*2,1)   = -(1-xsi(i))/4
      dNr(i*2,2)   = -(1+xsi(i))/4
      dNr(i*2,3)   =  (1+xsi(i))/4
      dNr(i*2,4)   =  (1-xsi(i))/4
   end do
   JT = matmul(dNr,transpose(coord))

   do i=1,ngp
        indx = (/ 2*i-1, 2*i /)
	detJ = det2(JT(indx,:))
        JTinv = JT(indx,:)
	call inv2(JTinv)
	dNx = matmul(JTinv,dNr(indx,:))
	do j=1,4
	   Be(1,j*2-1) = dNx(1,j)
	   Be(2,j*2)   = dNx(2,j)
	   Be(3,j*2-1) = dNx(2,j)
	   Be(3,j*2)   = dNx(1,j)
	end do
	do j=1,4
	   N2(1,j*2-1) = N(i,j)
	   N2(1,j*2)   = N(i,j)
	end do
	Ke = Ke + wp(i)*t*detJ*matmul(transpose(Be),matmul(D(:,:,i),Be))
	fe = fe + wp(i)*t*detJ*matmul(transpose(N2),b)
  end do

  deallocate(gp,w,xsi,eta,N,dNr,JT,wp)
      
  return
      
end subroutine plani4e_KeFe


subroutine plani4s_eset(et,elcoord,ed)
   implicit none

   double precision, intent(in)    :: elcoord(:,:), ed(:)
   double precision, intent(inout) :: et(:,:)
   double precision                :: eci(4,2)
   
   call plani4s_eseteci(et,elcoord,ed,eci)
   
   return
      
end subroutine plani4s_eset

   
subroutine plani4s_eseteci(et,elcoord,ed,eci)
   implicit none

   double precision, intent(in)  :: elcoord(:,:), ed(:)
   double precision, intent(out) :: eci(:,:), et(:,:)
   double precision              :: Be(3,8), t, g1, g2, w1, w2, coord(2,4), detJ, N2(2,8), JTinv(2,2), dNx(2,4) 
   double precision, allocatable :: gp(:,:), w(:,:), xsi(:), eta(:), N(:,:), dNr(:,:), JT(:,:), wp(:)
   integer                       :: ptype, ir, ngp, r2, i, j, indx(2)

   t          = 1d0   ! Element thickness
   ir         = 2   ! Integration order (1, 2 and 3 implemented)
   ngp        = ir*ir   ! Number of Gauss points
   coord(1,:) = elcoord(1,:) ! x-coordinates of element nodes
   coord(2,:) = elcoord(2,:) ! y-coordinates of element nodes

   ! Gauss points
   ! ------------
      allocate(gp(4,2),w(4,2),xsi(4),eta(4),N(4,4),wp(4))
      g1 = 0.577350269189626D0
      w1 = 1D0
      gp(:,1) = (/ -g1, g1, -g1, g1 /)
      gp(:,2) = (/ -g1, -g1, g1, g1 /)
      w(:,1)  = (/ w1, w1, w1, w1 /)
      w(:,2)  = (/ w1, w1, w1, w1 /)
      wp  = w(:,1)*w(:,2)      
      xsi = gp(:,1)
      eta = gp(:,2)
   r2  = ngp*2
   allocate(dNr(r2,4),JT(r2,2))
   dNr = 0D0
   JT  = 0D0
   N   = 0D0

   ! Shape functions
   ! ---------------
   do i=1,size(N,1)
      N(i,1) = (1-xsi(i))*(1-eta(i))/4
      N(i,2) = (1+xsi(i))*(1-eta(i))/4
      N(i,3) = (1+xsi(i))*(1+eta(i))/4
      N(i,4) = (1-xsi(i))*(1+eta(i))/4
   end do
   do i=1,size(gp,1)
      dNr(i*2-1,1) = -(1-eta(i))/4
      dNr(i*2-1,2) =  (1-eta(i))/4
      dNr(i*2-1,3) =  (1+eta(i))/4
      dNr(i*2-1,4) = -(1+eta(i))/4
      dNr(i*2,1)   = -(1-xsi(i))/4
      dNr(i*2,2)   = -(1+xsi(i))/4
      dNr(i*2,3)   =  (1+xsi(i))/4
      dNr(i*2,4)   =  (1-xsi(i))/4
   end do

   JT = matmul(dNr,transpose(coord))
   eci = matmul(N,transpose(coord))

      do i=1,ngp
         indx = (/ 2*i-1, 2*i /)
	 detJ = det2(JT(indx,:))
         JTinv = JT(indx,:)
	 call inv2(JTinv)
	 dNx = matmul(JTinv,dNr(indx,:))
	 do j=1,4
	    Be(1,j*2-1) = dNx(1,j)
	    Be(2,j*2)   = dNx(2,j)
	    Be(3,j*2-1) = dNx(2,j)
	    Be(3,j*2)   = dNx(1,j)
	 end do
         et(:,i) = matmul(Be,ed)
      end do

   deallocate(gp,w,xsi,eta,N,dNr,JT,wp)
    	  
   return

end subroutine plani4s_eseteci


subroutine c2d4_f(ef,elcoord,ep,es)
   implicit none

   double precision, intent(in)  :: elcoord(:,:), ep, es(:,:)
   double precision, intent(out) :: ef(:)
   double precision              :: Be(3,8), t, g1, g2, w1, w2, coord(2,4), detJ, N2(2,8), JTinv(2,2), dNx(2,4), fint(8), stress(3), tmp(8,4)
   !double precision              :: ex(4), ey(4)
   double precision, allocatable :: gp(:,:), w(:,:), xsi(:), eta(:), N(:,:), dNr(:,:), JT(:,:), wp(:)
   integer                       :: ptype, ir, ngp, r2, i, j, indx(2)

   stop 'c2d4_f is not debugged'
   
   ptype      = 1   ! Analysis type (1=plane stress, 2=plane strain)
   t          = ep   ! Element thickness
   ir         = 2   ! Integration order (1, 2 and 3 implemented)
   ngp        = ir*ir   ! Number of Gauss points
   coord(1,:) = elcoord(:,1) ! x-coordinates of element nodes
   coord(2,:) = elcoord(:,2) ! y-coordinates of element nodes
   fint       = 0D0     ! Internal force vector

   ! Gauss points
   ! ------------
   if (ir==1) then
      allocate(gp(1,2),w(1,2),xsi(2),eta(2),N(1,4),wp(1))
      g1 = 0D0
      w1 = 2D0
      gp(:,1) = (/ g1, g2 /)
      w(:,1)  = (/ w1, w1 /)
      wp  = w1*w1
      xsi = g1
      eta = g2
   else if (ir==2) then
      allocate(gp(4,2),w(4,2),xsi(4),eta(4),N(4,4),wp(4))
      g1 = 0.577350269189626D0
      w1 = 1D0
      gp(:,1) = (/ -g1, g1, -g1, g1 /)
      gp(:,2) = (/ -g1, -g1, g1, g1 /)
      w(:,1)  = (/ w1, w1, w1, w1 /)
      w(:,2)  = (/ w1, w1, w1, w1 /)
      wp  = w(:,1)*w(:,2)      
      xsi = gp(:,1)
      eta = gp(:,2)
   else if (ir==3) then
      allocate(gp(9,2),w(9,2),xsi(9),eta(9),N(9,4),wp(9))
      g1 = 0.774596699241483D0
      g2 = 0D0
      w1 = 0.555555555555555D0
      w2 = 0.888888888888888D0
      gp(:,1) = (/ -g1, -g2, g1, -g1, g2, g1, -g1, g2, g1 /)
      gp(:,2) = (/ -g1, -g1, -g1, g2, g2, g2, g1, g1, g1  /)
      w(:,1)  = (/ w1, w2, w1, w1, w2, w1, w1, w2, w1 /)
      w(:,2)  = (/ w1, w1, w1, w2, w2, w2, w1, w1, w1 /)
      wp  = w(:,1)*w(:,2)      
      xsi = gp(:,1)
      eta = gp(:,2)
   end if
   r2  = ngp*2
   allocate(dNr(r2,4),JT(r2,2))
   dNr = 0D0
   JT  = 0D0
   N   = 0D0

   ! Shape functions
   ! ---------------
   do i=1,size(N,1)
      N(i,1) = (1-xsi(i))*(1-eta(i))/4
      N(i,2) = (1+xsi(i))*(1-eta(i))/4
      N(i,3) = (1+xsi(i))*(1+eta(i))/4
      N(i,4) = (1-xsi(i))*(1+eta(i))/4
   end do
   do i=1,size(gp,1)
      dNr(i*2-1,1) = -(1-eta(i))/4
      dNr(i*2-1,2) =  (1-eta(i))/4
      dNr(i*2-1,3) =  (1+eta(i))/4
      dNr(i*2-1,4) = -(1+eta(i))/4
      dNr(i*2,1)   = -(1-xsi(i))/4
      dNr(i*2,2)   = -(1+xsi(i))/4
      dNr(i*2,3)   =  (1+xsi(i))/4
      dNr(i*2,4)   =  (1-xsi(i))/4
   end do

   JT = matmul(dNr,transpose(coord))

   ! Plane stress
   ! ------------
   if (ptype==1) then
      do i=1,ngp
         indx = (/ 2*i-1, 2*i /)
	 detJ = det2(JT(indx,:))
         JTinv = JT(indx,:)
	 call inv2(JTinv)
	 dNx = matmul(JTinv,dNr(indx,:))
	 do j=1,4
	    Be(1,j*2-1) = dNx(1,j)
	    Be(2,j*2)   = dNx(2,j)
	    Be(3,j*2-1) = dNx(2,j)
	    Be(3,j*2)   = dNx(1,j)
	 end do
         fint = fint + matmul(transpose(Be),es(:,i))*wp(i)*detJ*t
      end do
      ef(:) = fint(:)
   else if (ptype==2) then
      ! Plane strain
      ! ------------
      do i=1,ngp
	 indx = (/ 2*i-1, 2*i /)
	 detJ = det2(JT(indx,:))
         JTinv = JT(indx,:)
	 call inv2(JTinv)
	 dNx = matmul(JTinv,dNr(indx,:))
	 do j=1,4
	    Be(1,j*2-1) = dNx(1,j)
            Be(2,j*2)   = dNx(2,j)
            Be(3,j*2-1) = dNx(2,j)
            Be(3,j*2)   = dNx(1,j)
         end do
         stress(:) = es((/1,2,4/),i)
         fint = fint + matmul(transpose(Be),stress)*wp(i)*detJ*t
      end do
      ef(:) = fint(:)
   end if

   deallocate(gp,w,xsi,eta,N,dNr,JT,wp)
      	  
   return

end subroutine c2d4_f

! help routines

function DET2(A)
  IMPLICIT NONE
  DOUBLE PRECISION                :: DET2
  DOUBLE PRECISION                :: A(:,:)

  DET2=A(1,1)*A(2,2)-A(2,1)*A(1,2) 

RETURN
END function det2


! REWRITE THIS TO INVERSE OF 2x2 MATRIX
SUBROUTINE INV2(A)
   IMPLICIT NONE
   INTEGER                        :: N, INFO1, INFO2, LWORK
   DOUBLE PRECISION               :: A(:,:)
   INTEGER                        :: IPIV(size(A,1))
   DOUBLE PRECISION               :: WORK(size(A,1))
   N=size(A,1)
   LWORK=N
   CALL DGETRF(N,N,A,N,IPIV,INFO1)
   CALL DGETRI(N,A,N,IPIV,WORK,LWORK,INFO2)


RETURN
END subroutine inv2


end module elem_small_cont_2d
