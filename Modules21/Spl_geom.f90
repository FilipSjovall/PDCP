module Spl_geom

implicit none


contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	     Spline functions			!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Spl(Bspl,u,n,ip,Ua)
implicit none
! u = coordinates 
! n = nbr control vertices - 1
! ip = spline degree
! Ua = knot vector

! Input
integer :: i,m,n,ip,p
double precision :: u,Ua(2*(ip+1)) 

! Output
double precision :: Bspl(size(Ua)-1,size(Ua)/2)



m = SIZE(Ua)-1

do i=0,(m-1)
	if(u.GE.Ua(i+1).AND.u.LT.Ua(i+2)) then
		Bspl(i+1,1) = 1D0
   else
		Bspl(i+1,1) = 0D0
   end if
end do

do p=1,ip
	do i=0,n
		Bspl(i+1,p+1) = 0D0
		if((ABS(Ua(i+p+1)-Ua(i+1))).gt.0D0) then
			Bspl(i+1,p+1) = Bspl(i+1,p+1) + ( (u-Ua(i+1))/(Ua(i+p+1)-Ua(i+1)) )*Bspl(i+1,p-1+1)
		end if
		if((ABS(Ua(i+p+2)-Ua(i+2))).gt.0D0) then	
			Bspl(i+1,p+1) = Bspl(i+1,p+1) + ( (Ua(i+p+1+1)-u)/(Ua(i+p+1+1)-Ua(i+1+1)) )*Bspl(i+1+1,p-1+1) 
		end if 
	end do
end do

if (u.EQ.1) then
  Bspl(m-ip-1+1,:)=1D0
end if

	
end subroutine Spl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!+-------------------------------------------------+!
!| Create geomtetry from defined control vertices  |!
!+-------------------------------------------------+!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wash_geometry(coord,coorduv,dof,edof,edofn,bc,V0,V,DV,Uk,Vk,nelu,nelv,alpha,xmax,ymax,H0,ri,L)
implicit none
! Input
integer, intent(in)		 			 :: nelu,nelv
double precision, intent(in) :: alpha(:),xmax,ymax,H0,ri,L


! Output
integer          						  :: edof(nelu*nelv,5),edofn(nelu*nelv,9),dof((nelv+1)*(nelu+1),2)
double precision              :: coord((nelv+1)*(nelu+1),2),coorduv((nelv+1)*(nelu+1),2),V0(size(alpha)/4,2,2),V(size(alpha)/4,2,2),DV(size(alpha),2,size(alpha),2)
double precision 							:: Uk(size(alpha)/2),Vk(4),bc(nelv+1+2,2)!,bc(2,2)
double precision,dimension(size(alpha)/4) :: alphavl,alphahl,alphavu,alphahu	
! Temporaries
integer                       :: ndvars,p,q,i,j,k,ll,m,n,rr,s,Ukl,Vkl,Uvk,Vvk,nel,pt,elnr,ii
double precision 							:: u,vv,idx,R(2,(nelu+1)*(nelv+1)),Lx,tmp1,tmp2
double precision						  :: Mu(size(Uk)-1,size(Uk)/2),Mv(size(Vk)-1,size(Vk)/2)!!,Bspl(:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 							!
! 							!
! 							!
! 							!	
!   Description: input, output etc		!
! 							!
! 							!
!							  ! 
! 							!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(*,*) 'post subroutine call'
!write(*,*) 'nelu,nelv,alpha',nelu,nelv,alpha
!write(*,*) 'xmax,ymax,H0,ri,L',xmax,ymax,H0,ri,L

ndvars = size(alpha)/4	


Lx    = xmax/(2d0*ndvars) 
! Divide alpha into upper/lower , horizontal/vertical

alphavl = alpha(1:ndvars) 
alphavu = alpha(ndvars+1:2*ndvars)
alphahl = alpha(2*ndvars+1:3*ndvars) 
alphahu = alpha(3*ndvars+1:4*ndvars)

V0 = 0D0
V  = 0D0
! Assign origin state of control vertices
do i=1,ndvars
	V0(i,1,1) = ri+(dble(i)-1D0)*xmax/(dble(ndvars)-1D0)-Lx/2D0 	
	V0(i,2,1) = (H0-ymax)-(H0)*(dble(i)-1D0)/(dble(ndvars)-1D0)+(0.75D0)*ymax
	V0(i,1,2) = ri+(dble(i)-1D0)*xmax/(dble(ndvars)-1D0)-Lx/2D0
	V0(i,2,2) = H0-(H0)*(dble(i)-1D0)/(dble(ndvars)-1D0)+(0.75D0)*ymax
end do
  V0(1,1,1) = V0(1,1,1) + Lx/2D0 
! 
q=1
p=ndvars-1
! nbr control vertices horizontal and vertical 
m = 1
n = ndvars-1
! Spline degrees
rr = p+n+1
s = q+m+1
! Knot-vector length
Ukl = rr+1
Vkl = s+1


! Fill Knot-vectors
do i=1,Ukl
	Uk(i) = 0d0
	if(i.gt.Ukl-p-1) then
		Uk(i) = 1d0
	end if
end do

do i=1,Vkl
	Vk(i) = 0d0
	if(i.gt.Vkl-q-1) then
		Vk(i) = 1d0
	end if
end do

Vk(Vkl/2+1:Vkl) = 1D0
Uk(Ukl/2+1:Ukl) = 1D0		

! Element distribution
Uvk = Ukl-  2*(p+1)
Vvk = Vkl - 2*(q+1)
do i=1,Uvk
	Uk(p+1+i)=dble(i)/(dble(Uvk)+1D0)
end do	
do i=1,Vvk
	Vk(q+1+i)=dble(i)/(dble(Vvk)+1D0)
end do

! Assign current coordinates of control vertices
do i=1,ndvars
 V(i,1,1) = V0(i,1,1) + Lx*alphahl(i) 
 V(i,2,1) = V0(i,2,1) + L*alphavl(i)
 V(i,1,2) = V0(i,1,2) + Lx*alphahu(i)
 V(i,2,2) = V0(i,2,2) + L*alphavu(i)
end do


! Sensitivity of control vertices w.r.t alpha
do i=1,ndvars*4
	if(i.le.ndvars) then
		DV(i,1,i,1) = 0D0
		DV(i,2,i,1) = L
		DV(i,1,i,2) = 0D0
		DV(i,2,i,2) = 0D0
	else if(i.gt.ndvars.and.i.le.ndvars*2) then
		DV(i,1,i,1) = 0D0
		DV(i,2,i,1) = 0D0
		DV(i,1,i,2) = 0D0
		DV(i,2,i,2) = L
	else if(i.gt.ndvars*2.and.i.le.ndvars*3) then
		DV(i,1,i,1) = Lx
		DV(i,2,i,1) = 0D0
		DV(i,1,i,2) = 0D0
		DV(i,2,i,2) = 0D0
	else
		DV(i,1,i,1) = 0D0
		DV(i,2,i,1) = 0D0
		DV(i,1,i,2) = Lx	
		DV(i,2,i,2) = 0D0
	end if
end do

!! Modify vertices such that bc-nodes won't move    
! Lower left
V(ndvars,1,1) = V0(ndvars,1,1) + 0.5*Lx
V(ndvars,2,1) = V0(ndvars,2,1) + 0.5*L

DV(ndvars,1,ndvars,1) = 0D0
DV(ndvars,2,ndvars,1) = 0D0

DV(ndvars*3,1,ndvars*3,1) = 0D0
DV(ndvars*3,2,ndvars*3,1) = 0D0

V(1,1,2) = V0(1,1,2) + 0.5*Lx
V(1,2,2) = V0(1,2,2) + alphavu(1)*L

DV(ndvars*3+1,1,ndvars*3+1,2) = 0D0
DV(ndvars*3+1,2,ndvars*3+1,2) = 0D0


! Create geometry
nel = nelu*nelv
pt  = 0
R   = 0d0
Mu  = 0d0
Mv  = 0d0
coord= 0d0
coorduv= 0d0
u=0d0
vv=0d0
!write(*,*) 'u,n,p,Uk',u,n,p,Uk
!write(*,*) 'vv,m,q,Vk',vv,m,q,Vk

do j=1,nelv+1
	do i=1,nelu+1
		pt = pt + 1
		u   = (i*1d0-1D0)/(nelu*1d0) 
		vv  = (j*1d0-1D0)/(nelv*1d0)
		Mu  = 0d0
		Mv  = 0d0
		call Spl(Mu,u,n,p,Uk)
		call Spl(Mv,vv,m,q,Vk)
			do k=0,m
				do ll=0,n							
						 R(1,pt) = R(1,pt) + V(ll+1,1,k+1)*Mu(ll+1,p+1)*Mv(k+1,q+1) 
						 R(2,pt) = R(2,pt) + V(ll+1,2,k+1)*Mu(ll+1,p+1)*Mv(k+1,q+1) 																	                 	 					
				end do
			end do
			coord(pt,:)   = R(:,pt) 
			coorduv(pt,1) = u 
			coorduv(pt,2) = vv
	end do
end do 
!write(*,*) 'coord:' , coord

!write(*,*) ,'R',R
! Define element topology
elnr = 0
do i=1,nelu
	do j=1,nelv
		elnr = elnr+1
		edof(elnr,1) = elnr    
		edof(elnr,2) = 1+(j-1)*(nelu+1)+(i-1)
		edof(elnr,3) = 1+(j-1)*(nelu+1)+(i-1)+1
		edof(elnr,4) = 1+(j-1)*(nelu+1)+(i-1)+1+(1+nelu)
		edof(elnr,5) = 1+(j-1)*(nelu+1)+(i-1)+(1+nelu)
	end do
end do

! edofn will be used 'as edof'
do i=1,nel
	edofn(i,1) = edof(i,1)  
	edofn(i,2) = 2*edof(i,2)-1
	edofn(i,3) = 2*edof(i,2)
	edofn(i,4) = 2*edof(i,3)-1
	edofn(i,5) = 2*edof(i,3)
	edofn(i,6) = 2*edof(i,4)-1
	edofn(i,7) = 2*edof(i,4)
	edofn(i,8) = 2*edof(i,5)-1
	edofn(i,9) = 2*edof(i,5)
end do

! Degrees of freedom
do i=1,((nelv+1)*(nelu+1))
	dof(i,1) = 2*i-1
	dof(i,2) = 2*i
end do

!===== Boundary conditions =======
! tmp1, tmp2 can be replaced by coord(idx,:)
! Lower right node
tmp1=coord(1,1)
tmp2=coord(1,2)
!idx = 1
do i=1,size(coord,1)
		if(tmp1.le.coord(i,1)) then	
			tmp1= coord(i,1)	
			if(tmp2.ge.coord(i,2)) then
			idx = i
			tmp1= coord(i,1)
			tmp2= coord(i,2)
			end if
		end if
end do

bc(1,1) = idx*2d0 ! vertical dof
bc(1,2) = 0d0


! Upper left node (load)
idx = 1
tmp1=coord(1,1)
tmp2=coord(1,2)
do i=1,size(coord,1)
		if(tmp1.ge.coord(i,1)) then
			if(tmp2.ge.coord(i,2)) then
			idx = i
			tmp1= coord(i,1)
			tmp2= coord(i,2)
			end if
		end if
end do

! New bc set-up
tmp1=coord(1,1)
idx = 2
do i=1,size(coord,1)
		if(tmp1.ge.coord(i,1)) then
			idx = idx+1
			bc(idx,1) = 2d0*i-1d0
			bc(idx,2) = 0d0
		end if
end do

bc(2,1) = 2d0*idx
bc(2,2) = 0D0 ! vertical dof

!write(*,*) 'bc:' , bc

end subroutine wash_geometry


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          Update geometry						      											!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine update_geom(coord,V,nelu,nelv,alpha,V0,Uk,Vk,L,xmax)
implicit none
! Input
integer, intent(in)		 			 :: nelu,nelv
double precision, intent(in) :: alpha(:),xmax,L,V0(:,:,:),Uk(:),Vk(:)

! Output
double precision              :: coord((nelv+1)*(nelu+1),2),coorduv((nelv+1)*(nelu+1),2), V(size(alpha)/4,2,2)

! Temporaries
integer                       :: ndvars,p,q,i,j,k,ll,m,n,Ukl,Vkl,Uvk,Vvk,nel,pt,elnr,ii
double precision 							:: u,vv,idx,R(2,(nelu+1)*(nelv+1)),Lx,tmp1,tmp2
double precision						  :: Mu(size(Uk)-1,size(Uk)/2),Mv(size(Vk)-1,size(Vk)/2)!!,Bspl(:,:)
double precision,dimension(size(alpha)/4) :: alphavl,alphahl,alphavu,alphahu	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Horizontal-knot-length parameter as well as nbr of control-vertices along each edge
ndvars = size(alpha)/4
Lx    = xmax/(2d0*ndvars) 
q=1
p=ndvars-1
! nbr control vertices horizontal and vertical 
m = 1
n = ndvars-1
! Divide alpha into upper/lower , horizontal/vertical
alphavl = alpha(1:ndvars) 
alphavu = alpha(ndvars+1:2*ndvars)
alphahl = alpha(2*ndvars+1:3*ndvars) 
alphahu = alpha(3*ndvars+1:4*ndvars)

! Assign current coordinates of control vertices
do i=1,ndvars
 V(i,1,1) = V0(i,1,1) + Lx*alphahl(i) 
 V(i,2,1) = V0(i,2,1) + L*alphavl(i)
 V(i,1,2) = V0(i,1,2) + Lx*alphahu(i)
 V(i,2,2) = V0(i,2,2) + L*alphavu(i)
end do

!! Modify vertices such that bc-nodes won't move    
! Lower left
V(ndvars,1,1) = V0(ndvars,1,1) + 0.5*Lx
V(ndvars,2,1) = V0(ndvars,2,1) + 0.5*L

V(1,1,2) = V0(1,1,2) + 0.5*Lx
V(1,2,2) = V0(1,2,2) + alphavu(1)*L

! Create geometry
nel = nelu*nelv

pt  = 0
R   = 0d0
Mu  = 0d0
Mv  = 0d0
coord= 0d0
coorduv= 0d0
u=0d0
vv=0d0
do j=1,nelv+1
	do i=1,nelu+1
		pt = pt + 1
		u   = (i*1d0-1D0)/(nelu*1d0) 
		vv  = (j*1d0-1D0)/(nelv*1d0)
		Mu  = 0d0
		Mv  = 0d0
		call Spl(Mu,u,n,p,Uk)
		call Spl(Mv,vv,m,q,Vk)
			do k=0,m
				do ll=0,n							
						 R(1,pt) = R(1,pt) + V(ll+1,1,k+1)*Mu(ll+1,p+1)*Mv(k+1,q+1) 
						 R(2,pt) = R(2,pt) + V(ll+1,2,k+1)*Mu(ll+1,p+1)*Mv(k+1,q+1) 																	                 	 					
				end do
			end do
			coord(pt,:)   = R(:,pt) 
			coorduv(pt,1) = u 
				coorduv(pt,2) = vv
		end do
end do 

end subroutine update_geom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Sensitivity of nodes w.r.t alpha		    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Nodesens(dXdalpha,eu,ev,Uk,Vk,DV,npar)
implicit none
! Input
integer :: npar
double precision, intent(in) :: eu(4),ev(4),Uk(2*npar),Vk(4),DV(4*npar,2,4*npar,2)
! Output
double precision, intent(out) :: dXdalpha(4*npar,4,2)
! Temporaries
integer :: r,s,n,m,p,q,ncorner,ipar,icorner,i,j
double precision :: Mu(2*npar-1,npar),Mv(3,2),u,v
r = size(Uk)-1
s = size(Vk)-1
n = size(DV,1)/4-1
m = size(DV,2)-1
p = r-n-1
q = s-m-1
dXdalpha = 0d0
Mu=0d0
Mv=0d0
ncorner=size(eu)

do ipar=1,4*npar
	do icorner=1,ncorner
		u = eu(icorner)
		v = ev(icorner)
		call Spl(Mu,u,n,p,Uk)
		call Spl(Mv,v,m,q,Vk)
		dXdalpha(ipar,icorner,1) = 0d0
		dXdalpha(ipar,icorner,2) = 0d0	
		do j=0,m 
			do i=0,n
				if(ipar.le.npar) then
				dXdalpha(ipar,icorner,1) = dXdalpha(ipar,icorner,1) + DV(i+1,1,ipar,j+1)*Mu(i+1,p+1)*Mv(j+1,q+1)
				dXdalpha(ipar,icorner,2) = dXdalpha(ipar,icorner,2) + DV(i+1,2,ipar,j+1)*Mu(i+1,p+1)*Mv(j+1,q+1) 
				else if(ipar.gt.npar.and.ipar.le.npar*2) then
				dXdalpha(ipar,icorner,1) = dXdalpha(ipar,icorner,1) + DV(n+1+i+1,1,ipar,j+1)*Mu(i+1,p+1)*Mv(j+1,q+1)
				dXdalpha(ipar,icorner,2) = dXdalpha(ipar,icorner,2) + DV(n+1+i+1,2,ipar,j+1)*Mu(i+1,p+1)*Mv(j+1,q+1) 
				else if(ipar.gt.npar*2.and.ipar.le.3*npar) then
				dXdalpha(ipar,icorner,1) = dXdalpha(ipar,icorner,1) + DV(2*(n+1)+i+1,1,ipar,j+1)*Mu(i+1,p+1)*Mv(j+1,q+1)
				dXdalpha(ipar,icorner,2) = dXdalpha(ipar,icorner,2) + DV(2*(n+1)+i+1,2,ipar,j+1)*Mu(i+1,p+1)*Mv(j+1,q+1) 
				else
				dXdalpha(ipar,icorner,1) = dXdalpha(ipar,icorner,1) + DV(3*(n+1)+i+1,1,ipar,j+1)*Mu(i+1,p+1)*Mv(j+1,q+1)
				dXdalpha(ipar,icorner,2) = dXdalpha(ipar,icorner,2) + DV(3*(n+1)+i+1,2,ipar,j+1)*Mu(i+1,p+1)*Mv(j+1,q+1) 
				end if
			end do
		end do
	end do
end do

end subroutine Nodesens



! Använda any(res==coord(i)) för att ta bort duplicates

end module Spl_geom
