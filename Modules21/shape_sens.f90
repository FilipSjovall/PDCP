module shape_sens

implicit none

contains

subroutine plani4se_nl(DVe,DFe,D_el,ngp,npar,dXdalpha,ed,es,coord)
implicit none

INTEGER                                       :: r2, ngp, i,  ipar, ll, npar
DOUBLE PRECISION                              :: g1, g2, w1, w2,l11,l12,l22,l21,l33,detJ,xbar,D_J_D_alpha
DOUBLE PRECISION, PARAMETER                   :: pi=3.14159265358979D0
DOUBLE PRECISION, DIMENSION(ngp,4)            :: N
DOUBLE PRECISION, DIMENSION(ngp)              :: eta, xsi,wp, w_1, w_2, Phi
DOUBLE PRECISION, DIMENSION(2*ngp,4),TARGET   :: dNr
DOUBLE PRECISION, DIMENSION(4,8)              :: BL0, BL1, BL, D_B_L0_D_alpha, D_B_L1_D_alpha, D_B_L_D_alpha
DOUBLE PRECISION, DIMENSION(5,8)              :: BNL,H,D_H_D_alpha
DOUBLE PRECISION, DIMENSION(2*ngp,2)          :: JT
DOUBLE PRECISION, DIMENSION(2,4)              :: dNx,D_G_D_alpha
DOUBLE PRECISION, DIMENSION(4,2)              :: coord
DOUBLE PRECISION, DIMENSION(5,5)              :: S
DOUBLE PRECISION, POINTER                     :: D(:,:),Cp(:,:)
DOUBLE PRECISION, DIMENSION(ngp,4,4), TARGET  :: D_el
DOUBLE PRECISION, POINTER                     :: dNr_p(:,:)
DOUBLE PRECISION, DIMENSION(2,2)              :: JTinv, D_J_D_alphad
DOUBLE PRECISION, DIMENSION(8)                :: ed
DOUBLE PRECISION, DIMENSION(ngp,3,3)          :: es
DOUBLE PRECISION, DIMENSION(4*npar,1)         :: DVe
DOUBLE PRECISION, DIMENSION(4*npar,8,1)       :: DFe
DOUBLE PRECISION, DIMENSION(4,5)		          :: L,D_L_D_alpha
DOUBLE PRECISION, DIMENSION(4,1) 		          :: D_S_H
DOUBLE PRECISION, DIMENSION(4*npar,4,2)       :: dXdalpha
DOUBLE PRECISION, DIMENSION(4,2)		          :: D_X_D_alpha 


IF (ngp==1) THEN
	 g1=0.0D0; w1=2.0D0
	 xsi=(/g1/)
	 eta=(/g1/)
	 w_1=(/w1 /)
	 w_2=(/w1 /)
ELSE IF (ngp==4) THEN
	 g1=0.577350269189626D0; w1=1D0
    	 xsi=(/-g1, g1,-g1, g1/)
       eta=(/-g1, -g1,g1, g1/)
       w_1=(/w1, w1, w1, w1/)
       w_2=w_1
ELSE IF (ngp==9) THEN ! ????????????
       g1=0.774596699241483D0; g2=0.D0;
       w1=0.555555555555555D0; w2=0.888888888888888D0;
       xsi=(/-g1,-g2, g1,-g1, g2, g1,-g1, g2, g1/)
       eta=(/-g1,-g1,-g1, g2, g2, g2, g1, g1, g1/)
       w_1=(/w1, w2, w1, w1, w2, w1, w1, w2, w1/)
       w_2=(/w1, w1, w1, w2, w2, w2, w1, w1, w1/)
ELSE
    WRITE(*,*) 'Used number of integration points not implemented'
END IF
       wp=w_1*w_2


!--------- shape functions -----------------------------------
  N(:,1)=(1D0-xsi)*(1D0-eta)/4.D0;  N(:,2)=(1D0+xsi)*(1D0-eta)/4.D0
  N(:,3)=(1D0+xsi)*(1D0+eta)/4.D0;  N(:,4)=(1D0-xsi)*(1D0+eta)/4.D0

  dNr(1:2*ngp-1:2,1)=-(1D0-eta)/4.D0;     dNr(1:2*ngp-1:2,2)= (1D0-eta)/4.D0
  dNr(1:2*ngp-1:2,3)= (1D0+eta)/4.D0;     dNr(1:2*ngp-1:2,4)=-(1D0+eta)/4.D0
  dNr(2:2*ngp:2,1)=-(1D0-xsi)/4.D0;       dNr(2:2*ngp:2,2)=-(1D0+xsi)/4.D0
  dNr(2:2*ngp:2,3)= (1D0+xsi)/4.D0;       dNr(2:2*ngp:2,4)= (1D0-xsi)/4.D0

  
  
  JT=MATMUL(dNr,coord)

!DO i=1,npar
	!DFe(ipar,:) = 0d0
!END DO 


!--------- axisymmetric conditions -----------------------------
DVe = 0d0
DFe = 0d0
DO i=1,ngp
 
  detJ=JT(2*i-1,1)*JT(2*i,2)-JT(2*i-1,2)*JT(2*i,1)
  IF (detJ<0) THEN
    WRITE(*,*)'Jacobideterminant equal or less than zero!'
  END IF
  JTinv=1D0/detJ*RESHAPE((/JT(2*i,2), -JT(2*i,1),-JT(2*i-1,2),JT(2*i-1,1)/),(/2,2/))

  dNr_p=>dNr(2*i-1:2*i,:)

  dNx=MATMUL(JTinv,dNr_P)

  xbar=SUM(N(i,:)*coord(:,1))
  BL0(1,1:8:2)=dNx(1,:)
  BL0(2,2:8:2)=dNx(2,:)
  BL0(3,1:8:2)=dNx(2,:)  
  BL0(3,2:8:2)=dNx(1,:)
  BL0(4,1:8:2)=N(i,:)/xbar 	
  
  l11=SUM(dNx(1,:)*ed(1:8:2))
  l22=SUM(dNx(2,:)*ed(2:8:2))
  l12=SUM(dNx(2,:)*ed(1:8:2))
  l21=SUM(dNx(1,:)*ed(2:8:2))
  l33=SUM(N(i,:)*ed(1:8:2)/xbar)
		
  BL1(1,:) = (/ dNx(1,1)*l11, dNx(1,1)*l21, dNx(1,2)*l11, dNx(1,2)*l21,&
              	dNx(1,3)*l11, dNx(1,3)*l21, dNx(1,4)*l11, dNx(1,4)*l21/)
  BL1(2,:) = (/ dNx(2,1)*l12, dNx(2,1)*l22, dNx(2,2)*l12, dNx(2,2)*l22,&
              	dNx(2,3)*l12, dNx(2,3)*l22, dNx(2,4)*l12, dNx(2,4)*l22/)
  BL1(3,:)=(/(l11*dNx(2,1)+l12*dNx(1,1)),  (l21*dNx(2,1)+l22*dNx(1,1)),&
            (l11*dNx(2,2)+l12*dNx(1,2)),  (l21*dNx(2,2)+l22*dNx(1,2)),& 
            (l11*dNx(2,3)+l12*dNx(1,3)),  (l21*dNx(2,3)+l22*dNx(1,3)),& 
            (l11*dNx(2,4)+l12*dNx(1,4)),  (l21*dNx(2,4)+l22*dNx(1,4))/)                
  BL1(4,:) =   l33*BL0(4,:); 

  S(1:2,1:2)=es(i,1:2,1:2)
  S(3:5,3:5)=es(i,:,:)
  BL=BL0+BL1		
 
  D=>D_el(1:4,1:4,i)
 

	DO ipar=1,4*npar
		D_X_D_alpha = dXdalpha(ipar,:,:)		
		D_G_D_alpha = MATMUL(MATMUL(-dNx,D_X_D_alpha),dNx)
		D_J_D_alphad= MATMUL(dNx,D_X_D_alpha)
		D_J_D_alpha = 0d0

		! Trace
		do ll=1,size(D_J_D_alphad,2)
			D_J_D_alpha = D_J_D_alpha + detJ*D_J_D_alphad(ll,ll) 
		end do
		H = 0d0
		H(1,1:7:2) = dNx(1,:)
		H(2,1:7:2) = dNx(2,:)
		H(3,2:8:2) = dNx(1,:)
		H(4,2:8:2) = dNx(2,:)
		H(5,1:7:2) = N(i,:)
	
		L = 0d0
		L = reshape( (/ l11, 0d0, l21, 0d0, 0d0,&
		       					0d0, l12, 0d0, l22, 0d0,&
		       					l12, l11, l22, l21, 0d0,&
			 		 					0d0, 0d0, 0d0, 0d0, l33/xbar /), shape(D_L_D_alpha),ORDER=(/2,1/))	

		D_H_D_alpha          = 0d0
		D_H_D_alpha(1,1:7:2) = D_G_D_alpha(1,:)
		D_H_D_alpha(2,1:7:2) = D_G_D_alpha(2,:)
		D_H_D_alpha(3,2:8:2) = D_G_D_alpha(1,:)
		D_H_D_alpha(4,2:8:2) = D_G_D_alpha(2,:)
    D_H_D_alpha(5,:)     = 0d0

		D_L_D_alpha          = 0d0
		D_L_D_alpha(1,1:3:2) = (/ sum(D_G_D_alpha(1,:)*ed(1:7:2)), sum(D_G_D_alpha(1,:)*ed(2:8:2))  /)
		D_L_D_alpha(2,2:4:2) = (/ sum(D_G_D_alpha(2,:)*ed(1:7:2)), sum(D_G_D_alpha(2,:)*ed(2:8:2))  /)		
		D_L_D_alpha(3,:)     = (/ sum(D_G_D_alpha(2,:)*ed(1:7:2)), sum(D_G_D_alpha(1,:)*ed(1:7:2)), &
						  								sum(D_G_D_alpha(2,:)*ed(2:8:2)), sum(D_G_D_alpha(1,:)*ed(2:8:2)), 0d0 /)
		D_L_D_alpha(4,:)     = (/ 0d0, 0d0, 0d0, 0d0, sum(N(i,:)*D_X_D_alpha(:,1))*(-2d0*l33/(xbar**2)) /)		

		D_B_L0_D_alpha          = 0 
		D_B_L0_D_alpha(1,1:7:2) = D_G_D_alpha(1,:)
		D_B_L0_D_alpha(2,2:8:2) = D_G_D_alpha(2,:)
		D_B_L0_D_alpha(3,1:7:2) = D_G_D_alpha(2,:)
		D_B_L0_D_alpha(3,2:8:2) = D_G_D_alpha(1,:)
		D_B_L0_D_alpha(4,1:7:2) = (/ -N(i,1)*sum(N(i,:)*D_X_D_alpha(:,1))/(xbar**2), -N(i,2)*sum(N(i,:)*D_X_D_alpha(:,1))/(xbar**2), &
						     								 -N(i,3)*sum(N(i,:)*D_X_D_alpha(:,1))/(xbar**2), -N(i,4)*sum(N(i,:)*D_X_D_alpha(:,1))/(xbar**2) /)
	
		D_B_L1_D_alpha          = MATMUL(D_L_D_alpha,H) + MATMUL(L,D_H_D_alpha)

		D_B_L_D_alpha 	        = D_B_L0_D_alpha + D_B_L1_D_alpha

		D_S_H                   = MATMUL(D_B_L0_D_alpha,reshape(ed,(/8,1/))) + MATMUL(0.5*D_B_L1_D_alpha,reshape(ed,(/8,1/)))
		D_S_H                   = MATMUL(D,D_S_H)


		DFe(ipar,:,:)           = DFe(ipar,:,:) + ( MATMUL(TRANSPOSE(D_B_L_D_alpha),RESHAPE((/es(i,1,1),es(i,2,2),es(i,1,2),es(i,3,3)/),(/4,1/)))*detJ &
																		+ MATMUL(TRANSPOSE(BL),D_S_H)*detJ &
																		+ MATMUL(TRANSPOSE(BL),RESHAPE((/es(i,1,1),es(i,2,2),es(i,1,2),es(i,3,3)/),(/4,1/)))*D_J_D_alpha ) *2d0*pi*xbar*wp(i) & 
				                        		+ MATMUL(TRANSPOSE(BL),RESHAPE((/es(i,1,1),es(i,2,2),es(i,1,2),es(i,3,3)/),(/4,1/)))*sum(N(i,:)*D_X_D_alpha(:,1))*2d0*pi*detJ*wp(i)

    DVe(ipar,1) = DVe(ipar,1) + (D_J_D_alpha*xbar + detJ*sum(N(i,:)*D_X_D_alpha(:,1)))*2d0*pi*wp(i)

															
		
	END DO
	 !write(*,*) 'H ' , H 
	 !write(*,*) 'L ' , L
	 !write(*,*) 'D_L_D_alpha ' , D_L_D_alpha 
	 !write(*,*) 'D_H_D_alpha ' , D_H_D_alpha 
   !write(*,*) 'D_B_L1_D_alpha ' , D_B_L1_D_alpha 
END DO


end subroutine plani4se_nl

end module shape_sens
