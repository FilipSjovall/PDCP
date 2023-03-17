module elem_large_cont_2d6
! Element routines for 6-node triangular element
use matrix_util, only: inv2, det2, inv3, det3
use mater_large, only: getF
!
implicit none
!
double precision xsi(3,3), w(3), P0(2,3)
integer          index(3,3),ngp
!
parameter       (ngp  = 3)
parameter       (xsi  = (/ (/2.d0/3.d0,1.d0/3.d0,1.d0/3.d0/),(/1.d0/3.d0,2.d0/3.d0,1.d0/3.d0/),(/1.d0/3.d0,1.d0/3.d0,2.d0/3.d0/) /) )
!parameter       (xsi(2,:) = )
!parameter       (xsi(3,:) = )
parameter       (w    = [1.d0/3.d0,1.d0/3.d0,1.d0/3.d0])
parameter       (index= [[1,4,7],[2,5,8],[3,6,9]])
parameter       (P0   = (/ [0.d0,0.d0],[1.d0,0.d0],[0.d0,1.d0] /))
private xsi, w, index, P0

double precision N(3,6) ! Thress gauss points
data            (N(:,1) = xsi(1,:)*(2.d0*xsi(1,:)-1))
data            (N(:,2) = xsi(2,:)*(2.d0*xsi(2,:)-1))
data            (N(:,3) = xsi(3,:)*(2.d0*xsi(3,:)-1))
data            (N(:,4) = 4.d0*xsi(1,:)*xsi(2,:))
data            (N(:,5) = 4.d0*xsi(2,:)*xsi(3,:))
data            (N(:,6) = 4.d0*xsi(1,:)*xsi(3,:))
private N

double precision dNr(6,9)
! Gauss pt 1
data            (dNr(:,1) = (/4.d0*xsi(1,1)-1.d0, 0.d0, 0.d0, 4.d0*xsi(1,2), 0.d0, 4.d0*xsi(1,3)/))
data            (dNr(:,2) = (/0.d0, 4.d0*xsi(1,2)-1.d0, 0.d0, 4.d0*xsi(1,1), 4.d0*xsi(1,3), 0.d0/))
data            (dNr(:,3) = (/0.d0, 0.d0, 4.d0*xsi(1,3)-1.d0, 0.d0, 4.d0*xsi(1,2), 4.d0*xsi(1,1)/))
! Gauss pt 2
data            (dNr(:,4) = (/4.d0*xsi(2,1)-1.d0, 0.d0, 0.d0, 4.d0*xsi(2,2), 0.d0, 4.d0*xsi(2,3)/))
data            (dNr(:,5) = (/0.d0, 4.d0*xsi(2,2)-1.d0, 0.d0, 4.d0*xsi(2,1), 4.d0*xsi(2,3), 0.d0/))
data            (dNr(:,6) = (/0.d0, 0.d0, 4.d0*xsi(2,3)-1.d0, 0.d0, 4.d0*xsi(2,2), 4.d0*xsi(2,1)/))
! Gauss pt 3
data            (dNr(:,7) = (/4.d0*xsi(3,1)-1.d0, 0.d0, 0.d0, 4.d0*xsi(3,2), 0.d0, 4.d0*xsi(3,3)/))
data            (dNr(:,8) = (/0.d0, 4.d0*xsi(3,2)-1.d0, 0.d0, 4.d0*xsi(3,1), 4.d0*xsi(3,3), 0.d0/))
data            (dNr(:,9) = (/0.d0, 0.d0, 4.d0*xsi(3,3)-1.d0, 0.d0, 4.d0*xsi(3,2), 4.d0*xsi(3,1)/))
! Loop: do i=1,9 - index (3*i-2,3*i-1,3*i) in accessing dNr(:,i)
private dNr

contains 
subroutine size_check()
    !print * , xsi(:,1)
    print * , ""
    !print * , xsi(:,2)
    print * , ""
    !print * , xsi(:,3)
    print * , "Rank check: ", rank(P0)
    return
end subroutine size_check

! Deformation gradient
subroutine c2tl6_e(ke,coord,t,D,ed,es)
    double precision :: ke(:,:), coord(6,2), D(:,:,:)
    double precision :: ed(:), es(:,:)
    double precision :: t
    
    double precision :: JT(3,3), JTinv(3,3)
    double precision :: gradxy(2,6)
    double precision :: B0(3,12)
    double precision :: A(3,4)
    double precision :: A_temp(4)
    double precision :: Bl0(3,12)
    double precision :: H0(4,12)
    double precision :: R(4,4) 
    double precision :: dNx(2,6)
    double precision :: stress(2,2)
    double precision :: detJ
    integer          :: gp
    !
    JT(:,1)   = (/1.d0,1.d0,1.d0/)
    
    ke = 0.d0
    do gp = 1,ngp
        JT(:,2:3) = MATMUL( transpose(dNr(:,index(gp,:)) ),coord)
        call inv3(JTinv,JT)

        detJ      = det3(JT)

        dNx       = MATMUL( MATMUL(P0,JTinv), transpose(dNr(:,index(gp,:)) ) )

        Bl0           = 0.d0
        Bl0(1,1:12:2) = dNx(1,:)
        Bl0(2,2:12:2) = dNx(2,:)
        Bl0(3,1:12:2) = dNx(2,:)
        Bl0(3,2:12:2) = dNx(1,:)

        H0            = 0.d0
        H0(1,1:12:2)  = dNx(1,:)
        H0(2,1:12:2)  = dNx(2,:)
        H0(3,2:12:2)  = dNx(1,:)
        H0(4,2:12:2)  = dNx(2,:)

        A_temp        = matmul(H0,ed)
        A(1,:)        = (/A_temp(1), 0.d0, A_temp(3), 0.d0/)
        A(2,:)        = (/0.d0, A_temp(2), 0.d0, A_temp(4)/)
        A(3,:)        = (/A_temp(2), A_temp(1), A_temp(4), A_temp(3)/)
        B0            = 0.d0
        B0            = Bl0 + matmul(A,H0)
        
        stress(1,:)   = (/ es(1,gp), es(3,gp) /)
        stress(2,:)   = (/ es(3,gp), es(2,gp) /)

        R             = 0.d0
        R(1:2,1:2)    = stress
        R(3:4,3:4)    = stress

        Ke = Ke + ( MATMUL( transpose(B0), MATMUL( D(:,:,gp), B0) ) &
                + MATMUL( MATMUL(transpose(H0),R),H0) ) * detJ * t * w(gp)/2

    end do
        !print * , "R: " , R
        !print * , " "
end subroutine c2tl6_e

subroutine c2tl6_f(ef,coord,t,es,ed)
    double precision :: ef(:), coord(6,2)
    double precision :: es(:,:), ed(:)
    double precision :: t
    
    double precision :: JT(3,3), JTinv(3,3)
    double precision :: gradxy(2,6)
    double precision :: Bl0(3,12)
    double precision :: B0(3,12)
    double precision :: A(3,4)
    double precision :: A_temp(4)
    double precision :: H0(4,12)
    double precision :: S(3) 
    double precision :: dNx(2,6)
    double precision :: detJ
    integer          :: gp


    JT(:,1)   = (/1.d0,1.d0,1.d0/)
    ef        = 0.d0
    do gp = 1,ngp
        !
        JT(:,2:3) = MATMUL( transpose(dNr(:,index(gp,:)) ),coord)
        call inv3(JTinv,JT)
        !
        detJ          = det3(JT)
        !
        dNx           = MATMUL( MATMUL(P0,JTinv), transpose(dNr(:,index(gp,:)) ) )
        !
        Bl0           = 0.d0
        Bl0(1,1:12:2) = dNx(1,:)
        Bl0(2,2:12:2) = dNx(2,:)
        Bl0(3,1:12:2) = dNx(2,:)
        Bl0(3,2:12:2) = dNx(1,:)
        !
        H0            = 0.d0
        H0(1,1:12:2)  = dNx(1,:)
        H0(2,1:12:2)  = dNx(2,:)
        H0(3,2:12:2)  = dNx(1,:)
        H0(4,2:12:2)  = dNx(2,:)
        !
        A_temp        = matmul(H0,ed)
        A(1,:)        = (/A_temp(1), 0.d0, A_temp(3), 0.d0/)
        A(2,:)        = (/0.d0, A_temp(2), 0.d0, A_temp(4)/)
        A(3,:)        = (/A_temp(2), A_temp(1), A_temp(4), A_temp(3)/)
        !
        B0            = 0.d0
        B0            = Bl0 + matmul(A,H0)
        !
        S             = (/ es(1,gp), es(2,gp), es(3,gp) /)
        !
        ef            = ef + matmul(transpose(B0),S)*detJ*t*w(gp)/2
    end do
    !
end subroutine c2tl6_f

subroutine c2tl6_d(defgrad,ed,coord)
    double precision :: defgrad(3,2,2), ed(:)    
    integer          :: gp
    double precision :: coord(6,2)
    
    double precision :: JT(3,3), JTinv(3,3)
    double precision :: H0(4,12)
    double precision :: dNx(2,6)
    double precision :: detJ
    double precision :: temp(4)
    
    JT(:,1)         = (/1.d0,1.d0,1.d0/)
    defgrad(:,:,:)  = 0.d0
    do gp = 1,ngp
        JT(:,2:3) = MATMUL( transpose(dNr(:,index(gp,:)) ),coord)
        !print * , JT
        !print * , " "
        call inv3(JTinv,JT)
    
        detJ            = det3(JT)
        
        dNx             = MATMUL( MATMUL(P0,JTinv), transpose(dNr(:,index(gp,:)) ) )
        
        H0              = 0.d0
        H0(1,1:12:2)    = dNx(1,:)
        H0(2,1:12:2)    = dNx(2,:)
        H0(3,2:12:2)    = dNx(1,:)
        H0(4,2:12:2)    = dNx(2,:)
        temp            = MATMUL(H0,ed) 
        defgrad(gp,1,1) = temp(1) + 1.d0
        defgrad(gp,1,2) = temp(2)
        defgrad(gp,2,1) = temp(3)
        defgrad(gp,2,2) = temp(4) + 1.d0
    end do
end subroutine c2tl6_d

end module elem_large_cont_2d6