module mater_polyfcc

! Many of the routines use explicit shapes of arrays
! This should not be used as it require the compiler to
! allocate a new variable in the subroutine, instead
! use eg.  double precision :: F(:,:) in the declaration
! for the subroutine call

! Last modified
! -------------
! M. Ristinmaa 2011-05-15
!    - Initial implementation
!
! M. Ristinmaa 2011-05-17
!    - Included all polycrystal files in module
!
! H. Hallberg 2011-06-17
!    - Tidied up the code and added comments
!    - Added the slip resistance gg0 as a material parameter in mp(17)
!    - Added the constant parameter deg2rad to convert from degrees to radians
!    - Added the constant parameters pihalf, sqrt2inv, sqrt3inv
!    - Removed several "crap" variables
!    - Changed the subroutine "norm" such that the input parameter n is no longer needed and should be removed
!
! Y. Mellbin Later
!    - Tidied furter
!    - Removed kinematic hardening
!    - Removed temperature dependence
!    - Added recrystalization

use mater_large
use matrix_util
use fem_system, only: solveq

    implicit none

    type grains
        integer :: number_gra, curr_array_length, steps
        double precision, dimension(:,:),   pointer :: sizes
        double precision, dimension(:,:,:), pointer :: Fps
        double precision, dimension(:,:,:), pointer :: gs
        double precision, dimension(:,:,:), pointer :: Psis
        double precision, dimension(:,:),   pointer :: pa
    end type grains

    ! Define constant parameters
    ! --------------------------
    double precision, parameter :: onethird = 0.333333333333333D0 ! 1/3
    double precision, parameter :: twothird = 0.666666666666667D0 ! 2/3
    double precision, parameter :: onehalf  = 0.5D0               ! 1/2
    double precision, parameter :: onesixth = 0.166666666666667D0 ! 1/6
    double precision, parameter :: onepi    = 3.141592653589793D0 ! Pi
    double precision, parameter :: twopi    = 6.283185307179586D0 ! 2*Pi
    double precision, parameter :: pihalf   = 1.570796326794897D0 ! Pi/2
    double precision, parameter :: deg2rad  = 0.017453292519943D0 ! Pi/180
    double precision, parameter :: sqrt2inv = 0.707106781186547D0 ! 1/sqrt(2)
    double precision, parameter :: sqrt3inv = 0.577350269189626D0 ! 1/sqrt(3)

    ! Number of slip systems
    ! ----------------------
    integer, parameter            :: nsys=12
    private nsys

    ! Material parameters
    ! -------------------
    double precision :: mu, kappa, G0, gamma0, mm,qp, HH, RR, gg0
    double precision :: QQ, BB, theta0, p, cc, rho, alpha, xsi
    private mu, kappa, G0, gamma0, mm,qp, HH, RR
    private QQ, BB, theta0, p, cc, rho, alpha, xsi

    ! Problem size
    ! ------------
    integer :: number_ele, number_gau
    private number_ele, number_gau

	! Slip planes and slip directions
    double precision, dimension(nsys,3,1) :: Nr,Mr
    private Nr,Mr

    ! Internal variables
    ! ------------------
    integer                       :: new, old
    type(grains), allocatable     :: gp_store(:,:)
    private new, old

    !Stiffness matrix
    double precision, allocatable :: Dmatrix(:,:,:,:)

    double precision :: angle(1000)

!    integer, parameter            :: debugg=1
    integer                       :: debugg
    integer                       :: maxstep, debugg_stress

contains

! ==============================================================================

subroutine polycryfcc_init(mp,nelm,ngp,ngr)

    implicit none

    double precision :: mp(:)
    integer          :: nelm, ngp, ngr

    double precision :: R(3,3), R1(3,3), R2(3,3)
    double precision :: ty, phi1, phi2, phi3, sizesum
    integer          :: ierr, elnr, gp, gr
    integer          :: seed(1), i

    double precision :: phi(3), xsi, psi

   number_ele = nelm ! Number of elements
    number_gau = ngp  ! Number of Gauss points

    ! Material constants
    ! ==================
    mu      = mp(1)
    kappa   = mp(2)
    G0      = mp(3)
    gamma0  = mp(4)
    mm      = mp(5)
    qp      = mp(6)
    HH      = mp(7)
    RR      = mp(8)
    QQ      = mp(9)
    BB      = mp(10)
    theta0  = mp(11)
    p       = mp(12)
    cc      = mp(13)
    rho     = mp(14)
    alpha   = mp(15)
    xsi     = mp(16)
    gg0     = mp(17)


    ! Slip plane normals and slip directions, cf. Anand (2004)
    ! --------------------------------------------------------
    Nr(1,:,1)  = sqrt3inv*(/1D0,1D0,1D0/);    Mr(1,:,1)  = sqrt2inv*(/1D0,-1D0,0D0/)
    Nr(2,:,1)  = sqrt3inv*(/1D0,1D0,1D0/);    Mr(2,:,1)  = sqrt2inv*(/-1D0,0D0,1D0/)
    Nr(3,:,1)  = sqrt3inv*(/1D0,1D0,1D0/);    Mr(3,:,1)  = sqrt2inv*(/0D0,1D0,-1D0/)
    Nr(4,:,1)  = sqrt3inv*(/-1D0,1D0,1D0/);   Mr(4,:,1)  = sqrt2inv*(/1D0,0D0,1D0/)
    Nr(5,:,1)  = sqrt3inv*(/-1D0,1D0,1D0/);   Mr(5,:,1)  = sqrt2inv*(/-1D0,-1D0,0D0/)
    Nr(6,:,1)  = sqrt3inv*(/-1D0,1D0,1D0/);   Mr(6,:,1)  = sqrt2inv*(/0D0,1D0,-1D0/)
    Nr(7,:,1)  = sqrt3inv*(/1D0,-1D0,1D0/);   Mr(7,:,1)  = sqrt2inv*(/-1D0,0D0,1D0/)
    Nr(8,:,1)  = sqrt3inv*(/1D0,-1D0,1D0/);   Mr(8,:,1)  = sqrt2inv*(/0D0,-1D0,-1D0/)
    Nr(9,:,1)  = sqrt3inv*(/1D0,-1D0,1D0/);   Mr(9,:,1)  = sqrt2inv*(/1D0,1D0,0D0/)
    Nr(10,:,1) = sqrt3inv*(/-1D0,-1D0,1D0/);  Mr(10,:,1) = sqrt2inv*(/-1D0,1D0,0D0/)
    Nr(11,:,1) = sqrt3inv*(/-1D0,-1D0,1D0/);  Mr(11,:,1) = sqrt2inv*(/1D0,0D0,1D0/)
    Nr(12,:,1) = sqrt3inv*(/-1D0,-1D0,1D0/);  Mr(12,:,1) = sqrt2inv*(/0D0,-1D0,-1D0/)

    ! Allocate memory for internal variables.
    ! The last entry is used for keeping track of
    ! if it new or old vaiables that are stored
    ! ===========================================
    new = 1
    old = 2

    allocate(gp_store(number_gau,number_ele),stat=ierr)

    do elnr=1,nelm
        do gp=1,ngp
            gp_store(gp,elnr)%number_gra = ngr
            gp_store(gp,elnr)%curr_array_length = ngr*2
            gp_store(gp,elnr)%steps = 2
            allocate(gp_store(gp,elnr)%sizes(ngr*2,2),stat=ierr)
            allocate(gp_store(gp,elnr)%Fps(9,ngr*2,2),stat=ierr)
            allocate(gp_store(gp,elnr)%gs(nsys,ngr*2,2),stat=ierr)
            allocate(gp_store(gp,elnr)%Psis(6,ngr*2,2),stat=ierr)
        end do
    end do

    do elnr=1,nelm
        do gp=1,ngp
            sizesum = 0D0
            do gr=1,ngr
                call random_number(ty)
                gp_store(gp,elnr)%sizes(gr,1) = ty + 0.5
                !size_store(gr, gp, elnr,1) = 1.0
                sizesum = sizesum + gp_store(gp,elnr)%sizes(gr,1)
            end do
            do gr=1,ngr
                gp_store(gp,elnr)%sizes(gr,1) = gp_store(gp,elnr)%sizes(gr,1)/sizesum
                gp_store(gp,elnr)%sizes(gr,2) = gp_store(gp,elnr)%sizes(gr,1)
            end do
        end do
    end do

    ! Allocate memory for initial orientation distribution, slip plane normals and slip directions
    ! --------------------------------------------------------------------------------------------


    open(unit=12,access='sequential',form='unformatted', status='unknown')
    read(12)(angle(i),i=1,1000) ! Read angles from file
    close(12)

    !call euler2r(R_goss,0D0*pi/180D0,45D0*pi/180D0,0D0*pi/180D0)

    do elnr=1,nelm
        do gp=1,ngp
            do gr=1,ngr
                do i=1,nsys
                    gp_store(gp,elnr)%gs(i,gr,old) = gg0
                end do
                gp_store(gp,elnr)%Fps((/1,5,9/),gr,old) = 1D0
                gp_store(gp,elnr)%Fps((/2,3,4,6,7,8/),gr,old) = 0D0

                call random_number(ty);  xsi = ty*onepi
                call random_number(ty);  psi = ty*onepi
                call random_number(ty)
                phi = (/ dcos(xsi)*dsin(psi), &
                         dsin(xsi)*dsin(psi), &
                         dcos(psi) /) * angle(floor(ty*1000)+1)

                call random_number(phi1); phi1 = phi1*360D0
                call random_number(phi2); phi2 = phi2*180D0
                call random_number(phi3); phi3 = phi3*360D0
                phi = 0D0

                if (mod(gr,2).eq.0) then
                    gp_store(gp,elnr)%Psis(1,gr,old) = phi1*deg2rad         ! z
                    gp_store(gp,elnr)%Psis(2,gr,old) = phi2*deg2rad         ! x
                    gp_store(gp,elnr)%Psis(3,gr,old) = phi3*deg2rad         ! z
                else
                    gp_store(gp,elnr)%Psis(1,gr,old) = (360D0-phi1)*deg2rad ! z
                    gp_store(gp,elnr)%Psis(2,gr,old) = phi2*deg2rad         ! x
                    gp_store(gp,elnr)%Psis(3,gr,old) = (360D0-phi3)*deg2rad ! z
                end if

                gp_store(gp,elnr)%Psis(4,gr,old) = phi(1)
                gp_store(gp,elnr)%Psis(5,gr,old) = phi(2)
                gp_store(gp,elnr)%Psis(6,gr,old) = phi(3)

            end do


            if (ngr.eq.1) then
                gp_store(gp,elnr)%Psis(1,1,old) = 0
                gp_store(gp,elnr)%Psis(2,1,old) = 0
                gp_store(gp,elnr)%Psis(3,1,old) = 0

                phi(1) =  onehalf
                phi(2) = -onehalf-sqrt2inv
                phi(3) =  onehalf+sqrt2inv

                ty = -dacos(sqrt2inv-1d0/4d0)
                ty = ty/(2d0*dsin(ty))


                phi(1) = 0
                phi(2) = -pihalf/2d0
                phi(3) = 0
                call rotmat(R1,Phi)
                phi(2) = 0
                phi(3) = pihalf/2d0
                call rotmat(R2,Phi)
                R=matmul(R2, R1)
                call r2euler(phi1,phi2,phi3,R)
                gp_store(gp,elnr)%Psis(1,1,old) = phi1
                gp_store(gp,elnr)%Psis(2,1,old) = phi2
                gp_store(gp,elnr)%Psis(3,1,old) = phi3

                gp_store(gp,elnr)%Psis(4,1,old) = 0
                gp_store(gp,elnr)%Psis(5,1,old) = 0
                gp_store(gp,elnr)%Psis(6,1,old) = 0
            end if


            gp_store(gp,elnr)%Fps(:,:,new)      = gp_store(gp,elnr)%Fps(:,:,old)
            gp_store(gp,elnr)%gs(:,:,new)       = gp_store(gp,elnr)%gs(:,:,old)
            gp_store(gp,elnr)%Psis(:,:,new)     = gp_store(gp,elnr)%Psis(:,:,old)


        end do
    end do


! Stiffness matrix
   allocate (Dmatrix(6,6,number_gau,number_ele),stat=ierr)


    return

end subroutine polycryfcc_init

! ==============================================================================

subroutine polycryfcc_accept()

    implicit none

    integer          :: elnr, gp

    do elnr=1,number_ele
        do gp=1,number_gau
            gp_store(gp,elnr)%sizes(:,old)      = gp_store(gp,elnr)%sizes(:,new)
            gp_store(gp,elnr)%Fps(:,:,old)      = gp_store(gp,elnr)%Fps(:,:,new)
            gp_store(gp,elnr)%gs(:,:,old)       = gp_store(gp,elnr)%gs(:,:,new)
            gp_store(gp,elnr)%Psis(:,:,old)     = gp_store(gp,elnr)%Psis(:,:,new)
        end do
    end do
    if (new.eq.1) then
        new = 2
        old = 1
    else
        old = 2
        new = 1
    end if

    return

end subroutine polycryfcc_accept

! ================================================================================================================

subroutine polycryfcc(stype,se,elnr,dg,dgn,theta,dt, firstrun)
! ---------------------------------------------------------------------------------
!
! Purpose:
!    Calculate the updated stress tensor and stiffness tensor components for a
!    polycrystalline FCC material.
!
! Input:
!    stype       Stress type (currently 2nd Piola-Kirchhoff)
!    elnr        Element number
!    dg          Deformation gradient (previous) [9 x ngp x nel]
!    dgn         Deformation gradient (current)  [9 x ngp x nel]
!    theta       Temperature    THIS CAN POSSIBLY BE REMOVED IN AN ADIABATIC ANALYSIS
!    dt          Time increment
!
! Output:
!    se          Stress tensor components
!    Dm          Stiffness tensor
!
!----------------------------------------------------------------------------------

    implicit none

    character(len=*) :: stype(*)
    integer          :: elnr, ii
    double precision :: dg(:,:,:), dgn(:,:,:), se(:,:)
    double precision :: theta(:,:), dt
    double precision :: stress(8,3,3), F(8,3,3), Fold(8,3,3), Dmx(8,6,6)
    logical firstrun

    maxstep=0

    do ii=1,8
        stress(ii,1,:)=(/se(1,ii),se(4,ii),se(5,ii) /)
        stress(ii,2,:)=(/se(4,ii),se(2,ii),se(6,ii) /)
        stress(ii,3,:)=(/se(5,ii),se(6,ii),se(3,ii) /)
        F(ii,1,:)=(/dg(1,ii,elnr),dg(2,ii,elnr),dg(3,ii,elnr) /)
        F(ii,2,:)=(/dg(4,ii,elnr),dg(5,ii,elnr),dg(6,ii,elnr) /)
        F(ii,3,:)=(/dg(7,ii,elnr),dg(8,ii,elnr),dg(9,ii,elnr) /)
        Fold(ii,1,:)=(/dgn(1,ii,elnr),dgn(2,ii,elnr),dgn(3,ii,elnr) /)
        Fold(ii,2,:)=(/dgn(4,ii,elnr),dgn(5,ii,elnr),dgn(6,ii,elnr) /)
        Fold(ii,3,:)=(/dgn(7,ii,elnr),dgn(8,ii,elnr),dgn(9,ii,elnr) /)
    end do


	!call poly_upd_euler(Dmx,stress,F,Fold,dt,elnr)
    !call poly_upd(Dmx,stress,F,Fold,dt,elnr)
    !call poly_upd_gpu(Dmx,stress,F,Fold,dt,elnr)
    call poly_upd_gpu_stepsize(Dmx,stress,F,Fold,dt,elnr, firstrun)
    !call poly_upd_gpu_e(Dmx,stress,F,Fold,dt,elnr)

    if (debugg_stress.eq.1) write(*,*)'poly stress ',stress
    do ii=1,8
        se(:,ii)=(/stress(ii,1,1), stress(ii,2,2), stress(ii,3,3), &
                   stress(ii,1,2), stress(ii,1,3), stress(ii,2,3) /)
        Dmatrix(:,:,ii,elnr) = Dmx(ii,:,:)
    end do

!    write(*,"(A,I4)")'max number of steps in RK ',maxstep

    return

end subroutine polycryfcc

! ================================================================================================================

subroutine dpolycryfcc(D,elnr)
  implicit none
  double precision     :: D(:,:,:)
  integer              :: elnr

  D=Dmatrix(:,:,:,elnr)

end subroutine dpolycryfcc

! ==============================================================================

subroutine euler2r(R,phi1,theta,phi2)
! ---------------------------------------------------------------------------------
!
! Purpose:
!    Convert the three Euler angles (phi1, theta, phi2) to a rotation matrix R.
!
! Input:
!    phi1        Euler angle 1 [double precision scalar]
!    theta       Euler angle 2 [double precision scalar]
!    phi2        Euler angle 3 [double precision scalar]
!
! Output:
!    R           Rotation matrix [3 x 3]
!
!----------------------------------------------------------------------------------

    implicit none

    double precision :: phi1, phi2, theta
    double precision :: R(:,:)

    R(1,:)=(/ dcos(phi1)*dcos(phi2)-dsin(phi1)*dsin(phi2)*dcos(theta), &
              dsin(phi1)*dcos(phi2)+dcos(phi1)*dsin(phi2)*dcos(theta), &
              dsin(phi2)*dsin(theta) /)
    R(2,:)=(/-dcos(phi1)*dsin(phi2)-dsin(phi1)*dcos(phi2)*dcos(theta), &
             -dsin(phi1)*dsin(phi2)+dcos(phi1)*dcos(theta)*dcos(phi2), &
              dcos(phi2)*dsin(theta) /)
    R(3,:)=(/ dsin(phi1)*dsin(theta), &
             -dcos(phi1)*dsin(theta), &
              dcos(theta) /)

    return

end subroutine euler2r

! ================================================================================================================

subroutine r2euler(phi1,phi2,phi3,R)
! ---------------------------------------------------------------------------------
!
! Purpose:
!    Convert a rotation matrix R to the three Euler angles (phi1, theta, phi2).
!
! Input:
!    R           Rotation matrix [3 x 3]
!
! Output:
!    phi1        Euler angle 1 [double precision scalar]
!    theta       Euler angle 2 [double precision scalar]
!    phi2        Euler angle 3 [double precision scalar]
!
!----------------------------------------------------------------------------------

    implicit none

    double precision                 :: phi1,phi2,phi3
    double precision, dimension(3,3) :: R

    if (R(3,3).eq.0) then
        phi3 = pihalf
        phi1 = dasin(R(3,1))
        phi2 = dasin(R(1,3))
    else
        phi3 = dacos(R(3,3))
        if (dsin(phi3).eq.0D0) then
            phi1 = acos(R(1,1))*0.5D0
            phi2 = phi1
        else
            phi1 = dacos(-R(3,2)/dsin(phi3))
            phi2 = dacos(R(2,3)/dsin(phi3))
        end if
    end if

    return

end subroutine r2euler

! ================================================================================================================

subroutine rotmat(R,Phi)
! ---------------------------------------------------------------------------------
!
! Purpose:
!    Calculate the 3x3 rotation matrix R from the pseudovector Phi whos magnitude
!    is the rotation magnitude and whos direction is the axis of rotation.
!
! Input:
!    Phi         Pseudo-vector [3 x 1]
!
! Output:
!    R           Rotation matrix [3 x 3]
!
!----------------------------------------------------------------------------------

    implicit none

    double precision, intent(out)    :: R(3,3)
    double precision, intent(in)     :: Phi(3,1)
    double precision, dimension(3,1) :: n
    double precision, dimension(3,3) :: nxn, mcross, id
    double precision                 :: phi_v

    call eye(id,3)
    !phi_v = dsqrt(sum(Phi*Phi))
    call norm(Phi(:,1),phi_v,3)


    if (phi_v.eq.0D0) then
        R = id
    else
        n      = Phi/phi_v
        nxn    = matmul(n,transpose(n))
        mcross = reshape((/ 0D0,-n(3,1),n(2,1),n(3,1),0D0,-n(1,1),-n(2,1),n(1,1),0D0/),(/3,3/),ORDER=(/2,1/))
        R      = dcos(phi_v)*id + dsin(phi_v)*mcross + (1D0-dcos(phi_v))*nxn
    end if

    return

end subroutine rotmat

! ================================================================================================================

subroutine poly_upd_euler(DOUT,PK_H,F_n,F_o,dt,elnr)

    implicit none

    integer                                          :: ngr,a,elnr,gr,gp
    double precision                                 :: phi1,phi2,theta,dt,ty

    double precision, dimension(nsys,3,1)            :: N,M
    double precision, dimension(6,6)                 :: D
    double precision, dimension(3,3)                 :: stress
    double precision, dimension(number_gau,6,6)      :: Dout

    double precision, dimension(9)                   :: Fp_old, Fp_new
    double precision, dimension(nsys)                :: g_old, g_new, tau_old, tau_new
    double precision, dimension(number_gau,3,3)      :: PK_H, F_n, F_o
    double precision, dimension(3,3)                 :: R, R_scatt, R_pref
    double precision, dimension(3,1)                 :: phi

    type(grains)     :: eggrains


    if (debugg.ge.1) write(*,*)'start poly_upd'


    PK_H    = 0D0
    Dout    = 0D0

    do gp=1,number_gau
        eggrains = gp_store(gp, elnr)
        ngr = eggrains%number_gra
        do gr=1,ngr
            !write (*,*) 'gp', gp, 'grain', gr
            phi1     = eggrains%Psis(1,gr,old)
            theta    = eggrains%Psis(2,gr,old)
            phi2     = eggrains%Psis(3,gr,old)
            phi(1,1) = eggrains%Psis(4,gr,old)
            phi(2,1) = eggrains%Psis(5,gr,old)
            phi(3,1) = eggrains%Psis(6,gr,old)

            call euler2r(R_pref,phi1,theta,phi2)

            call rotmat(R_scatt,Phi)
            R = matmul(R_scatt,transpose(R_pref))

            do a=1,nsys
                N(a,:,:) = matmul(R,Nr(a,:,:))
                M(a,:,:) = matmul(R,Mr(a,:,:))
            end do

            Fp_old(:)    = eggrains%Fps(:,gr,old)
            g_old(:)     = eggrains%gs(:,gr,old)
            tau_old(:)   = 0d0!eggrains%taus(:,gr,old)

if (debugg_stress.eq.1) write(*,*)'ploy up gauss point',gp
if ((debugg_stress.eq.1).and.(gp.gt.2)) debugg=12

            call update(D,stress, Fp_new, g_new, tau_new,F_n(gp,:,:),F_o(gp,:,:),Fp_old,g_old, tau_old,N,M,nsys,dt)

if (debugg_stress.eq.1) write(*,*)'ploy up ',gp, stress

            eggrains%Fps(:,gr,new)    = Fp_new(:)
            eggrains%gs(:,gr,new)     = g_new(:)
            !eggrains%taus(:,gr,new)   = tau_new(:)
            eggrains%Psis(1,gr,new)   = phi1
            eggrains%Psis(2,gr,new)   = theta
            eggrains%Psis(3,gr,new)   = phi2
            eggrains%Psis(4,gr,new)   = phi(1,1)
            eggrains%Psis(5,gr,new)   = phi(2,1)
            eggrains%Psis(6,gr,new)   = phi(3,1)

            PK_H(gp,:,:) = PK_H(gp,:,:) + stress(:,:)*eggrains%sizes(gr, old)
            Dout(gp,:,:) = Dout(gp,:,:) + D(:,:)*eggrains%sizes(gr, old)

        end do
    end do

    if (debugg.ge.1) write(*,*)'end poly_upd'

    return

end subroutine poly_upd_euler

! ==============================================================================

subroutine update(Dout,PK, Fpv_new, g_new, tau_new,F_n,F_o,Fpv_old,g_old, tau_old,N,M,nsys,dt)

    implicit none

    double precision                           :: Dout(:,:), PK(:,:)
    double precision                           :: Fpv_old(:), Fpv_new(:), g_old(:), g_new(:), tau_old(:), tau_new(:)
    double precision                           :: F_n(:,:), F_o(:,:)
    double precision                           :: N(:,:,:), M(:,:,:)
    integer                                    :: nsys
    double precision                           :: dt




    integer                                    :: a, b, iter
    double precision                           :: tol, J, norm_dy, alpha, dgamma0
    double precision, dimension(3,3)           :: Fp, IFp_o, F_new, F_new_i, Fp_old, Fe_trial, Fp_new, IFp, F_old, F_old_i, Fe
    double precision, dimension(3,3)           :: Ce_trial_hat, Ce_i, Ap, ICe, LP, Ce
    double precision, dimension(3,3)           :: kr, tmp33
    double precision, dimension(nsys)          :: dg, g, dgamma, Gs, tau
    double precision, dimension(nsys*2,nsys*2) :: Jac, Jac_inv
    double precision, dimension(nsys,nsys)     :: q
    double precision, dimension(nsys*2)        :: Y, Res, DY

    ! ATS quantities
    double precision                           :: dRadtau, dRidtau, a1, a2, a3, trCe
    double precision, dimension(3,3)           :: C, IC, Se, Ce_trial, Apinv, Fe_new
    double precision, dimension(6,6)           :: diag_unit, DSeDCe, push_fw, DCe_trial_hat_DC, DRCeDCe_trial_hat
    double precision, dimension(6,6)           :: DCe_DCe_trial_hat, DCe_DC, tmp66, tmp66b, Kr_dyad_Kr
    double precision, dimension(6,6)           :: CeI_KR, DCe_hat_DCe
    double precision, dimension(6,1)           :: tmp61, Ce_M, C_IM, Ce_trialM, CE_IM, KR_M
    double precision, dimension(9,1)           :: tmp91
    double precision, dimension(1,6)           :: tmp16, Dgamma_DCe_trial_hat, tmp16b, Dg_DCe_trial_hat
    double precision, dimension(9,6)           :: DApinv_DC
    double precision, dimension(6,9)           :: BC
    double precision, dimension(nsys,6)        :: dtaudCe_trial
    double precision, dimension(nsys*2,6)      :: DRDCe_trial_hat, DYDCe_trial
    double precision, dimension(nsys,3,3)      :: dApIdgamma

    tol     = 1D-12

    call eye(kr,3)

    dgamma0 = gamma0*dt

    F_new=F_n(:,:)
    J = det3(F_new)
    F_new_i=J**(-onethird)*F_new
    F_old=F_o(:,:)
    F_old_i=(det3(F_old))**(-onethird)*F_old
    Fp_old=reshape((/Fpv_old(1:9)/),(/3,3/),order=(/2,1/))

    call inv3(IFp_o,Fp_old)

    Fe_trial     = matmul(F_new_i,IFp_o)
    Ce_trial_hat = matmul(transpose(Fe_trial),Fe_trial)

    q=qp
    do a=1,nsys
        q(a,a)     = 1D0
        g(a)       = g_old(a)
    end do
    Gs    = 0D0
    do a=1,nsys
        do b=1,nsys
            Gs(a)=Gs(a)+QQ*q(a,b)*(g(b))
        end do
    end do


    tau = 0D0
    do a=1,nsys
        !tau(a) = tau_old(a)
    end do
    call inv3(IFp,Fp_old)
    Fe   = matmul(F_old_i,IFp)
    Ce_i = matmul(transpose(Fe),Fe)
    do a=1,nsys
        tau(a) = mu*dot_product(matmul(M(a,:,1), Ce_i),N(a,:,1))
    end do





    do a=1,nsys
        Y(a*2-1) = dgamma0*(dabs(tau(a))/(G0+Gs(a)))**mm*sign(1D0,tau(a))
        Y(a*2)   = dabs(tau(a))/(G0+Gs(a))*(1D0-BB*g(a))*dabs(Y(a*2-1))
    end do

    norm_DY = 99D0
    iter    = 0
    DY      = 0D0
    do while ((norm_DY.gt.tol).and.(iter.lt.20))
        iter = iter + 1
        Y    = Y-DY
        call jacobian(RES,Jac,Y,Ce_trial_hat,N,M,g,nsys,dt,LP,Gs,tau,Apinv,dApIdgamma,J)

        dy=res
        call solveq(Jac,dy)
        call norm(res,norm_DY,nsys*2)
        if (iter.ge.60) then
            write(*,*)'Problem in subroutine "update". Program execution stopped!'
            stop
        end if

    end do

    call jacobian(RES,Jac,Y,Ce_trial_hat,N,M,g,nsys,dt,LP,Gs,tau,Apinv,dApIdgamma,J)
    call inverse(Jac,Jac_inv,nsys*2)




    do a=1,nsys
        dgamma(a) = Y(a*2-1)
        dg(a)     = Y(a*2)
    end do

    call inv3(AP,APinv)
    Fp_new = matmul(Ap,Fp_old)
    call inv3(IFp,Fp_new)
    Fe   = matmul(F_new_i,IFp)
    Ce_i = matmul(transpose(Fe),Fe)

    Ce   = (J**twothird) * Ce_i
    call inv3(ICe,Ce)

    Fpv_new(1:3) = Fp_new(1,1:3)
    Fpv_new(4:6) = Fp_new(2,1:3)
    Fpv_new(7:9) = Fp_new(3,1:3)

    do a=1,nsys
        g_new(a) = g(a) + Dg(a)
        tau_new(a) = tau(a)
    end do


    trCe = Ce(1,1)+Ce(2,2)+Ce(3,3)
    PK(:,:) = matmul(matmul(IFp,(kappa/2.0*(J**2-1)*ICe+mu*J**(-twothird)*(kr-trCe/3.0*ICe))),transpose(IFp))

!!!!!!!!!!!!!! ATS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    C = matmul(transpose(F_new),F_new)
    call inv3(IC,C)

    !SE = kappa*Dlog(J)*ICe + mu*MATMUL(ICe, logCe_i)
    Se = (kappa/2.0*(J**2-1)*ICe+mu*J**(-twothird)*(kr-trCe/3.0*ICe))

    CALL EYE(diag_unit,6)
    diag_unit(4,4)=.5D0
    diag_unit(5,5)=.5D0
    diag_unit(6,6)=.5D0
    Kr_dyad_Kr=0D0
    Kr_dyad_Kr(1:3,1:3)=1D0

    Fe_new = (J**onethird)*Fe


    Ce_IM  =RESHAPE((/ICe(1,1), ICe(2,2), ICe(3,3), ICe(1,2), ICe(1,3), ICe(2,3)/),(/6,1/))
    Ce_M   =RESHAPE((/ Ce(1,1),  Ce(2,2),  Ce(3,3),  Ce(1,2),  Ce(1,3),  Ce(2,3)/),(/6,1/))
    C_IM   =RESHAPE((/ IC(1,1),  IC(2,2),  IC(3,3),  IC(1,2),  IC(1,3),  IC(2,3)/),(/6,1/))
    Kr_M   =RESHAPE((/ 1,  1,  1,  0,  0,  0/),(/6,1/))



    CALL BAR_DYAD(CeI_Kr,ICe,Kr)

    DCe_hat_DCe=J**(-twothird)*(diag_unit-onethird*MATMUL(Ce_M,TRANSPOSE(Ce_IM)))





   call BAR_DYAD(tmp66,ICe, ICe)
    tmp66(:,4:6) = onehalf*tmp66(:,4:6)
    tmp66b = matmul(Kr_M,transpose(Ce_IM))+matmul(Ce_IM,transpose(Kr_M))

    a1 = (kappa/2.0)*(J**2.0)+(mu/9.0)*(J**(-twothird))*trCe
    a2 = (mu/3.0)*(J**(-twothird))
    a3 = (-kappa/4.0)*((J**2.0)-1.0)+(mu/6.0)*(J**(-twothird))*trCe

    DSeDCe = a1*matmul(Ce_IM,transpose(Ce_IM)) - a2*tmp66b + 2.0*a3*tmp66





!!!!!!!!! DCe_trial_hat_DC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call BAR_DYAD(push_fw, IFp, IFp)

    Ce_trial  = (J**twothird)*Ce_trial_hat
    Ce_trialM = reshape((/Ce_trial(1,1), Ce_trial(2,2), Ce_trial(3,3),Ce_trial(1,2), Ce_trial(1,3), Ce_trial(2,3)/),(/6,1/))

    call BAR_DYAD(tmp66,transpose(iFP_o),transpose(iFP_o))
    tmp66(:,4:6) = onehalf*tmp66(:,4:6)

    DCe_trial_hat_DC = (J**(-twothird))*(tmp66 - onethird*matmul(Ce_trialM, transpose(C_IM)))
    DCe_trial_hat_DC(4:6,:) = 2D0*DCe_trial_hat_DC(4:6,:)

!!!!!!!!!! DtauDCe_trial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call BAR_DYAD(DRCeDCe_trial_hat,transpose(Apinv),transpose(Apinv))
    DRCeDCe_trial_hat(:,4:6) = onehalf*DRCeDCe_trial_hat(:,4:6)

    do a=1,nsys
        tmp16b(1,:) = (/ N(a,1,1)*M(a,1,1),N(a,2,1)*M(a,2,1),N(a,3,1)*M(a,3,1),N(a,1,1)*M(a,2,1)+N(a,2,1)*M(a,1,1),&
                         N(a,1,1)*M(a,3,1)+N(a,3,1)*M(a,1,1),N(a,2,1)*M(a,3,1)+N(a,3,1)*M(a,2,1) /)
        tmp16      = mu*tmp16b          !dtaudCe
        tmp16b = matmul(tmp16,DRCeDCe_trial_hat )
        dtaudCe_trial(a,:)=(/tmp16b(1,1),tmp16b(1,2),tmp16b(1,3),&
                              tmp16b(1,4),tmp16b(1,5),tmp16b(1,6)/)
    end do

!!!!!!!!!! DYDCe_trial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do a=1,nsys
        dRadtau   = -Dgamma0*mm*(dabs(tau(a))/(G0+Gs(a)))**(mm-1D0)*1D0/(G0+Gs(a))
        dRidtau   = -sign(1D0,tau(a))/(G0+gs(a))*(1D0-BB*(g(a)+dg(a)))*dabs(dgamma(a))
        DRDCe_trial_hat(a*2-1,:) = dRadtau*dtaudCe_trial(a,:)
        DRDCe_trial_hat(a*2,:)   = dRidtau*dtaudCe_trial(a,:)
    end do
    DYDCe_trial = -matmul(Jac_inv,DRDCe_trial_hat)

!!!!!!!!!!! DApin_DC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    DApinv_DC=0

    do a=1,nsys
        Dgamma_DCe_trial_hat(1,:) = DYDCe_trial(a*2-1,1:6)   ! [1x6]

        tmp91 = reshape((/dApIdgamma(a,1,1),dApIdgamma(a,2,2),dApIdgamma(a,3,3),&
                          dApIdgamma(a,1,2),dApIdgamma(a,2,1),dApIdgamma(a,1,3),&
                          dApIdgamma(a,3,1),dApIdgamma(a,2,3),dApIdgamma(a,3,2)/),(/9,1/))

        DApinv_DC = DApinv_DC + matmul(matmul(tmp91,Dgamma_DCe_trial_hat),DCe_trial_hat_DC)
    end do

!!!!!!!!! DCe_DC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DCe_DCe_trial_hat = DRCeDCe_trial_hat
    tmp66 = matmul(DRCeDCe_trial_hat,DCe_trial_hat_DC)
    do a=1,nsys
        Dgamma_DCe_trial_hat(1,:) = DYDCe_trial(a*2-1,1:6)   ! [1x6]

        tmp33 = matmul(transpose(Apinv),matmul(Ce_trial_hat, dApIdgamma(a,:,:)))
        tmp61 = reshape((/2D0*tmp33(1,1),2D0*tmp33(2,2),2D0*tmp33(3,3),&
                         tmp33(1,2)+tmp33(2,1),tmp33(1,3)+tmp33(3,1),tmp33(2,3)+tmp33(3,2)/),(/6,1/))
        DCe_DCe_trial_hat = DCe_DCe_trial_hat + matmul(tmp61,Dgamma_DCe_trial_hat)

    end do

    DCe_DC = J**(twothird) * matmul(DCe_DCe_trial_hat,DCe_trial_hat_DC) + onethird*MATMUL(Ce_M,TRANSPOSE(C_IM))
    DCe_DC(4:6,:) = 2D0*DCE_DC(4:6,:)

    tmp33 = matmul(Se,transpose(iFp))

    BC = reshape((/&
         2*iFP_o(1,1)*tmp33(1,1), 2*iFP_o(1,2)*tmp33(2,1),2*iFP_o(1,3)*tmp33(3,1),&
         2*iFP_o(1,1)*tmp33(2,1), 2*iFP_o(1,2)*tmp33(1,1),2*iFP_o(1,1)*tmp33(3,1),&
         2*iFP_o(1,3)*tmp33(1,1), 2*iFP_o(1,2)*tmp33(3,1),2*iFP_o(1,3)*tmp33(2,1),&
         !
         2*iFP_o(2,1)*tmp33(1,2), 2*iFP_o(2,2)*tmp33(2,2),2*iFP_o(2,3)*tmp33(3,2),&
         2*iFP_o(2,1)*tmp33(2,2), 2*iFP_o(2,2)*tmp33(1,2),2*iFP_o(2,1)*tmp33(3,2),&
         2*iFP_o(2,3)*tmp33(1,2), 2*iFP_o(2,2)*tmp33(3,2),2*iFP_o(2,3)*tmp33(2,2),&
         !
         2*iFP_o(3,1)*tmp33(1,3), 2*iFP_o(3,2)*tmp33(2,3),2*iFP_o(3,3)*tmp33(3,3),&
         2*iFP_o(3,1)*tmp33(2,3), 2*iFP_o(3,2)*tmp33(1,3),2*iFP_o(3,1)*tmp33(3,3),&
         2*iFP_o(3,3)*tmp33(1,3), 2*iFP_o(3,2)*tmp33(3,3),2*iFP_o(3,3)*tmp33(2,3),&
         !
         iFP_o(1,1)*tmp33(1,2)+iFP_o(2,1)*tmp33(1,1), iFP_o(1,2)*tmp33(2,2)+iFP_o(2,2)*tmp33(2,1),iFP_o(1,3)*tmp33(3,2)+iFP_o(2,3)*tmp33(3,1),&
         iFP_o(1,1)*tmp33(2,2)+iFP_o(2,1)*tmp33(2,1), iFP_o(1,2)*tmp33(1,2)+iFP_o(2,2)*tmp33(1,1),iFP_o(1,1)*tmp33(3,2)+iFP_o(2,1)*tmp33(3,1),&
         iFP_o(1,3)*tmp33(1,2)+iFP_o(2,3)*tmp33(1,1), iFP_o(1,2)*tmp33(3,2)+iFP_o(2,2)*tmp33(3,1),iFP_o(1,3)*tmp33(2,2)+iFP_o(2,3)*tmp33(2,1),&
         !
         iFP_o(1,1)*tmp33(1,3)+iFP_o(3,1)*tmp33(1,1), iFP_o(1,2)*tmp33(2,3)+iFP_o(3,2)*tmp33(2,1),iFP_o(1,3)*tmp33(3,3)+iFP_o(3,3)*tmp33(3,1),&
         iFP_o(1,1)*tmp33(2,3)+iFP_o(3,1)*tmp33(2,1), iFP_o(1,2)*tmp33(1,3)+iFP_o(3,2)*tmp33(1,1),iFP_o(1,1)*tmp33(3,3)+iFP_o(3,1)*tmp33(3,1),&
         iFP_o(1,3)*tmp33(1,3)+iFP_o(3,3)*tmp33(1,1), iFP_o(1,2)*tmp33(3,3)+iFP_o(3,2)*tmp33(3,1),iFP_o(1,3)*tmp33(2,3)+iFP_o(3,3)*tmp33(2,1),&
         !
         iFP_o(2,1)*tmp33(1,3)+iFP_o(3,1)*tmp33(1,2), iFP_o(2,2)*tmp33(2,3)+iFP_o(3,2)*tmp33(2,2),iFP_o(2,3)*tmp33(3,3)+iFP_o(3,3)*tmp33(3,2),&
         iFP_o(2,1)*tmp33(2,3)+iFP_o(3,1)*tmp33(2,2), iFP_o(2,2)*tmp33(1,3)+iFP_o(3,2)*tmp33(1,2),iFP_o(2,1)*tmp33(3,3)+iFP_o(3,1)*tmp33(3,2),&
         iFP_o(2,3)*tmp33(1,3)+iFP_o(3,3)*tmp33(1,2), iFP_o(2,2)*tmp33(3,3)+iFP_o(3,2)*tmp33(3,2),iFP_o(2,3)*tmp33(2,3)+iFP_o(3,3)*tmp33(2,2)&
                  /),(/6,9/), order=(/2,1/))

    DOUT(:,:) = 2D0*(matmul(matmul(push_fw,DSeDCe),DCe_DC) + matmul(BC,DApinv_DC))

    return

end subroutine update

! ============================================================================================================================

subroutine jacobian(RES,JAC,Y,Ce_trial,N,M,gin,nsys,dt,LP,G,tau,Apinv,dApIdgamma,J)

    implicit none

    integer                                    :: a,b, i,nsys,alpha,beta
    double precision                           :: h, tmp, dt, dgamma0, J
    double precision                           :: dRadtau, dRadG, dRidtau, dRidG
    double precision, dimension(nsys,nsys)     :: q, kr_ab, dtaudgamma
    double precision, dimension(nsys*2,nsys*2) :: JAC
    double precision, dimension(nsys*2)        :: Y, RES
    double precision, dimension(3,3)           :: Ce_trial, Ce_i
    double precision, dimension(3,3)           :: Apinv, Lp
    double precision, dimension(nsys,3,3)      :: dApIdgamma
    double precision, dimension(3,3)           :: TMP33, TMP33b, kr3, dapinv
    double precision, dimension(nsys)          :: gin, g, dgamma, tau, dg
    double precision, dimension(nsys,3,1)      :: N, M
    double precision, dimension(1,6)           :: tmp16, tmp16b
    double precision, dimension(6,1)           :: tmp61

    dgamma0 = gamma0*dt

    call eye(kr3,3)

    kr_ab = 0D0
    q     = qp
    do a=1,nsys
        q(a,a)     = 1D0
        dgamma(a)  = Y(a*2-1)
        dg(a)      = Y(a*2)
        kr_ab(a,a) = 1D0
    end do

    G = 0D0
    do a=1,nsys
        do b=1,nsys
            G(a) = G(a) + QQ*q(a,b)*(gin(b) + dg(b))
        end do
    end do


    LP = 0D0
    do a=1,nsys
        LP = LP + dgamma(a)*matmul(M(a,:,:),transpose(N(a,:,:)))
    end do
    call padeexp(-LP,Apinv)


    TMP33=kr3+onehalf*lp
    call inv3(TMP33b,TMP33)
    tmp33=matmul(TMP33b,(kr3-onehalf*lp))+kr3

    dapinv=0d0
    dApIdgamma=0d0

    do a=1,3
        do b=1,3
            do alpha=1,3
                do beta=1,3
                    dapinv(alpha,beta)=onehalf*TMP33b(a,alpha)*(tmp33(beta,b))
                end do
            end do
            do i=1,nsys
                dApIdgamma(i,a,b)= - dot_product(matmul(M(i,:,1),dapinv),N(i,:,1))
            end do
        end do
    end do


    Ce_i = matmul(transpose(Apinv),matmul(Ce_trial,Apinv))

    tau = 0D0
    do a=1,nsys
        tau(a) = mu*dot_product(matmul(M(a,:,1), Ce_i),N(a,:,1))
    end do

! ============ Residual ========================================================================================

    do a=1,nsys
        res(a*2-1) = dgamma(a) - dgamma0*(dabs(tau(a))/(G0+G(a)))**mm*sign(1D0,tau(a))
        res(a*2)   = dg(a) - dabs(tau(a))/(G0+G(a))*(1D0-BB*(gin(a)+dg(a)))*dabs(dgamma(a))
    end do

! ============ Jacobian ========================================================================================



    do a=1,nsys
        do b=1,nsys
            ! ---------- dtau_dCe
            tmp16b(1,:) = (/ N(a,1,1)*M(a,1,1),N(a,2,1)*M(a,2,1),N(a,3,1)*M(a,3,1),N(a,1,1)*M(a,2,1)+N(a,2,1)*M(a,1,1),&
                             N(a,1,1)*M(a,3,1)+N(a,3,1)*M(a,1,1),N(a,2,1)*M(a,3,1)+N(a,3,1)*M(a,2,1) /)
            tmp16 = mu*tmp16b
            ! ---------- dCe_ddgamma
            tmp33 = matmul(transpose(dApIdgamma(b,:,:)),matmul(Ce_trial,APinv)) + matmul(transpose(Apinv),matmul(Ce_trial,dApIdgamma(b,:,:)))

            tmp61(:,1)      = (/ tmp33(1,1),tmp33(2,2),tmp33(3,3),tmp33(1,2),tmp33(1,3),tmp33(2,3) /)
            dtaudgamma(a,b) = dot_product(tmp16(1,:),tmp61(:,1))
        end do
    end do


    Jac = 0D0
    do a=1,nsys

        dRadtau   = -Dgamma0*mm*(dabs(tau(a))/(G0+g(a)))**(mm-1D0)*1D0/(G0+g(a))
        dRadG     =  Dgamma0*mm*(dabs(tau(a))/(G0+g(a)))**(mm)*1D0/(G0+g(a))*sign(1D0,tau(a))

        dRidtau   = -sign(1D0,tau(a))/(G0+g(a))*(1D0-BB*(gin(a)+dg(a)))*dabs(dgamma(a))
        dRidG     =  dabs(tau(a))/(G0+g(a))**2D0*(1D0-BB*(gin(a)+dg(a)))*dabs(dgamma(a))

        do b=1,nsys
            ! ---------- dRa_dgamma -----------------------------------------------------------------------------------
            Jac(a*2-1,b*2-1) = kr_ab(a,b) + dRadtau*dtaudgamma(a,b)
            ! ---------- dRa_dg ---------------------------------------------------------------------------------------
            Jac(a*2-1,b*2) = dRadG*QQ*q(a,b)
            ! -------   dRi_dgamma ------------------------------------------------------------------------------------
            Jac(a*2,b*2-1) = -dabs(tau(a))/(G0+g(a))*(1D0-BB*(gin(a)+dg(a)))*sign(1D0,dgamma(b))*kr_ab(a,b)+dRidtau*dtaudgamma(a,b)
            ! -------   dRi_dg ---------------------------------------------------------------------------------------
            Jac(a*2,b*2) = kr_ab(a,b) + dabs(tau(a))/(G0+g(a))*BB*abs(dgamma(a))*kr_ab(a,b)+dRidG*QQ*q(a,b)
        end do
    end do

    return
end subroutine jacobian

! ============================================================================================================================


subroutine poly_upd(DOUT,PK_H,F_n_gp,F_o_gp,dt,elnr)

    implicit none

    integer                                          :: ngr,a,elnr,gr,gp
    double precision                                 :: phi1,phi2,theta,dt
    double precision, dimension(nsys,3,1)            :: N,M
    double precision, dimension(6,6)                 :: D
    double precision, dimension(3,3)                 :: stress
    double precision, dimension(number_gau,6,6)      :: Dout

    double precision, dimension(9)                   :: Fpv_old, Fpv_new
    double precision, dimension(nsys)                :: g_old, g_new
    double precision, dimension(number_gau,3,3)      :: PK_H, F_n_gp, F_o_gp
    double precision, dimension(3,3)   			     :: F_n, F_o
    double precision, dimension(3,3)                 :: R, R_scatt, R_pref
    double precision, dimension(3,1)                 :: phi

    type(grains)     :: eggrains

    integer                                    :: b, i, k, steps
    double precision                           :: J, dgamma0, trCe, a1, a2, a3
    double precision, dimension(nsys+9*7)      :: Y0, Yend
    double precision, dimension(3,3)           :: F_new, F_new_i, Fp_new, IFp, F_old, Fe, Se
    double precision, dimension(3,3)           :: Ce_i, ICe, Ce, C_new, C_old, kr, tmp33, tmp33b
    double precision, dimension(nsys)          :: dg, g, tau
    double precision, dimension(nsys,nsys)     :: q
    double precision, dimension(9,6)           :: dFpdC
    double precision, dimension(6,6)           :: diag_unit, DSeDCe, push_fw, DCe_DC, tmp66, tmp66b
    double precision, dimension(6,1)           :: Kr_M, Ce_IM
    double precision, dimension(6,9)           :: tmp69a, tmp69b

	integer        								  :: ierr
    double precision, dimension(:,:),     pointer :: Y0_v, Yend_v, g_v
    double precision, dimension(:,:,:,:), pointer :: N_v, M_v


    if (debugg.ge.1) write(*,*)'start poly_upd_gpu'

    PK_H    = 0D0
    Dout    = 0D0

    call eye(kr,3)
    dgamma0 = gamma0*dt
    q=qp
    do a=1,nsys
        q(a,a)     = 1D0
    end do


    do gp=1,number_gau
        eggrains = gp_store(gp, elnr)
        ngr = eggrains%number_gra
        F_n = F_n_gp(gp,:,:)
        F_o = F_o_gp(gp,:,:)

	    F_new=F_n(:,:)
	    J = det3(F_new)
	    F_new_i=J**(-onethird)*F_new
	    F_old=F_o(:,:)

	    C_old = matmul(transpose(F_old),F_old)
	    C_new = matmul(transpose(F_new),F_new)


        allocate(Y0_v(nsys+9*7,ngr),stat=ierr)
        allocate(Yend_v(nsys+9*7,ngr),stat=ierr)
        allocate(g_v(nsys,ngr),stat=ierr)
        allocate(N_v(nsys,3,1,ngr),stat=ierr)
        allocate(M_v(nsys,3,1,ngr),stat=ierr)

		Yend_v = 0D0

        do gr=1,ngr
            !write (*,*) 'gp', gp, 'grain', gr
            phi1     = eggrains%Psis(1,gr,old)
            theta    = eggrains%Psis(2,gr,old)
            phi2     = eggrains%Psis(3,gr,old)
            phi(1,1) = eggrains%Psis(4,gr,old)
            phi(2,1) = eggrains%Psis(5,gr,old)
            phi(3,1) = eggrains%Psis(6,gr,old)
            eggrains%Psis(1,gr,new)   = eggrains%Psis(1,gr,old)
            eggrains%Psis(2,gr,new)   = eggrains%Psis(2,gr,old)
            eggrains%Psis(3,gr,new)   = eggrains%Psis(3,gr,old)
            eggrains%Psis(4,gr,new)   = eggrains%Psis(4,gr,old)
            eggrains%Psis(5,gr,new)   = eggrains%Psis(5,gr,old)
            eggrains%Psis(6,gr,new)   = eggrains%Psis(6,gr,old)

            call euler2r(R_pref,phi1,theta,phi2)

            call rotmat(R_scatt,Phi)
            R = matmul(R_scatt,transpose(R_pref))

            do a=1,nsys
                N(a,:,:) = matmul(R,Nr(a,:,:))
                M(a,:,:) = matmul(R,Mr(a,:,:))
            end do

            Fpv_old(:)    = eggrains%Fps(:,gr,old)
            g_old(:)     = eggrains%gs(:,gr,old)

			if (debugg_stress.eq.1) write(*,*)'ploy up gauss point',gp
			if ((debugg_stress.eq.1).and.(gp.gt.2)) debugg=12


		    if (debugg.ge.1) write(*,*)'start update'

		    do a=1,nsys
		        g(a)       = g_old(a)
		    end do


		    Y0(1:nsys) = 0D0
		    Y0((nsys+1):(nsys+9)) = Fpv_old(:)
		    Y0((nsys+10):(nsys+7*9)) = 0D0


        	Y0_v(:,gr) = Y0
        	g_v(:,gr) = g
        	N_v(:,:,:,gr) = N
        	M_v(:,:,:,gr) = M

		end do

	    steps=40
	    if (debugg.ge.1) write(*,*)'-----------steps----',steps
	    if (steps.gt.maxstep) maxstep=steps


		call odesplit(Y0_v, Yend_v, C_old, C_new,N_v,M_v,g_v,nsys, dt,steps, ngr)


		do gr=1,ngr

			if (isnan(1D0)) then
				return
			end if


			Yend(1:(nsys+9*7)) = Yend_v(1:(nsys+9*7),gr)
        	g = g_v(:,gr)
        	N = N_v(:,:,:,gr)
        	M = M_v(:,:,:,gr)

		    dg(:)     = Yend(1:nsys)
		    Fpv_new(:) = Yend((nsys+1):(nsys+9))

		    do a=1,9
		        dFpdC(a,:) = Yend((nsys+4+6*a):(nsys+9+6*a))
		    end do


		    Fp_new=reshape((/Fpv_new(1:9)/),(/3,3/),order=(/2,1/))
		    call inv3(IFp,Fp_new)
		    Fe   = matmul(F_new_i,IFp)
		    Ce_i = matmul(transpose(Fe),Fe)

		    tau = 0D0
		    do a=1,nsys
		        tau(a) = dot_product(matmul(M(a,:,1),mu*Ce_i),N(a,:,1))
		    end do


		    Ce   = (J**twothird) * Ce_i
		    call inv3(ICe,Ce)

		    do a=1,nsys
		        g_new(a) = g(a) + dg(a)
		    end do
		    trCe = Ce(1,1)+Ce(2,2)+Ce(3,3)
		    Se = (kappa/2.0*(J**2-1)*ICe+mu*J**(-twothird)*(kr-trCe/3.0*ICe))
		    stress(:,:) = matmul(matmul(IFp,(kappa/2.0*(J**2-1)*ICe+mu*J**(-twothird)*(kr-trCe/3.0*ICe))),transpose(IFp))



			!!!!!!!!!!!!!! ATS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!! DSeDCe !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		    Ce_IM  =RESHAPE((/ICe(1,1), ICe(2,2), ICe(3,3), ICe(1,2), ICe(1,3), ICe(2,3)/),(/6,1/))
		    Kr_M   =RESHAPE((/ 1,  1,  1,  0,  0,  0/),(/6,1/))

		    call BAR_DYAD(tmp66,ICe, ICe)
		    tmp66(:,4:6) = onehalf*tmp66(:,4:6)

		    tmp66b = matmul(Kr_M,transpose(Ce_IM))+matmul(Ce_IM,transpose(Kr_M))

		    a1 = (kappa/2.0)*(J**2.0)+(mu/9.0)*(J**(-twothird))*trCe
		    a2 = (mu/3.0)*(J**(-twothird))
		    a3 = (-kappa/4.0)*((J**2.0)-1.0)+(mu/6.0)*(J**(-twothird))*trCe

		    DSeDCe = a1*matmul(Ce_IM,transpose(Ce_IM)) - a2*tmp66b + 2.0*a3*tmp66

			!!!!!!!!! DCe_DC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		    CALL EYE(diag_unit,6)
		    call BAR_DYAD(push_fw, IFp, IFp)
		    push_fw(:,4:6) = onehalf*push_fw(:,4:6)

		    tmp33 = matmul(C_new,IFp)
		    tmp33b = matmul((IFp),Se)

		    tmp69a=0D0
		    tmp69b=0d0
		    do i=1,3
		        a=mod(i,3)
		        b=mod(i+1,3)
		        do k=1,3
		            tmp69b(i,i+3*(k-1)) = tmp33(i,k)
		            tmp69b(6-b, a+1+3*(k-1)) = tmp33(i,k)
		            tmp69b(6-a, b+1+3*(k-1)) = tmp33(i,k)
		            tmp69a(i,(i-1)*3+k) = tmp33b(k,i)
		            tmp69a(6-b, a*3+k) = tmp33b(k,i)
		            tmp69a(6-a, b*3+k) = tmp33b(k,i)
		        end do
		    end do

		    DCe_DC = matmul(transpose(push_fw), (diag_unit-2*matmul(tmp69b,dFpdC)))

			!!!!!!!!!!! DSDC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		    DSeDCe(:,4:6) = 2.0*DSeDCe(:,4:6)
		    tmp66 = matmul(DSeDCe,DCe_DC)
		    tmp66(4:6,:) = 2.0*tmp66(4:6,:)
		    tmp66b = matmul(tmp69a,dFpdC)
		    D(:,:) = 2D0*(matmul(push_fw,tmp66)-2D0*matmul(push_fw,tmp66b))


		    if (debugg.ge.1) write(*,*)'end update'

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			if (debugg_stress.eq.1) write(*,*)'poly up ',gp, stress

            eggrains%Fps(:,gr,new)    = Fpv_new(:)
            eggrains%gs(:,gr,new)     = g_new(:)

            PK_H(gp,:,:) = PK_H(gp,:,:) + stress(:,:)*eggrains%sizes(gr, old)
            Dout(gp,:,:) = Dout(gp,:,:) + D(:,:)*eggrains%sizes(gr, old)

        end do

        deallocate(Y0_v)
        deallocate(Yend_v)
        deallocate(g_v)
        deallocate(N_v)
        deallocate(M_v)
    end do

    if (debugg.ge.1) write(*,*)'end poly_upd'

    return


end subroutine poly_upd

! ==============================================================================


end module mater_polyfcc
