!############################################################
! The module contains variables needed to run opt_module.f90
!############################################################
! Created: 2020-09-30
! Anna Dalklint
! Lund University
!############################################################
module opt_variables

!$ use omp_lib
use abaqus_util
use matlab_util
use sparse_util

   ! Global variables
   implicit none
   integer*8                           :: matptr
   integer, allocatable                :: enod(:,:) , bcnod(:,:), flnod(:,:), fbnod(:,:)
   double precision, allocatable       :: coord(:,:), bcval(:), flval(:)
   integer, allocatable                :: fldof(:), fbdof(:), bcdof(:)
   character(len=40)                   :: filename
   integer                             :: ierr, startH
   integer                             :: nbc, dofnod, nnod, ndof, nelm, max_iter, C, ndes, ndes_orig, ndens
   integer                             :: ngp, nbr_eig, nbr_eigMax, iwait, nmax, imax, B_iter, p_iter, nbr_bound, nbr_k, nbr_sample
   integer                             :: nbr_eig_U, nbr_eig_L, band
   double precision                    :: Arclength, th
   double precision                    :: bound_U, bound_L, min_gconnect
   integer                             :: eldof, elnod, els, eldef, elm, middle_dof, nbr_bound_vars, nbr_disp_vars, nbr_extra_vars
  
   double precision                    :: Emin, Emax, pagg_k, pagg_eig, pAggFinal_k, paggFinal_eig, stepLength
   integer                             :: L, L0, opt_finished, loadstep_bandgap

   real(16), parameter                 :: pi = 4*atan(1.0_16)
   
   double precision                    :: start, finish
   character(80)                       :: aggregation, optimizer, FE_dim, penalization, initial_guess
   character(2)                        :: Lagrangian_type
   
   double precision, allocatable       :: ed(:,:)
   double precision, allocatable       :: a(:), ahist(:,:), energyKHist(:,:), energyMHist(:,:), ares_hist(:,:)
   double precision, allocatable       :: dg(:,:,:), Fvec(:)

   double precision, allocatable       :: Dsmall(:,:)
   double precision                    :: mp(2)

   double precision                    :: v, p, peig, omega, pold, Bold, pIncr, pFinal, E1, E2, massrho1, massrho2
   double precision, allocatable       :: E(:), massRho(:)
   integer                             :: nbr_phase
   double precision                    :: Vtot, Volfrac, Ve, DensTot, Mtot
   double precision, allocatable       :: tmp2(:), f0valhist(:), gvalhist(:,:), history(:,:)
   double precision, allocatable       :: Pattern(:)
   double precision, allocatable       :: R(:,:)
   double precision                    :: pert
   integer                             :: ipert, iterIncr, save_itr, ipert_sym
   double precision                    :: TOL, designTOL, LambdaNom

   double precision                    :: Scaleg0, SafetyMarginal, LambdaMax, scaleK
   double precision, allocatable       :: EVALS_TOT(:,:), EVECS_TOT(:,:), evals_U(:,:), evals_L(:,:), evals_U_inv(:,:)
   double precision, allocatable       :: kvec_mat(:,:)

   type(sparse)                        :: K, Kstar, M, Mstar
   type(sparse)                        :: K_Re_Im, M_Re_Im, Kbloch, Mbloch
   double precision, allocatable       :: Kfull(:,:), Mfull(:,:), Kstarfull(:,:), Kfilterfull(:,:)
   double precision, allocatable       :: Kfilterstarfull(:,:), Kfilterstarsymfull(:,:), Kfiltersymfull(:,:)
   double precision, allocatable       :: K_Re_Imfull(:,:), M_Re_Imfull(:,:), Kstar_Re_Imfull(:,:), Mstar_Re_Imfull(:,:)

   double precision                    :: Z, geps, gmove_topo, gmove_bound, penalty, gmove_disp
   double precision, allocatable       :: rho(:), gmove(:)
   double precision, allocatable       :: dGvol(:), dGeig(:,:), dGconnect(:)
   double precision, allocatable       :: df(:), df0(:), dfmma(:)
   double precision, allocatable       :: energyK(:), energyM(:)

   ! IPOPT
   double precision                    :: CTOL
   integer, allocatable                :: ACON1(:), AVAR1(:)
   double precision, allocatable       :: GLOW(:), GUPP(:), AVEC(:)
     
   ! PDE Filter 
   double precision, allocatable       :: Tmatrix(:,:), rhotilde(:,:)
   double precision, allocatable       :: ed_rho(:,:,:), gammagp(:,:)
   type(sparse)                        :: Kfilter, Kfilterstar
   double precision                    :: radius, r2, times_el, el_len

   ! Different on/off parameters  
   integer                             :: Hprojection, Gscaling, diffCheck, update_done, standardRun, printout, enable_filter
   integer                             :: diffCheck_iteration, irestart, restart_iteration, enpose_design_sym, design_sym
    
   ! Heaviside & gamma
   double precision                    :: B1, rho0, delta0_E, delta0_rho
   double precision                    :: B, eta, Binit, Bfinal  

   ! Restart
   double precision, allocatable       :: rho_restart(:), coord_restart(:,:), m_el_vec_restart(:), enod_restart(:,:)
   integer                             :: el_m_restart, refinement_iter
 
   ! Design symmetry
   type(sparse)                        :: Qproj
   double precision, allocatable       :: QProjFull(:,:)
   integer                             :: el_m, el_s
   double precision, allocatable       :: df0_sym(:), dfmma_sym(:), df_c(:), df_sym(:)
   integer, allocatable                :: m_el_vec(:)

   ! Periodic boundary
   type(sparse)                        :: Proj, Projf, ProjBloch
   double precision, allocatable       :: ProjFull(:,:), ProjfFull(:,:), ProjBlochFull(:,:)
   integer                             :: nod_m, nod_s, dof_m, dof_s
   integer, allocatable                :: nod_s_vec(:)
   double precision                    :: width, height
   double precision, allocatable       :: UX(:), UY(:), VX(:), VY(:)
   integer                             :: ux_control, uy_control, vx_control, vy_control, nbr_free
   double precision                    :: ux_load, uy_load, vx_load, vy_load


end module opt_variables








