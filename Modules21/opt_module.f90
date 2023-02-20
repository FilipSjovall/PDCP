!#####################################################
! The module contains subroutines needed to run MMA
!#####################################################
! Created: 2019-12-13
! Anna Dalklint
! Lund University
!#####################################################
module opt_module

!$ use omp_lib
use opt_variables
use abaqus_util
use matlab_util
use mater_hyperel
use sparse_util
use fem_system
use FE_element

contains

!=======================================================================

!            Allocate global variables and initiate problem

!=======================================================================
subroutine opt_Init()
   implicit none
   double precision                    :: tmpreal
   double precision, allocatable       :: Dgp(:,:,:)
   integer                             :: ie, je
   character(120)                      :: xstring, iter_str

   call matinit(matptr)
   
   call abareadinp(filename,coord,enod,bcnod,bcval,flnod,flval) 

   if (size(coord,1).eq.3) then
      FE_dim = '3D'
   
      Write(*,*) '##################################' 
      write(*,*) '#           3D model             #'
      Write(*,*) '##################################' 
   elseif (size(coord,1).eq.2) then
      FE_dim = '2D'
   
      Write(*,*) '##################################' 
      write(*,*) '#           2D model             #'
      Write(*,*) '##################################' 
   endif
  


   ! 8 node 3D element
   if (FE_dim.eq.'3D') then
      ! Number of Gauss points
      ngp = 8
      ! Degrees of freedom per node
      dofnod = 3
      ! Number of element nodes
      elnod = 8
      ! Number of stress components
      els = 6
      ! Number of material stiffness components
      elm = 6
      ! Number of deformation gradient components
      eldef = 9
      ! So that the correct element functions are used
      th = 0d0
      
      write(*,*) 'Bloch implementation not ready for 3D'
      stop

   ! 4 node 2D element
   elseif (FE_dim.eq.'2D') then
      ngp = 4
      dofnod = 2
      elnod = 4
      ! Only plane strain available
      els = 4
      elm = 3
      eldef = 4
      ! Number of Brillouin zone boundaries
      nbr_bound = 4
      ! Start and stop in k=[0,0]
      nbr_k = nbr_bound*(nbr_sample+1)
      !nbr_k = 1
      ! Domain dimensions 
      width = maxval(coord(1,:))
      height = maxval(coord(2,:))
   endif
   
   nbc = size(bcval)
   nnod = maxval(enod)
   nelm = size(enod,2)
   ndof = nnod*dofnod
   eldof = dofnod*elnod
   
   allocate(m_el_vec(nelm), stat=ierr) 

   ! The number of design variables
   nbr_bound_vars = 0
   nbr_disp_vars = 0
   
   ! If you want to utilize design symmetry
   if (enpose_design_sym.eq.1) then
      if (design_sym.eq.45) then
         call sym_design_45deg()
      elseif (design_sym.eq.90) then
         call sym_design_90deg()
      else
         write(*,*) 'Error imposing design symmetry: Symmetry choice has not been implemented'
         stop
      endif
      
      nbr_extra_vars = nbr_bound_vars + nbr_disp_vars
      ndes = el_m*nbr_phase + nbr_extra_vars
      ndens = nelm*nbr_phase
      ndes_orig = ndens + nbr_extra_vars 
      
      allocate(df0_sym(ndes), stat=ierr)
      allocate(DF_C(ndes_orig), stat=ierr)
      allocate(DFMMA_sym(ndes*C), stat=ierr)
      allocate(DF_sym(ndes), stat=ierr)
   else
      el_m = nelm
      ndens = nelm*nbr_phase
      ! For initial design
      do ie = 1,nelm
         m_el_vec(ie) = ie
      enddo

      nbr_extra_vars = nbr_bound_vars + nbr_disp_vars
      ndes = ndens + nbr_extra_vars
      ndes_orig = ndes 
   endif

   write(*,*) 'Number of material phases     :', nbr_phase
   write(*,*) 'Number of design variables    :', ndes
   write(*,*) 'Number of nods                :', nnod
   write(*,*) 'Number of dofs                :', ndof
   write(*,*) 'Number of elements            :', nelm
   write(*,*) 'Number of boundary conditions :', nbc
   write(*,*) 'Number of concentrated loads  :', size(flval)
   
   ! FE-quantities
   allocate(dg(eldef,ngp,nelm), stat=ierr)
   allocate(ed(eldof,nelm), stat=ierr)
   allocate(a(ndof), stat=ierr)
   allocate(Fvec(ndof), stat=ierr)
   ! Check
   allocate(Dsmall(elm,elm), stat=ierr)
   allocate(Dgp(elm,elm,ngp), stat=ierr)
   allocate(ahist(nmax,ndof), stat=ierr)
   allocate(ares_hist(imax,ndof), stat=ierr)
   allocate(energyKHist(nbr_eigmax,Max_iter), stat=ierr)
   allocate(energyMhist(nbr_eigmax,Max_iter), stat=ierr)
   allocate(Pattern(ndof), stat=ierr) 
   allocate(R(dofnod,dofnod), stat=ierr) 
   allocate(gammagp(ngp,nelm), stat=ierr)
   allocate(fldof(size(flval)), stat=ierr)
   allocate(bcdof(size(bcval)), stat=ierr)
   allocate(history(3,ndof), stat=ierr)
   allocate(EnergyK(nbr_eigmax),stat=ierr)
   allocate(EnergyM(nbr_eigmax),stat=ierr)
   allocate(kvec_mat(dofnod,nbr_k),stat=ierr)
   allocate(E(nbr_phase),stat=ierr)
   allocate(massRho(nbr_phase),stat=ierr)
     
   ! Optimization
   allocate(tmp2(9), stat=ierr)
   allocate(f0valhist(MAX_ITER), stat=ierr)
   allocate(gvalhist(C,MAX_ITER), stat=ierr)
   allocate(glow(C),stat=ierr)
   allocate(gupp(C),stat=ierr)
   allocate(rho(ndes), stat=ierr)
   allocate(gmove(ndes), stat=ierr)
   allocate(df0(ndes_orig), stat=ierr)
   allocate(df(ndes_orig*C), stat=ierr)   
   allocate(DFMMA(ndes_orig*C), stat=ierr)
   allocate(dGvol(ndes_orig), stat=ierr)   
   allocate(dGeig(2,ndes_orig), stat=ierr)   
   allocate(dGconnect(ndes_orig), stat=ierr)   
   if (enable_filter.eq.1) then
      allocate(rhotilde(nnod,nbr_phase), stat=ierr)
   else
      allocate(rhotilde(nelm,nbr_phase), stat=ierr)
   endif
   allocate(ed_rho(elnod,nelm,nbr_phase), stat=ierr)
   allocate(Tmatrix(elnod,nelm), stat=ierr) 

   gmove(1:ndes-nbr_extra_vars) = gmove_topo
   !gmove(ndes-3) = gmove_bound
   !gmove(ndes-2) = gmove_bound
   !gmove(ndes-1) = gmove_disp
   !gmove(ndes) = gmove_disp


   if (enable_filter.eq.1) then
      radius = times_el*el_len/(2d0*dsqrt(3d0))
      r2 = radius**2d0
      write(*,*) '----------------------------------------------------------------------------'
      write(*,*) 'Filter length scale: ', radius 
      write(*,*) 'Element length', el_len, ' times: ', times_el
      write(*,*) '----------------------------------------------------------------------------'
   endif
  
   ahist = 0d0
   f0valhist = 0d0
   gvalhist = 0d0
   energyKHist = 0d0
   energyMHist = 0d0 

   ! Find dofs for boundary and loads
   Pattern=0D0
   call finddof(fldof,flnod,dofnod)
   call finddof(bcdof,bcnod,dofnod)
   Pattern(fldof)=flval

   ! dg = [F11 F12 F13 F21 F22 F23 F31 F32 F33]
   ! dg = [F11 F12 F21 F22]

   ed=0d0
   a=0d0
   
   if (FE_dim.eq.'3D') then
      dg=0d0
      dg(1,:,:)=1d0
      dg(5,:,:)=1d0
      dg(9,:,:)=1d0
   elseif (FE_dim.eq.'2D') then
      dg=0d0
      dg(1,:,:)=1d0
      dg(4,:,:)=1d0
   endif


   ! Sparse structures
   !call spatopdef(K,enod,dofnod)
   !call spatopdef(M,enod,dofnod)
   if (enable_filter.eq.1) then
      call spatopdef_sym(Kfilter,enod)  
   endif
   call spatopdef_sym(K,enod,dofnod)
   call spatopdef_sym(M,enod,dofnod)

   ! Material parameters
   E(1) = E1
   massrho(1) = massrho1
   if (nbr_phase.ge.2) then
      E(2) = E2
      massrho(2) = massrho2
      ! Multiply by E to get true bulk and shear modulus
      mp(1) = 1d0/3d0/(1d0-2d0*v) ! K
      mp(2) = 1d0/2d0/(1d0+v) ! G
   else
      ! Multiply by E to get true bulk and shear modulus
      mp(1) = E(1)/3d0/(1d0-2d0*v) ! K
      mp(2) = E(1)/2d0/(1d0+v) ! G
   endif

   ! Obtains K and G
   call neohooke_init(mp)
     
   ! Material tangent stiffness for pure displacement
   call dneohooke('tl',Dgp,dg(:,:,1))  
  
   ! Dsmall is the constant linear stiffness
   Dsmall = Dgp(:,:,1) 

   ! Calculate element volumes
   Vtot=0D0
   do ie=1,nelm
      call elem_vol(ve,coord(:,enod(:,ie)),th)
      Vtot = Vtot + Ve
   enddo
   Vtot = Volfrac*Vtot 
   Mtot = Vtot*DensTot
   
   ! Initial primitive lattice vectors
   R(:,1) = (/width, 0d0/)
   R(:,2) = (/0d0, height /)

   ! Update R due to deformation
   if (lagrangian_type.eq.'ul') then
      call update_lattice_vecs()
   endif
   
   ! Compute periodic boundary conditions
   call periodic_bc()

   ! Compute sampled wave vectors
   call compute_wave_vecs()
   
   ! Construct system matrices for real+imaginary part
   call spa_real_imaginary(K_re_im, M_re_im, K, M, ndof)

   ! IF We want the upper triangle of K as sparse pattern to obtain more from paradiso
   !call upper_triangular_matrix()
   
   ! Initiate PDE-filter matrices
   if (enable_filter.eq.1) then
      call PDE_filter_init()
   endif
   
   write( iter_str, '(i10)' )  refinement_iter
   xstring = 'data_common_'
   xstring=xstring(1:12)//iter_str(10:10)//'.mat'
   call matwrt2f(xstring,enod,'enod','w') 
   call matwrt2f(xstring,coord,'coord','u')  
    !return
end subroutine opt_Init

!=======================================================================

!                     Deallocate global variables

!=======================================================================
subroutine opt_destroy()
   implicit none
   integer :: ierr

   allocate(rho_restart(ndes),stat=ierr)
   allocate(coord_restart(dofnod,nnod),stat=ierr)
   allocate(m_el_vec_restart(nelm),stat=ierr)
   allocate(enod_restart(elnod,nelm),stat=ierr)

   rho_restart = rho
   coord_restart = coord
   m_el_vec_restart = m_el_vec
   el_m_restart = el_m
   enod_restart = enod

   deallocate(enod, stat=ierr)
   deallocate(bcnod, stat=ierr)
   deallocate(flnod, stat=ierr)
   deallocate(fbnod, stat=ierr)
   deallocate(coord, stat=ierr)
   deallocate(bcval, stat=ierr)
   deallocate(flval, stat=ierr)
   deallocate(fldof, stat=ierr)
   deallocate(fbdof, stat=ierr)
   deallocate(bcdof, stat=ierr)
   deallocate(df0_sym, stat=ierr)
   deallocate(DF_C, stat=ierr)
   deallocate(DFMMA_sym, stat=ierr)
   deallocate(DF_sym, stat=ierr)
   deallocate(dg, stat=ierr)
   deallocate(ed, stat=ierr)
   deallocate(a, stat=ierr)
   deallocate(Fvec, stat=ierr)
   deallocate(Dsmall, stat=ierr)
   deallocate(ahist, stat=ierr)
   deallocate(ares_hist, stat=ierr)
   deallocate(energyKHist, stat=ierr)
   deallocate(energyMhist, stat=ierr)
   deallocate(Pattern, stat=ierr) 
   deallocate(R, stat=ierr) 
   deallocate(gammagp, stat=ierr)
   deallocate(history, stat=ierr)
   deallocate(EnergyK,stat=ierr)
   deallocate(EnergyM,stat=ierr)
   deallocate(kvec_mat,stat=ierr)
   deallocate(E,stat=ierr)
   deallocate(massRho,stat=ierr)
   deallocate(tmp2, stat=ierr)
   deallocate(f0valhist, stat=ierr)
   deallocate(gvalhist, stat=ierr)
   deallocate(glow,stat=ierr)
   deallocate(gupp,stat=ierr)
   deallocate(rho, stat=ierr)
   deallocate(gmove, stat=ierr)
   deallocate(df0, stat=ierr)
   deallocate(df, stat=ierr)   
   deallocate(DFMMA, stat=ierr)
   deallocate(dGvol, stat=ierr)   
   deallocate(dGeig, stat=ierr)   
   deallocate(dGconnect, stat=ierr)   
   deallocate(rhotilde, stat=ierr)
   deallocate(ed_rho, stat=ierr)
   deallocate(Tmatrix, stat=ierr) 
   deallocate(EVALS_TOT, stat=ierr)
   deallocate(EVECS_TOT, stat=ierr)
   deallocate(evals_L, stat=ierr)
   deallocate(evals_U, stat=ierr)
   deallocate(evals_U_inv, stat=ierr)
   deallocate(AVAR1, stat=ierr)
   deallocate(ACON1, stat=ierr) 
   deallocate(AVEC, stat=ierr) 
   deallocate(m_el_vec, stat=ierr) 
   deallocate(nod_s_vec, stat=ierr)
   deallocate(UX, stat=ierr)
   deallocate(UY, stat=ierr)
   deallocate(VX, stat=ierr)
   deallocate(VY, stat=ierr)

   return

end subroutine opt_destroy


!=======================================================================

!                             MMA

!=======================================================================

subroutine MMA()
   implicit none
   double precision         :: Z, geps, tmpreal
   double precision         :: QMMA(NDES*C), P0MMA(NDES), Q0(NDES), UU(C), PMMA(NDES*C)
   integer                        :: IYFREE(C)
   double precision               :: GRADF(C), DSRCH(C), HESSF(C*(C+1)/2)
   double precision,allocatable   :: XVAL(:), XVAL_SYM(:), XLOW(:), XUP(:), XMIN(:), XMAX(:)
   double precision               :: ALFA(NDES), BETA(NDES)
   double precision               :: AMMA(C), BMMA(C), CMMA(C), YMMA(C), ULAM(C)
   double precision               :: FVAL(C), F0VAL
   double precision               :: DAT(1), rhomin, rhomax
   integer                        :: IDAT(2), ie, je, ki, i, IMMA, NEW_X, res_iter, je_res, ie_res
   double precision               :: xc, yc
   double precision               :: FMAX(C)
   ! Not used but required in iter_cb
   integer                        :: ALG_MODE, LS_TRIAL, ISTOP, NZ, TASK, ACON(ndes*C), AVAR(ndes*C)
   double precision               :: OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE
   double precision               :: ALPHA_DU, ALPHA_PR
   double precision               :: designRes(ndes), designResidual
   double precision,allocatable   :: xold1(:), xold2(:)
   integer                        :: indx1, ic
   character(80)                  :: s_plot
   character(120)                 :: xstring

   
   allocate(AVAR1(ndes_orig*C), stat=ierr)
   allocate(ACON1(ndes_orig*C), stat=ierr) 
   allocate(AVEC(ndes_orig*C), stat=ierr) 

   allocate(XVAL(ndes_orig), stat=ierr) 
   allocate(XVAL_sym(ndes), stat=ierr) 
   allocate(XLOW(ndes), stat=ierr) 
   allocate(XUP(ndes), stat=ierr) 
   allocate(XMIN(ndes), stat=ierr) 
   allocate(XMAX(ndes), stat=ierr) 
   allocate(XOLD1(ndes), stat=ierr) 
   allocate(XOLD2(ndes), stat=ierr) 
   
   GEPS =1d-9 ! MMA Tolerance parameter for the constraints.
   NZ = ndes_orig*C
   designTol = 1d-3
   IMMA = 0

   ! -----------------------------------------------------------------------
   !                         Derivative test
   ! -----------------------------------------------------------------------      
   if (diffCheck.eq.1) then

      if (enpose_design_sym.eq.1) then
         if (ipert.le.nelm) then
            ipert_sym = Qproj%ja(ipert)
         else
            ipert_sym = (ipert - nelm) + el_m
         endif
      endif

      if (restart_iteration.gt.1) then
         xstring = 'load data'
         indx1 = 9
         res_iter = restart_iteration - 1
         
         write(s_plot, '(i10)')  res_iter
         
         if (res_iter.lt.10) xstring=xstring(1:indx1)//'_000' //s_plot(10:10) //''' )'
         if (res_iter.ge.10) xstring=xstring(1:indx1)//'_00' //s_plot(9:10) //''' )'
         if (res_iter.ge.100) xstring=xstring(1:indx1)//'_0' //s_plot(8:10) //''' )'
         if (res_iter.ge.1000) xstring=xstring(1:indx1)//'_'//s_plot(7:10) //''' )'
         indx1=indx1+5
         xstring=xstring(1:indx1)//'.mat'
         indx1=indx1+4
         write(*,*) '##################################################################'
         write(*,*) 
         write(*,*) 'Loading initial density distribution in iteration', res_iter
         call matcommand(matptr,xstring(1:indx1))
         write(*,*) 
         write(*,*) '##################################################################'
         
         call matGetVal(matptr,'rho',rho)
         IMMA = res_iter
         write(*,*) '------------------------------------------------------------------'
         write(*,*) 'Performing a pertubation at optimization-step,', IMMA  
         write(*,*) 'Size of perturbation: ', pert
         write(*,*) 'Design variable perturbated: ', ipert
         write(*,*) '------------------------------------------------------------------'
         if (enpose_design_sym.eq.1) then
            rho(ipert_sym) = rho(ipert_sym)+pert
         else
            rho(ipert) = rho(ipert)+pert
         endif

         XVAL_SYM = rho
         !bound_U = rho(ndes-3)
         !bound_L = rho(ndes-2)
         !ux_load = rho(ndes-1)
         !vy_load = rho(ndes)
         
         call matGetVal(matptr,'tmp2',tmp2)
         p=tmp2(1)
         B=tmp2(2)
         iwait = tmp2(5)
         
         update_done = 0
      else         
         write(*,*) '------------------------------------------------------------------'
         write(*,*) 'Performing a pertubation at optimization-step,', IMMA  
         write(*,*) 'Size of perturbation: ', pert
         write(*,*) 'Design variable perturbated (sym): ', ipert_sym
         ! Initial guess
         do ie=1,el_m
            je = m_el_vec(ie)
            if ((sum(coord(1,enod(:,je)))/4d0.gt.300d0).and.(sum(coord(1,enod(:,je)))/4d0.lt.700d0)) then
               if ((sum(coord(2,enod(:,je)))/4d0.gt.300d0).and.(sum(coord(2,enod(:,je)))/4d0.lt.700d0)) then
                  ! Material 1
                  XVAL_SYM(ie) = 0.05d0
                  ! Material 2
                  if (nbr_phase.ge.2) then
                     XVAL_SYM(ie+el_m) = 0.01d0
                  endif
               else
                  XVAL_SYM(ie) = 0.5d0
                  if (nbr_phase.ge.2) then
                     XVAL_SYM(ie+el_m) = 0.7d0
                  endif
               endif
            else
               XVAL_SYM(ie) = 0.5d0
               if (nbr_phase.ge.2) then
                  XVAL_SYM(ie+el_m) = 0.7d0
               endif
            endif
         enddo

         ! Bound variables
         !XVAL_sym(ndes-3) = 2.1d0
         !XVAL_sym(ndes-2) = 2d0
         ! Displacement variables
         !XVAL_sym(ndes-1) = 1d0
         !XVAL_sym(ndes) = 1d0

         if (enpose_design_sym.eq.1) then
            XVAL_SYM(ipert_sym) = XVAL_SYM(ipert_sym)+pert
            write(*,*) 'Value of perturbated design variable: ', xval_sym(ipert_sym)
         else
            XVAL_SYM(ipert) = XVAL_SYM(ipert)+pert
            write(*,*) 'Value of perturbated design variable: ', xval_sym(ipert)
         endif
         write(*,*) '------------------------------------------------------------------'

         !bound_U = xval_sym(ndes-3)
         !bound_L = xval_sym(ndes-2)
         !ux_load = xval_sym(ndes-1)
         !vy_load = xval_sym(ndes)
      endif
   endif

   if ((irestart.eq.1).and.(diffCheck.ne.1)) then
      xstring = 'load data'
      indx1 = 9
      
      res_iter = restart_iteration - 1
         
      write(s_plot, '(i10)')  res_iter
         
      if (res_iter.lt.10) xstring=xstring(1:indx1)//'_000' //s_plot(10:10) //''' )'
      if (res_iter.ge.10) xstring=xstring(1:indx1)//'_00' //s_plot(9:10) //''' )'
      if (res_iter.ge.100) xstring=xstring(1:indx1)//'_0' //s_plot(8:10) //''' )'
      if (res_iter.ge.1000) xstring=xstring(1:indx1)//'_'//s_plot(7:10) //''' )'
      indx1=indx1+5
      xstring=xstring(1:indx1)//'.mat'
      indx1=indx1+4
      write(*,*) '##################################################################'
      write(*,*) 
      write(*,*) 'Loading initial density distribution in iteration', res_iter
      call matcommand(matptr,xstring(1:indx1))
      write(*,*) 
      write(*,*) '##################################################################'   
      call matGetVal(matptr,'rho',rho)

      XVAL_SYM = rho
      !bound_U = rho(ndes-3)
      !bound_L = rho(ndes-2)
      !ux_load = rho(ndes-1)
      !vy_load = rho(ndes)
           
      call matGetVal(matptr,'tmp2',tmp2)
      p=tmp2(1)
      B=tmp2(2)
      iwait = tmp2(5)
      IMMA = res_iter
   endif
   
   ! Regular start
   if ((irestart.eq.0).and.(diffCheck.eq.0).and.((refinement_iter.eq.0).or.(standardRun.eq.1))) then
      ! Initial guess
      if (initial_guess.eq.'random') then
         call random_seed()
         do ie=1,ndes-nbr_extra_vars
            CALL RANDOM_NUMBER(tmpreal)
            XVAL_SYM(ie) = tmpreal
         enddo
         !bound_L = 4.0d0
         !bound_U = 1.0d0
         !ux_load = 1d0
         !vy_load = 1d0
      elseif (initial_guess.eq.'circle') then
         ! Circle of softer material in middle
         do ie=1,el_m
            je = m_el_vec(ie)
            if (dsqrt((sum(coord(1,enod(:,je)))/4d0 - width/2d0)**2d0 + (sum(coord(2,enod(:,je)))/4d0 - height/2d0)**2d0).lt.350d0) then
               XVAL_SYM(ie) = 0.1d0
         !      XVAL_SYM(ie+el_m) = 0.1d0
            else
               XVAL_SYM(ie) = 0.5d0
         !      XVAL_SYM(ie+el_m) = 0.7d0
            endif
         enddo
         !bound_L = 4.0d0
         !bound_U = 1.0d0
         !ux_load = 1d0
         !vy_load = 1d0
      elseif (initial_guess.eq.'restart') then
         xstring = 'load initial_design.mat'
         call matcommand(matptr,xstring(1:23))
         call matGetVal(matptr,'Xval_sym',XVAL_SYM)
         !bound_L = xval_sym(ndes-3)
         !bound_U = xval_sym(ndes-2)
         !ux_load = xval_sym(ndes-1)
         !vy_load = xval_sym(ndes)
      endif


      !XVAL_SYM(ndes-3) = bound_U
      !XVAL_SYM(ndes-2) = bound_L
      !XVAL_SYM(ndes-1) = ux_load
      !XVAL_SYM(ndes) = vy_load
   endif

   if ((standardRun.eq.0).and.(refinement_iter.ne.0)) then
      do ie_res=1,el_m_restart
         je_res = m_el_vec_restart(ie_res)
         do ie=1,el_m
            je = m_el_vec(ie)
            xc = sum(coord(1,enod(:,je)))/(1d0*elnod)
            yc = sum(coord(2,enod(:,je)))/(1d0*elnod)
            if ((xc.lt.maxval(coord_restart(1,enod_restart(:,je_res)))).and.(yc.lt.maxval(coord_restart(2,enod_restart(:,je_res))))) then
               if ((xc.gt.minval(coord_restart(1,enod_restart(:,je_res)))).and.(yc.gt.minval(coord_restart(2,enod_restart(:,je_res))))) then
                  XVAL_SYM(ie) = rho_restart(ie_res)
               endif
            endif 
         enddo
      enddo
      !bound_U = rho_restart(size(rho_restart,1)-3)
      !bound_L = rho_restart(size(rho_restart,1)-2)
      !ux_load = rho_restart(size(rho_restart,1)-1)
      !vy_load = rho_restart(size(rho_restart,1))
      !XVAL_SYM(ndes-3) = bound_U 
      !XVAL_SYM(ndes-2) = bound_L 
      !XVAL_SYM(ndes-1) = ux_load
      !XVAL_SYM(ndes) = vy_load 
   endif

   if (allocated(rho_restart)) then
      deallocate(rho_restart)
      deallocate(coord_restart)
      deallocate(m_el_vec_restart)
      deallocate(enod_restart)
   endif

   if (enpose_design_sym.eq.1) then
      ! XVAL[dim]=nelm+2
      ! XVAL_SYM[dim]=ndes
      do je=1,nbr_phase
         XVAL(nelm*(je-1)+1:nelm*je) = matmul(Qproj, XVAL_SYM(el_m*(je-1)+1:el_m*je))
      enddo
   else
      ! No symmetry imposed
      XVAL(1:ndes_orig-nbr_extra_vars) = XVAL_SYM(1:ndes_orig-nbr_extra_vars)
   endif

   !XVAL(ndes_orig-3) = bound_U
   !XVAL(ndes_orig-2) = bound_L
   !XVAL(ndes_orig-1) = ux_load
   !XVAL(ndes_orig) = vy_load

   if (refinement_iter.eq.0) then
      call matwrt2f('initial_design_coarse.mat',Xval,'Xval','w')  
      call matwrt2f('initial_design_coarse.mat',Xval_sym,'Xval_sym','u')  
   else
      call matwrt2f('initial_design.mat',Xval,'Xval','w')  
      call matwrt2f('initial_design.mat',Xval_sym,'Xval_sym','u')  
   endif
   
   ! rho[dim]=ndes
   rho = XVAL_SYM
   XOLD1 = XVAL_SYM
   XOLD2 = XVAL_SYM

   IYFREE = 1
   
   ! Min/max design variables
   if (nbr_phase.gt.1) then
      rhomin = 1d-4
   else
      rhomin = 0.00000D0
   endif
   rhomax = 1D0
   
   XMIN(1:ndes-nbr_extra_vars) = rhomin  
   XMAX(1:ndes-nbr_extra_vars) = rhomax 
   ! Bound variables
   !xmin(ndes-3) = 1d-9 ! Close to zero
   !xmax(ndes-3) = 1d2 ! infinity
   !xmin(ndes-2) = 1d-9 ! Close to zero
   !xmax(ndes-2) = 1d2 ! infinity
   ! Displacement variables
   !xmin(ndes-1) = 0.0000d0 
   !xmax(ndes-1) = 1d3
   !xmin(ndes) = 0.0000d0 
   !xmax(ndes) = 1d3 

   AMMA = 0d0
   CMMA = 1D3
   
   FMAX = GUPP
  
   XLOW = 0d0
   XUP = 0D0
 
   FVAL = 0d0
   F0VAL = 0d0
   
   TASK = 1

   designResidual = 1
   NEW_X = 1
   IDAT = 0
   
   opt_finished = 0

   ! Run optimization loop
   do while ((designResidual.gt.designTOL).and.(IMMA.lt.max_iter))
      start = omp_get_wtime()
      
      IMMA = IMMA + 1
      ! IDAT(1) is the iteration counter
      ! IDAT(2) is the checker for new solve of eigenvalueproblem
      IDAT(1) = IMMA   
         
      ! Obtain value of objective
      call eval_F(NDENS, XVAL, NEW_X, F0VAL, IDAT, DAT, IERR)
      if (IERR.ne.0) write(*,*) 'Error in eval_F'

      ! Now everything has been updated
      NEW_X = 0
      
      ! Obtain sensitivities of objective
      call eval_grad_F(NDENS, XVAL, NEW_X, df0, IDAT, DAT, IERR)
      if (IERR.ne.0) write(*,*) 'Error in eval_grad_F'
      
      ! Obtain values of constraints
      call eval_G(NDENS, XVAL, NEW_X, C, FVAL, IDAT, DAT, IERR)
      if (IERR.ne.0) write(*,*) 'Error in eval_G'

      ! Obtain sensitivities of constraints
      call eval_jac_G(TASK, NDENS, XVAL, NEW_X, C, NZ, ACON, AVAR, AVEC, IDAT, DAT, IERR)
      
      do i = 1,C
         do ie = 1,ndes_orig
            ki = (ie-1)*C + i
            DFMMA(ki) = df((i-1)*ndes_orig+ie)
         enddo
      enddo
      if (IERR.ne.0) write(*,*) 'Error in eval_jac_G'

      if (enpose_design_sym.eq.1) then
         do je=1,nbr_phase
            df0_sym(el_m*(je-1)+1:el_m*je) = matmultrans(Qproj, df0(nelm*(je-1)+1:nelm*je))
         enddo
         ! Accomodate bound formulation
         df0_sym(ndes-(nbr_extra_vars-1):ndes) = df0(ndes_orig-(nbr_extra_vars-1):ndes_orig)

         do i = 1,C
            do ie = 1,ndes_orig
               ki = (ie-1)*C + i
               DF_C(ie) = DFMMA(KI)
            enddo

            do je=1,nbr_phase
               DF_sym(el_m*(je-1)+1:el_m*je) = matmultrans(Qproj,DF_C(nelm*(je-1)+1:nelm*je))
            enddo

            ! Accomodate bound formulation
            do je = 0,nbr_extra_vars-1
               DF_sym(ndes-je) = DFMMA(((ndes_orig-je)-1)*C + i)
            enddo

            do ie = 1,ndes
               ki = (ie-1)*C + i
               DFMMA_sym(ki) = DF_SYM(ie)
            enddo
         enddo

      else
         df0_sym = df0
         dfmma_sym = dfmma
      endif
      
      
      if ((irestart.eq.0).and.(diffCheck.eq.0)) then
         Write(*,*) 'Starting MMA update'
         ! Solve the subproblem
         !call MMASUB(IMMA,C,NDES,GEPS,GMOVE,IYFREE,XVAL_SYM,RHO, & 
         call MMASUB(IMMA,C,NDES,GEPS,IYFREE,XVAL_SYM,RHO, & 
                    XMIN,XMAX,XOLD1,XOLD2,XLOW,XUP, &
                    ALFA,BETA,AMMA,BMMA,CMMA,YMMA,Z,ULAM, &
                    F0VAL,FVAL,FMAX,DF0_SYM,DFMMA_SYM, &
                    PMMA,QMMA,P0MMA,Q0,UU,GRADF,DSRCH,HESSF)
         write(*,*) 'MMA update completed'
      endif

      XOLD2 = XOLD1 ! Design variables two iterations ago
      XOLD1 = XVAL_SYM  ! Design variables one iteration ago
      XVAL_SYM = RHO
     
      if (enpose_design_sym.eq.1) then
         ! XVAL[dim]=nelm+2
         ! XVAL_SYM[dim]=ndes
         do je=1,nbr_phase
            XVAL(nelm*(je-1)+1:nelm*je) = matmul(Qproj, XVAL_SYM(el_m*(je-1)+1:el_m*je))
         enddo
         do je = 0,nbr_extra_vars-1
            XVAL(ndes_orig-je) = XVAL_sym(ndes-je)
         enddo
      else
         ! No symmetry imposed
         XVAL = XVAL_SYM
      endif

      !bound_U = XVAL_SYM(ndes-3)
      !bound_L = XVAL_SYM(ndes-2)
      !ux_load = XVAL_SYM(ndes-1)
      !vy_load = XVAL_SYM(ndes)
      
      NEW_X = 1
      
      finish = omp_get_wtime()
   
      if (printout.eq.1) then
         write(*,*) 'Time taken: ', finish-start
      endif
      
      ! Calculate the design residual
      designRes = XVAL_SYM - XOLD1 
      
      if (printout.eq.1) then
         write(*,*) 'Maximum change in design (original) variable: ', maxval(dabs(designRes(1:ndes-nbr_extra_vars)))
      endif
      
      designResidual = dsqrt(dot_product(designRes,designRes))
      
      !call update_PDE_filter(XVAL, NEW_X)

      if (printout.eq.1) then
         write(*,*) '##################################################################'            
         write(*,*)
         write(*,*) '            Update of design in iteration', IMMA
         write(*,*)
         write(*,*) '##################################################################'
         write(*,*) 'Objective function: ', F0VAL
         write(*,*) 'Number of constraints: ', C
         write(*,*) '----------------------------------------------------'
         do i = 1,C
            write(*,*) 'Constraint', i, ': ', FVAL(i)
         enddo
         write(*,*) '----------------------------------------------------'
         write(*,*)
         write(*,*) 'Design residual: ', designResidual   
         write(*,*)   
         write(*,*) 'Min/max original design variable: '
         write(*,*) minval(rho(1:ndes-nbr_extra_vars)), -minval(-rho(1:ndes-nbr_extra_vars))
         !write(*,*) 'Bound design variables: '
         !write(*,*) rho(ndes-3), rho(ndes-2)
         !write(*,*) 'Displacement design variables: '
         !write(*,*) rho(ndes-1), rho(ndes)
         
         if (enable_filter.eq.1) then
            write(*,*) 'Min/max filtered original design variable: '
            do je=1,nbr_phase
               write(*,*) 'Material phase: ', je
               write(*,*) minval(rhotilde(:,je)), -minval(-rhotilde(:,je))
            enddo
         endif
         write(*,*)
      endif
      
      if (diffCheck.eq.1) then
         write(*,*) '----------------------------------------------------'   
         if (enpose_design_sym.eq.1) then
            write(*,*) 'df0(ipert): ', df0_sym(ipert_sym)
            do ic=1,C
               write(*,*) 'Constraint',i, ': df(ipert): ', dfmma_sym((ipert_sym-1)*C+ic)
            enddo
         else
            write(*,*) 'df0(ipert): ', df0(ipert)
            do ic=1,C
               write(*,*) 'Constraint', i, ': df(ipert): ', df(ndes_orig*(ic-1)+ipert)
            enddo
         endif
         write(*,*) '----------------------------------------------------'
      endif
      
      write(*,*) '##################################################################'

      ! Callback
      call iter_CB(ALG_MODE, IMMA, F0VAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT, DAT, ISTOP)
      if (IERR.ne.0) write(*,*) 'Error in iter_CB'
     
      if ((irestart.eq.1).or.(diffCheck.eq.1)) then
         return 
      endif
   enddo
   

   ! Output:
   if (IERR.eq.0) then
      write(*,*)
      write(*,*) 'The solution was found.'
      write(*,*)
      write(*,*) 'Objective function value ', F0VAL
      write(*,*)
   else
      write(*,*)
      write(*,*) 'An error occoured.'
      write(*,*) 'The error code is ',IERR
      write(*,*)
   endif

   opt_finished = 1
   
   loadstep_bandgap = 0
   call iter_CB(ALG_MODE, IMMA, F0VAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT, DAT, ISTOP)

   ! Perform final computations of bandgap evolution over the loadsteps
   call bandgap_evol(XVAL, DAT, IDAT)
       
   return
end subroutine MMA 






!=======================================================================

!                  Create symmetric density pattern
!             (Implemented for 45 degree design symmetry)
!               (Implemented for a square design domain)
!                      TODO: Implement for 3d

!=======================================================================
subroutine sym_design_45deg()
   implicit none
   double precision     :: xc, yc, xc_s, yc_s, c_s, c_m, xc_m, yc_m
   integer              :: i, j, ii, master, slave, n, im, is
   double precision     :: ncoord(dofnod,elnod), mcoord(dofnod,elnod)
   integer,allocatable  :: cellQ(:,:), s_el_vec(:)

   write(*,*)
   write(*,*) 'Start to create symmetric 45 degree design pattern'

   allocate(cellQ(nelm,2), stat=ierr) 
   allocate(s_el_vec(nelm), stat=ierr) 

   el_m = 0
   el_s = 0
   m_el_vec = 0
   s_el_vec = 0

   ! Locate "master"/"slave" elements
   do i=1,nelm
      ncoord = coord(:,enod(:,i))
      master = 0
      slave = 0
      xc = sum(ncoord(1,:))/(1d0*elnod) - width/2d0
      yc = sum(ncoord(2,:))/(1d0*elnod) - height/2d0

      if ((xc.gt.0d0).and.(yc.gt.0d0).and.((yc-xc).le.1d-3)) then
         master = 1
      endif

      if (master.eq.1) then
         el_m = el_m + 1
         m_el_vec(el_m) = i
      else
         el_s = el_s + 1
         s_el_vec(el_s) = i
      endif
   enddo
   
   ! "Master" elements
   do i=1,el_m
      im = m_el_vec(i)
      cellQ(im,:) = (/im,i/)
   enddo

   ! "Slave" elements
   !$OMP PARALLEL DO DEFAULT(NONE) &
   !$OMP SHARED(cellQ,s_el_vec,m_el_vec,coord,enod,width,height,elnod,el_s,el_m) & 
   !$OMP PRIVATE(i,j,is,im,ncoord,mcoord,xc_s,xc_m,yc_s,yc_m,c_s,c_m) 
   do i=1,el_s
      is = s_el_vec(i)
      ncoord = coord(:,enod(:,is))
      ! Translate coordinate system such that the center is in (0,0)
      xc_s = sum(ncoord(1,:))/(1d0*elnod) - width/2d0
      yc_s = sum(ncoord(2,:))/(1d0*elnod) - height/2d0
      c_s = dsqrt(xc_s**2d0 + yc_s**2d0)

      do j=1,el_m
         im = m_el_vec(j)
         mcoord = coord(:,enod(:,im))
         xc_m = sum(mcoord(1,:))/(1d0*elnod) - width/2d0 
         yc_m = sum(mcoord(2,:))/(1d0*elnod) - height/2d0 
         c_m = dsqrt(xc_m**2d0 + yc_m**2d0)

         if ((dabs(c_s-c_m).le.1d-4).and.(dabs(maxval((/dabs(xc_m), dabs(yc_m)/))-maxval((/dabs(xc_s), dabs(yc_s)/))).lt.1d-4)) then
            cellQ(is,:) = (/is,j/)
         endif
      enddo
   enddo
   !$OMP END PARALLEL DO

   call spaCellDef(Qproj,cellQ,nelm,el_m)
   Qproj%a=1D0
   
   !allocate(QProjfull(nelm,el_m), stat=ierr)  
   !call spa2full(Qprojfull,Qproj)
   !call matwrt2f('banan.mat',Qprojfull,'Qproj','w')  

   write(*,*) 'Symmetric design pattern created'
   write(*,*)
    
   return
end subroutine sym_design_45deg



!=======================================================================

!                  Create symmetric density pattern
!             (Implemented for 90 degree design symmetry)
!               (Implemented for a square design domain)
!                      TODO: Implement for 3d

!=======================================================================
subroutine sym_design_90deg()
   implicit none
   double precision     :: xc, yc, xc_s, yc_s, c_s, c_m, xc_m, yc_m
   integer              :: i, j, ii, master, slave, n, im, is
   double precision     :: ncoord(dofnod,elnod), mcoord(dofnod,elnod)
   integer,allocatable  :: cellQ(:,:), s_el_vec(:)

   write(*,*)
   write(*,*) 'Start to create symmetric 90 degree design pattern'

   allocate(cellQ(nelm,2), stat=ierr) 
   allocate(s_el_vec(nelm), stat=ierr) 

   el_m = 0
   el_s = 0
   m_el_vec = 0
   s_el_vec = 0

   ! Locate "master"/"slave" elements
   do i=1,nelm
      ncoord = coord(:,enod(:,i))
      master = 0
      slave = 0
      xc = sum(ncoord(1,:))/(1d0*elnod) - width/2d0
      yc = sum(ncoord(2,:))/(1d0*elnod) - height/2d0

      if ((xc.gt.0d0).and.(yc.gt.0d0)) then
         master = 1
      endif

      if (master.eq.1) then
         el_m = el_m + 1
         m_el_vec(el_m) = i
      else
         el_s = el_s + 1
         s_el_vec(el_s) = i
      endif
   enddo
   
   ! "Master" elements
   do i=1,el_m
      im = m_el_vec(i)
      cellQ(im,:) = (/im,i/)
   enddo

   ! "Slave" elements
   !$OMP PARALLEL DO DEFAULT(NONE) &
   !$OMP SHARED(cellQ,s_el_vec,m_el_vec,coord,enod,width,height,elnod,el_s,el_m) & 
   !$OMP PRIVATE(i,j,is,im,ncoord,mcoord,xc_s,xc_m,yc_s,yc_m,c_s,c_m) 
   do i=1,el_s
      is = s_el_vec(i)
      ncoord = coord(:,enod(:,is))
      ! Translate coordinate system such that the center is in (0,0)
      xc_s = sum(ncoord(1,:))/(1d0*elnod) - width/2d0
      yc_s = sum(ncoord(2,:))/(1d0*elnod) - height/2d0
      c_s = dsqrt(xc_s**2d0 + yc_s**2d0)

      do j=1,el_m
         im = m_el_vec(j)
         mcoord = coord(:,enod(:,im))
         xc_m = sum(mcoord(1,:))/(1d0*elnod) - width/2d0 
         yc_m = sum(mcoord(2,:))/(1d0*elnod) - height/2d0 
         c_m = dsqrt(xc_m**2d0 + yc_m**2d0)

         if ((dabs(c_s-c_m).le.1d-4).and.((dabs(dabs(xc_s)-xc_m)).lt.1d-4).and.((dabs(dabs(yc_s)-yc_m)).lt.1d-4)) then
            cellQ(is,:) = (/is,j/)
         endif
      enddo
   enddo
   !$OMP END PARALLEL DO

   call spaCellDef(Qproj,cellQ,nelm,el_m)
   Qproj%a=1D0
   
   !allocate(QProjfull(nelm,el_m), stat=ierr)  
   !call spa2full(Qprojfull,Qproj)
   !call matwrt2f('banan.mat',Qprojfull,'Qproj','w')  
   !stop

   write(*,*) 'Symmetric design pattern created'
   write(*,*)
    
   return
end subroutine sym_design_90deg




!=======================================================================

!             Create periodic boundary conditions matrices
!               (Implemented for a square design domain)
!                      TODO: Implement for 3d

!=======================================================================
subroutine periodic_bc()
   implicit none
   integer                             :: i, j, ii, Xmaster, Ymaster, Xslave, Yslave, n, im, is, corner, m_corner_nod
   integer                             :: m_dof(dofnod), s_dof(dofnod), m_corner_dof(dofnod)
   double precision                    :: ncoord(dofnod), mcoord(dofnod)
   double precision                    :: corner_points(elnod-1,dofnod)
   integer, allocatable                :: cellPBlochTMP(:,:), cellP(:,:),cellPf(:,:), cellPBloch(:,:)
   integer,allocatable                 :: nod_m_vec(:)
   integer                             :: side_nodes, extra_terms, counter, counter_col

   nod_m = 0
   nod_s = 0
   
   allocate(nod_s_vec(nnod), stat=ierr)
   allocate(nod_m_vec(nnod), stat=ierr)

   allocate(cellPf(nnod,2), stat=ierr)
   allocate(UX(ndof), stat=ierr)
   allocate(UY(ndof), stat=ierr)
   allocate(VX(ndof), stat=ierr)
   allocate(VY(ndof), stat=ierr)
   UX = 0
   UY = 0
   VX = 0
   VY = 0

   write(*,*)
   write(*,*) 'Start to create periodic boundary pattern'
   nod_m_vec = 0
   nod_s_vec = 0
   
   do n=1,nnod
      ncoord = coord(:,n)
      Xmaster = 0
      Ymaster = 0
      Xslave = 0
      Yslave = 0

      if (dabs(ncoord(1)-width).lt.1d-5) then
         Xslave = 1
         UX(n*dofnod-1) = 1
         UY(n*dofnod) = 1
      endif

      if (dabs(ncoord(2)-height).lt.1d-5) then
         Yslave = 1
         VX(n*dofnod-1) = 1
         VY(n*dofnod) = 1
      endif

      if ((Xslave.eq.0).and.(Yslave.eq.0)) then
         nod_m = nod_m + 1
         nod_m_vec(nod_m) = n
      else
         nod_s = nod_s + 1
         nod_s_vec(nod_s) = n
      endif
   enddo

   ! Number of nodes on one side of the square domain
   side_nodes = sum(UX)
   ! Assuming norm(u_scalar) = 1 or 0
   nbr_free = abs(ux_control-1 + uy_control-1 + vx_control-1 + vy_control-1)
   extra_terms = nbr_free*side_nodes

   allocate(cellP(ndof+extra_terms,2), stat=ierr)

   dof_m = nod_m*dofnod
   dof_s = nod_s*dofnod
   
   !allocate(Projfull(ndof,dof_m+nbr_free), stat=ierr)  
   !allocate(Projffull(nnod,nod_m), stat=ierr)  
   !allocate(Kstarfull(dof_m,dof_m), stat=ierr)  
   !allocate(Kfilterstarfull(nod_m,nod_m), stat=ierr)  

   !allocate(Mfull(ndof,ndof), stat=ierr)  
   !allocate(Kfull(ndof,ndof), stat=ierr)  
   !allocate(ProjBlochFull(ndof*2,dof_m*2), stat=ierr)
    
   ! Fix allocations of global variables
   
   allocate(cellPBloch(2*ndof+2*dof_s,2), stat=ierr) ! [Real; Imaginary]
   
   allocate(cellPBlochTMP(2*dof_s,2), stat=ierr) ! [Real; Imaginary]

   ! Master dofs
   !$OMP PARALLEL DO DEFAULT(NONE) &
   !$OMP SHARED(cellP,cellPf,nod_m,nod_m_vec,coord,dofnod,m_corner_dof,m_corner_nod) & 
   !$OMP PRIVATE(im,ncoord,m_dof,i,j) 
   do i=1,nod_m
      im = nod_m_vec(i)
      ncoord = coord(:,im)
      m_dof = (/im*dofnod-1,im*dofnod/)

      do j=1,dofnod
         cellP(m_dof(j),:) = (/m_dof(j), (i-1)*dofnod + j/)
      enddo
      cellPf(im,:) = (/im, i/)

      ! Find master corner node
      if ((dabs(ncoord(1)).lt.1d-6).and.(dabs(ncoord(2)).lt.1d-6)) then
         m_corner_dof = (/im*dofnod-1,im*dofnod/)
         m_corner_nod = im 
      endif
   enddo
   !$OMP END PARALLEL DO

   corner_points(1,:) = (/width, 0d0/)
   corner_points(2,:) = (/width, height/)
   corner_points(3,:) = (/0d0, height/)

   ! Slave dofs
   !$OMP PARALLEL DO DEFAULT(NONE) &
   !$OMP SHARED(cellP,cellPf,cellPBlochtmp,nod_m,nod_m_vec,coord,dofnod,m_corner_dof) & 
   !$OMP SHARED(corner_points,nod_s,nod_s_vec,m_corner_nod,dof_m,ndof,dof_s) & 
   !$OMP PRIVATE(is,im,ncoord,mcoord,s_dof,i,j,ii,corner) 
   do i=1,nod_s
      is = nod_s_vec(i)
      ncoord = coord(:,is)
      s_dof = (/is*dofnod-1,is*dofnod/)
      corner = 0

      ! Check if corner point
      do j=1,size(corner_points,1)
         if ((dabs(ncoord(1)-corner_points(j,1)).lt.1d-5).and.(dabs(ncoord(2)-corner_points(j,2)).lt.1d-5)) then
            corner = 1
         endif
      enddo

      if (corner.eq.1) then
         do j=1,dofnod
            cellP(s_dof(j),:) = (/s_dof(j), m_corner_dof(j)/)
            ! Imaginary and real part
            cellPBlochTMP((i-1)*dofnod+j,:) = (/s_dof(j), m_corner_dof(j)+dof_m/) ! Upper right
            cellPBlochTMP((i-1)*dofnod+j+dof_s,:) = (/s_dof(j)+ndof, m_corner_dof(j)/) !Lower left
         enddo
         cellPf(is,:) = (/is, m_corner_nod/)
      else
         do j=1,nod_m
            im = nod_m_vec(j)
            mcoord = coord(:,im)
            
            if ((dabs(mcoord(1)).lt.1d-6).and.(dabs(mcoord(2)-ncoord(2)).lt.1d-5)) then
               do ii=1,dofnod
                  cellP(s_dof(ii),:) = (/s_dof(ii), (j-1)*dofnod + ii/)
                  ! Imaginary and real part
                  cellPBlochTMP((i-1)*dofnod+ii,:) = (/s_dof(ii), (j-1)*dofnod + ii + dof_m/) ! Upper right
                  cellPBlochTMP((i-1)*dofnod+ii+dof_s,:) = (/s_dof(ii)+ndof, (j-1)*dofnod + ii/) !Lower left
               enddo
               cellPf(is,:) = (/is, j/)
            elseif ((dabs(mcoord(2)).lt.1d-6).and.(dabs(mcoord(1)-ncoord(1)).lt.1d-5)) then
               do ii=1,dofnod
                  cellP(s_dof(ii),:) = (/s_dof(ii), (j-1)*dofnod + ii/)
                  ! Imaginary and real part
                  cellPBlochTMP((i-1)*dofnod+ii,:) = (/s_dof(ii), (j-1)*dofnod + ii + dof_m/) ! Upper right
                  cellPBlochTMP((i-1)*dofnod+ii+dof_s,:) = (/s_dof(ii)+ndof, (j-1)*dofnod + ii/) !Lower left
               enddo
               cellPf(is,:) = (/is, j/)
            endif 
         enddo
      endif
   enddo
   !$OMP END PARALLEL DO

   counter = 1
   counter_col = 1
   if (ux_control.eq.0) then
      do i=1,nod_s
         is = nod_s_vec(i)
         ncoord = coord(:,is)
         s_dof = (/is*dofnod-1,is*dofnod/)

         ! Right side
         if (dabs(ncoord(1)-width).lt.1d-5) then
            cellP(ndof+counter,:) = (/s_dof(1), dof_m+counter_col/)
            counter = counter + 1
         endif
      enddo
      counter_col = counter_col+1
   endif

   if (uy_control.eq.0) then
      do i=1,nod_s
         is = nod_s_vec(i)
         ncoord = coord(:,is)
         s_dof = (/is*dofnod-1,is*dofnod/)

         ! Right side
         if (dabs(ncoord(1)-width).lt.1d-5) then
            cellP(ndof+counter,:) = (/s_dof(2), dof_m+counter_col/)
            counter = counter + 1
         endif
      enddo
      counter_col = counter_col+1
   endif
   
   if (vx_control.eq.0) then
      do i=1,nod_s
         is = nod_s_vec(i)
         ncoord = coord(:,is)
         s_dof = (/is*dofnod-1,is*dofnod/)

         ! Top side
         if (dabs(ncoord(2)-height).lt.1d-5) then
            cellP(ndof+counter,:) = (/s_dof(1), dof_m+counter_col/)
            counter = counter + 1
         endif
      enddo
      counter_col = counter_col+1
   endif

   if (vy_control.eq.0) then
      do i=1,nod_s
         is = nod_s_vec(i)
         ncoord = coord(:,is)
         s_dof = (/is*dofnod-1,is*dofnod/)

         ! Top side
         if (dabs(ncoord(2)-height).lt.1d-5) then
            cellP(ndof+counter,:) = (/s_dof(2), dof_m+counter_col/)
            counter = counter + 1
         endif
      enddo
      counter_col = counter_col+1
   endif

   cellPBloch(1:ndof,:) = cellP(1:ndof,:)
   cellPBloch(ndof+1:2*ndof,1) = cellP(1:ndof,1)+ndof
   cellPBloch(ndof+1:2*ndof,2) = cellP(1:ndof,2)+dof_m ! [Real; Imaginary]
   cellPBloch(2*ndof+1:2*ndof+2*dof_s,:) = cellPBlochTMP(:,:) 

   call spaCellDef(Proj,cellP,ndof,dof_m+nbr_free)
   Proj%a=1D0
   !call spa2full(Projfull,Proj)
   !call matwrt2f('banan.mat',Projfull,'Projfull','u')  

    
   call spaCellDef(ProjBloch,cellPBloch,2*ndof,2*dof_m)
   ProjBloch%a=1D0

   allocate(EVALS_TOT(nbr_eigMax,nbr_k), stat=ierr)
   allocate(EVECS_TOT(nbr_eigMax,(ndof*2)*nbr_k), stat=ierr)
   allocate(evals_L(nbr_eig_l,nbr_k), stat=ierr)
   allocate(evals_U(nbr_eig_U,nbr_k), stat=ierr)
   allocate(evals_U_inv(nbr_eig_U,nbr_k), stat=ierr)
   
   call spaCellDef(Projf,cellPf,nnod,nod_m)
   Projf%a=1D0
    
   write(*,*) 'Periodic boundary pattern created'
   write(*,*)
  

   return
end subroutine periodic_bc


!=======================================================================

!   Update the lattice vectors according to the macroscopic deformation
!                               gradient 
!               (Implemented for a square design domain)
!                      TODO: Implement for 3d

!=======================================================================
subroutine update_lattice_vecs()
   implicit none
   double precision     :: F_M(2,2)

   F_M = 0d0
   F_M(1,1) = 1d0
   F_M(2,2) = 1d0

   F_M(1,1) = F_M(1,1) + ux_load/width 
   F_M(2,2) = F_M(2,2) + vy_load/height
   F_M(1,2) = F_M(1,2) + uy_load/width 
   F_M(2,1) = F_M(2,1) + vx_load/height

   R = matmul(F_M, R)

   write(*,*)
   write(*,*) 'Lattice vectors updated with F^M'
   write(*,*)

   return
end subroutine update_lattice_vecs


!=======================================================================

!         Update Bloch periodic boundary conditions matrices 
!               (Implemented for a square design domain)
!                      TODO: Implement for 3d

!=======================================================================
!subroutine periodic_bc_bloch(kvec, ProjMat, conjugate)
subroutine periodic_bc_bloch(kvec, conjugate)
   implicit none
   !type(Sparse)         :: ProjMat
   double precision     :: kvec(dofnod), arg1, arg2, arg12, arg, ncoord(dofnod)
   integer              :: i, j, rows, is, s_dof(dofnod), indx_re, indx_im, conjugate

   if (printout.eq.1) then
      write(*,*) 'Start to create periodic Bloch boundary pattern for kvec: ', kvec
   endif
   
   if (conjugate.eq.1) then
      arg1 = dot_product(kvec,R(:,1))
      arg2 = dot_product(kvec,R(:,2))
      arg12 = dot_product(kvec,R(:,1)+R(:,2))
   elseif (conjugate.eq.-1) then
      arg1 = -dot_product(kvec,R(:,1))
      arg2 = -dot_product(kvec,R(:,2))
      arg12 = -dot_product(kvec,R(:,1)+R(:,2))
   else
      write(*,*) 'Periodic_bc_bloch: Wrong input, stopping'
      stop
   endif

   do i = 1,nod_s
      is = nod_s_vec(i)
      ncoord = coord(:,is)
      s_dof = (/is*dofnod-1,is*dofnod/)
      
      !if ((dabs(ncoord(1).eq.width).and.(ncoord(2).ne.height)) then
      if ((dabs(ncoord(1)-width).lt.1d-5).and.(dabs(ncoord(2)-height).gt.1d-2)) then
         arg = arg1
      !elseif ((ncoord(1).ne.width).and.(ncoord(2).eq.height)) then
      elseif ((dabs(ncoord(1)-width).gt.1d-2).and.(dabs(ncoord(2)-height).lt.1d-5)) then
         arg = arg2
      else
         arg = arg12
      endif

      do j = 1,size(s_dof)
         ! Index of value in sparse matrix
         indx_re = ProjBloch%ia(s_dof(j)) 
         ! We know that there are two entries on each row
         indx_im = ProjBloch%ia(s_dof(j))+1 
         
         ! Top right quadrant
         ProjBloch%a(indx_re) = dcos(arg)
         ProjBloch%a(indx_im) = -dsin(arg)

         ! Index of value in sparse matrix
         indx_re = ProjBloch%ia(s_dof(j)+ndof) 
         ! We know that there are two entries on each row
         indx_im = ProjBloch%ia(s_dof(j)+ndof)+1 

         ! Bottom left quadrant
         ProjBloch%a(indx_re) = dsin(arg)
         ProjBloch%a(indx_im) = dcos(arg)
      enddo
   enddo

   if (printout.eq.1) then
      write(*,*) 'Bloch periodic boundary pattern created'
   endif
  

   return
end subroutine periodic_bc_bloch






!=======================================================================

!            Compute the wave vectors in the Brillouin zone
!               (Implemented for a square design domain)
!                      TODO: Implement for 3d

!=======================================================================
subroutine compute_wave_vecs()
   implicit none
   double precision     :: kx_start, ky_start, kx_end, ky_end
   double precision     :: kx, ky
   integer              :: i, incr, j, jstart, jend
   double precision     :: R_3D(3,3), G_3D(3,3), G(dofnod,dofnod), norm_a1xa2

   write(*,*)
   write(*,*) 'Computing kvec_mat where each boundary is sampled: ', nbr_sample, ' times'
   
   ! nbr_sample = number of interior sampling points on each boundary of the irreductible
   ! Brillouin zone
   incr = 0

   R_3D = 0d0
   R_3D(1:2,1) = R(:,1)
   R_3D(1:2,2) = R(:,2)

   ! Compute out of a_1,a_2-plane vector
   call cross_product(R_3D(:,1),R_3D(:,2),R_3D(:,3))
   norm_a1xa2 = dsqrt(dot_product(R_3D(:,3),R_3D(:,3)))

   ! Compute reciprocal lattice vectors
   call cross_product(R_3D(:,2),R_3D(:,3)/norm_a1xa2,G_3D(:,1))
   call cross_product(R_3D(:,3)/norm_a1xa2,R_3D(:,1),G_3D(:,2))

   G = 0d0
   G(:,1) = (2d0*pi/norm_a1xa2)*G_3D(1:2,1)
   G(:,2) = (2d0*pi/norm_a1xa2)*G_3D(1:2,2)

   !write(*,*) 'G1', G(:,1)
   !write(*,*) 'G2', G(:,2)
   !write(*,*) 'a1*b2', dot_product(R(:,1),G(:,2))
   !write(*,*) 'a1*b1', dot_product(R(:,1),G(:,1))

   do i=1,nbr_Bound
      if (i.eq.1) then
         ky_start = 0d0
         ky_end = 0d0
         kx_start = 0d0
         kx_end = G(1,1)/2d0 !pi/R(1,1)
      elseif (i.eq.2) then
         ky_start = 0d0
         ky_end = G(2,2)/2d0 !pi/R(2,2)
         kx_start = G(1,1)/2d0 !pi/R(1,1)
         kx_end = G(1,1)/2d0 !pi/R(1,1)
      elseif (i.eq.3) then
         !ky_start = G(2,2)/2d0 !pi/R(2,2)
         !ky_end = 0d0
         !kx_start = G(1,1)/2d0 !pi/R(1,1)
         !kx_end = 0d0
         ky_start = G(2,2)/2d0
         ky_end = G(2,2)/2d0
         kx_start = G(1,1)/2d0
         kx_end = 0d0
      elseif (i.eq.4) then
         ky_start = G(2,2)/2d0
         ky_end = 0d0
         kx_start = 0d0
         kx_end = 0d0
      endif
      
      do j=1,nbr_sample+1
         incr = incr + 1 

         kx = kx_start + (j-1)*(kx_end-kx_start)/(1d0*(nbr_sample+1d0))
         ky = ky_start + (j-1)*(ky_end-ky_start)/(1d0*(nbr_sample+1d0))

         kvec_mat(:,incr) = (/kx, ky/)
      enddo
   enddo

   ! Don't compute for k = [0,0] again...
   !kvec_mat(:,incr+1) = (/0d0,0d0/) !(/1d-3, 1d-3/)

   write(*,*) 'kvec_mat created'
   write(*,*)
  

   return
end subroutine compute_wave_vecs


!=======================================================================

!                    Initiate PDE-filter matrices

!=======================================================================
subroutine PDE_filter_init()
   implicit none
   integer            :: ie
   double precision  :: ken(elnod,elnod), Men(elnod,elnod), fen(elnod)
   
   ! Calculate T, M_filter and K_filter
   ! M_filter and K_filter are added together for simplification     
   Kfilter%a=0d0
     
   do ie=1,nelm
      call elem_flow_stiff(ken,coord(:,enod(:,ie)),th,r2)  
      call elem_flow_mass(Men,coord(:,enod(:,ie)),th,1d0)

      ken = ken + Men
      call elem_flow_bf(fen,coord(:,enod(:,ie)),th,1d0)

      Tmatrix(:,ie) = fen
      call assem(Kfilter,ken,enod(:,ie),1) 
   enddo  

   return
end subroutine PDE_filter_init
 

!=======================================================================

!              Update filtered design variables

!=======================================================================
subroutine update_PDE_filter(XVAL, NEW_X)
   implicit none
   double precision               :: XVAL(:)
   integer                        :: ie, NegativeE, NEW_X, i
   double precision,allocatable   :: ffilterstar(:), ffilter(:)
   double precision               :: g, g_pert, gp
   
   allocate(ffilterstar(nod_m), stat=ierr) 
   allocate(ffilter(nnod), stat=ierr) 

   if (NEW_X.eq.1) then
      if (printout.eq.1) then
         write(*,*) '### Update PDE filter ###'
      endif
      !rho = XVAL 
      
      do i=1,nbr_phase
         ffilter = 0d0

         ! Denotes T*XVAL as ffilter
         do ie=1,nelm
            ffilter(enod(:,ie)) = ffilter(enod(:,ie)) + Tmatrix(:,ie)*XVAL(nelm*(i-1)+ie) 
         enddo
         NegativeE = 0d0
      
         ! Solve to obtain rho_tilde
         call spamul_symprod(Kfilterstar,Projf,Kfilter) !P'*K*P
         
         ffilterstar = matmulTrans(Projf,ffilter) !P'*F
         
         call solveq(Kfilterstar,ffilterstar,NegativeE) 

         rhotilde(:,i)=matmul(Projf,ffilterstar) 

         ! Obtain element rho_tilde (ed_rho)
         call extract(ed_rho(:,:,i),rhotilde(:,i),enod,1)
      enddo
   endif

   deallocate(ffilterstar)
   deallocate(ffilter)

   return
end subroutine update_PDE_filter


!=======================================================================

!                       Update gamma

!=======================================================================
subroutine update_gamma(XVAL)
   implicit none
   integer            :: ie, igp, je
   double precision   :: XVAL(:), rhogp(nbr_phase,ngp)

   if (printout.eq.1) then
      write(*,*) '### Update gamma-field ###'
   endif
   
   gammagp = 0d0

   !$OMP PARALLEL DO DEFAULT(NONE) &
   !$OMP SHARED(nelm,XVAL,ed_rho,gammagp,ngp,p,nbr_phase,enable_filter) & 
   !$OMP PRIVATE(rhogp,ie,igp,je) 
   do ie=1,nelm
      do je=1,nbr_phase
         if (enable_filter.eq.1) then
            call elem_gp(rhogp(je,:),ed_rho(:,ie,je)) 
         else
            rhogp(je,:) = XVAL((je-1)*nelm + ie) 
         endif
      enddo
      do igp=1,ngp      
         ! Only use the design variable controlling the void/solid phase
         call gammafun(gammagp(igp,ie),rhogp(:,igp),p)
         if ((gammagp(igp,ie)).lt.0D0)  then
            write(*,*) 'This xval value yields a negative gammagp', rhogp(:,igp)          
         endif   
      enddo
     enddo 
   !$OMP end parallel do

   return
end subroutine update_gamma


!=======================================================================

!                     Choose the optimizer

!=======================================================================

subroutine opt(choice)
   implicit none
   character(80) :: choice
   
   if (choice.eq.'MMA') then
      call MMA()
   elseif (choice.eq.'IPOPT') then
      write(*,*) 'IPOPT not implemented for usage on lunarc!'
      stop
   else
      stop 'This optimizer has not been implemented!'
   endif


end subroutine opt



!=======================================================================

!                   Computation of objective function

!=======================================================================
subroutine eval_F(NELM, XVAL, NEW_X, F, IDAT, DAT, IERR)
   implicit none
   integer                              :: NELM, NEW_X
   double precision                     :: F, F0VAL, XVAL(:)
   double precision                     :: DAT(:)
   integer                              :: IDAT(:)
   integer                              :: IERR, je, ie
   
   if (printout.eq.1) then
      write(*,*) '### Compute objective function ###'
   endif

   if (enable_filter.eq.1) then
      call update_PDE_filter(XVAL, NEW_X)
   else
      if (NEW_X.eq.1) then
         do je=1,nbr_phase
            do ie=1,nelm
               rhotilde(ie,je) = XVAL((je-1)*nelm + ie)
            enddo
         enddo
      endif
   endif
      
   call update_params(XVAL, NEW_X, IDAT)
   
   !call eval_Newton_dg(XVAL, NMAX, 1d3, NEW_X, DAT, IDAT)
   
   call eval_BANDGAP(XVAL, NEW_X, F0VAL, IDAT, DAT, IERR)
   !call eval_BANDGAP_BOUND(XVAL, NEW_X, F0VAL, IDAT, DAT, IERR)
   !call eval_BANDGAP_VOL(XVAL, NEW_X, F0VAL, IDAT, DAT, IERR)
      
   F = F0VAL/Scaleg0
   
   f0valhist(IDAT(1)) = F
      
   IERR = 0
   return
end subroutine eval_F


!=======================================================================

!                   Computation of objective function
!                             (Mid-bandgap)

!=======================================================================
subroutine eval_BANDGAP(XVAL, NEW_X, F0VAL, IDAT, DAT, IERR)
   implicit none
   integer                              :: NEW_X
   double precision                     :: F0VAL, F0val_undef, F0val_def
   double precision                     :: DAT(:), XVAL(:)
   integer                              :: IDAT(:)
   integer                              :: IERR
   double precision                     :: val_L(nbr_k), val_U(nbr_k), g_L, g_U, g_L_0, g_U_0
   integer                              :: i

   ! No deformation
   call eval_Newton_dg(XVAL, 1, 1d15, NEW_X, DAT, IDAT)
   call eval_EigenProb(NEW_X, DAT, IDAT)
      
   do i=1,nbr_k
      call aggfun(aggregation,val_L(i),evals_L(:,i),nbr_eig_L,pAgg_eig)
      call aggfun(aggregation,val_U(i),evals_U_inv(:,i),nbr_eig_U,pAgg_eig)
   enddo

   call aggfun(aggregation,g_L_0,val_L,nbr_k,pAgg_k)
   call aggfun(aggregation,g_U_0,val_U,nbr_k,pAgg_k)

   F0val_undef = (g_L_0 - 1d0/g_U_0)/(g_L_0 + 1d0/g_U_0)

   val_L = 0d0
   val_U = 0d0
   
   call eval_Newton_dg(XVAL, NMAX, 1d3, NEW_X, DAT, IDAT)
   call eval_EigenProb(NEW_X, DAT, IDAT)
      
   do i=1,nbr_k
      call aggfun(aggregation,val_L(i),evals_L(:,i),nbr_eig_L,pAgg_eig)
      call aggfun(aggregation,val_U(i),evals_U_inv(:,i),nbr_eig_U,pAgg_eig)
   enddo

   call aggfun(aggregation,g_L,val_L,nbr_k,pAgg_k)
   call aggfun(aggregation,g_U,val_U,nbr_k,pAgg_k)

   F0val_def = (g_L - 1d0/g_U)/(g_L + 1d0/g_U)

   F0val = F0val_undef/F0val_def
      
   IERR = 0
   return
end subroutine eval_BANDGAP




!=======================================================================

!                   Computation of objective function
!                   (Mid-bandgap, bound formulation)

!=======================================================================
subroutine eval_BANDGAP_bound(XVAL, NEW_X, F0VAL, IDAT, DAT, IERR)
   implicit none
   integer                              :: NEW_X
   double precision                     :: F0VAL
   double precision                     :: DAT(:), XVAL(:)
   integer                              :: IDAT(:)
   integer                              :: IERR

   ! Minus since max(g0) = -min(-g0)
   F0val = -(bound_U - bound_L)/(bound_U + bound_L)
      
   IERR = 0
   return
end subroutine eval_BANDGAP_bound

!=======================================================================

!                   Computation of objective function
!                        (Mid-bandgap + Volume)

!=======================================================================
subroutine eval_BANDGAP_VOL(XVAL, NEW_X, F0VAL, IDAT, DAT, IERR)
   implicit none
   integer                              :: NEW_X
   double precision                     :: F0VAL
   double precision                     :: DAT(:), XVAL(:)
   integer                              :: IDAT(:)
   integer                              :: IERR
   integer                              :: ie, igp, NegativeE, je
   double precision                     :: rhogp(nbr_phase,ngp), Scalar(ngp)
   double precision                     :: tmp, volCurr,g
 
   volCurr = 0d0
   scalar = 0d0

   !$OMP PARALLEL DO DEFAULT(NONE) &
   !$OMP SHARED(nelm,coord,ed_rho,ngp,enod,Mtot,th,Hprojection,nbr_phase) & 
   !$OMP SHARED(penalization,p,massrho,XVAL,enable_filter) & 
   !$OMP PRIVATE(ie,rhogp,igp,Scalar,tmp,je,g) &
   !$OMP REDUCTION(+:volCurr)
   do ie=1,nelm   
      do je=1,nbr_phase
         if (enable_filter.eq.1) then
            call elem_gp(rhogp(je,:),ed_rho(:,ie,je)) 
         else
            rhogp(je,:) = XVAL((je-1)*nelm + ie) 
         endif
      enddo
      do igp=1,ngp
         call gfun(penalization,g,rhogp(:,igp),1d0,massRho,'rho')
         scalar(igp) = g
         !if (Hprojection.eq.1) then
         !   call Hfun(Scalar(igp),rhogp(:,igp))
         !else
         !   Scalar(igp) = rhogp(:,igp)
         !endif    
      enddo 
      call elem_vol(tmp,coord(:,enod(:,ie)),th,scalar)
      volCurr = volCurr + tmp/Mtot
   enddo
   !$OMP end parallel do

   ! Minus since max(g0) = -min(-g0)
   F0val = -(bound_U - bound_L)/(bound_U + bound_L) + penalty*volcurr
      
   IERR = 0
   return
end subroutine eval_BANDGAP_VOL



!=======================================================================

!               Computation of gradient of objective function

!=======================================================================
subroutine eval_grad_F(NELM, XVAL, NEW_X, GRAD, IDAT, DAT, IERR)
   implicit none
   integer                              :: NELM, NEW_X
   double precision                     :: GRAD(:), XVAL(:)
   double precision                     :: DAT(1)
   integer                              :: IDAT(:)
   integer                              :: IERR, ie, je
   
   if (printout.eq.1) then
      write(*,*) '### Compute sensitivity of objective function ###'
   endif
        
   if (enable_filter.eq.1) then
      call update_PDE_filter(XVAL, NEW_X)
   else
      if (NEW_X.eq.1) then
         do je=1,nbr_phase
            do ie=1,nelm
               rhotilde(ie,je) = XVAL((je-1)*nelm + ie)
            enddo
         enddo
      endif
   endif
         
   call update_params(XVAL, NEW_X, IDAT)

   call eval_Newton_dg(XVAL, NMAX, 1d3, NEW_X, DAT, IDAT)
       
   call eval_grad_BANDGAP(XVAL, NEW_X, GRAD, IDAT, DAT, IERR)
   !call eval_grad_BANDGAP_VOL(XVAL, NEW_X, GRAD, IDAT, DAT, IERR)

   GRAD = GRAD/scaleg0
   
   df0 = GRAD

   IERR = 0
   return
end subroutine eval_grad_F






!=======================================================================

!               Computation of gradient of objective function
!                             (Mid-bandgap)

!=======================================================================
subroutine eval_grad_BANDGAP(XVAL, NEW_X, GRADF0, IDAT, DAT, IERR)
   implicit none
   integer                              :: NEW_X
   double precision                     :: GRADF0(:)
   double precision                     :: DAT(:), XVAL(:)
   integer                              :: IDAT(:)
   integer                              :: IERR

   GRADF0 = 0d0
   GRADF0(ndes_orig-3) = (-1d0)*2d0*bound_L/((bound_U + bound_L)**2d0)
   GRADF0(ndes_orig-2) = (-1d0)*(-2d0)*bound_U/((bound_U + bound_L)**2d0)

   IERR = 0
   return
end subroutine eval_grad_BANDGAP


!=======================================================================

!               Computation of gradient of objective function
!                             (Mid-bandgap + Volume)

!=======================================================================
subroutine eval_grad_BANDGAP_VOL(XVAL, NEW_X, GRADF0, IDAT, DAT, IERR)
   implicit none
   integer                              :: NEW_X
   double precision                     :: GRADF0(:)
   double precision                     :: DAT(:), XVAL(:)
   integer                              :: IDAT(:)
   integer                              :: IERR
   integer                              :: ie, igp, NegativeE, je
   double precision                     :: rhogp(nbr_phase,ngp), Scalar(ngp), fen(ngp), gp
   double precision,allocatable         :: ffilter(:,:), ffilterstar(:,:)
   
   if (enable_filter.eq.1) then
      allocate(ffilter(nnod,nbr_phase), stat=ierr) 
      allocate(ffilterstar(nod_m,nbr_phase), stat=ierr) 
   else
      allocate(ffilter(nelm,nbr_phase), stat=ierr) 
   endif

   scalar = 0d0 
   ffilter=0D0
   GRADF0 = 0d0
   
   !$OMP PARALLEL DO DEFAULT(NONE) &
   !$OMP SHARED(nelm,coord,ed_rho,ngp,enod,th,Hprojection,nbr_phase) & 
   !$OMP SHARED(penalization,p,massrho,XVAL,enable_filter) &
   !$OMP PRIVATE(ie,rhogp,igp,Scalar,fen,je,gp) &
   !$OMP reduction(+:ffilter)
   do ie=1,nelm  
      do je=1,nbr_phase
         if (enable_filter.eq.1) then
            call elem_gp(rhogp(je,:),ed_rho(:,ie,je)) 
         else
            rhogp(je,:) = XVAL((je-1)*nelm + ie) 
         endif

         do igp=1,ngp
            call gfunprim(penalization,gp,rhogp(:,igp),1d0,massRho,'rho',je)
            scalar(igp) = gp
            !if (Hprojection.eq.1) then
            !   call Hprimfun(Scalar(igp),rhogp(:,igp))
            !else
            !   Scalar(igp)=1d0
            !endif   
         enddo

         call elem_flow_bf(fen,coord(:,enod(:,ie)),th,Scalar) 
         if (enable_filter.eq.1) then
            ffilter(enod(:,ie),je) = ffilter(enod(:,ie),je) + fen     
         else
            ffilter(ie,je) = ffilter(ie,je) + sum(fen)
         endif
      enddo
   enddo
   !$OMP end parallel do

   if (enable_filter.eq.1) then
      NegativeE = 0d0 
      do je=1,nbr_phase
         ffilterstar(:,je) = matmulTrans(Projf,ffilter(:,je))
      enddo
      
      ! First obtain af from Kf af = Ff
      call solveq(Kfilterstar,ffilterstar,negativeE,nbr_phase)
      
      do je=1,nbr_phase
         ffilter(:,je) = matmul(Projf,ffilterstar(:,je)) 
      enddo
   endif
  
   do ie=1,nelm
      do je=1,nbr_phase
         if (enable_filter.eq.1) then
            GRADF0(ie+(je-1)*nelm) = penalty*dot_product(Tmatrix(:,ie),ffilter(enod(:,ie),je))/Vtot
         else
            GRADF0(ie+(je-1)*nelm) = penalty*ffilter(ie,je)/Vtot
         endif
      enddo
   enddo 

   GRADF0(ndes_orig-3) = (-1d0)*2d0*bound_L/((bound_U + bound_L)**2d0)
   GRADF0(ndes_orig-2) = (-1d0)*(-2d0)*bound_U/((bound_U + bound_L)**2d0)
   
   deallocate(ffilter)
   if (allocated(ffilterstar)) then
      deallocate(ffilterstar)
   endif

   IERR = 0
   return
end subroutine eval_grad_BANDGAP_VOL


!=======================================================================

!                     Computation of constraints

!=======================================================================

subroutine eval_G(NELM, XVAL, NEW_X, CONST, FVAL, IDAT, DAT, IERR)
   implicit none
   integer                              :: NELM, NEW_X, CONST, I
   double precision                     :: FVAL(:), XVAL(:)
   double precision                     :: DAT(:)
   integer                              :: IDAT(:)
   integer                              :: IERR, ie, je
   double precision                     :: Gvol, Geig(2), Gconnect

   if (printout.eq.1) then
      write(*,*) '### Compute constraint(s) ###'
   endif

   if (enable_filter.eq.1) then
      call update_PDE_filter(XVAL, NEW_X)
   else
      if (NEW_X.eq.1) then
         do je=1,nbr_phase
            do ie=1,nelm
               rhotilde(ie,je) = XVAL((je-1)*nelm + ie)
            enddo
         enddo
      endif
   endif
   
   call update_params(XVAL,NEW_X, IDAT)

   call eval_Newton_dg(XVAL, NMAX, 1d3, NEW_X, DAT, IDAT)

   FVAL = 0d0

   if (C.ge.2) then
      call eval_GEIG(XVAL, NEW_X, Geig, IDAT, DAT, IERR)
      FVAL(1) = Geig(1)/bound_L - 1d0
      FVAL(2) = 1d0 - Geig(2)/bound_U
   endif

   ! sqrt(S11*S22) >= min_gconnect
   if (C.ge.3) then
      call eval_Gconnect(XVAL, NEW_X, Gconnect, IDAT, DAT, IERR)
      FVAL(3) = 1d0 - Gconnect/min_gconnect
   endif
   
   if (C.ge.4) then
      call eval_Gvol(XVAL, NEW_X, Gvol, IDAT, DAT, IERR)
      FVAL(4) = Gvol
   endif

   do i=1,C
      gvalhist(i,IDAT(1)) = FVAL(i)
   enddo

   IERR = 0
   
   return
end subroutine eval_G


!=======================================================================

!                Computation of volume/mass-constraint

!=======================================================================
subroutine eval_GVOL(XVAL, NEW_X, GvolCurr, IDAT, DAT, IERR)
   implicit none
   integer                              :: NEW_X
   double precision                     :: DAT(:), XVAL(:)
   integer                              :: IDAT(:)
   integer                              :: IERR
   integer                              :: ie, igp, NegativeE, je
   double precision                     :: rhogp(nbr_phase,ngp), Scalar(ngp)
   double precision                     :: tmp, GvolCurr, g
 
   if (nbr_phase.gt.1) then
      GvolCurr = -Mtot/Mtot!-Vtot/Vtot
   else
      GvolCurr = -Vtot/Vtot
   endif
   scalar = 0d0

   !$OMP PARALLEL DO DEFAULT(NONE) &
   !$OMP SHARED(nelm,coord,ed_rho,ngp,enod,Mtot,th,Hprojection,nbr_phase) & 
   !$OMP SHARED(penalization,p,massrho,Vtot,XVAL,enable_filter) &
   !$OMP PRIVATE(ie,rhogp,igp,Scalar,tmp,je,g) &
   !$OMP REDUCTION(+:GvolCurr)
   do ie=1,nelm   
      do je=1,nbr_phase
         if (enable_filter.eq.1) then
            call elem_gp(rhogp(je,:),ed_rho(:,ie,je)) 
         else
            rhogp(je,:) = XVAL((je-1)*nelm + ie) 
         endif
      enddo
      do igp=1,ngp
         call gfun(penalization,g,rhogp(:,igp),1d0,massRho,'rho')
         scalar(igp) = g
      enddo 
      call elem_vol(tmp,coord(:,enod(:,ie)),th,scalar)
      if (nbr_phase.gt.1) then
         GvolCurr = GvolCurr + tmp/Mtot
      else
         GvolCurr = GvolCurr + tmp/Vtot
      endif
   enddo
   !$OMP end parallel do

   IERR = 0
   return
end subroutine eval_GVOL


!=======================================================================

!                Computation of eigenvalue-constraint

!=======================================================================
subroutine eval_GEIG(XVAL, NEW_X, GeigCurr, IDAT, DAT, IERR)
   implicit none
   integer                              :: NEW_X
   double precision                     :: DAT(:), XVAL(:)
   integer                              :: IDAT(:)
   integer                              :: IERR
   double precision                     :: GeigCurr(2)
   double precision                     :: val_L(nbr_k), val_U(nbr_k), g_L, g_U
   integer                              :: i

   call eval_EigenProb(NEW_X, DAT, IDAT)
      
   GeigCurr = 0d0

   do i=1,nbr_k
      call aggfun(aggregation,val_L(i),evals_L(:,i),nbr_eig_L,pAgg_eig)
      call aggfun(aggregation,val_U(i),evals_U_inv(:,i),nbr_eig_U,pAgg_eig)
   enddo

   call aggfun(aggregation,g_L,val_L,nbr_k,pAgg_k)
   call aggfun(aggregation,g_U,val_U,nbr_k,pAgg_k)
      
   GeigCurr(1) = g_L!*scaleK
   GeigCurr(2) = 1d0/(g_U)!*scaleK

   IERR = 0
   return
end subroutine eval_GEIG



!=======================================================================

!               Computation of connectivity-constraint

!=======================================================================
subroutine eval_GCONNECT(XVAL, NEW_X, GconnectCurr, IDAT, DAT, IERR)
   implicit none
   integer                              :: NEW_X
   double precision                     :: DAT(:), XVAL(:)
   integer                              :: IDAT(:)
   integer                              :: IERR
   double precision                     :: GconnectCurr

   GconnectCurr = 1d0

   if ((ux_control.eq.1).and.(ux_load.ne.0d0)) then
      GconnectCurr = GconnectCurr*dot_product(Fvec,UX)
   endif
   !if ((uy_control.eq.1).and.(uy_load.ne.0d0)) then
   !   GconnectCurr = GconnectCurr*dot_product(Fvec,UY)
   !endif
   !if ((vx_control.eq.1).and.(vx_load.ne.0d0)) then
   !   GconnectCurr = GconnectCurr*dot_product(Fvec,VX)
   !endif
   if ((vy_control.eq.1).and.(vy_load.ne.0d0)) then
      GconnectCurr = GconnectCurr*dot_product(Fvec,VY)
   endif

   GconnectCurr = dsqrt(GconnectCurr)

   IERR = 0
   return
end subroutine eval_GCONNECT




!=======================================================================

!               Computation of Jacobian of constraints

!=======================================================================

subroutine eval_jac_G(TASK, NELM, XVAL, NEW_X, CONST, NZ, ACON, AVAR, AVEC, IDAT, DAT, IERR)
   implicit none
   integer                              :: TASK, NELM, NEW_X, CONST, NZ, I, J
   double precision                     :: XVAL(:), AVEC(NZ)
   integer                              :: ACON(NZ), AVAR(NZ)
   double precision                     :: DAT(:)
   integer                              :: IDAT(:)
   integer                              :: IERR, ie, je

   if (printout.eq.1) then
      write(*,*) '### Compute sensitivites of constraints ###'
   endif
   
   if (TASK.eq.0) then
      if (optimizer.eq.'IPOPT') write(*,*) 'Task = 0'
      do I = 1, NZ
         ! Col (design variable)
         AVAR(I) = AVAR1(I)
         ! Row (constraint)
         ACON(I) = ACON1(I)
      enddo
   else
      if (optimizer.eq.'IPOPT') write(*,*) 'Task != 0'
      
      if (enable_filter.eq.1) then
         call update_PDE_filter(XVAL, NEW_X)
      else
         if (NEW_X.eq.1) then
            do je=1,nbr_phase
               do ie=1,nelm
                  rhotilde(ie,je) = XVAL((je-1)*nelm + ie)
               enddo
            enddo
         endif
      endif
      
      call update_params(XVAL, NEW_X, IDAT)
     
      call eval_Newton_dg(XVAL, NMAX, 1d3, NEW_X, DAT, IDAT) 
      
      if (C.ge.2) then
         call eval_grad_GEIG(XVAL, NEW_X, IDAT, DAT, IERR)
      endif

      if (C.ge.3) then
         call eval_grad_Gconnect(XVAL, NEW_X, IDAT, DAT, IERR)
      endif
       
      if (C.ge.4) then
         call eval_grad_Gvol(XVAL, NEW_X, IDAT, DAT, IERR)
      endif
      
      df = 0d0
      AVEC = 0d0

      do J = 1, ndes_orig
         do I = 1, 2
            df((I-1)*ndes_orig+J) = dGeig(I,J)
            AVEC((I-1)*ndes_orig+J) = dGeig(I,J)
         enddo
         if (C.ge.3) then
            df((3-1)*ndes_orig+J) = dGconnect(J)
            AVEC((3-1)*ndes_orig+J) = dGconnect(J)
         endif
         if (C.ge.4) then
            df((4-1)*ndes_orig+J) = dGvol(J)
            AVEC((4-1)*ndes_orig+J) = dGvol(J)
         endif
      enddo

   endif
   
   IERR = 0
   return
end subroutine eval_jac_G

      
!=======================================================================

!                Gradient of volume-constraint

!=======================================================================
subroutine eval_grad_GVOL(XVAL, NEW_X, IDAT, DAT, IERR)
   implicit none
   integer                              :: NEW_X
   double precision                     :: DAT(:), XVAL(:)
   integer                              :: IDAT(:)
   integer                              :: IERR, je
   integer                              :: ie, igp, NegativeE
   double precision                     :: rhogp(nbr_phase,ngp), Scalar(ngp), fen(ngp), gp
   double precision,allocatable         :: ffilter(:,:), ffilterstar(:,:)

   if (enable_filter.eq.1) then
      allocate(ffilter(nnod,nbr_phase), stat=ierr) 
      allocate(ffilterstar(nod_m,nbr_phase), stat=ierr) 
   else
      allocate(ffilter(nelm,nbr_phase), stat=ierr) 
   endif

   scalar = 0d0 
   dGvol = 0d0
   ffilter=0D0
   
   !$OMP PARALLEL DO DEFAULT(NONE) &
   !$OMP SHARED(nelm,coord,ed_rho,ngp,enod,th,Hprojection,nbr_phase) & 
   !$OMP SHARED(penalization,p,massrho,XVAL,enable_filter) &
   !$OMP PRIVATE(ie,rhogp,igp,Scalar,fen,je,gp) &
   !$OMP reduction(+:ffilter)
   do ie=1,nelm  
      do je=1,nbr_phase
         if (enable_filter.eq.1) then
            call elem_gp(rhogp(je,:),ed_rho(:,ie,je)) 
         else
            rhogp(je,:) = XVAL((je-1)*nelm + ie) 
         endif
      enddo
      do je=1,nbr_phase
         do igp=1,ngp
            call gfunprim(penalization,gp,rhogp(:,igp),1d0,massRho,'rho',je)
            scalar(igp) = gp
         enddo

         call elem_flow_bf(fen,coord(:,enod(:,ie)),th,Scalar) 
         if (enable_filter.eq.1) then
            ffilter(enod(:,ie),je) = ffilter(enod(:,ie),je) + fen     
         else
            ffilter(ie,je) = ffilter(ie,je) + sum(fen)
         endif
      enddo
   enddo
   !$OMP end parallel do
  
   if (enable_filter.eq.1) then
      NegativeE = 0d0 
      do je=1,nbr_phase
         ffilterstar(:,je) = matmulTrans(Projf,ffilter(:,je))
      enddo
      
      ! First obtain af from Kf af = Ff
      call solveq(Kfilterstar,ffilterstar,negativeE,nbr_phase)
      
      do je=1,nbr_phase
         ffilter(:,je) = matmul(Projf,ffilterstar(:,je)) 
      enddo
   endif
   
   do ie=1,nelm
      if (nbr_phase.gt.1) then
         do je=1,nbr_phase
            if (enable_filter.eq.1) then
               dGvol(ie+(je-1)*nelm) = dot_product(Tmatrix(:,ie),ffilter(enod(:,ie),je))/Mtot
            else 
               dGvol(ie+(je-1)*nelm) = ffilter(ie,je)/Mtot
            endif
         enddo
      else
         if (enable_filter.eq.1) then
            dGvol(ie) = dot_product(Tmatrix(:,ie),ffilter(enod(:,ie),1))/Vtot
         else
            dGvol(ie) = ffilter(ie,1)/Vtot
         endif
      endif
   enddo 

   deallocate(ffilter)
   if (allocated(ffilterstar)) then
      deallocate(ffilterstar)
   endif
 
   IERR = 0
   return
end subroutine eval_grad_GVOL


!=======================================================================

!                    Gradient of eigenvalue

!=======================================================================
subroutine eval_grad_GEIG(XVAL, NEW_X, IDAT, DAT, IERR)
   implicit none
   integer                              :: NEW_X
   double precision                     :: DAT(:), XVAL(:)
   integer                              :: IDAT(:)
   integer                              :: IERR
   integer                              :: ie, NegativeE, i, j, ii, je, jee
   double precision,allocatable         :: ffilter_mult(:,:), RHSstar(:), RHS_re_im(:), RHS_save(:)
   double precision                     :: eig
   double precision,allocatable         :: ffilter(:,:), ffilter_re_im(:,:), ffilterTMP(:,:), RHS(:)
   double precision,allocatable         :: ed_evec(:,:), ed_adjoint(:,:)
   double precision                     :: dfdarg_U(nbr_eig_U), dfdarg_L(nbr_eig_L)
   double precision                     :: dfdarg_tmp(nbr_k)
   double precision,allocatable         :: deig_U(:,:), deig_L(:,:), deig_U_du(:,:), deig_L_du(:,:)
   double precision,allocatable         :: dAgg_L(:,:), dAgg_U(:,:), dAgg_U_du(:,:), dAgg_L_du(:,:)
   double precision,allocatable         :: dg_L(:), dg_U(:), eigvec(:), dg_L_du(:), dg_U_du(:)
   double precision                     :: val_L(nbr_k), val_U(nbr_k), g_U, g_L

   if (enable_filter.eq.1) then
      allocate(ffilter_mult(nod_m,2*nbr_phase),stat=ierr)
      allocate(ffilter(nnod,nbr_phase),stat=ierr)
      allocate(ffilter_re_im(nnod,nbr_phase),stat=ierr)
      allocate(ffilterTMP(nnod,nbr_phase),stat=ierr)

      allocate(deig_U(nbr_eig_U,nbr_phase*nnod),stat=ierr)
      allocate(deig_L(nbr_eig_L,nbr_phase*nnod),stat=ierr)
      allocate(dAgg_L(nbr_k,nbr_phase*nnod),stat=ierr)
      allocate(dAgg_U(nbr_k,nbr_phase*nnod),stat=ierr)
      allocate(dg_L(nbr_phase*nnod),stat=ierr)
      allocate(dg_U(nbr_phase*nnod),stat=ierr)
   else
      allocate(ffilter_mult(nelm,2*nbr_phase),stat=ierr)
      allocate(ffilter(nelm,nbr_phase),stat=ierr)
      allocate(ffilter_re_im(nelm,nbr_phase),stat=ierr)
      allocate(ffilterTMP(nelm,nbr_phase),stat=ierr)
      
      allocate(deig_U(nbr_eig_U,nbr_phase*nelm),stat=ierr)
      allocate(deig_L(nbr_eig_L,nbr_phase*nelm),stat=ierr)
      allocate(dAgg_L(nbr_k,nbr_phase*nelm),stat=ierr)
      allocate(dAgg_U(nbr_k,nbr_phase*nelm),stat=ierr)
      allocate(dg_L(nbr_phase*nelm),stat=ierr)
      allocate(dg_U(nbr_phase*nelm),stat=ierr)
   endif

   allocate(RHS(ndof),stat=ierr)
   allocate(RHS_save(ndof),stat=ierr)
   allocate(RHSstar(dof_m+nbr_free),stat=ierr)
   allocate(RHS_re_im(ndof),stat=ierr)

   allocate(ed_evec(eldof,nelm),stat=ierr)
   allocate(ed_adjoint(eldof,nelm),stat=ierr)
   allocate(eigvec(2*ndof),stat=ierr)

   allocate(deig_U_du(nbr_eig_U,nbr_disp_vars),stat=ierr)
   allocate(deig_L_du(nbr_eig_L,nbr_disp_vars),stat=ierr)
   allocate(dAgg_U_du(nbr_k,nbr_disp_vars),stat=ierr)
   allocate(dAgg_L_du(nbr_k,nbr_disp_vars),stat=ierr)
   allocate(dg_L_du(nbr_disp_vars),stat=ierr)
   allocate(dg_U_du(nbr_disp_vars),stat=ierr)

   call eval_EigenProb(NEW_X, DAT, IDAT)
     
   call extract(ed,a,enod,dofnod)         

   NegativeE = 0d0 

   ! TODO: Consider orphan parallelization?
   do i=1,nbr_k
      if (printout.eq.1) then
         write(*,*) 'Perform sensitivity analysis for k-vector: ', kvec_mat(:,i)
      endif
      do j=1,nbr_eigMax
         !write(*,*) 'Perform sensitivity analysis for eigenvalue: ', j
         
         ! The eigenvalue
         eig = EVALS_TOT(j,i)*scaleK
         ! The eigenvector
         eigvec = EVECS_TOT(j,(ndof*2)*(i-1)+1:(ndof*2)*i)
               
         ! Calculate right hand side of adjoint equation
         RHS = 0d0
         do ii = 1,2
            ! Extract eigenvector v = [Re, Im]
            call extract(ed_evec,eigvec((ii-1)*ndof+1:ii*ndof),enod,dofnod)
            call adjoint(xval,RHS_re_im,ed_evec,p)
            RHS = RHS + RHS_re_im  
         enddo

         ! For sensitivity w.r.t u
         RHS_save = RHS

         ! Solve linear system
         RHSstar = matmulTrans(Proj,RHS)
         call solveq(Kstar,RHSstar,bcnod,bcval,dofnod,NegativeE)  

         ! Obtain true solution
         RHS = matmul(Proj,RHSstar)
         call extract(ed_adjoint,RHS,enod,dofnod)

         ! Now calculate the sensitivity of the eigenvalue
         ffilter = 0d0
         do ii = 1,2
            ! Extract eigenvector v = [Re, Im]
            call extract(ed_evec,eigvec((ii-1)*ndof+1:ii*ndof),enod,dofnod)
            call deigdrho(xval,ffilter_re_im,eig,ed_evec,p) 
            ffilter = ffilter + ffilter_re_im  
         enddo

         ! Calculate the last contribution from dFintdrho
         call dFintdrho(xval,ffilterTMP,ed_adjoint,p)   

         if (j.le.band) then
            do je=1,nbr_phase
               if (enable_filter.eq.1) then
                  dEig_L(j,(je-1)*nnod+1:je*nnod) = ffilter(:,je) - ffilterTMP(:,je)
               else
                  dEig_L(j,(je-1)*nelm+1:je*nelm) = ffilter(:,je) - ffilterTMP(:,je)
               endif
            enddo
            dfdarg_L(j) = 1d0
            ! The sensitivity with respect to u
            dEig_L_du(j,1) = dot_product(UX,RHS_save) - dot_product(RHS,matmul_sym(K,UX))
            dEig_L_du(j,2) = dot_product(VY,RHS_save) - dot_product(RHS,matmul_sym(K,VY))
         elseif (j.gt.band) then
            do je=1,nbr_phase
               if (enable_filter.eq.1) then
                  dEig_U(j-band,(je-1)*nnod+1:je*nnod) = ffilter(:,je) - ffilterTMP(:,je)
               else
                  dEig_U(j-band,(je-1)*nelm+1:je*nelm) = ffilter(:,je) - ffilterTMP(:,je)
               endif
            enddo
            dfdarg_U(j-band) = -(1d0/(evals_U(j-band,i)**2d0))
            ! The sensitivity with respect to u
            dEig_U_du(j-band,1) = dot_product(UX,RHS_save) - dot_product(RHS,matmul_sym(K,UX))
            dEig_U_du(j-band,2) = dot_product(VY,RHS_save) - dot_product(RHS,matmul_sym(K,VY))
         endif
      enddo
      ! Obtain first part of chainrule for sensitivity
      ! Call twice to obtain for both constraint functions
      ! L
      call Aggprim(aggregation,dAgg_L(i,:),evals_L(:,i),dfdarg_L,dEig_L,nbr_eig_L,pAgg_eig)
      ! U
      call Aggprim(aggregation,dAgg_U(i,:),evals_U_inv(:,i),dfdarg_U,dEig_U,nbr_eig_U,pAgg_eig)
      ! Displacement
      ! L
      call Aggprim(aggregation,dAgg_L_du(i,:),evals_L(:,i),dfdarg_L,dEig_L_du,nbr_eig_L,pAgg_eig)
      ! U
      call Aggprim(aggregation,dAgg_U_du(i,:),evals_U_inv(:,i),dfdarg_U,dEig_U_du,nbr_eig_U,pAgg_eig)
   enddo
  
   ! Second part of chain rule
   do i=1,nbr_k
      call aggfun(aggregation,val_L(i),evals_L(:,i),nbr_eig_L,pAgg_eig)
      call aggfun(aggregation,val_U(i),evals_U_inv(:,i),nbr_eig_U,pAgg_eig)
   enddo
   ! L
   dfdarg_tmp = 1d0
   call Aggprim(aggregation,dg_L,val_L,dfdarg_tmp,dAgg_L,nbr_k,pAgg_k)
   ! U
   call Aggprim(aggregation,dg_U,val_U,dfdarg_tmp,dAgg_U,nbr_k,pAgg_k)
   ! Displacement
   ! L
   dfdarg_tmp = 1d0
   call Aggprim(aggregation,dg_L_du,val_L,dfdarg_tmp,dAgg_L_du,nbr_k,pAgg_k)
   ! U
   call Aggprim(aggregation,dg_U_du,val_U,dfdarg_tmp,dAgg_U_du,nbr_k,pAgg_k)
   
   call aggfun(aggregation,g_L,val_L,nbr_k,pAgg_k)
   call aggfun(aggregation,g_U,val_U,nbr_k,pAgg_k)

   jee = 0
   do je=1,2*nbr_phase
      if (je.le.nbr_phase) then
         if (enable_filter.eq.1) then
            ffilter_mult(:,je) = matmulTrans(Projf,dg_L((je-1)*nnod+1:je*nnod))
         else
            ffilter_mult(:,je) = dg_L((je-1)*nelm+1:je*nelm)
         endif
      else
         jee = jee + 1
         if (enable_filter.eq.1) then
            ffilter_mult(:,je) = matmulTrans(Projf,dg_U((jee-1)*nnod+1:jee*nnod))
         else
            ffilter_mult(:,je) = dg_U((jee-1)*nelm+1:jee*nelm)
         endif
      endif
   enddo

   if (enable_filter.eq.1) then
      ! First obtain af from Kf af = Ff
      call solveq(Kfilterstar,ffilter_mult,NegativeE,2*nbr_phase) 
   endif

   ! Then obtain deigdz by scalarproduct with T  
   do jee=1,2
      do je=1,nbr_phase
         if (enable_filter.eq.1) then
            ffilter(:,je) = matmul(Projf,ffilter_mult(:,je+(jee-1)*nbr_phase)) 
         else
            ffilter(:,je) = ffilter_mult(:,je+(jee-1)*nbr_phase)
         endif
         do ie=1,nelm
            if (jee.eq.1) then
               if (enable_filter.eq.1) then
                  dGeig(1,ie+(je-1)*nelm) = (1d0/bound_L/scaleK)*dot_product(Tmatrix(:,ie),ffilter(enod(:,ie),je))
               else
                  dGeig(1,ie+(je-1)*nelm) = (1d0/bound_L/scaleK)*ffilter(ie,je)
               endif
            else
               if (enable_filter.eq.1) then
                  dGeig(2,ie+(je-1)*nelm) = ((1d0/(g_U**2d0))/bound_U/scaleK)*dot_product(Tmatrix(:,ie),ffilter(enod(:,ie),je))
               else
                  dGeig(2,ie+(je-1)*nelm) = ((1d0/(g_U**2d0))/bound_U/scaleK)*ffilter(ie,je)
               endif
            endif
         enddo
      enddo
   enddo

   ! Bound_U
   !dGeig(1,ndes_orig-3) = 0d0
   !dGeig(2,ndes_orig-3) = (1d0/g_U)*(1d0/(bound_U)**2d0)!*scaleK!1d0   
   ! Bound_L
   !dGeig(1,ndes_orig-2) = g_L*(-1d0/(bound_L)**2d0)!*scaleK!-1d0 
   !dGeig(2,ndes_orig-2) = 0d0
   ! Displacement variables
   ! ux_load
   !dGeig(1,ndes_orig-1) = (1d0/bound_L/scaleK)*dg_L_du(1)
   !dGeig(2,ndes_orig-1) = ((1d0/(g_U**2d0))/bound_U/scaleK)*dg_U_du(1)
   ! vy_load
   !dGeig(1,ndes_orig) = (1d0/bound_L/scaleK)*dg_L_du(2)
   !dGeig(2,ndes_orig) = ((1d0/(g_U**2d0))/bound_U/scaleK)*dg_U_du(2)
   
   IERR = 0
   return
end subroutine eval_grad_GEIG


!=======================================================================

!           Computation of gradient of connectivity constraint

!=======================================================================
subroutine eval_grad_GCONNECT(XVAL, NEW_X, IDAT, DAT, IERR)
   implicit none
   integer                              :: NEW_X
   double precision                     :: DAT(:), XVAL(:)
   integer                              :: IDAT(:)
   integer                              :: IERR
   integer                              :: ie, NegativeE, je
   double precision,allocatable         :: ffilter(:,:), ffilterstar(:,:), RHSstar(:), RHS(:), RHS_save(:)
   double precision,allocatable         :: ed_adjoint(:,:)
   double precision                     :: S11, S22

   if (enable_filter.eq.1) then
      allocate(ffilter(nnod,nbr_phase),stat=ierr)
      allocate(ffilterstar(nod_m,nbr_phase),stat=ierr)
   else
      allocate(ffilter(nelm,nbr_phase),stat=ierr)
   endif
   allocate(RHSstar(dof_m+nbr_free),stat=ierr)
   allocate(RHS(ndof),stat=ierr)
   allocate(RHS_save(ndof),stat=ierr)
   allocate(ed_adjoint(eldof,nelm),stat=ierr)

   call eval_Newton_dg(XVAL, NMAX, 1d3, NEW_X, dat, idat)
   
   call extract(ed,a,enod,dofnod)

   S11 = 0d0
   S22 = 0d0
 
   if ((ux_control.eq.1).and.(ux_load.ne.0d0)) then
      S11 = dot_product(Fvec,UX)
   endif
   if ((vy_control.eq.1).and.(vy_load.ne.0d0)) then
      S22 = dot_product(Fvec,VY)
   endif

   ! Obtain adjoint
   RHSstar = (-1d0/(2d0*dsqrt(S11*S22)))*(S11*matmulTrans(Proj,matmul_sym(K,VY)) + S22*matmulTrans(Proj,matmul_sym(K,UX)))
   call solveq(Kstar,RHSstar,bcnod,bcval,dofnod,NegativeE)  

   ! Obtain true solution
   RHS = matmul(Proj,RHSstar)
   RHS_save = RHS
   RHS = RHS + (1d0/(2d0*dsqrt(S11*S22)))*(S11*VY + S22*UX)
   call extract(ed_adjoint,RHS,enod,dofnod)  
   
   ffilter = 0d0
   call dFintdrho(xval,ffilter,ed_adjoint,p)   

   if (enable_filter.eq.1) then
      do je=1,nbr_phase
         ffilterstar(:,je) = matmulTrans(Projf,ffilter(:,je))
      enddo
      
      ! First obtain af from Kf af = Ff
      call solveq(Kfilterstar,ffilterstar,negativeE,nbr_phase)
      
      do je=1,nbr_phase
         ffilter(:,je) = matmul(Projf,ffilterstar(:,je)) 
      enddo
   endif
   
   ! Scalarproduct with T    
   dGconnect = 0d0
   do ie=1,nelm
      do je=1,nbr_phase
         if (enable_filter.eq.1) then
            dGconnect(ie+(je-1)*nelm) =-dot_product(Tmatrix(:,ie),ffilter(enod(:,ie),je))/min_gconnect
         else
            dGconnect(ie+(je-1)*nelm) =-ffilter(ie,je)/min_gconnect
         endif
      enddo
   enddo

   ! ux_load
   !dGconnect(ndes_orig-1) = -((1d0/(2d0*dsqrt(S11*S22)))*(S11*dot_product(VY,matmul_sym(K,UX))+S22*dot_product(UX,matmul_sym(K,UX))) &
   !   + dot_product(RHS_save,matmul_sym(K,UX)))/min_gconnect
   ! vy_load
   !dGconnect(ndes_orig) = -((1d0/(2d0*dsqrt(S11*S22)))*(S11*dot_product(VY,matmul_sym(K,VY))+S22*dot_product(UX,matmul_sym(K,VY))) &
   !   + dot_product(RHS_save,matmul_sym(K,VY)))/min_gconnect


   IERR = 0
   return
end subroutine eval_grad_Gconnect




!=======================================================================

!               Computation of Hessian of Lagrangian

!=======================================================================

subroutine eval_HESS(TASK, NELM, XVAL, NEW_X, OBJFACT, CONST, LAM, NEW_LAM, NNZH, IRNH, ICNH, HESS, IDAT, DAT, IERR)
   implicit none
   integer TASK, NELM, NEW_X, CONST, NEW_LAM, NNZH, i
   double precision XVAL(:), OBJFACT, LAM(CONST), HESS(NNZH)
   integer IRNH(NNZH), ICNH(NNZH)
   double precision DAT(:)
   integer IDAT(:)
   integer IERR
   
   write(*,*) 'Hessians are not supported at this point'
  
   IERR = 1
   return
end subroutine eval_HESS


!=======================================================================

!                  Callback method called once per iteration

!=======================================================================

subroutine iter_CB(ALG_MODE, ITER_COUNT, OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT, DAT, ISTOP)
   implicit none
   integer                   :: ALG_MODE, ITER_COUNT, LS_TRIAL
   double precision          :: OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE
   double precision          :: ALPHA_DU, ALPHA_PR
   double precision          :: DAT(:)
   integer                   :: IDAT(:)
   integer                   :: ISTOP
   integer                   :: indx1
   character(80)             :: s_plot, load_str
   character(120)            :: xstring
   
   update_done = 0

   if (optimizer.ne.'MMA') then
      ! Update iteration count
      IDAT(1) = ITER_COUNT+1
      write(*,*) 'Callback in iteration ', ITER_COUNT
   endif
   xstring='data'
   indx1=4

   if ((irestart.eq.1).or.(diffCheck.eq.1)) then
      xstring='restart'
      indx1=7
   endif

   tmp2(3) = dble(IDAT(1))
   tmp2(4) = dble(C)
   tmp2(6) = nbr_phase
   tmp2(8) = dble(nbr_k)
  
   if (refinement_iter.eq.0) then
      if (opt_finished.eq.1) then
         xstring='data_coarse_final'
         indx1=17
      else
         xstring='data_coarse'
         indx1=11
      endif
   else
      if (opt_finished.eq.1) then
         if (loadstep_bandgap.eq.0) then
            xstring='data_final'
            indx1=10
         else
            write(load_str, '(i10)' )  loadstep_bandgap
            xstring='data_loadstep_'
            indx1=14
            if (loadstep_bandgap.lt.10) then
               xstring=xstring(1:indx1)//load_str(10:10) //''' )'
               indx1 = 15
            endif
            if (loadstep_bandgap.ge.10) then
               xstring=xstring(1:indx1)//load_str(9:10) //''' )'
               indx1 = 16
            endif
         endif
      endif 
   endif


   if ((mod(IDAT(1),save_itr).eq.0).or.(opt_finished.eq.1)) then 
      write( s_plot, '(i10)' )  ITER_COUNT
              
      if (ITER_COUNT.lt.10) xstring=xstring(1:indx1)//'_000' //s_plot(10:10) //''' )'
      if (ITER_COUNT.ge.10) xstring=xstring(1:indx1)//'_00' //s_plot(9:10) //''' )'
      if (ITER_COUNT.ge.100) xstring=xstring(1:indx1)//'_0' //s_plot(8:10) //''' )'
      if (ITER_COUNT.ge.1000) xstring=xstring(1:indx1)//'_'//s_plot(7:10) //''' )'
      indx1=indx1+5
      xstring=xstring(1:indx1)//'.mat'
      indx1=indx1+4
 
      call matwrt2f(xstring(1:indx1),a,'a','w')  
      call matwrt2f(xstring(1:indx1),EVECS_TOT,'evecs','u')  
      call matwrt2f(xstring(1:indx1),EVALS_TOT,'evals','u')  
      call matwrt2f(xstring(1:indx1),df0,'df0dx','u')  
      call matwrt2f(xstring(1:indx1),df,'dfdx','u')
      call matwrt2f(xstring(1:indx1),rhotilde,'rhotilde','u') 
      call matwrt2f(xstring(1:indx1),rho,'rho','u')  
      call matwrt2f(xstring(1:indx1),tmp2,'tmp2','u') 
      call matwrt2f(xstring(1:indx1),kvec_mat,'kvec_mat','u') 
              
      
      xstring='hist'
      indx1=4
      if (refinement_iter.eq.0) then
         if (opt_finished.eq.1) then
            xstring='hist_coarse_final'
            indx1=17
         else
            xstring='hist_coarse'
            indx1=11
         endif
      else
         if (opt_finished.eq.1) then
            xstring='hist_final'
            indx1=10
         endif 
      endif

      if ((irestart.eq.1).or.(diffCheck.eq.1)) then
         xstring='restart_hist'
         indx1=12
      endif
              
      write( s_plot, '(i10)' )  ITER_COUNT
              
      if (ITER_COUNT.lt.10) xstring=xstring(1:indx1)//'_000' //s_plot(10:10) //''' )'
      if (ITER_COUNT.ge.10) xstring=xstring(1:indx1)//'_00' //s_plot(9:10) //''' )'
      if (ITER_COUNT.ge.100) xstring=xstring(1:indx1)//'_0' //s_plot(8:10) //''' )'
      if (ITER_COUNT.ge.1000) xstring=xstring(1:indx1)//'_'//s_plot(7:10) //''' )'
      indx1=indx1+5
      xstring=xstring(1:indx1)//'.mat'
      indx1=indx1+4

      call matwrt2f(xstring(1:indx1),history,'history','w')
      call matwrt2f(xstring(1:indx1),ahist,'ahist','u')  
      call matwrt2f(xstring(1:indx1),ares_hist,'ares_hist','u')  
      call matwrt2f(xstring(1:indx1),f0valhist,'f0val','u') 
      call matwrt2f(xstring(1:indx1),gvalhist,'gval','u')  
   endif   
   return
end subroutine iter_CB


!=======================================================================

!                       Update parameters

!=======================================================================
subroutine update_params(XVAL,NEW_X,IDAT)
   implicit none
   integer            :: NEW_X, IDAT(:), iter
   double precision   :: XVAL(:)

   iter = IDAT(1)

   pold = p
   Bold = B

   if ((update_done.eq.0).and.(NEW_X.eq.1)) then
      if (printout.eq.1) then
         write(*,*) '### Update parameters ###'
      endif
      update_done = 1
      
      ! Update penalization parameter
      if (mod(iter,p_iter).eq.0) then
         p = p + pIncr   
      endif
     
      if (p.gt.pfinal) p=pfinal
      ! Update beta in Heaviside
      if (Hprojection.eq.1) then
          ! Wait a few iterations
          if (dabs(p-pfinal).lt.1D-8) then
            iwait=iwait+1
          endif
          ! Update sharpness of Heaviside filter
          if (iwait.ge.StartH) then
             if ((mod(iter,B_iter).eq.0)) B = 1.3d0*B !B+1D0
          endif
          if (B.gt.Bfinal) B=Bfinal
          if (p.lt.pfinal) B=Binit
     endif
     
     tmp2(1)=p
     tmp2(2)=B
     tmp2(5)=iwait

     tmp2(8)=nbr_k
     tmp2(9)=scaleK
     call update_gamma(XVAL)

     if (printout.eq.1) then
        write(*,*) '####################################################'            
        write(*,*)
        write(*,*) '                 Parameter values                   '
        write(*,*)
        write(*,*) '####################################################'      
        write(*,*) 'Penalty exponent', p
        write(*,*) 'Heaviside coefficient', B 
        write(*,*) '####################################################'
     endif
   endif

   
   return
end subroutine update_params


!=======================================================================

!                  Calculates dfint/drho

!=======================================================================
subroutine dFintdrho(XVAL,ffilterCurr,ed_adjointCurr,pcurr)   
   implicit none
   double precision                 :: ffilterCurr(:,:), ed_adjointCurr(:,:), XVAL(:)
   double precision                 :: CRAP98(eldef,ngp), Ones(ngp), Dgp(elm,elm,ngp), gp
   double precision                 :: esscaled(els,ngp), Bconst(elm,ngp), esscaledsmall(els,ngp)
   double precision                 :: fke(elnod), g, ffe(elnod)  
   double precision                 :: fen(elnod), dggamma(eldef,ngp), esTMP(els,ngp)
   integer                          :: ie, igp, je, ierr
   double precision                 :: gammap, rhogp(nbr_phase,ngp), pcurr
   double precision,allocatable     :: es(:,:,:)
  
   allocate(es(els,ngp,nelm),stat=ierr)

   Ones=1D0   
   ffiltercurr = 0d0
   
   !$OMP PARALLEL DO DEFAULT(NONE) &
   !$OMP SHARED(nelm,ed_rho,ed,enod,coord,gammagp,nbr_phase,E,es,pcurr,ones,lagrangian_type) & 
   !$OMP SHARED(ngp,dg,ed_adjointCurr,dofnod,th,elm,penalization,els,Dsmall,XVAL,enable_filter) & 
   !$OMP PRIVATE(fke,ffe,fen,dggamma,crap98,Dgp,je,esscaled,gp) &
   !$OMP PRIVATE(ie,igp,g,esTMP,esscaledsmall,rhogp,Bconst,gammap) &
   !$OMP REDUCTION(+:ffilterCurr) 
   do ie=1,nelm
      ! Extract gauss point values from nodal quantities
      do je=1,nbr_phase
         if (enable_filter.eq.1) then
            call elem_gp(rhogp(je,:),ed_rho(:,ie,je)) 
         else
            rhogp(je,:) = XVAL((je-1)*nelm + ie) 
         endif
      enddo

      ! The deformation gradient 
      call elem_dg(dg(:,:,ie),coord(:,enod(:,ie)),ed(:,ie))
                 
      ! The gamma-scaled deformation gradient 
      call calcdggamma(dggamma,dg(:,:,ie),gammagp(:,ie))
     
      ! Calculate the tangent stiffness with gamma
      call dneohooke('tl',Dgp,dggamma)
            
      ! Calculate the 2nd Piola Kirchhoff stress with gamma
      call neohooke('2ndPiola',es(:,:,ie),dggamma)

      do je=1,nbr_phase
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !                   NONLINEAR STIFFNESS MATRIX CONTRIBUTION
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! The material behaviour
         do igp=1,ngp
            call gfun(penalization,g,rhogp(:,igp),pcurr,E,'E')
            call gammaprimfun(gammap,rhogp(:,igp),pcurr,je)
            
            esscaled(:,igp) = g*es(:,igp,ie)*gammap*gammagp(igp,ie)
            Dgp(:,:,igp) = g*Dgp(:,:,igp)*gammap*gammagp(igp,ie)
         enddo
                 
         ! Calculate the large def. element stiffness matrix contribution
         call elem_stiff_gamma_Sens('tl',fke,coord(:,enod(:,ie)),th,Dgp,ed(:,ie),ed_adjointCurr(:,ie),ed(:,ie),esscaled,gammagp(:,ie))   

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !                    NONLINEAR RESIDUAL CONTRIBUTION
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Thresholding
         do igp=1,ngp
            call gfun(penalization,g,rhogp(:,igp),pcurr,E,'E')  
            call gfunprim(penalization,gp,rhogp(:,igp),pcurr,E,'E',je)
            call gammaprimfun(gammap,rhogp(:,igp),pcurr,je)
            ! Division by g done
            esscaled(:,igp)=es(:,igp,ie)*(gammagp(igp,ie)*gp + g*gammap)    
         enddo
         
         ! Calculate the nonlinear element forces
         call elem_force_gamma_sens('tl',ffe,coord(:,enod(:,ie)),th,ed(:,ie),ed_adjointCurr(:,ie),esscaled,gammagp(:,ie)) 

         ! Add the contributions
         fen = fke + ffe
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !                      LINEAR RESIDUAL CONTRIBUTION
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculate Bconst for small strain stress
         call elem_Etilde1(coord(:,enod(:,ie)),0D0*ed(:,ie),ed(:,ie),Bconst,Crap98,Ones) 
                    
         ! Calculate small strain contribution
         esscaledsmall = 0d0
         esTMP = 0d0
         do igp=1,ngp
            call gfun(penalization,g,rhogp(:,igp),pcurr,E,'E')
            call gfunprim(penalization,gp,rhogp(:,igp),pcurr,E,'E',je)
            call gammaprimfun(gammap,rhogp(:,igp),pcurr,je)

            ! Calculate sigma (Cauchy stress)
            if (els.eq.4) then
               esscaledsmall([1,2,4],igp) = matmul(Dsmall,Bconst(:,igp))
               ! Division by g done
               esTMP([1,2,4],igp) = ((1D0-gammagp(igp,ie)**2D0)*gp)*esscaledsmall([1,2,4],igp)
               esTMP([1,2,4],igp) = esTMP([1,2,4],igp) + g*(-2d0*(gammagp(igp,ie))*gammap)*esscaledsmall([1,2,4],igp)
            else
               esscaledsmall(:,igp) = matmul(Dsmall,Bconst(:,igp))
               esTMP(:,igp) = ((1D0-gammagp(igp,ie)**2D0)*gp)*esscaledsmall(:,igp)
               esTMP(:,igp) = esTMP(:,igp) + g*(-2d0*(gammagp(igp,ie))*gammap)*esscaledsmall(:,igp)
            endif 
         enddo

         esscaledsmall = esTmp

         ! Calculate the linear element forces 
         call elem_force_gamma_sens('tl',ffe,coord(:,enod(:,ie)),th,0D0*ed(:,ie),ed_adjointCurr(:,ie),esscaledsmall,ones) 
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         ! Add the contributions
         fen = fen + ffe

         ! Assemble. ffilter is dfdrhotilde
         if (enable_filter.eq.1) then
            ffilterCurr(enod(:,ie),je) = ffilterCurr(enod(:,ie),je) + fen           
         else
            ffilterCurr(ie,je) = ffilterCurr(ie,je) + sum(fen)
         endif
      enddo
   enddo
   !$OMP end parallel do

   return
end subroutine dFintdrho





!=======================================================================

!     Calculate adjoint in sensitivity for stability constraint
!                !!! ed_evec = [evec_re; evec_im] !!!
!                       !!! K = [K 0; 0 K]; !!!

!=======================================================================
subroutine adjoint(XVAL,RHSCurr,ed_evec,pcurr)
   implicit none
   double precision                 :: RHSCurr(:), ed_evec(:,:), XVAL(:)
   double precision                 :: rhogp(nbr_phase,ngp), dggamma(eldef,ngp), pcurr
   double precision                 :: CRAP98(eldef,ngp), GradF(elm,ngp)
   double precision                 :: Dgp(elm,elm,ngp), Bl(elm,ngp), arg1(elm,ngp)
   double precision                 :: arg2(elm,ngp), g, gtilde, rhse(eldof), arg(eldof)
   double precision                 :: E_a_evec(elm,ngp), E_evec_evec(elm,ngp), E_0_evec(elm,ngp), arg1TMP(els,ngp)
   double precision                 :: arg2TMP(els,ngp)
   integer                          :: ie, igp, ierr, je

   RHScurr = 0d0

   !$OMP PARALLEL DO DEFAULT(NONE) &
   !$OMP SHARED(nelm,ed_rho,ed,enod,coord,gammagp,nbr_phase,E,lagrangian_type) & 
   !$OMP SHARED(ngp,dg,pcurr,ed_evec,dofnod,th,elm,penalization,xval,enable_filter) & 
   !$OMP PRIVATE(E_a_evec,E_0_evec,E_evec_evec,GRADF,Bl,arg2,arg1) &
   !$OMP PRIVATE(dggamma,crap98,Dgp,arg,je) &
   !$OMP PRIVATE(ie,igp,g,rhse,rhogp,arg1TMP,arg2TMP) &
   !$OMP REDUCTION(+:RHScurr) 
   ! Calculate the adjoint force vector
   do ie=1,nelm   
      ! Extract gauss point values from nodal quantities
      do je=1,nbr_phase
         if (enable_filter.eq.1) then
            call elem_gp(rhogp(je,:),ed_rho(:,ie,je)) 
         else
            rhogp(je,:) = XVAL((je-1)*nelm + ie) 
         endif
      enddo
      ! Calculate the different E-tilde
      ! E(a,evec) 
      call elem_Etilde1(coord(:,enod(:,ie)),ed(:,ie),ed_evec(:,ie),E_a_evec,crap98,gammagp(:,ie)) 
      ! E(0,evec)
      call elem_Etilde1(coord(:,enod(:,ie)),0D0*ed(:,ie),ed_evec(:,ie),E_0_evec,crap98,gammagp(:,ie))
      ! E(evec,evec)
      call elem_Etilde1(coord(:,enod(:,ie)),ed_evec(:,ie),ed_evec(:,ie),E_evec_evec,crap98,gammagp(:,ie))  

      ! The deformation gradient 
      call elem_dg(dg(:,:,ie),coord(:,enod(:,ie)),ed(:,ie))
      
      ! The gamma-scaled deformation gradient
      call calcdggamma(dggamma,dg(:,:,ie),gammagp(:,ie)) 

      ! Calculate nabla(f_AT)
      call dfat(GradF,dggamma,E_a_evec,E_a_evec)            
      
      ! Calculate the tangent 
      call dneohooke('tl',Dgp,dggamma)

      ! Introduce quantity
      Bl=E_evec_evec-E_0_evec  

      ! Thresholding               
      do igp=1,ngp
          call gfun(penalization,g,rhogp(:,igp),pcurr,E,'E') 
          ! Note that 2D0 is replaced by 1D0 since the factor 2D0 is inside sub. dfat   
          arg1(:,igp) = (matmul(Dgp(:,:,igp),Bl(:,igp)) + GradF(:,igp))*gammagp(igp,ie)*g
          arg2(:,igp) = 2d0*matmul(Dgp(:,:,igp),E_a_evec(:,igp))*gammagp(igp,ie)*g
      enddo

      ! First term in right hand side of adjoint linear system
      if (elm.eq.3) then
         arg1TMP = 0d0
         arg1TMP([1,2,4],:) = arg1
         call elem_force_gamma('tl',rhse,coord(:,enod(:,ie)),th,ed(:,ie),arg1TMP,gammagp(:,ie))
      else
         call elem_force_gamma('tl',rhse,coord(:,enod(:,ie)),th,ed(:,ie),arg1,gammagp(:,ie))
      endif

      ! The second term is split in two and then added to he
      if (elm.eq.3) then
         arg2TMP = 0d0
         arg2TMP([1,2,4],:) = arg2
         call elem_force_gamma('tl',arg,coord(:,enod(:,ie)),th,ed_evec(:,ie),arg2TMP,gammagp(:,ie)) 
      else
         call elem_force_gamma('tl',arg,coord(:,enod(:,ie)),th,ed_evec(:,ie),arg2,gammagp(:,ie)) 
      endif

      rhse=rhse+arg

      if (elm.eq.3) then
         arg2TMP = 0d0
         arg2TMP([1,2,4],:) = arg2
         call elem_force_gamma('tl',arg,coord(:,enod(:,ie)),th,0D0*ed_evec(:,ie),arg2TMP,gammagp(:,ie)) 
      else
         call elem_force_gamma('tl',arg,coord(:,enod(:,ie)),th,0D0*ed_evec(:,ie),arg2,gammagp(:,ie)) 
      endif

      rhse=rhse-arg

      ! Instert in global h
      call insert(RHScurr,rhse,enod(:,ie),dofnod)  
   enddo  
   !$OMP end parallel do

   return
end subroutine adjoint










!=======================================================================

!                    Calculates deig/drho

!=======================================================================
subroutine deigdrho(XVAL,ffilterCurr,eigcurr,ed_evec,pcurr)   
   implicit none
   double precision                 :: ffilterCurr(:,:), eigcurr, ed_evec(:,:), XVAL(:)
   double precision                 :: rhogp(nbr_phase,ngp), pcurr
   double precision                 :: fke(elnod), fme(elnod), fkesmall(elnod), Bl(elm,ngp), Bl_evec_ed(elm,ngp)
   double precision                 :: valf, dggamma(eldef,ngp), g, gtilde 
   double precision                 :: CRAP98(eldef,ngp), Ones(ngp), Dgp(elm,elm,ngp), gp, gptilde
   double precision                 :: esscaled(els,ngp), GradF(elm,ngp), fen(elnod)
   double precision                 :: E_ed_evec(elm,ngp), E_evec_evec(elm,ngp), AT(elm,ngp), E_ed_ed(elm,ngp)
   double precision                 :: Scalar(ngp), E_0_evec(elm,ngp)
   integer                          :: ie, igp, ierr, je
   double precision                 :: gammap, E_0_ed(elm,ngp), E_evec_ed(elm,ngp), E_ed_evecgamma(elm,ngp)
   double precision,allocatable     :: es(:,:,:)
  
   allocate(es(els,ngp,nelm),stat=ierr)

   Ones=1D0
   ffiltercurr = 0d0

   !$OMP PARALLEL DO DEFAULT(NONE) &
   !$OMP SHARED(nelm,ed_rho,ed,enod,coord,gammagp,th,penalization,hprojection,lagrangian_type) & 
   !$OMP SHARED(ngp,dg,es,pcurr,eigcurr,ed_evec,Dsmall,Ones,massrho,nbr_phase,E,XVAL,enable_filter) & 
   !$OMP PRIVATE(Bl,Bl_evec_ed,E_ed_evec,E_evec_evec,E_ed_ed,Scalar,E_ed_evecgamma,GradF) &
   !$OMP PRIVATE(E_0_ed,E_evec_ed,E_0_evec,dggamma,crap98) &
   !$OMP PRIVATE(ie,igp,g,gp,gammap,Dgp,esscaled) &
   !$OMP PRIVATE(rhogp,fke,fme,fkesmall,fen,je) &
   !$OMP REDUCTION(+:ffiltercurr) 
   do ie=1,nelm
      ! Extract gauss point values from nodal quantities
      do je=1,nbr_phase
         if (enable_filter.eq.1) then
            call elem_gp(rhogp(je,:),ed_rho(:,ie,je)) 
         else
            rhogp(je,:) = XVAL((je-1)*nelm + ie) 
         endif
      enddo

      ! The deformation gradient 
      call elem_dg(dg(:,:,ie),coord(:,enod(:,ie)),ed(:,ie))
      
      ! The gamma-scaled deformation gradient
      call calcdggamma(dggamma,dg(:,:,ie),gammagp(:,ie)) 

      ! Calculate the tangent stiffness 
      call dneohooke('tl',Dgp,dggamma) 

      ! Calculate the 2nd Piola Kirchhoff stress
      call neohooke('2ndPiola',es(:,:,ie),dggamma) 
      
      do je=1,nbr_phase
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !              NONLINEAR STIFFNESS MATRIX CONTRIBUTION 
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! The material behaviour
         do igp=1,ngp
            call gfun(penalization,g,rhogp(:,igp),pcurr,E,'E')
            call gfunprim(penalization,gp,rhogp(:,igp),pcurr,E,'E',je) 
            call gammaprimfun(gammap,rhogp(:,igp),pcurr,je)
            
            ! Division by g done
            esscaled(:,igp)= es(:,igp,ie)*(g*2d0*gammap*gammagp(igp,ie) + gammagp(igp,ie)**2d0*gp)  
            Dgp(:,:,igp)   = Dgp(:,:,igp)*(g*2d0*gammap*gammagp(igp,ie) + gammagp(igp,ie)**2d0*gp)
         enddo
                    
         call elem_stiff_gamma_sens('tl',fke,coord(:,enod(:,ie)),th,Dgp,ed(:,ie),ed_evec(:,ie),ed_evec(:,ie),esscaled,gammagp(:,ie))   
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !                        MASS MATRIX CONTRIBUTION
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! The material behaviour
         do igp=1,ngp 
            call gfun(penalization,g,rhogp(:,igp),1d0,massrho,'rho')
            call gfunprim(penalization,gp,rhogp(:,igp),1d0,massrho,'rho',je) 
            !call gammaprimfun(gammap,rhogp(:,igp),pcurr,je)

            !scalar(igp) = g*2d0*gammagp(igp,ie)*gammap + gp*gammagp(igp,ie)**2d0
            scalar(igp) = gp
         enddo
                    
         if (nbr_phase.gt.1) then
            call elem_mass_sens(fme,coord(:,enod(:,ie)),th,1d0,Scalar,ed_evec(:,ie),ed_evec(:,ie))   
         else
            call elem_mass_sens(fme,coord(:,enod(:,ie)),th,massrho(1),Scalar,ed_evec(:,ie),ed_evec(:,ie))   
         endif
         
         fme = -1d0*eigcurr*fme
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !                   LINEAR STIFFNESS MATRIX CONTRIBUTION
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         Dgp = 0d0
         ! The material behaviour
         do igp=1,ngp 
              call gfun(penalization,g,rhogp(:,igp),pcurr,E,'E')
              call gfunprim(penalization,gp,rhogp(:,igp),pcurr,E,'E',je)
              call gammaprimfun(gammap,rhogp(:,igp),pcurr,je)
             
              ! Division by g done 
              Dgp(:,:,igp)=((1D0-gammagp(igp,ie)**2D0)*gp - g*2d0*(gammagp(igp,ie))*gammap)*Dsmall
         enddo
                    
         call elem_stiff_gamma_sens('tl',fkesmall,coord(:,enod(:,ie)),th,Dgp,0d0*ed(:,ie),ed_evec(:,ie),ed_evec(:,ie),0d0*esscaled,ones)   
            
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !                     h CONTRIBUTION (SHOULD FIX SUBROUTINE)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculate the different E-tilde
         ! E(gamma a,v)
         call elem_Etilde2(coord(:,enod(:,ie)),ed(:,ie),ed_evec(:,ie),E_ed_evec,crap98,gammagp(:,ie),ones) 
         ! E(gamma a,gamma v)
         call elem_Etilde1(coord(:,enod(:,ie)),ed(:,ie),ed_evec(:,ie),E_ed_evecgamma,crap98,gammagp(:,ie)) 
         ! E(0,gamma v)
         call elem_Etilde1(coord(:,enod(:,ie)),0D0*ed(:,ie),ed_evec(:,ie),E_0_evec,crap98,gammagp(:,ie))
         
         ! E(gamma v,gamma phi)
         call elem_Etilde1(coord(:,enod(:,ie)),ed_evec(:,ie),ed_evec(:,ie),E_evec_evec,crap98,gammagp(:,ie))
         ! E(gamma a, a)
         call elem_Etilde2(coord(:,enod(:,ie)),ed(:,ie),ed(:,ie),E_ed_ed,crap98,gammagp(:,ie),ones)  
         
         !!!!!!!!!!!!!!!!!!!! DIFFERENT FROM MATHIAS WORK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! E(gamma phi,gamma a)
         call elem_Etilde1(coord(:,enod(:,ie)),ed_evec(:,ie),ed(:,ie),E_evec_ed,crap98,gammagp(:,ie))
         ! E(0,gamma a)
         call elem_Etilde1(coord(:,enod(:,ie)),0d0*ed(:,ie),ed(:,ie),E_0_ed,crap98,gammagp(:,ie)) 

         ! Calculate nabla(f_AT)
         call dfat(GradF,dggamma,E_ed_evecgamma,E_ed_evecgamma)            

         ! Recalculate the tangent stiffness 
         call dneohooke('tl',Dgp,dggamma)

         !!!!!!!!!!!!!!!!!!!! DIFFERENT FROM MATHIAS WORK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Introduce quantity
         Bl = E_evec_evec-E_0_evec  
         Bl_evec_ed = E_evec_ed-E_0_ed
         
         ! Thresholding               
         do igp=1,ngp
             call gfun(penalization,g,rhogp(:,igp),pcurr,E,'E')
             call gammaprimfun(gammap,rhogp(:,igp),pcurr,je)
             ! Note that 2D0 is replaced by 1D0 since the factor 2D0 is inside sub. dfat   
             Scalar(igp) = 0D0
             Scalar(igp) = dot_product(E_ed_ed(:,igp),(matmul(Dgp(:,:,igp),Bl(:,igp))*g*gammap + GradF(:,igp)*g*gammap))
             Scalar(igp) = Scalar(igp) + 2d0*dot_product(Bl_evec_ed(:,igp),matmul(Dgp(:,:,igp),E_ed_evec(:,igp)))*g*gammap
         enddo
         
         ! Calculate as an body force
         call elem_flow_bf(fen,coord(:,enod(:,ie)),th,Scalar)
        
         fen = fen + fke + fme + fkesmall

         if (enable_filter.eq.1) then
            ffiltercurr(enod(:,ie),je) = ffiltercurr(enod(:,ie),je) + fen
         else
            ffiltercurr(ie,je) = ffiltercurr(ie,je) + sum(fen)
         endif
      enddo
   enddo
   !$OMP end parallel do
   
   return
end subroutine deigdrho


!=======================================================================

!             Compute the band-diagram in each load step

!=======================================================================
subroutine bandgap_evol(XVAL, DAT, IDAT)
   implicit none
   integer           :: iload, NEW_X, IDAT(:)
   double precision  :: DAT(:), XVAL(:)
   ! Not used but required in iter_cb
   integer           :: ALG_MODE, LS_TRIAL, ISTOP
   double precision  :: INF_PR, INF_DU, MU, DNORM, REGU_SIZE
   double precision  :: ALPHA_DU, ALPHA_PR

   do iload=1,nmax
      write(*,*) 'Compute dispersion diagram in loadstep: ', iload
      NEW_X = 1
      call eval_Newton_dg(XVAL, iload, 1d3, NEW_X, DAT, IDAT)
      call eval_EigenProb(NEW_X, DAT, IDAT)
      loadstep_bandgap = iload
      call iter_CB(ALG_MODE, IDAT(1), 0d0, INF_PR, INF_DU, MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT, DAT, ISTOP)
   enddo

   return
end subroutine bandgap_evol


!=======================================================================

!                  Solves the equilibrium iterations

!=======================================================================
subroutine eval_Newton_dg(XVAL, NMAX_loc, scale0, NEW_X, DAT, IDAT)
   implicit none
   double precision                     :: DAT(:), fe(eldof), fesmall(eldof), feb(eldof), kesmall(eldof,eldof), ke(eldof,eldof)
   double precision                     :: XVAL(:), scale0
   double precision                     :: dggamma(eldef,ngp), Dgp(elm,elm,ngp), me(eldof,eldof), scalar(ngp)
   double precision                     :: esscaled(els,ngp), esscaledsmall(els,ngp), dg_square(dofnod,dofnod)
   double precision                     :: rhogp(nbr_phase,ngp), CRAP98(eldef,ngp), Ones(ngp), Bconst(elm,ngp)
   double precision                     :: residual, g
   integer                              :: cDone, NEW_X, n, i, igp, ie, NegativeE, IDAT(:)
   integer                              :: ii, jj, ic, n_iter, NEW_EIG, nmax_loc
   integer                              :: Lnew, info, je
   double precision,allocatable         :: deltaastar(:), resstar(:), deltaa(:)
   double precision                     :: loadScale, loadIncr, elcoord(dofnod,elnod)
   double precision,allocatable         :: es(:,:,:), aold(:), edinc(:,:)
  
   allocate(es(els,ngp,nelm),stat=ierr)
   allocate(deltaastar(dof_m+nbr_free),stat=ierr)
   allocate(deltaa(ndof),stat=ierr)
   allocate(aold(ndof),stat=ierr)
   allocate(resstar(dof_m+nbr_free),stat=ierr)
   allocate(edinc(eldof,nelm),stat=ierr)

   if (NEW_X.eq.1) then
      Ones=1D0
      
      NEW_EIG = 1
      IDAT(2) = NEW_Eig
      
      ed=0d0
      edinc = 0d0
      es=0d0
      a=0d0
      if (FE_dim.eq.'3D') then
         dg=0d0
         dg(1,:,:)=1d0
         dg(5,:,:)=1d0
         dg(9,:,:)=1d0
      elseif (FE_dim.eq.'2D') then
         dg=0d0
         dg(1,:,:)=1d0
         dg(4,:,:)=1d0
      endif
      
      Deltaa=0D0
      NegativeE = 0d0
      loadScale = 0d0

      n_iter = nmax_loc
      loadIncr = 1d0/n_iter

      n = 0
      
      if (printout.eq.1) then
         write(*,*) '###################################'
         write(*,*) '       Start of NR iterations  '    
         write(*,*) '###################################' 
      endif
      
      ! Starting the for-loop (NR)
      ! Equilibrium loop starts
      Do while (n.lt.n_iter)  
         ! Variable checking which load step we're on
         n = n + 1

         aold = a

         !if ((n.eq.1).and.(idat(1).eq.1)) then
         if ((n.eq.1)) then
            loadIncr = 1d0/scale0
         elseif ((n.eq.2)) then
            loadIncr = (1d0-loadIncr)/(nmax-1d0)
         endif

         ! Update scaling of load
         loadScale = loadScale + loadIncr
 
         ! Fix loading
         if (n.eq.1) then
            !if (IDAT(1).eq.1) then
               if ((ux_control.eq.1).and.(ux_load.ne.0d0)) then
                  a = a + loadScale*ux_load*UX
               endif
               if ((uy_control.eq.1).and.(uy_load.ne.0d0)) then
                  a = a + loadScale*uy_load*UY
               endif
               if ((vx_control.eq.1).and.(vx_load.ne.0d0)) then
                  a = a + loadScale*vx_load*VX
               endif
               if ((vy_control.eq.1).and.(vy_load.ne.0d0)) then
                  a = a + loadScale*vy_load*VY
               endif
            !else
               !a = ahist(1,:)
            !endif
         else      
            a = a*(loadScale/(loadScale - loadIncr))
         endif

         call extract(ed,a,enod,dofnod)
         
         !call matwrt2f('banan.mat',a,'ainit','u')  
         !stop
         !if (n.eq.1) then
         !   call matwrt2f('banan.mat',UX,'UX','u')  
         !   call matwrt2f('banan.mat',VY,'VY','u')  
         !   stop
         !endif

         if (lagrangian_type.eq.'ul') then
            call extract(edinc,a-aold,enod,dofnod)
         endif
         
         resstar = 0d0
         resstar(bcdof)=0d0
         
         if (printout.eq.1) then
            write(*,*) 'Starting equlibrium iteration at loadstep, scale of load:', n, loadScale
         endif
       
         residual=1d0
         i=0

         do while ((i.lt.(imax)).and.(residual.gt.TOL))
            ! Some variable that checks which iteration we're on
            i=i+1
            
            K%a = 0d0
            Fvec = 0d0
            if ((n.eq.1).and.(i.eq.1)) then
               M%a = 0d0
            endif

            !$OMP PARALLEL DO DEFAULT(NONE) &
            !$OMP SHARED(nelm,coord,ed,enod,dg,es,gammagp,i,Hprojection,n,massrho,E,ones,lagrangian_type) & 
            !$OMP SHARED(K,M,dofnod,ed_rho,Dsmall,ngp,p,th,penalization,nbr_phase,els,edinc) & 
            !$OMP SHARED(XVAL,enable_filter) &
            !$OMP PRIVATE(ke,me,kesmall,Dgp,esscaled,ie,igp,rhogp,dggamma,g,scalar,je) & 
            !$OMP PRIVATE(crap98,fe,fesmall,esscaledsmall,Bconst,elcoord) &
            !$OMP REDUCTION(+:Fvec) 
            ! Calculate stiffness matrix
            do ie=1,nelm  
               ! Extract gauss point values from nodal quantities
               do je=1,nbr_phase
                  if (enable_filter.eq.1) then
                     call elem_gp(rhogp(je,:),ed_rho(:,ie,je)) 
                  else
                     rhogp(je,:) = XVAL((je-1)*nelm + ie) 
                  endif
               enddo

               ! The deformation gradient 
               call elem_dg(dg(:,:,ie),coord(:,enod(:,ie)),ed(:,ie))

               ! The gamma-scaled deformation gradient 
               call calcdggamma(dggamma,dg(:,:,ie),gammagp(:,ie))

               ! Calculate the tangent stiffness with gamma
               !call dneohooke(lagrangian_type,Dgp,dg(:,:,ie),dggamma)
               call dneohooke(lagrangian_type,Dgp,dg(:,:,ie))
             
               if (lagrangian_type.eq.'tl') then
                  ! Calculate the 2nd Piola Kirchhoff stress with gamma
                  !call neohooke('2ndPiola',es(:,:,ie),dg(:,:,ie),dggamma)
                  call neohooke('2ndPiola',es(:,:,ie),dg(:,:,ie))
               elseif (lagrangian_type.eq.'ul') then
                  !call neohooke('Cauchy',es(:,:,ie),dg(:,:,ie),dggamma)
                  call neohooke('Cauchy',es(:,:,ie),dg(:,:,ie))
               endif

               ! Now determine the material behaviour
               do igp=1,ngp
                  call gfun(penalization,g,rhogp(:,igp),p,E,'E')
                  esscaled(:,igp)= g*es(:,igp,ie)*gammagp(igp,ie)**2D0
                  Dgp(:,:,igp)   = g*Dgp(:,:,igp)*gammagp(igp,ie)**2D0
                  
                  call gfun(penalization,g,rhogp(:,igp),1d0,massrho,'rho')

                  !scalar(igp) = g*gammagp(igp,ie)**2d0
                  scalar(igp) = g
               enddo
               
               ! Calculate the large def. element stiffness matrix with gamma
               call elem_stiff_gamma(lagrangian_type,ke,coord(:,enod(:,ie)),th,Dgp,ed(:,ie),esscaled,gammagp(:,ie))

               ! Small strain contribution
               do igp=1,ngp 
                  call gfun(penalization,g,rhogp(:,igp),p,E,'E')
                  Dgp(:,:,igp)=g*(1D0-gammagp(igp,ie)**2D0)*Dsmall 
               enddo
              
               ! Calculate the small def. element stiffness matrix with gamma contribution
               call elem_stiff(lagrangian_type,kesmall,coord(:,enod(:,ie)),th,Dgp,0D0*ed(:,ie),0D0*esscaled)
               
               ! Assembling
               ! ke = gamma^2 (B'DB + G'YG), ksmall = (1-gamma^2) Bc'DBc
               !$OMP CRITICAL
               call assem(K,ke+kesmall,enod(:,ie),dofnod)
               !$OMP END CRITICAL

               ! May use initial mass matrix for both TL and UL
               if ((n.eq.1).and.(i.eq.1)) then
                  ! Mass matrix, what about gamma?
                  if (nbr_phase.gt.1) then
                     call elem_mass(me,coord(:,enod(:,ie)),th,1d0,Scalar)
                  else
                     call elem_mass(me,coord(:,enod(:,ie)),th,massrho(1),Scalar)  
                  endif

                  !$OMP CRITICAL
                  call assem(M,me,enod(:,ie),dofnod)
                  !$OMP END CRITICAL
               endif
              
               ! Thresholding
               do igp=1,ngp
                  call gfun(penalization,g,rhogp(:,igp),p,E,'E')  
                  esscaled(:,igp) = g*es(:,igp,ie)*gammagp(igp,ie)    
               enddo
                  
               ! Calculate Bconst
               call elem_Etilde1(coord(:,enod(:,ie)),0D0*ed(:,ie),ed(:,ie),Bconst,Crap98,Ones) 
                  
               ! Calculate small strain contribution
               esscaledsmall = 0d0
               do igp=1,ngp
                  call gfun(penalization,g,rhogp(:,igp),p,E,'E')
                  ! Calculate sigma (Cauchy stress)
                  if (els.eq.4) then
                     esscaledsmall([1,2,4],igp) = matmul(Dsmall,Bconst(:,igp))
                     esscaledsmall([1,2,4],igp) = g*(1D0-gammagp(igp,ie)**2D0)*esscaledsmall([1,2,4],igp)
                  else
                     esscaledsmall(:,igp) = matmul(Dsmall,Bconst(:,igp))
                     esscaledsmall(:,igp) = g*(1D0-gammagp(igp,ie)**2D0)*esscaledsmall(:,igp)
                  endif
               enddo

               ! Calculate element forces (linear + non-linear)
               call elem_force(lagrangian_type,fesmall,coord(:,enod(:,ie)),th,0D0*ed(:,ie),esscaledsmall)           
               
               call elem_force_gamma(lagrangian_type,fe,coord(:,enod(:,ie)),th,ed(:,ie),esscaled,gammagp(:,ie)) 
            
               ! Assembling
               call insert(Fvec,fe+fesmall,enod(:,ie),dofnod)
            enddo
            !$OMP END PARALLEL DO

            resstar = matmulTrans(Proj,Fvec) 
            resstar(bcdof) = 0d0
            deltaastar = -resstar

            call spamul_symprod(Kstar,Proj,K) !P'*K*P

            call solveq(Kstar,deltaastar,bcnod,bcval,dofnod,NegativeE) 
           
            if (NegativeE.gt.0) then
               write(*,*) 'Stopping due to the existence of negative eigenvalues!'
               !stop
            endif
            
            deltaa = matmul(Proj,deltaastar)
          
            a = a + Deltaa
            
            ares_hist(i,:) = a
            
            call extract(ed,a,enod,dofnod)
            call extract(edinc,deltaa,enod,dofnod)
           
            !write(*,*) 'resstar', resstar 
            
            residual=dsqrt(dot_product(resstar,resstar))                                  

            if (printout.eq.1) then
               write(*,*) 'Residual:', residual                                       
            endif
         enddo ! End of equlibrium loop
            
         if (residual.gt.tol) then
            call matwrt2f('data_error_newton.mat',ares_hist,'ares_hist','w') 
            call matwrt2f('data_error_newton.mat',ahist,'ahist','u') 
            write(*,*) 'Error in eval_Newton_dg: Newton iterations not converging.'
            stop
         endif

         ahist(n,:) = a

         if (printout.eq.1) then
            write(*,*) 'Horizontal macroscopic stress', dot_product(UX,Fvec), dot_product(UY,Fvec)
            write(*,*) 'Vertical macroscopic stress', dot_product(VX,Fvec), dot_product(VY,Fvec)
         endif

      enddo ! End of load steps 
  
      if (isnan(residual)) then
          stop
      endif

      if (printout.eq.1) then
         write(*,*) 'End of Newton-loop'
      endif
   endif ! End if(NEW_X.eq.1)

   return
end subroutine eval_Newton_dg






!=======================================================================

!         Solves the eigenvalue problem at the last load step

!=======================================================================
subroutine eval_EigenProb(NEW_X, DAT, IDAT)
   implicit none
   double precision                     :: DAT(:), kvec(dofnod), f_time, s_time, max_L, min_U
   integer                              :: NEW_X, NEW_EIG, IDAT(:),  i, info, j 
   integer                              :: restart, indx, ierr
   double precision,allocatable         :: EVECS(:,:), EVALS(:), EVECSOLD(:,:)
   character(1)                         :: which
   type(sparse)                         :: Kstar_Re_Im, Mstar_Re_Im

   allocate(EVECS(dof_m*2,dof_m*2),stat=ierr)
   allocate(EVALS(dof_m*2),stat=ierr)
   allocate(EVECSOLD(dof_m*2,dof_m*2),stat=ierr)
    
   NEW_EIG = IDAT(2)
         
   ! Search for smallest eigenvalues
   which = 'S'
 
   if (NEW_EIG.eq.1) then
      
      if (printout.eq.1) then
         write(*,*) 'Eigenvalue sensitivity analysis; set up eigenvalue problem'
      endif

      ! Copy over values
      K_Re_Im%a(1:size(K%a)) = K%a/scaleK
      K_Re_Im%a(size(K%a)+1:size(K_Re_Im%a)) = K%a/scaleK

      M_Re_Im%a(1:size(M%a)) = M%a
      M_Re_Im%a(size(M%a)+1:size(M_Re_Im%a)) = M%a
      
      ! Reset
      evecs_tot = 0d0
      evals_tot = 0d0
      evals_L = 0d0
      evals_U = 0d0
      evals_U_inv = 0d0

      !!$ s_time = omp_get_wtime()

      do i=1,nbr_k
         if (printout.eq.1) then
            write(*,*) '##############################################################################'
            write(*,*) 'Computing eigenvalues for k-vector nbr ', i, ' out of ', nbr_k
            write(*,*) '##############################################################################'
         endif
         kvec = kvec_mat(:,i)
         
         ! Update boundary condition matrices 
         call periodic_bc_bloch(kvec, 1)
         
         call spamul_symprod(Kstar_Re_Im,ProjBloch,K_Re_Im) !P'*K*P
         call spamul_symprod(Mstar_Re_Im,ProjBloch,M_Re_Im) !P'*K*P

         if (i.gt.1) then
            EVECSOLD = EVECS
         endif
         
         restart = 1

         !if (i.eq.1) then
         !   Emax = 1d3
         !   L0 = 100
         !elseif (i.eq.2) then
         !   Emax = 1d10
         !   L0 = 2000
         !endif

         if (L0.gt.dof_m*2) then
            L0 = dof_m*2
         endif

         do while (restart.eq.1)
            info = 99
            call eigenvalue(Kstar_Re_Im,Mstar_Re_Im,dof_m*2,L0,L,EVALS,EVECS,EVECSOLD,Emin,Emax,i,info,printout)
            
            if (printout.eq.1) then
               write(*,*) 'Number of eigenvalues found: ', L
            endif
            !write(*,*) 'eigenvalues',evals(1:L)
            
            if (info.eq.0) then
               if (L.lt.nbr_eigMax) then
                  write(*,*) '##############################################################################'
                  write(*,*) 'eval_EigProb: Too few unique eigenvalues were found; increase search interval.'
                  restart = 1
                  Emax = 1.2d0*Emax 
                  L0 = 1.5*L0
                  if (L0.gt.dof_m*2) then
                     L0 = dof_m*2   
                  endif
                  write(*,*) 'New Emax: ', Emax
                  write(*,*) 'New L0: ', L0
                  write(*,*) '##############################################################################'
               else
                  restart = 0
                  !if (i.ne.1) then
                     L0 = L*1.5d0
                  !endif
               endif
            elseif (info.eq.1) then
               write(*,*) '##############################################################################'
               write(*,*) 'eval_EigProb: No eigenvalues were found; increase search interval.'
               restart = 1
               Emax = 1.5d0*Emax
               L0 = 1.5*L0
               if (L0.gt.dof_m*2) then
                  L0 = dof_m*2   
               endif
               write(*,*) 'New Emax: ', Emax
               write(*,*) 'New L0: ', L0
               write(*,*) '##############################################################################'
            elseif (info.eq.3) then
               write(*,*) '##############################################################################'
               write(*,*) 'eval_EigProb: Size of subspace L0 too small; decrease search interval.'
               restart = 1
               Emax = 0.7*Emax
               L0 = 0.8*L0
               write(*,*) 'New Emax: ', Emax
               write(*,*) 'New L0: ', L0
               write(*,*) '##############################################################################'
            elseif (info.le.-100) then
               write(*,*) 'eval_EigProb: Error in FEAST! Stopping.'
               stop
            endif
         enddo

         do j=1,nbr_eigMax
            EVALS_TOT(j,i) = dabs(EVALS(j))!*scaleK
            EVECS_TOT(j,(ndof*2)*(i-1)+1:(ndof*2)*i) = matmul(ProjBloch,EVECS(:,j))
            ! Obtain sets above and under bandgap
            if (j.le.band) then
               evals_L(j,i) = evals(j)
            elseif (j.gt.band) then
               evals_U(j-band,i) = evals(j)
               !evals_U_inv(j-band,i) = bound_U/evals(j)
               evals_U_inv(j-band,i) = 1d0/evals(j)
            endif
         enddo
 
         if (printout.eq.1) then
            write(*,*) 'Eigenvalues (unique) above bandgap: ', EVALS_U(1:nbr_eig_U:2,i)
            write(*,*) 'Eigenvalues (unique) below bandgap: ', EVALS_L(1:nbr_eig_L:2,i)
         endif
      enddo
 
      ! If on first design iteration determine bound variables
      !if (IDAT(1).eq.1) then
      !   max_L = 0d0
      !   min_U = 1d20
      !   do i=1,nbr_k
      !      if (evals_L(nbr_eig_l,i).gt.max_L) then
      !         max_L = evals_L(nbr_eig_l,i)
      !      endif
      !      if (evals_U(1,i).lt.min_U) then
      !         min_U = evals_U(1,i)
      !      endif
      !   enddo
      !   ! Set the variables to the middle distance
      !   bound_L = (min_U+max_L)/2d0
      !   bound_U = bound_L
      !   write(*,*) 'Bound variables:'
      !   write(*,*) 'Bound_L', bound_l
      !   write(*,*) 'Bound_U', bound_U
      !endif
      
      !!$ f_time = omp_get_wtime()
      !write(*,*) 'Time taken: ', (f_time-s_time), 'seconds'
      
      NEW_EIG = 0
      IDAT(2) = NEW_EIG 
   endif

   return
end subroutine eval_EigenProb


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

   if (size(dAggcurr,1).ne.size(dargdrho,2)) then
      write(*,*) 'Size incoherence in call to: Aggprim'
      stop
   endif

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
!                   Only controlled by void/solid phase

!=======================================================================
subroutine gammafun(val, rhogpCurr, pcurr)
   implicit none
   double precision         :: val2, rhogpCurr(:), val, pcurr

   val = 0d0
   val2 = 0d0

   if (Hprojection.eq.1) then
      call Hfun(val2,rhogpCurr(nbr_phase))
   else
      val2 = rhogpCurr(nbr_phase)
   endif

   if (val2.lt.0d0) then
      val2 = 0d0
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
!                   Only controled by void/solid phase

!=======================================================================
subroutine gammaprimfun(val,rhogpCurr,pcurr,sens_indx)
   implicit none
   double precision                 :: val2, val3, rhogpCurr(:), val, pcurr
   integer                          :: sens_indx

   if (sens_indx.eq.nbr_phase) then
      val = 0d0
      val2 = 0d0
      val3 = 0d0

      if (Hprojection.eq.1) then
         call Hfun(val2,rhogpCurr(sens_indx))
      else
         val2 = rhogpCurr(sens_indx)
      endif
      if (val2.lt.0d0) then
         val2 = 0d0
      endif

      val=0D0
      if (Gscaling.eq.1) then
         val = (0.1D1-dtanh(B1*(val2**pcurr-rho0))**2D0) *  & 
                          B1*pcurr*val2**(pcurr-1d0)/(dtanh(B1*rho0) + dtanh(B1 * (0.1D1 - rho0)))

         if (Hprojection.eq.1) then
            call Hprimfun(val3, rhogpCurr(sens_indx))
            val = val*val3
         endif
      else
         val=0D0
      endif
   else
      val = 0d0
   endif
   
   return
end subroutine gammaprimfun


!=======================================================================

!                 Choose penalization function

!=======================================================================
subroutine gfun(pentype,val,rhogpcurr,pcurr,Mater,MaterType)
   implicit none
   double precision                 :: val, rhogpCurr(:), pCurr, Mater(:)
   character(len=*)                 :: pentype, MaterType
   
   if (size(rhogpCurr,1).eq.1) then
      if (pentype.eq.'SIMP') then
         call simpfun(val,rhogpcurr(1),pcurr,MaterType)
      elseif (pentype.eq.'RAMP') then
         call rampfun(val,rhogpcurr(1),pcurr,MaterType)
      else
        stop 'This type of aggregation function has not been implemented!'
      endif
   else
      if (pentype.eq.'SIMP') then
         call simpfun_multi(val,rhogpcurr,pcurr,Mater,MaterType)
      else
        stop 'This type of aggregation function has not been implemented!'
      endif
   endif

   return
end subroutine gfun



!=======================================================================

!                    Calculates the Xi-function
!                            (SIMP)

!=======================================================================
subroutine simpfun(val,rhogpCurr,pcurr,MaterType)
   implicit none
   double precision                 ::  rhogpCurr, val2, val, pcurr, delta0
   character(len=*)                 ::  MaterType
   
   val2 = 0d0
   val = 0d0

   if (MaterType.eq.'E') then
      delta0 = delta0_E
   elseif (MaterType.eq.'rho') then
      delta0 = delta0_rho
   else
      Write(*,*) 'Material type is not implemented.'
      stop
   endif

   if (Hprojection.eq.1) then
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
!                         (SIMP-multiphase)

!=======================================================================
subroutine simpfun_multi(val,rhogpCurr,pcurr,Mater,MaterType)
   implicit none
   double precision           :: rhogpCurr(:), val2(nbr_phase), val, pcurr, delta0
   double precision           :: Mater(:)
   integer                    :: ie      
   character(len=*)           :: MaterType
   
   if (MaterType.eq.'E') then
      delta0 = delta0_E
   elseif (MaterType.eq.'rho') then
      delta0 = delta0_rho
   else
      Write(*,*) 'Material type is not implemented.'
      stop
   endif

   val2 = 0d0
   val = 0d0

   do ie = 1,nbr_phase
      if (Hprojection.eq.1) then
         call Hfun(val2(ie),rhogpCurr(ie))
      else
         val2(ie) = rhogpCurr(ie)
      endif
      if (val2(ie).lt.0d0) then
         val2(ie) = delta0
      endif
   enddo

   val = (val2(nbr_phase)**pcurr)*Mater_out(val2,pcurr,Mater,nbr_phase-1)

   return
end subroutine simpfun_multi



!=======================================================================

!      Obtain interpolation of material parameter using recursion
!                           (SIMP)

!=======================================================================
recursive function Mater_out(rhogpcurr,pcurr,Mater,i)
   implicit none
   double precision        :: Mater_out, rhogpcurr(:), Mater(:),pcurr
   integer                 :: i

   if (i.eq.1) then
      Mater_out = rhogpcurr(1)**pcurr*Mater(nbr_phase) + (1d0-rhogpcurr(1)**pcurr)*Mater(nbr_phase-1)
   else
      Mater_out = rhogpcurr(i)**pcurr*Mater_out(rhogpcurr,pcurr,Mater,i-1) + (1d0-rhogpcurr(i)**pcurr)*Mater(nbr_phase-i)
   endif

   return
end function Mater_out



!=======================================================================

!                    Calculates the Xi-function
!                            (RAMP)

!=======================================================================
subroutine rampfun(val,rhogpCurr,pcurr,MaterType)
   implicit none
   double precision                 ::  rhogpCurr, val2, val, pcurr, delta0
   character(len=*)                 ::  MaterType
   
   val2 = 0d0
   val = 0d0

   if (MaterType.eq.'E') then
      delta0 = delta0_E
   elseif (MaterType.eq.'rho') then
      delta0 = delta0_rho
   else
      Write(*,*) 'Material type is not implemented.'
      stop
   endif

   if (Hprojection.eq.1) then
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
subroutine gfunprim(pentype,val,rhogpCurr,pcurr,Mater,MaterType,sens_indx)
   implicit none
   double precision                 :: val, rhogpCurr(:), pcurr, Mater(:)
   integer                          :: sens_indx
   character(len=*)                 :: pentype, MaterType
   
   if (size(rhogpCurr,1).eq.1) then
      if (pentype.eq.'SIMP') then
        call simpfunprim(val,rhogpcurr(1),pcurr,MaterType)
      elseif (pentype.eq.'RAMP') then
        call rampfunprim(val,rhogpcurr(1),pcurr,MaterType)
      else
        stop 'This type of aggregation function has not been implemented!'
      endif
   else
      if (pentype.eq.'SIMP') then
         call simpfunprim_multi(val,rhogpcurr,pcurr,Mater,MaterType,sens_indx)
      else
        stop 'This type of aggregation function has not been implemented!'
      endif
   endif

   return
end subroutine gfunprim




!=======================================================================

!           Calculates the derivative of Xi-function
!                            (SIMP)

!=======================================================================
subroutine simpfunprim(val,rhogpCurr,pcurr,MaterType)
   implicit none
   double precision          :: rhogpCurr, val2, val3, val, pcurr, delta0
   character(len=*)          ::  MaterType

   val = 0d0
   val2 = 0d0
   val3 = 0d0
   
   if (MaterType.eq.'E') then
      delta0 = delta0_E
   elseif (MaterType.eq.'rho') then
      delta0 = delta0_rho
   else
      Write(*,*) 'Material type is not implemented.'
      stop
   endif

   if (Hprojection.eq.1) then
      call Hfun(val2, rhogpCurr)
      call Hprimfun(val3, rhogpCurr)
      val=pcurr*(1D0-delta0)*val2**(pcurr-1D0)*val3
   else
      val=pcurr*(1D0-delta0)*rhogpCurr**(pcurr-1D0)
   endif

   if (rhogpCurr.lt.0D0)  val=0D0

   return
end subroutine simpfunprim


!=======================================================================

!           Calculates the derivative of Xi-function
!                        (SIMP-multiphase)

!=======================================================================
subroutine simpfunprim_multi(val,rhogpCurr,pcurr,Mater,MaterType,sens_indx)
   implicit none
   double precision          :: rhogpCurr(:), val2(nbr_phase), val3(nbr_phase), val, pcurr
   double precision          :: Mater(:)
   integer                   :: ie, sens_indx
   character(len=*)          :: MaterType

   val = 0d0
   val2 = 0d0
   val3 = 0d0
   
   do ie = 1,nbr_phase
      if (Hprojection.eq.1) then
         call Hfun(val2(ie),rhogpCurr(ie))
      else
         val2(ie) = rhogpCurr(ie)
      endif
      if (val2(ie).lt.0d0) then
         val2(ie) = 0d0
      endif
   enddo

   val3 = 0d0
   if (Hprojection.eq.1) then
      call Hprimfun(val3(sens_indx), rhogpCurr(sens_indx))
   else
      val3(sens_indx) = 1d0
   endif

   if (sens_indx.ne.nbr_phase) then
      val = (val2(nbr_phase)**pcurr)*Mater_out_prim(val2,val3,pcurr,Mater,nbr_phase-1,sens_indx)
   else
      val = (pcurr*val2(nbr_phase)**(pcurr-1d0))*Mater_out(val2,pcurr,Mater,nbr_phase-1)*val3(nbr_phase)
   endif
   
   !if (rhogpCurr(sens_indx).lt.0D0) val=0D0


   return
end subroutine simpfunprim_multi



!=======================================================================

!       Obtain derivative of material parameter using recursion
!                           (SIMP)

!=======================================================================
recursive function Mater_out_prim(rhogpcurr,rhogpcurr_prim,pcurr,Mater,i,sens_i)
   implicit none
   double precision        :: Mater_out_prim, rhogpcurr(:), rhogpcurr_prim(:), Mater(:), pcurr
   integer                 :: i, sens_i

   if (i.eq.1) then
      if (i.ne.sens_i) then
         Mater_out_prim = rhogpcurr(1)**pcurr*Mater(nbr_phase) &
            + (1d0-rhogpcurr(1)**pcurr)*Mater(nbr_phase-1)
      else
         Mater_out_prim = pcurr*rhogpcurr(1)**(pcurr-1d0)*Mater(nbr_phase)*rhogpcurr_prim(1) & 
            + (-pcurr*rhogpcurr(1)**(pcurr-1d0))*Mater(nbr_phase-1)*rhogpcurr_prim(1)
      endif
   else
      if (i.ne.sens_i) then
         Mater_out_prim = rhogpcurr(i)**pcurr*Mater_out_prim(rhogpcurr,rhogpcurr_prim,pcurr,Mater,i-1,sens_i) &
            + (1d0-rhogpcurr(i)**pcurr)*Mater(nbr_phase-i)
      else
         Mater_out_prim = pcurr*rhogpcurr(i)**(pcurr-1d0)*Mater_out_prim(rhogpcurr,rhogpcurr_prim,pcurr,Mater,i-1,sens_i)*rhogpcurr_prim(i) & 
            + (-pcurr*rhogpcurr(i)**(pcurr-1d0))*Mater(nbr_phase-i)*rhogpcurr_prim(i)
      endif
   endif

end function Mater_out_prim



!=======================================================================

!           Calculates the derivative of Xi-function
!                            (RAMP)

!=======================================================================
subroutine rampfunprim(val,rhogpCurr,pcurr,MaterType)
   implicit none
   double precision          :: rhogpCurr, val2, val3, val, pcurr, delta0
   character(len=*)          :: MaterType

   val = 0d0
   val2 = 0d0
   val3 = 0d0
   
   if (MaterType.eq.'E') then
      delta0 = delta0_E
   elseif (MaterType.eq.'rho') then
      delta0 = delta0_rho
   else
      Write(*,*) 'Material type is not implemented.'
      stop
   endif

   if (Hprojection.eq.1) then
      call Hfun(val2, rhogpCurr)
      call Hprimfun(val3, rhogpCurr)
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
subroutine Hfun(val, rhogpCurr)
   implicit none
   double precision                 ::  rhogpCurr, val
      
   val = 0d0  
   val = (dtanh(B*eta)+dtanh(B*(rhogpCurr-eta)))/(dtanh(B*eta)+dtanh(B*(0.1D1-eta)))  
   
   return
end subroutine Hfun


    
!=======================================================================

!        Calculates the derivative of the Heaviside function

!=======================================================================   
subroutine Hprimfun(val, rhogpCurr)
   implicit none
   double precision                 ::  rhogpCurr, val
    
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




end module opt_module









