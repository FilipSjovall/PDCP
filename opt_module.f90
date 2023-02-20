!#####################################################

! The module contains subroutines needed to run either
! MMA or IPOPT

!#####################################################
! Created: 2019-12-13
! Anna Dalklint
! Lund University
!#####################################################
module opt_module

use abaqus_util
use matlab_util
use sparse_util
use mater_hyperel
use elem_large_cont_3d
use elem_flow_3d

   ! Global variables
   implicit none
	integer*8                           :: matptr
	integer, allocatable           		:: enod(:,:) , bcnod(:,:), flnod(:,:)
	double precision, allocatable			:: coord(:,:), bcval(:), flval(:)
	character(len=40)                   :: filename
	integer										:: numThread, ierr, restart
	integer										:: nbc, dofnod, nnod, ndof, nelm, max_iter, C
	integer										:: ngp, NcMax, iwait, nmax
  
	double precision                    :: Emin, Emax, pagg, pAggFinal, stepLength, scaleCrit
	integer									   :: L, linB, Ncrit

	double precision							:: start, finish
	character(80)								:: aggregation, optimizer
	
	double precision, allocatable       :: ed(:,:), es(:,:,:), esx_nod(:), esy_nod(:), esxy_nod(:), gamma_nod(:)
	double precision, allocatable       :: res(:), designRes(:), a(:), anom(:), acrit(:), ahist(:,:)
	double precision, allocatable       :: Deltaa(:), deltaaold(:), DeltaaSens(:)
	double precision, allocatable       :: dg(:,:,:), edvarphi(:,:), varphi(:)

	double precision                    :: DeltaLambda, DeltaLambdaSens, DeltaLambdaOld, cPhi
	double precision, allocatable       :: aold(:), dag(:), dap(:), dapSens(:), dapOld(:) 
	integer, allocatable                :: fldof(:), bcdof(:)
	double precision                    :: Dsmall(6,6), mp(2)
	double precision                    :: unitVec(9)
	double precision                    :: designResidual
	double precision, allocatable       :: h(:)
	double precision                    :: E, v, p, omega, pold
	double precision                    :: Vtot, Volfrac, Ve
	double precision, allocatable       :: tmp2(:), f0valhist(:), eigenValhist(:,:), resvalhist(:), f1valhist(:), f2valhist(:), history(:,:)
	double precision, allocatable       :: Pattern(:), Fvec(:), dSfac(:,:)
	double precision                    :: pert
	integer										:: ipert
	double precision                    :: TOL, designTOL, LambdaNom, LambdaCrit, LambdaCritVal

	double precision                    :: Scaleg0, SafetyMarginal, LambdaMax
	double precision, allocatable       :: phi(:,:), phiC(:), tmpPhi(:), sfactor(:), sfactorInv(:), Xold(:,:), X(:,:)
	double precision, allocatable       :: edphi(:,:), edmu(:,:)

	type(sparse)                        :: K, Ktmp, Ko, Kg
   double precision, allocatable       :: Kofull(:,:), Kgfull(:,:), Kfilterfull(:,:) 

	double precision              		:: Z, geps, GMOVE, Gvol, Gstab
	double precision, allocatable 		:: rho(:), xold1(:), xold2(:)
	double precision, allocatable			:: dGvol(:), dGstab(:)
	double precision, allocatable 		:: df(:), df0(:), dAgg(:)

  
	! IPOPT
	double precision						   :: CTOL
	integer, allocatable						:: ACON1(:), AVAR1(:)
	double precision, allocatable			:: GLOW(:), GUPP(:)
	  
	! PDE Filter 
	double precision, allocatable       :: Tmatrix(:,:),rhotilde(:),ffilter(:),ffilterTMP(:)
	double precision, allocatable       :: edrho(:,:), gammagp(:,:)
	type(sparse)                        :: Kfilter, KfilterTrans
	double precision                    :: radius, r2

	! Different on/off parameters  
	integer									   :: Hprojection, diffCheck, new_p, update_done
	integer										:: diffCheck_iteration, irestart, restart_iteration
    
	! Heaviside & gamma
	double precision                    :: B1, rho0, delta0 
	double precision                    :: B, eta, p1, Binit, Bfinal  

	! used when saving matlab file
	integer         						   :: indx1
	character(80)   						   :: s_plot
	character(120)  						   :: xstring

   ! For creating symmetric K
	type(sparse)                        :: KSYM
	integer, allocatable                :: cellSym(:,:)
	double precision, allocatable       :: valSym(:)   



contains


!=======================================================================

!                 Allocate global variables

!=======================================================================
subroutine opt_Init()
	implicit none
	double precision 		:: Dgp(6,6,8), tmpreal
	integer					:: ie, je

	open(unit=5, file='data.out', status='unknown')
	call matinit(matptr)
	
	call abareadinp(filename,coord,enod,bcnod,bcval,flnod,flval) 

	nbc = size(bcval)
	dofnod = 3
	nnod = maxval(enod)
	ndof = nnod*dofnod
	nelm = size(enod,2)
  
	write(*,*) 'Number of nods                :', nnod
	write(*,*) 'Number of dofs                :', ndof
	write(*,*) 'Number of elements            :', nelm
	write(*,*) 'Number of boundary conditions :', nbc
	write(*,*) 'Number of concentrated loads  :', size(flval)
	

	! FE-quantities
	allocate(Kofull(ndof,ndof), stat=ierr)  	
	allocate(Kgfull(ndof,ndof), stat=ierr)  	
	allocate(Kfilterfull(nnod,nnod), stat=ierr)  
	allocate(esx_nod(nnod), stat=ierr)
	allocate(esy_nod(nnod), stat=ierr)
	allocate(esxy_nod(nnod), stat=ierr)
	allocate(gamma_nod(nnod), stat=ierr)	
	allocate(ed(24,nelm), stat=ierr)
	allocate(es(6,ngp,nelm), stat=ierr)
	allocate(dg(9,8,nelm), stat=ierr)
	allocate(a(ndof), stat=ierr)
	allocate(ahist(1000,ndof), stat=ierr)
	allocate(anom(ndof), stat=ierr) 	
	allocate(aold(ndof), stat=ierr) 	
	allocate(acrit(ndof), stat=ierr)
	allocate(res(ndof), stat=ierr)
	allocate(designRes(ndof), stat=ierr)   
   allocate(h(ndof), stat=ierr)
	allocate(tmpPhi(ndof), stat=ierr)
	allocate(phiC(ndof), stat=ierr)		
	allocate(dag(ndof), stat=ierr)
	allocate(dap(ndof), stat=ierr)   
	allocate(Pattern(ndof), stat=ierr) 
	allocate(varphi(ndof), stat=ierr)    
	allocate(edvarphi(24,nelm), stat=ierr)
	allocate(Deltaaold(ndof), stat=ierr)
	allocate(Deltaa(ndof), stat=ierr)	
	allocate(sfactor(ndof), stat=ierr)
	allocate(sfactorInv(NcMax), stat=ierr)
	allocate(edphi(24,nelm), stat=ierr)
	allocate(edmu(24,nelm), stat=ierr)
	allocate(fldof(size(flval)), stat=ierr)
	allocate(bcdof(size(bcval)), stat=ierr)
	allocate(Fvec(ndof), stat=ierr)  
	allocate(eigenValhist(NcMax,MAX_ITER), stat=ierr)
	allocate(history(2,ndof), stat=ierr)
	allocate(dAgg(nelm), stat=ierr)
	allocate(phi(ndof,Ncmax), stat=ierr)
	
	allocate(dapSens(ndof), stat=ierr)
	allocate(dapOld(ndof), stat=ierr)
	allocate(DeltaaSens(ndof), stat=ierr)
		
		
	! Optimization
	allocate(dSfac(NcMax,nelm), stat=ierr)
	allocate(tmp2(4), stat=ierr)
	allocate(f0valhist(MAX_ITER), stat=ierr)
	allocate(f1valhist(MAX_ITER), stat=ierr)
	allocate(f2valhist(MAX_ITER), stat=ierr)
	allocate(resvalhist(MAX_ITER), stat=ierr)
	allocate(glow(C),stat=ierr)
	allocate(gupp(C),stat=ierr)

!	allocate(xval(nelm), stat=ierr)
	allocate(xold1(nelm), stat=ierr)
	allocate(xold2(nelm), stat=ierr)
	allocate(rho(nelm), stat=ierr)
	allocate(df0(nelm), stat=ierr)
	allocate(df(nelm*C), stat=ierr)	
	allocate(dGvol(nelm), stat=ierr)	
	allocate(dGstab(nelm), stat=ierr)	
		
	allocate(gammagp(8,nelm), stat=ierr)
	allocate(rhotilde(nnod), stat=ierr)
	allocate(ffilter(nnod), stat=ierr)
	allocate(ffilterTMP(nnod), stat=ierr)
	allocate(edrho(8,nelm), stat=ierr)
	allocate(Tmatrix(8,nelm), stat=ierr) 

  
   UnitVec   =0D0
	UnitVec(1)=1d0
	UnitVec(5)=1d0
	UnitVec(9)=1d0 
	
	ahist = 0d0
	eigenValhist =0d0
	f0valhist = 0d0
	f1valhist = 0d0
	f2valhist = 0d0
	resvalhist = 0d0
!	history = 0d0
	
	! Find dofs for boundary and loads
	Pattern=0D0
   call finddof(fldof,flnod,dofnod)
   call finddof(bcdof,bcnod,dofnod)
   bcval = 0d0
	Pattern(fldof)=flval
	
	ed=0d0
	es=0d0
	dg=0d0
	dg(1,:,:)=1d0
	dg(5,:,:)=1d0
	dg(9,:,:)=1d0
	a=0d0
	res=0d0
	
	! Sparse structures
	call spatopdef(K,enod,dofnod)
   	call spatopdef(Ko,enod,dofnod)
   	call spatopdef(Kg,enod,dofnod)
  	call spatopdef(Kfilter,enod)  
	
	p1 = p
	pold = p
   
   ! Material parameters
   mp(1)		= E/3d0/(1d0-2d0*v)  ! K
	mp(2)		= E/2d0/(1d0+v)      ! G
	! Obtains K and G
   call neohooke_init(mp)
  	
  	! Material tangent stiffness for pure displacement
   call dneohooke('tl',Dgp,dg(:,:,1))  
  
   ! Dsmall is the constant linear stiffness
   Dsmall = Dgp(:,:,1) 
   
   ! Calculate element volumes
	Vtot=0D0
	do ie=1,nelm
   	call elmvol(Ve,coord(:,enod(:,ie)))  
   	Vtot = Vtot + Ve
	enddo
	Vtot=Volfrac*Vtot 
	!write(*,*) 'Vtot, Volfrac', Vtot, Volfrac

   ! We want the upper triangle of K as sparse pattern to obtain more from paradiso
	call upper_triangular_matrix()

	! Initiate PDE-filter matrices
	call PDE_filter_init()
	
	new_p = 0
	
	call matwrt2f('data_common.mat',enod,'enod','w')  
   call matwrt2f('data_common.mat',coord,'coord','u')  

	return
end subroutine opt_Init


!=======================================================================

!           Create upper triangular sparse pattern for K
!	    Used in Pardiso to obtain more information

!=======================================================================
subroutine upper_triangular_matrix()
	implicit none
	integer				:: j, i, ii

	write(*,*)
   write(*,*) 'Start to create symmetric sparse pattern'
   ! Assign K = 1d0
   K%a=1d0
   j=0 

   do i=1,ndof
   	do ii=i,ndof
      	if ((spagetval(K,i,ii)).gt.0D0) then
      		j=j+1    
      	endif              
      enddo
   enddo     
   
   allocate(cellSym(j,2), stat=ierr)
   allocate(valSym(j), stat=ierr)  
   j=0

   do i=1,ndof
   	do ii=i,ndof
      	if ((spagetval(K,i,ii)).gt.0D0) then
         	j=j+1
         	cellSym(j,:)=(/i,ii/)    
         endif
      enddo
   enddo   

   call spacelldef(KSYM,cellSym,ndof,ndof)  

   K%a=0d0     

   write(*,*) 'Symmetric sparse pattern created'
	write(*,*)

	return
end subroutine upper_triangular_matrix


!=======================================================================

!                    Initiate PDE-filter matrices

!=======================================================================
subroutine PDE_filter_init()
	implicit none
	integer				:: ie
	double precision  :: ken(8,8), Men(8,8), fen(8)
	
	! Calculate T, M_filter and K_filter
	! M_filter and K_filter are added together for simplification 	 
  	Kfilter%a=0d0
  	
   do ie=1,nelm
      call fl3d8_e(ken,coord(:,enod(:,ie)),r2)  
      call fl3d8_m(Men,coord(:,enod(:,ie)),1d0)
      ken = ken + Men
      call fl3d8_bf(fen,coord(:,enod(:,ie)),1d0)
      
      call assem(Kfilter,ken,enod(:,ie),1) 
      Tmatrix(:,ie) = fen
   enddo  
	
	return
end subroutine PDE_filter_init
 

!=======================================================================

!              Update filtered design variables

!=======================================================================
subroutine update_PDE_filter(XVAL, NEW_X)
	implicit none
	double precision	:: XVAL(NELM)
	integer				:: ie, NegativeE, NEW_X

	if (NEW_X.eq.1) then
		ffilter  = 0D0
		
		rho = XVAL 
		
		! Denotes T*XVAL as ffilter
		do ie=1,nelm
		   ffilter(enod(:,ie)) = ffilter(enod(:,ie)) + Tmatrix(:,ie)*XVAL(ie) 
		enddo
		
		! Solve to obtain rho_tilde
		call solveq(Kfilter,ffilter,NegativeE) 

		rhotilde=ffilter 
		
		! Obtain element rho_tilde (edrho)
		call extract(edrho,ffilter,enod,1)
	endif
	
	return
end subroutine update_PDE_filter



!=======================================================================

!                     Choose the optimizer

!=======================================================================

subroutine opt(choice)
	implicit none
	character(80) :: choice
	
	if (choice.eq.'MMA') then
		call MMA()
	elseif (choice.eq.'IPOPT') then
		!call IPOPT()
		stop 'This optimizer has not been implemented!'
	else
		stop 'This optimizer has not been implemented!'
	endif


end subroutine opt


!=======================================================================

!                             MMA

!=======================================================================

subroutine MMA()
	implicit none
	double precision			:: Z, geps, GMOVE, tmpreal
	double precision			:: QMMA(NELM*C), P0MMA(NELM), Q0(NELM), UU(C), PMMA(NELM*C)
	integer						:: IYFREE(C)
	double precision			:: GRADF(C), DSRCH(C), HESSF(C*(C+1)/2)
	double precision			:: XVAL(NELM), XLOW(NELM), XUP(NELM), XMIN(NELM), XMAX(NELM)
	double precision			:: ALFA(NELM), BETA(NELM)
	double precision			:: AMMA(C), BMMA(C), CMMA(C), YMMA(C), ULAM(C)
	double precision 			:: FVAL(C), F0VAL
	double precision 			:: DAT(1), rhomin, rhomax
	integer 			  			:: IDAT(1), ie, je, ki, i, IMMA, NEW_X, res_iter
	double precision			:: DFMMA(NELM*C), FMAX(C)
	! Not used but required in iter_cb
	integer 						:: ALG_MODE, LS_TRIAL, ISTOP, NZ, TASK, ACON(nelm*C), AVAR(nelm*C)
	double precision			:: OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE, AJAC(nelm*C)
	double precision			:: ALPHA_DU, ALPHA_PR

	   
	GEPS =1d-9 ! MMA Tolerance parameter for the constraints.
   GMOVE = 1d0
   NZ = nelm*C
   designTol = 1d-3
	IMMA = 0

	! -----------------------------------------------------------------------
	! 								Derivative test
	! -----------------------------------------------------------------------	   
   if (diffCheck.eq.1) then
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
   	  	rho(ipert) = rho(ipert)+pert
	  		
	  		XVAL = rho
	  		
	  		call matGetVal(matptr,'tmp2',tmp2)
		  	p=tmp2(1)
		  	B=tmp2(2)
		  	p1 = p
		  	
		  	update_done = 0
		  	
     	else   		
   	  	write(*,*) 'Performing a pertubation at optimization-step,', IMMA  

   	  	write(*,*) 'Size of perturbation: ', pert
   	  	write(*,*) 'Design variable perturbated: ', ipert
     		! Initial guess
			do ie=1,nelm 
				!CALL RANDOM_NUMBER(tmpreal)
				!XVAL(ie)=(Volfrac-0.05D0)+0.1d0*tmpreal
				XVAL(ie) = Volfrac  
			enddo
			
			rho = XVAL
			rho(ipert) = rho(ipert)+pert
			XVAL = rho
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
		XVAL = rho
	  		
		call matGetVal(matptr,'tmp2',tmp2)
		p=tmp2(1)
		B=tmp2(2)
		p1 = p
   	IMMA = res_iter
   endif
   
   ! Regular start
   if ((irestart.eq.0).and.(diffCheck.eq.0)) then
		! Initial guess
		do ie=1,nelm 
			!CALL RANDOM_NUMBER(tmpreal)
			!XVAL(ie)=(Volfrac-0.05D0)+0.1d0*tmpreal
			XVAL(ie) = Volfrac   
		enddo
	endif
	
	IYFREE = 1

	rho = XVAL
	XOLD1 = XVAL
	XOLD2 = XVAL
	
	! Min/max design variables
	rhomin = 0.00000D0
	rhomax = 1D0
	
	XMIN = rhomin  
	XMAX = rhomax 

	AMMA = 0d0
	CMMA = 1D3
	
	FMAX = GUPP
	
	FVAL = 0d0
	F0VAL = 0d0
	
	TASK = 1

	designResidual = 1
	! If NEW_X=1 we recently updated the design variables, otherwise NEW_X=0
	NEW_X = 1

	! Run optimization loop
	do while ((designResidual.gt.designTOL).and.(IMMA.lt.max_iter))
		
		IMMA = IMMA + 1
		! Iteration counter
		IDAT(1) = IMMA	
			
		! Obtain value of objective
		call eval_F(NELM, XVAL, NEW_X, F0VAL, IDAT, DAT, IERR)
		if (IERR.ne.0) write(*,*) 'Error in eval_F'
		
		! Now everything has been updated for the new design X
		NEW_X = 0
		
		! Obtain sensitivities of objective
		call eval_grad_F(NELM, XVAL, NEW_X, df0, IDAT, DAT, IERR)
		if (IERR.ne.0) write(*,*) 'Error in eval_grad_F'
	
		! Obtain values of constraints
		call eval_G(NELM, XVAL, NEW_X, C, FVAL, IDAT, DAT, IERR)
		if (IERR.ne.0) write(*,*) 'Error in eval_G'
		
		! Obtain sensitivities of constraints
		call eval_jac_G(TASK, NELM, XVAL, NEW_X, C, NZ, ACON, AVAR, AJAC, IDAT, DAT, IERR)
		
		! MMA requires a special sorting of the sensitivities
		do i = 1,C
			do ie = i,nelm
				ki = (ie-1)*C + i
				DFMMA(ki) = df((i-1)*nelm+ie)
			enddo
		enddo
		
		if (IERR.ne.0) write(*,*) 'Error in eval_jac_G'

		if ((irestart.eq.0).and.(diffCheck.eq.0)) then
			! Solve the subproblem
			call MMASUB(IMMA,C,NELM,GEPS,GMOVE,IYFREE,XVAL,RHO, & 
              		XMIN,XMAX,XOLD1,XOLD2,XLOW,XUP, &
             		ALFA,BETA,AMMA,BMMA,CMMA,YMMA,Z,ULAM, &
              		F0VAL,FVAL,FMAX,DF0,DFMMA, &
              		PMMA,QMMA,P0MMA,Q0,UU,GRADF,DSRCH,HESSF)
		endif
		
		write(*,*) 'MMA update completed'
		XOLD2 = XOLD1 ! Design variables two iterations ago
		XOLD1 = XVAL  ! Design variables one iteration ago
		XVAL = RHO
		
		NEW_X = 1
		
		! Calculate the design residual
		designRes = XVAL - XOLD1 
		designResidual = dsqrt(dot_product(designRes,designRes))
		
		call update_PDE_filter(XVAL, NEW_X)
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
		write(*,*) 'Min/max design variable: '
		write(*,*) minval(rho), -minval(-rho)
		
		write(*,*) 'Min/max filtered design variable: '
		write(*,*) minval(rhotilde), -minval(-rhotilde)
		write(*,*)
		
		if (diffCheck.eq.1) then
			write(*,*) '----------------------------------------------------'	
			write(*,*) 'df0(ipert): ', df0(ipert)
			write(*,*) 'df1(ipert): ', df(ipert)
			if (C.eq.2) write(*,*) 'df2(ipert): ', df(nelm+ipert)
			write(*,*) '----------------------------------------------------'
		endif
		
		write(*,*) '##################################################################'

		
		! Callback
		call iter_CB(ALG_MODE, IMMA, F0VAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT, DAT, ISTOP)
		if (IERR.ne.0) write(*,*) 'Error in iter_CB'
		
						
		if ((irestart.eq.1).or.(diffCheck.eq.1)) then
			write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
			write(*,*) '                   Stopping due to restart!                       '
			write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
			stop
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
   	write(*,*) 'The error code is ', IERR
   	write(*,*)
	endif

 9000 continue

	stop

 9990 continue
	write(*,*) 'Error setting an option'
	goto 9000
      
	return
end subroutine MMA 




!=======================================================================

!                   Computation of objective function

!=======================================================================
subroutine eval_F(NELM, XVAL, NEW_X, F, IDAT, DAT, IERR)
	implicit none
   integer										:: NELM, NEW_X
   double precision							:: F, F0VAL, XVAL(NELM)
   double precision							:: DAT(1)
   integer										:: IDAT(1)
   integer										:: IERR
   
   write(*,*) 'Entering eval_F'
   	
   call update_PDE_filter(XVAL, NEW_X)
   	
   call update_params(NEW_X,IDAT)
      	    
   call eval_NR(NEW_X, DAT, IDAT)
      	
   call eval_F0(NEW_X, F0VAL, IDAT, DAT, IERR)
      
	F = F0VAL/Scaleg0
	
	f0valhist(IDAT(1)) = F
		
	IERR = 0
	return
end subroutine eval_F


!=======================================================================

!                   Computation of objective function
!						  			(End-compliance)

!=======================================================================
subroutine eval_F0(NEW_X, F0VAL, IDAT, DAT, IERR)
	implicit none
   integer										:: NEW_X
   double precision							:: F0VAL
   double precision							:: DAT(1)
   integer										:: IDAT(1)
   integer										:: IERR
      
	F0val = dot_product(Pattern,a)
		
	IERR = 0
	return
end subroutine eval_F0


!=======================================================================

!               Computation of gradient of objective function

!=======================================================================
subroutine eval_grad_F(NELM, XVAL, NEW_X, GRAD, IDAT, DAT, IERR) 
	implicit none
	integer										:: NELM, NEW_X, I
	double precision							:: GRAD(NELM), XVAL(NELM)
	double precision							:: DAT(1)
	integer										:: IDAT(1)
	integer										:: IERR
   
   write(*,*) 'Entering eval_grad_F'
	      
	call update_PDE_filter(XVAL, NEW_X)
			
	call update_params(NEW_X, IDAT)

	call eval_NR(NEW_X, DAT, IDAT)
	   	 
	call eval_grad_F0(NEW_X, GRAD, IDAT, DAT, IERR)

	GRAD = GRAD/scaleg0
	
	df0 = GRAD

	IERR = 0
	return
end subroutine eval_grad_F


!=======================================================================

!               Computation of gradient of objective function
!						  			(End-compliance)

!=======================================================================
subroutine eval_grad_F0(NEW_X, GRADF0, IDAT, DAT, IERR)
	implicit none
	integer										:: NEW_X, I
	double precision							:: GRADF0(NELM)
	double precision							:: DAT(1)
	integer										:: IDAT(1)
	integer										:: IERR
   integer										:: ie, NegativeE

	GRADF0 = 0d0
	
	! Obtain adjoint
	varphi = Pattern
	call solveq(K,varphi,bcnod,bcval,dofnod,NegativeE,KSYM,cellSym)
	call extract(edvarphi,varphi,enod,dofnod)  
	call extract(ed,a,enod,dofnod)
	
	ffilter = 0d0
	! Calculate dFint/drho
	call dFintdrho(ffilter,edvarphi)	

	! First obtain af from Kf af = Ff
	call solveq(Kfilter,ffilter,negativeE)
	
	! Then obtain df0dz by scalarproduct with T    
	do ie=1,nelm 
   	GRADF0(ie) = dot_product(Tmatrix(:,ie),ffilter(enod(:,ie)))
	enddo 

	IERR = 0
	return
end subroutine eval_grad_F0


!=======================================================================

!                     Computation of constraints

!=======================================================================

subroutine eval_G(NELM, XVAL, NEW_X, CONST, FVAL, IDAT, DAT, IERR)
	implicit none
	integer										:: NELM, NEW_X, CONST, I
	double precision							:: FVAL(CONST), XVAL(NELM)
	double precision							:: DAT(1)
	integer										:: IDAT(1)
	integer										:: IERR

	write(*,*) 'Entering eval_G'

	call update_PDE_filter(XVAL, NEW_X)
			
	call update_params(NEW_X, IDAT)
	
	call eval_NR(NEW_X, DAT, IDAT)
	
	FVAL = 0d0
	if (C.eq.1) then
		! Usually the first constraint is the volume constraint but this might of course differ
		call eval_GVOL(NEW_X, Gvol, IDAT, DAT, IERR)
		FVAL(1) = Gvol
		f1valhist(IDAT(1)) = Gvol
	elseif (C.eq.2) then
		! Here you constraint nbr 2 and so on...
	endif

   IERR = 0
   
	return
end subroutine eval_G


!=======================================================================

!                Computation of volume-constraint

!=======================================================================
subroutine eval_GVOL(NEW_X, GvolCurr, IDAT, DAT, IERR)
	implicit none
	integer										:: NEW_X
	double precision							:: DAT(1)
	integer										:: IDAT(1)
	integer										:: IERR
   integer										:: ie, igp, NegativeE
   double precision							:: rhogp(8), Scalar(8)
   double precision							:: g, tmp, GvolCurr	
 
	! Leaving this code...
	GvolCurr = -Vtot/Vtot
	scalar = 0d0

	do ie=1,nelm   
		call fl3d8_gp(rhogp,edrho(:,ie)) 
		do igp=1,ngp
			if (Hprojection.eq.1) then
		   	call Hfun(Scalar(igp),rhogp(igp))
			else
		   	Scalar(igp) = rhogp(igp)
			endif    
		enddo 
		call fl3d8_vi(tmp,coord(:,enod(:,ie)),scalar)
		GvolCurr = GvolCurr + tmp/Vtot
	enddo

   IERR = 0
	return
end subroutine eval_GVOL



!=======================================================================

!               Computation of Jacobian of constraints
!						(i.e. sensitivities of constraints)
! 	   	 The setup follows what's required for usage of IPOPT

!=======================================================================

subroutine eval_jac_G(TASK, NELM, XVAL, NEW_X, CONST, NZ, ACON, AVAR, AJAC, IDAT, DAT, IERR)
	implicit none
	integer										:: TASK, NELM, NEW_X, CONST, NZ, I, J
	double precision							:: XVAL(NELM), AJAC(NZ)
	integer										:: ACON(NZ), AVAR(NZ)
	double precision							:: DAT(1)
	integer										:: IDAT(1)
	integer										:: IERR

	write(*,*) 'Entering eval_jac_G'

	if (TASK.eq.0) then
		write(*,*) 'Task = 0'
		do I = 1, NZ
			! Col (design variable)
      	AVAR(I) = AVAR1(I)
			! Row (constraint)
      	ACON(I) = ACON1(I)
      enddo
	else
		write(*,*) 'Task != 0'
		
		call update_PDE_filter(XVAL, NEW_X)
		
		call update_params(NEW_X, IDAT)

		call eval_NR(NEW_X, DAT, IDAT)
		
		
		if (C.eq.1) then
			call eval_grad_GVOL(NEW_X, IDAT, DAT, IERR)
		elseif (C.eq.2) then
			call eval_grad_GVOL(NEW_X, IDAT, DAT, IERR)
			! Second constraint and so on...
		endif
		
		AJAC = 0d0
		df = 0d0
		do I = 1, C
   		do J = 1, nelm
   			! Should change this....
				if (I.eq.1) then
      			AJAC((I-1)*nelm+J) = dGvol(J)
      			df((I-1)*nelm+J) = dGvol(J)
      		elseif (I.eq.2) then
      			! Second constraint and so on...
      		endif
      	enddo
   	enddo

	endif
	
	IERR = 0
	return
end subroutine eval_jac_G

		
!=======================================================================

!                Gradient of volume-constraint

!=======================================================================
subroutine eval_grad_GVOL(NEW_X, IDAT, DAT, IERR)
	implicit none
	integer										:: NEW_X
	double precision							:: DAT(1)
	integer										:: IDAT(1)
	integer										:: IERR
   integer										:: ie, igp, NegativeE
   double precision							:: rhogp(8), Scalar(8), fen(8)
   double precision							:: g

	! Leaving this code...
	scalar = 0d0 
	dGvol = 0d0
	ffilter=0D0
	
	do ie=1,nelm  
		call fl3d8_gp(rhogp,edrho(:,ie)) 
		do igp=1,ngp
			if (Hprojection.eq.1) then
		   	call Hprimfun(Scalar(igp),rhogp(igp))
			else
		   	Scalar(igp)=1d0
		   endif   
		enddo
		call fl3d8_bf(fen,coord(:,enod(:,ie)),Scalar) 
		ffilter(enod(:,ie))=ffilter(enod(:,ie))+fen            
	enddo
	
	call solveq(Kfilter,ffilter,negativeE)

	do ie=1,nelm
		dGvol(ie) = dot_product(Tmatrix(:,ie),ffilter(enod(:,ie)))/Vtot
	enddo 
 
   IERR = 0
	return
end subroutine eval_grad_GVOL



!=======================================================================

!               Computation of Hessian of Lagrangian
!					 !!!			This is not used		 !!!

!=======================================================================

subroutine eval_HESS(TASK, NELM, XVAL, NEW_X, OBJFACT, CONST, LAM, NEW_LAM, NNZH, IRNH, ICNH, HESS, IDAT, DAT, IERR)
	implicit none
	integer TASK, NELM, NEW_X, CONST, NEW_LAM, NNZH, i
	double precision XVAL(NELM), OBJFACT, LAM(CONST), HESS(NNZH)
	integer IRNH(NNZH), ICNH(NNZH)
	double precision DAT(1)
	integer IDAT(1)
	integer IERR
   
	write(*,*) 'Hessians are not supported at this point'
  
   IERR = 1
	return
end subroutine eval_HESS


!=======================================================================

!                  Callback method called once per iteration
!							  E.g. Outputting data to Matlaba

!=======================================================================

subroutine iter_CB(ALG_MODE, ITER_COUNT, OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT, DAT, ISTOP)
	implicit none
	integer 						:: ALG_MODE, ITER_COUNT, LS_TRIAL
	double precision			:: OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE
	double precision			:: ALPHA_DU, ALPHA_PR
	double precision			:: DAT(1)
	integer 						:: IDAT(1)
	integer 						:: ISTOP
	
	update_done = 0

	! Only needed for IPOPT
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
		     
	write( s_plot, '(i10)' )  ITER_COUNT
		     
	if (ITER_COUNT.lt.10) xstring=xstring(1:indx1)//'_000' //s_plot(10:10) //''' )'
	if (ITER_COUNT.ge.10) xstring=xstring(1:indx1)//'_00' //s_plot(9:10) //''' )'
	if (ITER_COUNT.ge.100) xstring=xstring(1:indx1)//'_0' //s_plot(8:10) //''' )'
	if (ITER_COUNT.ge.1000) xstring=xstring(1:indx1)//'_'//s_plot(7:10) //''' )'
	indx1=indx1+5
	xstring=xstring(1:indx1)//'.mat'
	indx1=indx1+4

	call integ2nod(esx_nod,es(1,:,:),'brick8',enod)
	call integ2nod(esy_nod,es(2,:,:),'brick8',enod)
	call integ2nod(esxy_nod,es(3,:,:),'brick8',enod)
	call integ2nod(gamma_nod,gammagp(:,:),'brick8',enod)
	call matwrt2f(xstring(1:indx1),a,'a','w')  
	call matwrt2f(xstring(1:indx1),ahist,'ahist','u')  
	call matwrt2f(xstring(1:indx1),anom,'anom','u')  
	call matwrt2f(xstring(1:indx1),acrit,'acrit','u') 
	call matwrt2f(xstring(1:indx1),esx_nod,'sigx_nod','u') 
	call matwrt2f(xstring(1:indx1),esy_nod,'sigy_nod','u') 
	call matwrt2f(xstring(1:indx1),esxy_nod,'sigxy_nod','u') 
	call matwrt2f(xstring(1:indx1),df0,'df0dx','u')  
	call matwrt2f(xstring(1:indx1),df,'dfdx','u')
	call matwrt2f(xstring(1:indx1),rhotilde,'rhotilde','u') 
	call matwrt2f(xstring(1:indx1),rho,'rho','u')  
	call matwrt2f(xstring(1:indx1),edrho,'edrho','u')  
	call matwrt2f(xstring(1:indx1),f0valhist,'f0val','u') 
	call matwrt2f(xstring(1:indx1),f1valhist,'f1val','u')  
	call matwrt2f(xstring(1:indx1),f2valhist,'f2val','u')  
	call matwrt2f(xstring(1:indx1),eigenvalhist,'eigenvalhist','u') 
	call matwrt2f(xstring(1:indx1),resvalhist,'resvalhist','u') 
	call matwrt2f(xstring(1:indx1),tmp2,'tmp2','u') 
	call matwrt2f(xstring(1:indx1),gamma_nod,'gamma_nod','u')
		     
	xstring='history'
	indx1=4
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
      
	return
end subroutine iter_CB


!=======================================================================

!           			   Update parameters

!=======================================================================
subroutine update_params(NEW_X,IDAT)
	implicit none
	integer				:: NEW_X, IDAT(1), iter

	iter = IDAT(1)
	write(*,*) 'Update parameters'
	pold = p

	if ((update_done.eq.0).and.(NEW_X.eq.1)) then
		update_done = 1
		
		! Update penalization parameter
		if (mod(iter,5).eq.0) then
			p = p + 0.05D0	
		endif
		
		if (p.gt.3D0) p=3D0
		
		! Update beta in Heaviside
		if (Hprojection.eq.1) then
			 ! Wait a few iterations
			 if (dabs(p-3D0).lt.1D-8) then
		      iwait=iwait+1
		    endif
		  	 ! Update sharpness of Heaviside filter
		  	 if (iwait.gt.20) then
		       if ((mod(iter,5).eq.0)) B = B+1D0 !2*B!
		    endif
		    if (B.gt.Bfinal) B=Bfinal
		    if (p.lt.2.97D0) B=Binit
	  endif

	  ! Update aggregation parameter
	  if (iwait.gt.20) then
		  if ((mod(iwait,50).eq.0)) pAgg = 2d0*pAgg
	  endif
	  if (pAgg.gt.pAggFinal) pAgg=pAggFinal
	  
	  p1 = p     
	  
	  tmp2(1)=p
     tmp2(2)=B
     
     call update_gamma()
     
      write(*,*) '####################################################'				
		write(*,*)
		write(*,*) '                 Parameter values                   '
		write(*,*)
		write(*,*) '####################################################'		
   	write(*,*) 'Penalty exponent', p
   	write(*,*) 'Heaviside coefficient', B 
   	write(*,*) '####################################################'
	endif

   
	return
end subroutine update_params


!=======================================================================

!           			   Update gamma

!=======================================================================
subroutine update_gamma()
	implicit none
	integer				:: ie, igp
	double precision	:: rhogp(8)

	write(*,*) 'Update gamma'
	
	gammagp = 0d0

	!$OMP PARALLEL DO DEFAULT(NONE) &
	!$OMP SHARED(nelm,edrho,gammagp,ngp) & 
	!$OMP PRIVATE(rhogp,ie,igp) 
	do ie=1,nelm
      call fl3d8_gp(rhogp,edrho(:,ie))    
      do igp=1,ngp      
         call gammafun(gammagp(igp,ie), rhogp(igp))
         if ((gammagp(igp,ie)).lt.0D0)  then
         	write(*,*) 'This xval value yields a negative gammagp', rhogp(igp)          
      	endif	
      enddo
  	enddo 
	!$OMP end parallel do

	return
end subroutine update_gamma


!=======================================================================

!                  Calculates dfint/dz

!=======================================================================
subroutine  dFintdrho(ffilterCurr,edmuCurr)	
	implicit none
	double precision						:: ffilterCurr(:), edmuCurr(:,:)
	double precision                 :: CRAP98(9,8), Ones(8), Dgp(6,6,8), gp, esscaled(6,8), Bconst(6,8), esscaledsmall(6,8)
	double precision                 :: fke(8), g, ffe(8), sfegp(8), Scalar(8), fen(8), dggamma(9,8)
   integer                          :: ie, igp
   double precision 					   :: gammap, rhogp(8)
  
   Ones=1D0   
   scalar = 0d0
  
   dg=0d0
   dg(1,:,:)=1d0
   dg(5,:,:)=1d0
   dg(9,:,:)=1d0
   es=0d0

   do ie=1,nelm
         ! Extract gauss point values from nodal quantities
         call fl3d8_gp(rhogp,edrho(:,ie)) 

			! The deformation gradient 
			call c3dtl8_d(dg(:,:,ie),coord(:,enod(:,ie)),ed(:,ie))
				  		
			! The gamma-scaled deformation gradient 
			call calcdggamma(dggamma,dg(:,:,ie),gammagp(:,ie))
     	
     		! Calculate the tangent stiffness with gamma
			call dneohooke('tl',Dgp,dggamma)
					
			! Calculate the 2nd Piola Kirchhoff stress with gamma
			call neohooke('2ndPiola',es(:,:,ie),dggamma)
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!						 NONLINEAR STIFFNESS MATRIX CONTRIBUTION
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! The material behaviour
			do igp=1,ngp
				call gfun(g,rhogp(igp))
            call gammaprimfun(gammap,rhogp(igp))
            if (gammagp(igp,ie).gt.0D0) then
					esscaled(:,igp)= g*es(:,igp,ie)*gammap/gammagp(igp,ie)  
            	Dgp(:,:,igp)   = g*Dgp(:,:,igp)*gammap/gammagp(igp,ie)  
            else
            	esscaled(:,igp)= 0d0 
            	Dgp(:,:,igp)   = 0d0
            endif
         enddo
				  	
			! Calculate the large def. element stiffness matrix contribution
			call c3dtl8_egammaSens(fke,coord(:,enod(:,ie)),Dgp,ed(:,ie),edmuCurr(:,ie),ed(:,ie),esscaled,gammagp(:,ie))	
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!						  NONLINEAR RESIDUAL CONTRIBUTION
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Thresholding
			do igp=1,ngp
            call gfun(g,rhogp(igp))  
            call gfunprim(gp,rhogp(igp))
            call gammaprimfun(gammap,rhogp(igp))
            if (gammagp(igp,ie).gt.0D0) then
            	esscaled(:,igp)=g*es(:,igp,ie)*(gp/g + gammap/gammagp(igp,ie))    
            else
            	esscaled(:,igp)= 0d0 
            endif 
         enddo
         
         ! Calculate the nonlinear element forces
         call c3dtl8_fgammaSens(ffe,coord(:,enod(:,ie)),ed(:,ie),edmuCurr(:,ie),esscaled,gammagp(:,ie)) 
         ! Add the contributions
         fen = fke + ffe
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !							 LINEAR RESIDUAL CONTRIBUTION
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculate Bconst
         call c3dtl8_emu(coord(:,enod(:,ie)),0D0*ed(:,ie),ed(:,ie),Bconst,Crap98,Ones) 
              		
         ! Calculate small strain contribution
         do igp=1,ngp
            call gfun(g,rhogp(igp))
            call gfunprim(gp,rhogp(igp))
            call gammaprimfun(gammap,rhogp(igp))
            
            ! Calculate sigma (Cauchy stress)
            if (gammagp(igp,ie).gt.0D0) then
            	esscaledsmall(:,igp) = matmul(Dsmall,Bconst(:,igp))
            	esscaledsmall(:,igp) = g*((1D0-gammagp(igp,ie)**2D0)*gp/g-2d0*(gammagp(igp,ie)**2D0)*gammap/gammagp(igp,ie))*esscaledsmall(:,igp)
            else
            	esscaled(:,igp)= 0d0 
            endif 
         enddo

			! Calculate the linear element forces 
			call c3dtl8_fgammaSens(ffe,coord(:,enod(:,ie)),0D0*ed(:,ie),edmuCurr(:,ie),esscaledsmall,ones) 
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         ! Add the contributions
			fen = fen + ffe

			! Assemble. ffilter is dfdrhotilde
	 		ffilterCurr(enod(:,ie)) = ffilterCurr(enod(:,ie)) - fen			  
   enddo

   return
end subroutine dFintdrho






!=======================================================================

!                Solves the equilbrium iterations

!=======================================================================
subroutine eval_NR(NEW_X, DAT, IDAT)
   implicit none
   double precision                    :: DAT(1)
   double precision							:: fe(24), fesmall(24), kesmall(24,24), ke(24,24)
   double precision							:: dggamma(9,8), Dgp(6,6,8), esscaled(6,8), esscaledsmall(6,8)
   double precision							:: rhogp(8), CRAP98(9,8), Ones(8), Bconst(6,8)
   double precision							:: residual, g
   integer                             :: cDone, NEW_X, n, i, igp, imax, ie, NegativeE, IDAT(1)
   double precision			  	 			:: a1, a2, a3, a4, a5, sol1, sol2, dlambda, lambda, lambdaOld
   integer                             :: ii, jj, ic
   integer										:: Lnew, info
   
	imax = 8
	Ones=1D0
	
	restart = 0

	if (NEW_X.eq.1) then
	  
     edphi = 0d0
	  edmu=0d0
	  history = 0d0
	  ahist = 0d0
	     
     ed=0d0
     es=0d0
     dg=0d0
     dg(1,:,:)=1d0
     dg(5,:,:)=1d0
     dg(9,:,:)=1d0
     a=0d0
     res=0d0
	  
     Lambda=0D0
     LambdaOld=0D0
     
     Deltaa=0D0
     DeltaLambda=LambdaMax/nmax

	  n = 0
	  
	  write(*,*) '###################################'
	  write(*,*) '      Start of NR iterations  '	 
	  write(*,*) '          Load-control		' 
	  write(*,*) '###################################' 
	  
     ! Starting the for-loop (NR)
	  ! Equilibrium loop starts
     Do while (n.lt.nmax)  
       	! Variable checking which load step we're on
  		 	n=n+1

			! Update load
			Lambda = Lambda + DeltaLambda

			write(*,*) 'Starting equlibrium iter at loadstep, lambda:', n, lambda
		
			res = 0d0
       	residual=1d0
       	i=0

       	do while ((i.lt.(imax)).and.(residual.gt.TOL))
		      ! Some variable that checks which iteration we're one
            i=i+1
            
		      K%a = 0d0
		      Fvec = 0d0
		      
		      ! Calculate stiffness matrix
		      do ie=1,nelm  
		        		! Extract gauss point values from nodal quantities
		      	  	call fl3d8_gp(rhogp,edrho(:,ie)) 

					  	! The deformation gradient 
					  	call c3dtl8_d(dg(:,:,ie),coord(:,enod(:,ie)),ed(:,ie))
					  	
					  	! The gamma-scaled deformation gradient 
					  	call calcdggamma(dggamma,dg(:,:,ie),gammagp(:,ie))

					  	! Calculate the tangent stiffness with gamma
						call dneohooke('tl',Dgp,dggamma)
						
					  	! Calculate the 2nd Piola Kirchhoff stress with gamma
					  	call neohooke('2ndPiola',es(:,:,ie),dggamma)

					  	! Now determine the material behaviour
					  	do igp=1,ngp
							call gfun(g,rhogp(igp)) 
							esscaled(:,igp)= g*es(:,igp,ie)*gammagp(igp,ie)**2D0
		               Dgp(:,:,igp)   = g*Dgp(:,:,igp)*gammagp(igp,ie)**2D0
		            enddo
					  	
					  	! Calculate the large def. element stiffness matrix with gamma
						call c3dtl8_egamma(ke,coord(:,enod(:,ie)),Dgp,ed(:,ie),esscaled,gammagp(:,ie))

						! Small strain contribution
						do igp=1,ngp 
		             	call gfun(g,rhogp(igp))
		             	Dgp(:,:,igp)=g*(1D0-gammagp(igp,ie)**2D0)*Dsmall  
		          	enddo
		          	
		          	! Calculate the small def. element stiffness matrix with gamma contribution
		          	call c3dtl8_e(kesmall,coord(:,enod(:,ie)),Dgp,0D0*ed(:,ie),0D0*esscaled)
		          	
					  	! Assembling
					  	! ke = gamma^2 (B'DB + G'YG), ksmall = (1-gamma^2) Bc'DBc
					  	call assem(K,ke+kesmall,enod(:,ie),dofnod)
	  			enddo

		      deltaa=-res
		      
		      call solveq(K,deltaa,bcnod,bcval,dofnod,NegativeE,KSYM,cellSym) 
    
		      a = a + Deltaa
		      
		      call extract(ed,a,enod,dofnod)
		      
		      ! Calculate internal forces
				do ie=1,nelm
		        		! Extract gauss point values from nodal quantities
		      	  	call fl3d8_gp(rhogp,edrho(:,ie)) 

					  	! The deformation gradient 
					  	call c3dtl8_d(dg(:,:,ie),coord(:,enod(:,ie)),ed(:,ie))
					  		
					  	! The gamma-scaled deformation gradient 
					  	call calcdggamma(dggamma,dg(:,:,ie),gammagp(:,ie))

					  	! Calculate the 2nd Piola Kirchhoff stress
					 	call neohooke('2ndPiola',es(:,:,ie),dggamma)

					  	! Thresholding
					  	do igp=1,ngp
		               call gfun(g,rhogp(igp))
		               esscaled(:,igp)=g*es(:,igp,ie)*gammagp(igp,ie)    
		           	enddo

		           	! Calculate Bconst
		           	call c3dtl8_emu(coord(:,enod(:,ie)),0D0*ed(:,ie),ed(:,ie),Bconst,Crap98,Ones) 

		           	! Calculate small strain contribution
		           	do igp=1,ngp
		           		call gfun(g,rhogp(igp))
		           		! Calculate sigma (Cauchy stress)
		            	esscaledsmall(:,igp) = matmul(Dsmall,Bconst(:,igp))	
		            	esscaledsmall(:,igp) = g*(1D0-gammagp(igp,ie)**2D0)*esscaledsmall(:,igp)
		          	enddo

					  	! Calculate element forces (linear + non-linear)
					  	call c3dtl8_f(fesmall,coord(:,enod(:,ie)),0D0*ed(:,ie),esscaledsmall)           
		         	call c3dtl8_fgamma(fe,coord(:,enod(:,ie)),ed(:,ie),esscaled,gammagp(:,ie)) 
					  
					  	! Assembling
					  	call insert(Fvec,fe+fesmall,enod(:,ie),dofnod)
				enddo
			
		      res = Fvec - Lambda*Pattern
		      res(bcdof) = 0d0
		      
		      residual=dsqrt(dot_product(res,res))                  
		      write(*,*) 'Residual:', residual                                     
			enddo ! End of equlibrium loop

       	history(i,:) = a
     enddo ! End of load steps 

	  if (isnan(residual)) then
			stop
	  endif
	  
	  write(*,*) 'End of NR'
   
   
   endif ! End if(NEW_X.eq.1)
  
end subroutine eval_NR





!=======================================================================

!                    Calculates the Xi-function

!=======================================================================
subroutine gfun(val,rhogpCurr)
   implicit none
   double precision                 ::  rhogpCurr, val2, val
	
	val2 = 0d0
	val = 0d0

   if (Hprojection.eq.1) then
   	call Hfun(val2,rhogpCurr)
   else
   	val2 = rhogpCurr
   endif

   ! If Helmholtz filter is used xval can be negative
   if (val2.ge.0D0) then
   	val = (1D0-delta0)*val2**p+delta0
   else
   	val = delta0
   endif


   return
end subroutine gfun


  
!=======================================================================

!           Calculates the derivative of Xi-function

!=======================================================================
subroutine gfunprim(val,rhogpCurr)
   implicit none
   double precision          :: rhogpCurr, val2, val3, val

	val = 0d0
	val2 = 0d0
	val3 = 0d0

   if (Hprojection.eq.1) then
   	call Hfun(val2, rhogpCurr)
   	call Hprimfun(val3, rhogpCurr)
   	val=p*(1D0-delta0)*val2**(p-1D0)*val3
   else
   	val=p*(1D0-delta0)*rhogpCurr**(p-1D0)
   endif

   if (rhogpCurr.lt.0D0)  val=0D0

   return
end subroutine gfunprim



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
   double precision                    :: dggamma(9,8), dge(9,8), gammagpe(8)
   integer                             :: igp

   do igp=1,ngp   
     dggamma(:,igp)=UnitVec+(dge(:,igp)-UnitVec)*gammagpe(igp)
   enddo

end subroutine calcdggamma


!=======================================================================

!                 Calculates the gamma(=Heaviside) function

!=======================================================================
subroutine gammafun(val, rhogpCurr)
   implicit none
   double precision			:: val2, rhogpCurr, val

	val = 0d0
	val2 = 0d0

   if (Hprojection.eq.1) then
   	call Hfun(val2,rhogpCurr)
   else
   	val2 = rhogpCurr
   endif

   val = 0D0
   if (val2.gt.1D-12) then
   	val = (dtanh(B1*rho0)+dtanh(B1*(val2**p1-rho0)))/(dtanh(B1*rho0)+dtanh(B1*(0.1D1-rho0)))   
   endif

   !  val=1D0
  
   return
end subroutine gammafun



!=======================================================================

!      Calculates the derivative of the gamma(=Heaviside)-function

!=======================================================================
subroutine gammaprimfun(val,rhogpCurr)
   implicit none
   double precision                 :: val2, val3, rhogpCurr, val

	val = 0d0
	val2 = 0d0
	val3 = 0d0

   if (Hprojection.eq.1) then
   	call Hfun(val2,rhogpCurr)
   else
   	val2 = rhogpCurr
   endif

   val=0D0
   if (val2.gt.1d-12) then
      val = (0.1D1-dtanh(B1*(val2**p1-rho0))**2D0) *  & 
                    B1*val2**p1*p1/val2/(dtanh(B1*rho0) + dtanh(B1 * (0.1D1 - rho0)))
   endif

   if (Hprojection.eq.1) then
   	call Hprimfun(val3, rhogpCurr)
   	val = val*val3
   endif

   ! val=0D0

   return
end subroutine gammaprimfun



end module opt_module








