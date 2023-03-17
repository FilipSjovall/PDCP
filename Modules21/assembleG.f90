module assembleG

use omp_lib
use abaqus_util
use matlab_util
use sparse_util
use mater_hyperel
use elem_large_cont_3d
use elem_flow_3d
use plani4axi
use Spl_geom
use shape_sens

implicit none

contains

subroutine assemGlobal(K,F,coord,ed,ke,fe,dg,ef,J_OUT,edp,temp,edof,enod,Dgp,esh,es,temp2,nelm,ie,dofnod,ierr,ngp,ndof)
type(sparse), intent(inout)      :: K
double precision, intent(inout)  :: F(:)
double precision, intent(inout)  :: coord(:,:), ed(:,:) 
double precision, intent(inout)  :: ke(:,:), fe(:), dg(:,:,:,:), ef(:,:)
double precision, intent(inout)  :: J_OUT(:,:), edp(:,:), temp(:,:)
integer, intent(inout)			 :: edof(:,:), enod(:,:)


double precision, intent(inout)  :: Dgp(4,4,4),esh(4,4), es(4,3,3), temp2(8)
integer, intent(inout)           :: nelm, ie, dofnod, ierr, ngp, ndof
integer						     :: el,gp

! Deformation gradient for each element	
edp = 0d0		
edp = transpose(ed)			
K%a=0d0
F=0d0
call	DEFORMATION4I_AXI(dg,J_out,coord,edof,edp,ngp,nelm,ndof)
    ! Calculate stiffness and mass matrix
	Ke = 0d0
	Fe = 0d0
    do ie=1,nelm

	    ! use temp to avoid  using non-contiguous call of coord
	    temp  = coord(enod(:,ie),:)
	    temp2 = ed(:,ie)
	      ! == 
	    do gp=1,ngp
	        ef(:,gp) = (/ dg(ie,gp,1,1), dg(ie,gp,1,2), dg(ie,gp,2,1), dg(ie,gp,2,2), dg(ie,gp,3,3) /)
	    end do	
	    ! == 								
	    ! Calculate the tangent stiffness [nelm x ngp x 3 x 3]
	    call dneohooke('tl',Dgp,ef) 
        ! Calculate the 2nd Piola Kirchhoff stress
        call neohooke('2ndPiola',esh,ef)        

		! reshape stress tensor
	    do gp=1,ngp
	    	es(gp,1,:) = (/ esh(1,gp), esh(4,gp), 0d0  /)  
	    	es(gp,2,:) = (/ esh(4,gp), esh(2,gp), 0d0  /) 
	    	es(gp,3,:) = (/ 0d0     , 0d0     , esh(3,gp) /)									 
	    end do
		
        ! Calculate the element stiffness matrix
	    call PLANI4E_AXI(Ke,temp,Dgp,temp2,es,ngp) 
        call PLANI4F_AXI(fe,temp,temp2,es,ngp)

        ! Assembling
        call insert(F,fe,enod(:,ie),dofnod)	
	    call assem(K,ke,enod(:,ie),dofnod)

		!if(ie.eq.1) then
		!	print * , "ef", ef(:,:)
		!	print * , "coord", temp
		!	!print * , "D" , Dgp
		!	print * , "ed" , temp2
		!	print * , "es1" , esh(:,1)
		!	print * , "es2" , esh(:,2)
		!	print * , "es3" , esh(:,3)
		!	print * , "es4" , esh(:,4)
		!	!print * , "ngp", ngp
		!	print * , "Ke is assemGlobal " , Ke
		!endif

    end do ! End calculation of stiffness matrix
end subroutine assemGlobal

subroutine assemGlobalF(F,coord,ed,fe,dg,ef,J_OUT,edp,temp,edof,enod,Dgp,esh,es,temp2,nelm,ie,dofnod,ierr,ngp,ndof)
	
	double precision, intent(inout)  :: F(:)
	double precision, intent(inout)  :: coord(:,:), ed(:,:) 
	double precision, intent(inout)  :: fe(:), dg(:,:,:,:), ef(:,:)
	double precision, intent(inout)  :: J_OUT(:,:), edp(:,:), temp(:,:)
	integer, intent(inout)			 :: edof(:,:), enod(:,:)
	
	
	double precision, intent(inout)  :: Dgp(4,4,4),esh(4,6), es(4,3,3), temp2(8)
	integer, intent(inout)           :: nelm, ie, dofnod, ierr, ngp, ndof
	integer						     :: el,gp
	
	
	
	! Deformation gradient for each element	
	edp = 0d0		
	edp = transpose(ed)		
	F=0d0
	call	DEFORMATION4I_AXI(dg,J_out,coord,edof,edp,ngp,nelm,ndof)
		! Calculate stiffness and mass matrix
	
		do ie=1,nelm  
	
			! use temp to avoid  using non-contiguous call of coord
			temp  = coord(enod(:,ie),:)
			temp2 = ed(:,ie)
			  ! == 
			do gp=1,ngp
				ef(:,gp) = (/ dg(ie,gp,1,1), dg(ie,gp,1,2), dg(ie,gp,2,1), dg(ie,gp,2,2), dg(ie,gp,3,3) /)
			end do	
			! == 								
			! Calculate the tangent stiffness [nelm x ngp x 3 x 3]
			call dneohooke('tl',Dgp,ef) 
			! Calculate the 2nd Piola Kirchhoff stress
			call neohooke('2ndPiola',esh,ef)        
	
			! reshape stress tensor
			do gp=1,ngp
				es(gp,1,:) = (/ esh(1,gp), esh(4,gp), 0d0  /)  
				es(gp,2,:) = (/ esh(4,gp), esh(2,gp), 0d0  /) 
				es(gp,3,:) = (/ 0d0     , 0d0     , esh(3,gp) /)									 
			end do
			
			! Calculate the element stiffness matrix
			call PLANI4F_AXI(fe,temp,temp2,es,ngp)
	
			! Assembling
			call insert(F,fe,enod(:,ie),dofnod)	
	
		end do ! End calculation of stiffness matrix
end subroutine assemGlobalF

subroutine assemSens(K_hist,DFe,dVe,dXdalpha,nelm,temp,temp2,tempeu,tempev,ngp,dg,Dgp,ef,es,esh,Uk,Vk,DV,npar,dF,dVol,ke,enod,edof,eu,ev,coord,ed,n,dofnod,Vmax)
	double precision, intent(inout) :: K_hist(:,:,:), DFe(:,:,:), dVe(:,:), dXdalpha(:,:,:), temp(:,:), temp2(:), tempeu(:), tempev(:), ed(:,:), Vmax, Dgp(:,:,:)
	double precision, intent(inout) :: dg(:,:,:,:), ef(:,:), es(:,:,:), esh(:,:), Uk(:), Vk(:), DV(:,:,:,:), dF(:,:,:,:), dVol(:), ke(:,:), eu(:,:), ev(:,:), coord(:,:)

	integer, intent(inout)          :: nelm, ngp, npar, enod(:,:), edof(:,:), n, dofnod
	integer 						:: ie, gp, ipar
		
		K_hist(n,:,:) = 0d0
		DFe = 0d0
		dVol= 0d0
		dXdalpha=0d0
		do ie=1,nelm

			! use temp to avoid  using non-contiguous calls
			temp   = coord(enod(:,ie),:)
			temp2  = ed(:,ie)
			tempeu = eu(ie,:)  
			tempev = ev(ie,:) 	

			! Reshape deformation gradient 
			do gp=1,ngp
				ef(:,gp) = (/ dg(ie,gp,1,1), dg(ie,gp,2,1), dg(ie,gp,1,2), dg(ie,gp,2,2), dg(ie,gp,3,3) /)
			end do									
					
			! Calculate the tangent stiffness [nelm x ngp x 3 x 3]
			call dneohooke('tl',Dgp,ef)

			! Calculate the 2nd Piola Kirchhoff stress
			call neohooke('2ndPiola',esh,ef)     


		  ! reshape stress tensor
			do gp=1,ngp
				es(gp,1,:) = (/ esh(1,gp), esh(4,gp), 0d0       /)  
				es(gp,2,:) = (/ esh(4,gp), esh(2,gp), 0d0       /) 
				es(gp,3,:) = (/ 0d0      , 0d0      , esh(3,gp) /)									 
			end do						

			! Calculate sensitivity of nodes w.r.t to control vertix parameters
			call Nodesens(dXdalpha,tempeu,tempev,Uk,Vk,DV,npar)

			! Calculate sensitivities of volume and internal force
			call plani4se_nl(DVe,DFe,Dgp,ngp,npar,dXdalpha,temp2,es,temp)				
			
			do ipar=1,4*npar
			 !dF(ipar,:,:) = dF(ipar,:,:) + DFe(ipar)								
				call insert(dF(n,ipar,:,1),DFe(ipar,:,1),enod(:,ie),dofnod)									
				dVol(ipar)     = dVol(ipar) + DVe(ipar,1)/Vmax 						
			end do						
				
! Calculate the element stiffness matrix
			call PLANI4E_AXI(Ke,temp,Dgp,temp2,es,ngp) 								

			! Assemble
			call assem(K_hist(n,:,:),ke,enod(:,ie),dofnod)
		end do
end subroutine assemSens

subroutine adjoint(M,K_hist,It,f_hist,f_tar,df0dx,coeff,f0val,bcval,lam,ii,bcdof,nloadsteps,npar,ipar,dF)
	double precision, intent(inout)	:: M(:), K_hist(:,:,:), It(:), f_hist(:,:), f_tar(:), df0dx(:), coeff(:), f0val, bcval(:), lam(:), dF(:,:,:,:)
	Integer, intent(inout)	        :: ii, bcdof(:), nloadsteps, ipar, npar
	
		do ii = 1,nloadsteps
			M       = 2*matmul(K_hist(ii,:,:),It)*(sum(It*f_hist(ii,:))-f_tar(ii))
			lam     = -M
			call solveq(K_hist(ii,:,:),lam,bcdof,0*bcval)
			coeff   = lam + 2*It*(sum(It*f_hist(ii,:))-f_tar(ii)) 				
			f0val   = f0val + (sum(It*f_hist(ii,:))-f_tar(ii))**2
			! Loop over each design variable				
			do ipar = 1,4*npar
				df0dx(ipar) = df0dx(ipar) + sum(coeff*dF(ii,ipar,:,1))
			end do
		 end do
end subroutine adjoint

subroutine assemGlobalPar(K,F,coord,ed,ke,fe,dg,ef,J_OUT,edp,temp,edof,enod,Dgp,esh,es,temp2,nelm,ie,dofnod,ierr,ngp,ndof)
	type(sparse)      :: K
	double precision  :: F(:)
	double precision  :: coord(:,:), ed(:,:) 
	double precision  :: ke(:,:), fe(:), dg(:,:,:,:), ef(:,:)
	double precision  :: J_OUT(:,:), edp(:,:), temp(:,:)
	integer			  :: edof(:,:), enod(:,:)


	double precision  :: Dgp(4,4,4),esh(4,4), es(4,3,3), temp2(8)
	integer           :: nelm, ie, dofnod, ierr, ngp, ndof
	integer			  :: el,gp

	! Deformation gradient for each element	
	edp = 0d0		
	edp = transpose(ed)			
	K%a=0d0
	F=0d0
	call	DEFORMATION4I_AXI(dg,J_out,coord,edof,edp,ngp,nelm,ndof)

    ! Calculate stiffness and mass matrix´ !!collapse(1)
	!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
	!$OMP SHARED(nelm,ngp,dg,coord,enod,ed,K,dofnod) &
	!$OMP PRIVATE(ie,gp,ef,Dgp,esh,es,Ke,fe,temp,temp2) &
	!$OMP REDUCTION(+:F)
    do ie=1,nelm  

	    ! == 
		temp  = coord(enod(:,ie),:)
	    temp2 = ed(:,ie)

	    do gp=1,ngp
	        ef(:,gp) = (/ dg(ie,gp,1,1), dg(ie,gp,1,2), dg(ie,gp,2,1), dg(ie,gp,2,2), dg(ie,gp,3,3) /)	
	    enddo	

	    ! == 								
	    ! Calculate the tangent stiffness [nelm x ngp x 3 x 3]
	    call dneohooke('tl',Dgp,ef) 
	
        ! Calculate the 2nd Piola Kirchhoff stress
        call neohooke('2ndPiola',esh,ef)        

		! reshape stress tensor
	    do gp=1,ngp
	    	es(gp,1,:) = (/ esh(1,gp), esh(4,gp), 0d0  /)  
	    	es(gp,2,:) = (/ esh(4,gp), esh(2,gp), 0d0  /) 
	    	es(gp,3,:) = (/ 0d0     , 0d0     , esh(3,gp) /)									 
	    enddo

		
        ! Calculate the element stiffness matrix
		call PLANI4E_AXI(Ke,temp,Dgp,temp2,es,ngp) 
        call PLANI4F_AXI(fe,temp,temp2,es,ngp)


		if(ie.eq.1) then
			!print * , "Ke", Ke
			!print * , "Fe", fe
			!print * , "ef-3", ef(:,3)
			!print * , "ef-4", ef(:,4)
		endif
        ! Assembling
        call insert(F,fe,enod(:,ie),dofnod)	

		!$OMP CRITICAL
	    call assem(K,ke,enod(:,ie),dofnod)
		!$OMP END CRITICAL
			!!!! SHARED(nelm,ngp,dg,coord,enod,ed,K,F,dofnod) 
			!!!! PRIVATE(ie,gp,ef,Dgp,esh,es,Ke,fe,temp,temp2) 
		!print * , "ie"  , ie
		!print * , "gp", gp
		!print * , "ef" , ef(:,1)
		!print * , "D"  , Dgp(:,:,1)
		!print * , "esh", esh(:,1)
		!print * , "es" , es(1,:,:)

    enddo ! End calculation of stiffness matrix
	!$OMP END PARALLEL DO


end subroutine assemGlobalPar

subroutine assemGlobalParF(F,coord,ed,fe,dg,ef,J_OUT,edp,temp,edof,enod,Dgp,esh,es,temp2,nelm,ie,dofnod,ierr,ngp,ndof)
	
	double precision, intent(inout)  :: F(:)
	double precision, intent(inout)  :: coord(:,:), ed(:,:) 
	double precision, intent(inout)  :: fe(:), dg(:,:,:,:), ef(:,:)
	double precision, intent(inout)  :: J_OUT(:,:), edp(:,:), temp(:,:)
	integer, intent(inout)			 :: edof(:,:), enod(:,:)
	
	
	double precision, intent(inout)  :: Dgp(4,4,4),esh(4,6), es(4,3,3), temp2(8)
	integer, intent(inout)           :: nelm, ie, dofnod, ierr, ngp, ndof
	integer						     :: el,gp
	
	
	
	! Deformation gradient for each element	
	edp = 0d0		
	edp = transpose(ed)		
	F=0d0
	call	DEFORMATION4I_AXI(dg,J_out,coord,edof,edp,ngp,nelm,ndof)
		! Calculate stiffness and mass matrix
		!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
		!$OMP SHARED(nelm,ngp,dg,coord,enod,ed,dofnod) &
		!$OMP PRIVATE(ie,gp,ef,esh,es,fe,temp,temp2) &
		!$OMP REDUCTION(+:F)
		do ie=1,nelm  
	
			! use temp to avoid  using non-contiguous call of coord
			temp  = coord(enod(:,ie),:)
			temp2 = ed(:,ie)
			  ! == 
			do gp=1,ngp
				ef(:,gp) = (/ dg(ie,gp,1,1), dg(ie,gp,1,2), dg(ie,gp,2,1), dg(ie,gp,2,2), dg(ie,gp,3,3) /)
			end do	
			! == 	
			! Calculate the 2nd Piola Kirchhoff stress
			call neohooke('2ndPiola',esh,ef)        
	
			! reshape stress tensor
			do gp=1,ngp
				es(gp,1,:) = (/ esh(1,gp), esh(4,gp), 0d0  /)  
				es(gp,2,:) = (/ esh(4,gp), esh(2,gp), 0d0  /) 
				es(gp,3,:) = (/ 0d0     , 0d0     , esh(3,gp) /)									 
			end do
			
			! Calculate the element stiffness matrix
			call PLANI4F_AXI(fe,temp,temp2,es,ngp)
	
			! Assembling
			!!!$OMP CRITICAL
			call insert(F,fe,enod(:,ie),dofnod)	
		end do ! End calculation of force vector
		!$OMP END PARALLEL DO

end subroutine assemGlobalParF

end module assembleG