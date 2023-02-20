module fem_system


! last modified
! M. Ristinmaa 2011-04-12
! M. Ristinmaa 2011-04-27
!   - added node based subroutines extract, solveq, assem
!   - cleaning out old subrutines
! M. Ristinmaa 2011-04-28
!   - new call insert implemented	
! M. Ristinmaa 2011-05-02
!   - change in extract ed now contains columns with element values	
! M. Ristinmaa 2011-05-16
!   - private routines
!   - bugg in extract_mn caused by change of rows and columns in enod
! M. Ristinmaa 2011-10-27
!   - rewrote sparse equations solver _sparse_bcn to _sparse_bcnx 
!     large speed increase
! M. Ristinmaa 2012-01-18
!   - bugg in solveq when edof was used together with bcdof
! ------------------------------------------------------------------------------
! A. Dalklint 2018-05-07
!   - Added routine which sends out the three smallest eigenvectors aswell
! A. Dalklint 2018-09-12
!	 - Added routine provided by M. Wallin which solves using symmetric sparse pattern
! A. Dalklint 2018-10-23
!	 - Added routines which performs some sensitivity calculations


use sparse_util
! Added!!!
use mater_hyperel
use elem_large_cont_3d
use elem_flow_3d


implicit none
 
interface insert
   module procedure insert1
   module procedure insert2
end interface
private insert1, insert2

interface extract
   module procedure extract_s
   module procedure extract_m
   module procedure extract_sn
   module procedure extract_mn
end interface
private extract_s, extract_m, extract_sn, extract_mn

interface assem
   module procedure assem_full
   module procedure assem_full_f
   module procedure assem_fulln
   module procedure assem_fulln_f
   module procedure assem_sparse
   module procedure assem_sparse_f
   module procedure assem_sparsen
   module procedure assem_sparsen_f
end interface
private assem_full, assem_full_f, assem_fulln, assem_fulln_f, assem_sparse
private assem_sparse_f, assem_sparsen, assem_sparsen_f

interface solveq
   module procedure eq_solve
   module procedure eq_solve_rs
   module procedure eq_solve_bc
   module procedure eq_solve_bcn
   module procedure eq_solve_sparse_bc
   module procedure eq_solve_sparse_bcnx
   module procedure eq_solve_sparse2_bcnx
   module procedure eq_solve_sparse_sym_bcnx
   module procedure eq_solve_sparse_sym
   module procedure eq_solve_sparse2
   module procedure eq_solve_sparse
end interface
private eq_solve, eq_solve_rs, eq_solve_bc, eq_solve_bcn, eq_solve_sparse_bc
private eq_solve_sparse_bcnx,eq_solve_sparse2_bcnx, eq_solve_sparse2, eq_solve_sparse
private eq_solve_sparse_sym_bcnx, eq_solve_sparse_sym
private eq_solve_sparse_bcn  ! never used

interface eigenvalue
   module procedure eq_eigenvalue_sparse
   module procedure eq_geneigenvalue_sparse
   module procedure eq_eigenvalue_sparse_bcnx
   module procedure eq_geneigenvalue_sparse_bcnx
   module procedure eq_geneigenvalue3_sparse
   module procedure eq_geneigenvalue3_sparse_bcnx
end interface
private eq_eigenvalue_sparse, eq_eigenvalue_sparse_bcnx, eq_geneigenvalue_sparse, eq_geneigenvalue_sparse_bcnx, eq_geneigenvalue3_sparse, eq_geneigenvalue3_sparse_bcnx

interface adjointCal
	module procedure adjoint
	module procedure fsk
	module procedure fskinter
end interface
private adjoint, fsk, fskinter

integer, parameter :: debugg=0
private debugg

!-------------------------------------------------------------------------------

contains

subroutine insert1(f,fe,enod,dofnod)
  implicit none
  double precision                :: f(:), fe(:)
  integer                         :: enod(:), dofnod
  integer                         :: edof(size(enod)*dofnod)
  integer                         :: noel
  integer                         :: ii, tmp1, tmp2

  noel=size(enod)
  do ii=1,noel
    tmp1=dofnod*(ii-1)+1
    tmp2=dofnod*(enod(ii)-1)+1
    edof(tmp1:tmp1+dofnod-1)=(/tmp2:(tmp2+dofnod-1)/)
  enddo

  f(edof)=f(edof)+fe

  return
end subroutine insert1

subroutine insert2(f,fe,edof)
  implicit none
  double precision                :: f(:), fe(:)
  integer                         :: edof(:)

  f(edof)=f(edof)+fe

  return
end subroutine insert2

subroutine extract_s(ed,a,edof)
  implicit none
  integer                         :: edof(:)
  double precision                :: ed(:), a(:)
 
  ed=a(edof)

return
end subroutine extract_s


subroutine extract_m(ed,a,edof)
  implicit none
  integer                         :: edof(:,:)
  double precision                :: ed(:,:), a(:)
  integer                         :: nrows, ie

  nrows=size(edof,2)
  do ie=1,nrows
    ed(:,ie)=a(edof(:,ie))
  end do

return
end subroutine extract_m


subroutine extract_sn(ed,a,enod,dofnod)
  implicit none
  integer                         :: enod(:), dofnod
  double precision                :: ed(:), a(:)
  integer                         :: edof(size(enod)*dofnod)
  integer                         :: node, tmp1, tmp2, i
 
  do i=1,size(enod)
    node=enod(i)
    tmp1=dofnod*(i-1)+1
    tmp2=dofnod*(node-1)+1
    edof(tmp1:tmp1+dofnod-1)=(/tmp2:tmp2-dofnod-1/)
  enddo
  ed=a(edof)

return
end subroutine extract_sn


subroutine extract_mn(ed,a,enod,dofnod)
  implicit none
  integer                         :: enod(:,:), dofnod
  double precision                :: ed(:,:), a(:)
  integer                         :: nrows, ncols, ie
  integer                         :: edof(size(enod,1)*dofnod)
  integer                         :: node, tmp1, tmp2, i

  nrows=size(enod,2)
  ncols=size(enod,1)
  do ie=1,nrows
    do i=1,ncols
      node=enod(i,ie)
      tmp1=dofnod*(i-1)+1
      tmp2=dofnod*(node-1)+1
      edof(tmp1:tmp1+dofnod-1)=(/tmp2:tmp2+dofnod-1/)
    enddo
    ed(:,ie)=a(edof)
  end do

return
end subroutine extract_mn


subroutine assem_full(K,Ke,edof)
  implicit none
  double precision                :: K(:,:), Ke(:,:)
  integer                         :: edof(:)
  integer                         :: ndof, icol, irow

  ndof=size(Ke,1)

  do icol=1,ndof
    do irow=1,ndof
      K(edof(irow),edof(icol))= K(edof(irow),edof(icol))+Ke(irow,icol)
    end do
  end do

return
end subroutine assem_full

subroutine assem_full_f(K,Ke,f,fe,edof)
  implicit none
  double precision                :: K(:,:), Ke(:,:), f(:), fe(:)
  integer                         :: edof(:)
  integer                         :: ndof, icol, irow

  ndof=size(Ke,1)

  do icol=1,ndof
    do irow=1,ndof
      K(edof(irow),edof(icol))= K(edof(irow),edof(icol))+Ke(irow,icol)
    end do
  end do
  f(edof)=f(edof)+fe

return
end subroutine assem_full_f


subroutine assem_fulln(K,Ke,enod,dofnod)
  implicit none
  double precision                :: K(:,:), Ke(:,:)
  integer                         :: enod(:), dofnod
  integer                         :: edof(size(Ke,1))
  integer                         :: ndof, icol, irow
  integer                         :: i, tmp1, tmp2

  ndof=size(Ke,1)

  do i=1,size(enod)
    tmp1=dofnod*(i-1)+1
    tmp2=dofnod*(enod(i)-1)+1
    edof(tmp1:(tmp1+dofnod-1))=(/tmp2:tmp2+dofnod-1/)
  enddo
  

  do icol=1,ndof
    do irow=1,ndof
      K(edof(irow),edof(icol))= K(edof(irow),edof(icol))+Ke(irow,icol)
    end do
  end do

return
end subroutine assem_fulln

subroutine assem_fulln_f(K,Ke,f,fe,enod,dofnod)
  implicit none
  double precision                :: K(:,:), Ke(:,:), f(:), fe(:)
  integer                         :: enod(:), dofnod
  integer                         :: edof(size(Ke,1))
  integer                         :: ndof, icol, irow
  integer                         :: i, tmp1, tmp2

  ndof=size(Ke,1)

  do i=1,size(enod)
    tmp1=dofnod*(i-1)+1
    tmp2=dofnod*(enod(i)-1)+1
    edof(tmp1:(tmp1+dofnod-1))=(/tmp2:tmp2+dofnod-1/)
  enddo
  

  do icol=1,ndof
    do irow=1,ndof
      K(edof(irow),edof(icol))= K(edof(irow),edof(icol))+Ke(irow,icol)
    end do
  end do

  f(edof)=f(edof)+fe

return
end subroutine assem_fulln_f


subroutine assem_sparse(K,Ke,edof)
  implicit none
  type(sparse)                    :: K
  integer                         :: edof(:)
  double precision                :: Ke(:,:)

  integer                         :: n, dofe, rad, i, m

  dofe=size(edof)

  do n=1,dofe
    rad=K%ia(edof(n))
    do m=1,dofe
       do  i=K%ia(edof(n)),K%ia(edof(n)+1)-1
         if (K%ja(i)==edof(m)) then
            K%a(i)=K%a(i)+Ke(n,m)
         end if
       end do
    end do
  end do
 
return
end subroutine assem_sparse

subroutine assem_sparse_f(K,Ke,f,fe,edof)
  implicit none
  type(sparse)                    :: K
  integer                         :: edof(:)
  double precision                :: Ke(:,:), f(:), fe(:)

  integer                         :: n, dofe, rad, i, m

  dofe=size(edof)

  do n=1,dofe
    rad=K%ia(edof(n))
    do m=1,dofe
       do  i=K%ia(edof(n)),K%ia(edof(n)+1)-1
         if (K%ja(i)==edof(m)) then
            K%a(i)=K%a(i)+Ke(n,m)
         end if
       end do
    end do
  end do
 
  f(edof)=f(edof)+fe

return
end subroutine assem_sparse_f


subroutine assem_sparsen(K,Ke,enod,dofnod)
  implicit none
  type(sparse)                    :: K
  integer                         :: enod(:), dofnod
  double precision                :: Ke(:,:)
  integer                         :: edof(size(ke,1))

  integer                         :: n, dofe, rad, i, m
  integer                         :: tmp1, tmp2

  do i=1,size(enod)
    tmp1=dofnod*(i-1)+1
    tmp2=dofnod*(enod(i)-1)+1
    edof(tmp1:(tmp1+dofnod-1))=(/tmp2:tmp2+dofnod-1/)
  enddo
  
  dofe=size(edof)

  do n=1,dofe
    rad=K%ia(edof(n))
    do m=1,dofe
       do  i=K%ia(edof(n)),K%ia(edof(n)+1)-1
         if (K%ja(i)==edof(m)) then
            K%a(i)=K%a(i)+Ke(n,m)
         end if
       end do
    end do
  end do
 
return
end subroutine assem_sparsen

subroutine assem_sparsen_f(K,Ke,f,fe,enod,dofnod)
  implicit none
  type(sparse)                    :: K
  integer                         :: enod(:), dofnod
  double precision                :: Ke(:,:), f(:), fe(:)
  integer                         :: edof(size(ke,1))

  integer                         :: n, dofe, rad, i, m
  integer                         :: tmp1, tmp2

  do i=1,size(enod)
    tmp1=dofnod*(i-1)+1
    tmp2=dofnod*(enod(i)-1)+1
    edof(tmp1:(tmp1+dofnod-1))=(/tmp2:tmp2+dofnod-1/)
  enddo
  
  dofe=size(edof)

  do n=1,dofe
    rad=K%ia(edof(n))
    do m=1,dofe
       do  i=K%ia(edof(n)),K%ia(edof(n)+1)-1
         if (K%ja(i)==edof(m)) then
            K%a(i)=K%a(i)+Ke(n,m)
         end if
       end do
    end do
  end do
 
  f(edof)=f(edof)+fe

return
end subroutine assem_sparsen_f


subroutine EQ_SOLVE(K,F)
  implicit none
  double precision                :: K(:,:), F(:)
  integer                         :: INFO

  integer                         :: NRHS=1, ndof, status
  integer, allocatable            :: IPIV(:)

  ndof=size(K,1)
  allocate(IPIV(ndof),STAT=status)
       
  CALL DGESV( ndof, NRHS,K, ndof, IPIV, F, ndof, INFO )

  deallocate(IPIV)

RETURN
END subroutine eq_solve


subroutine EQ_SOLVE_RS(K,F)
  implicit none
  double precision                :: K(:,:), F(:,:)
  integer                         :: INFO

  integer                         :: NRHS, ndof, status
  integer, allocatable            :: IPIV(:)

  ndof=size(K,1)
  NRHS=size(F,2)
  allocate(IPIV(ndof),STAT=status)
       
  CALL DGESV( ndof, NRHS,K, ndof, IPIV, F, ndof, INFO )

  deallocate(IPIV)

RETURN
END subroutine eq_solve_RS


subroutine EQ_SOLVE_BC(K,F,BCdof,BCval)
  implicit none
  double precision                :: K(:,:), F(:), BCval(:)
  integer                         :: BCdof(:)
  integer                         :: INFO

  double precision, allocatable   :: REDUCED_FORCE(:), Ktmp(:,:)
  integer, allocatable            :: IPIV(:), fdof(:), TMP0(:)
  integer                         :: status
  integer                         :: ndof, nbc, na
  integer                         :: NRHS=1, i, j

  ndof=size(K,1)
  nbc=size(BCdof)
  na=ndof-nbc

  allocate(TMP0(ndof),STAT=status)
  allocate(REDUCED_FORCE(na),STAT=status)
  allocate(IPIV(na),STAT=status)
  allocate(fdof(na),STAT=status)
  allocate(Ktmp(na,na),STAT=status)
   
  DO i=1,ndof
    TMP0(i)=i
  END DO 
  DO i=1,nbc
    TMP0(BCdof(i))=0
  END DO
  j=0
  DO i=1,ndof
    if (TMP0(i)==0) THEN
      j=j+1
    ELSE
      fdof(i-j)=i
    END IF 
  END DO
 
  DO i=1,na
    REDUCED_FORCE(i)=F(fdof(i))
    DO j=1,nbc
      REDUCED_FORCE(i)=REDUCED_FORCE(i)-K(fdof(i),BCdof(j))*BCval(j)
    END DO
  END DO

  Ktmp=K(fdof,fdof)
  CALL DGESV( na, NRHS,Ktmp, na, IPIV, REDUCED_FORCE, na, INFO )
  F(BCdof)=BCval
  F(fdof)=REDUCED_FORCE

  deallocate(TMP0)
  deallocate(REDUCED_FORCE)
  deallocate(IPIV)
  deallocate(fdof)
  deallocate(Ktmp)

return
end subroutine eq_solve_bc

subroutine EQ_SOLVE_BCn(K,F,BCnod,BCval,dofnod)
  implicit none
  double precision                :: K(:,:), F(:), BCval(:)
  integer                         :: bcnod(:,:), BCdof(size(bcnod,1)), dofnod
  integer                         :: INFO

  double precision, allocatable   :: REDUCED_FORCE(:), Ktmp(:,:)
  integer, allocatable            :: IPIV(:), fdof(:), TMP0(:)
  integer                         :: status
  integer                         :: ndof, nbc, na, node
  integer                         :: NRHS=1, i, j

! find bcdof
  do i=1,size(bcnod,1)
    node=bcnod(i,1)
    bcdof(i)=dofnod*(node-1)+bcnod(i,2)
  enddo


  ndof=size(K,1)
  nbc=size(BCdof)
  na=ndof-nbc

  allocate(TMP0(ndof),STAT=status)
  allocate(REDUCED_FORCE(na),STAT=status)
  allocate(IPIV(na),STAT=status)
  allocate(fdof(na),STAT=status)
  allocate(Ktmp(na,na),STAT=status)
   
  DO i=1,ndof
    TMP0(i)=i
  END DO 
  DO i=1,nbc
    TMP0(BCdof(i))=0
  END DO
  j=0
  DO i=1,ndof
    if (TMP0(i)==0) THEN
      j=j+1
    ELSE
      fdof(i-j)=i
    END IF 
  END DO
 
  DO i=1,na
    REDUCED_FORCE(i)=F(fdof(i))
    DO j=1,nbc
      REDUCED_FORCE(i)=REDUCED_FORCE(i)-K(fdof(i),BCdof(j))*BCval(j)
    END DO
  END DO

  Ktmp=K(fdof,fdof)
  CALL DGESV( na, NRHS,Ktmp, na, IPIV, REDUCED_FORCE, na, INFO )
  F(BCdof)=BCval
  F(fdof)=REDUCED_FORCE

  deallocate(TMP0)
  deallocate(REDUCED_FORCE)
  deallocate(IPIV)
  deallocate(fdof)
  deallocate(Ktmp)

return
end subroutine eq_solve_bcn


subroutine eq_solve_sparse(K,F,NegativeE)
  implicit none
  type(sparse)                     :: K
  double precision                 :: f(:)
  integer                          :: iparm(64), pt(64), NegativeE
  integer                          :: maxfct, mnum, mtype, phase, error, msglvl, idum
  integer                          :: nrhs, n, ki
  double precision                 :: x(size(F,1)), ddum
  integer                          :: idumn(size(F,1))
  double precision, allocatable    :: Kfull(:,:)
  double precision                 ::  da(size(F,1)), df(size(F,1))



  nrhs=1
  n=size(F,1)

  iparm= 0
  iparm(1) = 1 !
  iparm(2) = 2 !
  iparm(3) = 1 !
  iparm(4) = 0 !
  iparm(5) = 0 !
  iparm(6) = 0 !
  iparm(7) = 0 !
  iparm(8) = 9 !
  iparm(9) = 0 !
  iparm(10) = 13
  iparm(11) = 1 
  iparm(12) = 0 
  iparm(13) = 0 
  iparm(14) = 0 
  iparm(15) = 0 
  iparm(16) = 0 
  iparm(17) = 0 
  iparm(18) = -1
  iparm(19) = -1
  iparm(20) = 0 
  iparm(56) = 1 

  pt = 0
  error=0
  msglvl=0
  mtype=11
  maxfct=1
  mnum=1
  ddum=0d0


  phase=13
  call pardiso(pt, maxfct, mnum, mtype, phase, n, K%a, K%ia, K%ja, &
               idumn, nrhs, iparm, msglvl, f,x, error)



  call pardiso_getdiag(pt, df, da, mnum, error)


  NegativeE=0
  do ki=1,n
    if (df(ki).lt.0D0) NegativeE=NegativeE+1
  enddo
!  write(*,*) 'Nbr of negaivte eigenvalues', NegativeE

  if (minval(df).lt.0D0) then
    NegativeE=1   
  else
    NegativeE=-1
  endif


 ! NegativeE=-1


  f=x
  phase = -1 ! release internal memory
  call pardiso(pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
               idumn, nrhs, iparm, msglvl, ddum, ddum, error)

return
end subroutine eq_solve_sparse


! Same as above only with symmetric sparse pattern
subroutine eq_solve_sparse_SYM(K,F,NegativeE,dum,bcdof,nbcNew)
  implicit none
  type(sparse)                     :: K
  double precision                 :: f(:)
  integer                          :: iparm(64), pt(64), NegativeE, ierr
  integer                          :: maxfct, mnum, mtype, phase, error, msglvl, idum, ii, jj, bcdof(:), nbrIslands
  integer                          :: nrhs, n, ki, dum, indx, nbrbc, nbcNew
  double precision                 :: x(size(F,1)), ddum
  integer                          :: idumn(size(F,1)), perm(size(F,1))
  double precision, allocatable    :: Kfull(:,:)
  double precision                 ::  da(size(F,1)), df(size(F,1))
  logical								  :: refactorize
  !integer, allocatable				  :: islandbc(:), tmpdof(:), tmp(:)
  nrhs=1
  n=size(F,1)
  refactorize = .False.
  !nbrbc = size(bcdof,1)

  !allocate(islandbc(n), stat=ierr)

  iparm= 0
  iparm(1) = 1 !
  iparm(2) = 2 !
  iparm(3) = 1 ! Shouldn't it be 0?
  iparm(4) = 0 !
  iparm(5) = 2 ! ! 0 before, 2 to obtain the permutation matric used in phase 1 in idumn
  ! IF set = 1 one can provide the permutation array as below

  !do jj = 1,n
  !	idumn(jj) = jj
  !enddo

  iparm(6) = 0 !
  iparm(7) = 0 ! OUTPUT: number of iterative refinements steps
  iparm(8) = 0 ! 0 or 9?
  iparm(9) = 0 !
  !iparm(10) = 0 ! 13 for nonsymmetric matrices. 8 for symmetric.
  iparm(11) = 0 ! 0 for symmetric matrices. 1 for nonsymmetric. Weighting of elements.
  iparm(12) = 0 
  iparm(13) = 0 
  iparm(14) = 0 ! OUTPUT: Number of perturbated pivots. 
  iparm(15) = 0 
  iparm(16) = 0 
  iparm(17) = 0 
  iparm(18) = -1 ! Does report the number of non-zero elements in factors.
  iparm(19) = -1
  iparm(20) = 0 
  iparm(21) = 1 ! Maybe 3 or 0,2? Or 1. If 3: iterative refinements are not used during pivoting.
  iparm(36) = 0 ! Controls phase 331, 332 and 333
  iparm(56) = 1 ! Enables use of getdiag and pivot.

  pt = 0
  error=0
  msglvl=0
  mtype=-2
  maxfct=1
  mnum=1
  ddum=0d0
  x = 0d0

  ! First perform the analysis
  !phase=11
  !call pardiso(pt, maxfct, mnum, mtype, phase, n, K%a, K%ia, K%ja, &
  !             idumn, nrhs, iparm, msglvl, f, x, error)


 ! write(*,*) 'Analysis step completed.'

  ! First perform the analysis and factorization step
  phase=12
  call pardiso(pt, maxfct, mnum, mtype, phase, n, K%a, K%ia, K%ja, &
               idumn, nrhs, iparm, msglvl, f, x, error)
  
 ! write(*,*) 'Numerical factorization step completed.'

  !write(*,*) 'P', idumn
 ! write(*,*) '--------------------------------------------------------------'

  call pardiso_getdiag(pt, df, da, mnum, error)

 ! write(*,*) '--------------------------------------------------------------'
  !write(*,*) 'df', df
  !write(*,*) '--------------------------------------------------------------'
 ! write(*,*) 'da', da
 ! write(*,*)  '--------------------------------------------------------------'

 ! write(*,*) 'Replacing zeros on diagonal with ones.'
 ! write(*,*) 'K', K%a
	
  ! Tänk på att om vi inte har user-permutation så kanske det inte blir nod 3 som "bildar en nolla" först i gauss eliminering. Frågan är om man ska ta bort kraften från alla noder på den lösa delen?

  !islandbc = 0
  nbrIslands = 0
  do ii=1,n
		indx = idumn(ii)
		if (df(ii).lt.1d-7) then
			if (.not.any(bcdof==indx)) then
			   write(*,*) '!!!!!! Structural island exists !!!!!!'
				refactorize = .true.
				!islandbc(indx) = 1
				nbrIslands = nbrIslands + 1
				do jj=1,n
					if (indx.ne.jj) then
						if (spaGetVal(K, indx, jj).ne.0d0) then
							call spaPutVal(K, indx, jj, 0d0)
						endif
						if (spaGetVal(K, jj, indx).ne.0d0) then
							call spaPutVal(K, jj, indx, 0d0)
						endif
					endif
				enddo
				call spaPutVal(K, indx, indx, 1d0)
				!F(indx) = 1d-8	
				F(indx) = 0d0
				bcdof(nbcNew+nbrIslands) = indx
			endif	
		endif
  enddo
  nbcNew = nbcNew + nbrIslands

  ! Restart factorization
  if (refactorize) then
  		!write(*,*) 'K', K%a
  		! Not sure if it uses the same permutation vector?? Save it
  		perm = idumn

  		phase = -1 ! release internal memory
  		call pardiso(pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
               	idumn, nrhs, iparm, msglvl, ddum, ddum, error)

  		!write(*,*) '...............Restarting.......................'

  		! Choose the previous permutation vector
  		!iparm(5) = 1
  		!idumn = perm

  		!allocate(tmpdof(nbrbc+nbrIslands), stat=ierr)
  		!allocate(tmp(nbrbc+nbrIslands), stat=ierr)
		!jj = 0
		!do ii = 1,n
		!	if (islandbc(ii).eq.1) then
		!		jj = jj + 1
		!		tmp(nbrbc+jj) = ii
		!	endif
		!enddo

		!do ii = 1,nbrbc+nbrIslands
		!	if (ii.lt.nbrbc+nbrIslands) then
		!		tmpdof(ii) = bcdof(ii)
		!	else
		!		tmpdof(ii) = tmp(ii)
		!	endif
		!enddo
  		phase = 13
  		call pardiso(pt, maxfct, mnum, mtype, phase, n, K%a, K%ia, K%ja, &
               	idumn, nrhs, iparm, msglvl, f, x, error)

  else
	   phase = 33
  		call pardiso(pt, maxfct, mnum, mtype, phase, n, K%a, K%ia, K%ja, &
               	idumn, nrhs, iparm, msglvl, f, x, error)

  endif

  call pardiso_getdiag(pt, df, da, mnum, error)

  !write(*,*) '--------------------------------------------------------------'
  !write(*,*) 'df', df
  !write(*,*) '--------------------------------------------------------------'
  !write(*,*) 'da', da
  !write(*,*)  '--------------------------------------------------------------'
  !write(*,*) 'P', idumn
  !write(*,*) '--------------------------------------------------------------'

  !write(*,*) 'Number of iterative refinements', iparm(7)
  !write(*,*) '--------------------------------------------------------------'

  !write(*,*) 'Number of zero or negative pivots', iparm(30)
  !write(*,*) '--------------------------------------------------------------'

  !write(*,*) 'Number of perturbated pivots', iparm(14)
  !write(*,*) '--------------------------------------------------------------'

  NegativeE=0
  do ki=1,n
    if (df(ki).lt.0D0) NegativeE=NegativeE+1
  enddo
 ! write(*,*) 'Nbr of negaivte eigenvalues', NegativeE

  if (minval(df).lt.0D0) then
    NegativeE=1   
  else
    NegativeE=-1
  endif

!  write(*,*) 'Solving for real and structurally symmetric K', minval(da), minval(df)
 ! NegativeE=-1


  if (error.ne.0) then
   write(*,*) 'Problem in pardiso_getdiag, stopping'
   stop
  endif


  f=x

  phase = -1 ! release internal memory
  call pardiso(pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
               idumn, nrhs, iparm, msglvl, ddum, ddum, error)

  !if (allocated(tmpdof)) then
	!	deallocate(tmpdof)
	!	deallocate(tmp)
  !endif

  !if (allocated(islandbc)) then
	!	deallocate(islandbc)
  !endif

return
end subroutine eq_solve_sparse_SYM

subroutine eq_solve_sparse_bc_old(K,F,BCdof,BCval)
  implicit none
  type(sparse)                    :: K
  double precision                :: F(:), BCval(:)
  integer                         :: BCdof(:)

  integer                         :: iparm(64), pt(64)
  integer                         :: maxfct, mnum, mtype, phase, error, msglvl, idum
  integer                         :: nrhs, n, nbc, na, j, j2, ndof, ne
  double precision                :: ddum

  double precision, allocatable   :: REDUCED_FORCE(:), x(:), dp(:), res(:)
  integer, allocatable            :: fdof(:), idumn(:), ia(:), ja(:), INX(:)
  integer                         :: i, status, p, ni, ival(1), nj

  nrhs=1
  n=size(F,1)
  nbc=size(bcval)
  ndof=size(f)
  ne=size(K%a)
  n=ndof-nbc

  allocate(dp(ndof),stat=status)
  allocate(res(ndof),stat=status)
  allocate(INX(ndof),stat=status)

  allocate(REDUCED_FORCE(ndof-nbc),stat=status)
  allocate(x(ndof-nbc),stat=status)
  allocate(fdof(ndof-nbc),stat=status)
  allocate(idumn(ndof-nbc),stat=status)

  fdof=0
  x=0d0
  REDUCED_FORCE=0d0

  dp=0D0
  INX=0
  INX(bcdof)=1

  na=0
  DO i=1,ndof
    if (INX(i)==1) then
      na=na+1
      INX(i)=-1
    ELSE
      INX(i)=na
      fdof(i-na)=i
    END IF    
  END DO

  dp(bcdof)=bcval
  f=f-matmul(K,dp)
  REDUCED_FORCE=F(fdof)

! reduce the sparse matrix
! do not change K%ia and K%ja

  nj=size(K%ja)
  allocate(ja(nj),stat=status)
  ja=K%ja

  do i=1,nbc
    do j=1,ne
      if (bcdof(i).eq.ja(j)) then
        ja(j)=-1
      endif
    enddo
    do j=K%ia(bcdof(i)),K%ia(bcdof(i)+1)-1
       ja(j)=-1
    enddo
  enddo

  j=0
  do i=1,ne
    if (ja(i).eq.-1) then
      j=j+1
    else
      K%a(i-j)=K%a(i)
    endif
  enddo
  j2=i-j
  
  ni=size(K%ia)
  allocate(ia(ni),stat=status)
  ia=0
  ia(1)=1
  i=1
  do p=1,ni-1
    do j=K%ia(p),K%ia(p+1)-1
      if (ja(j).ne.-1) then
        i=i+1
      else
        j2=j2-1
      endif   
      ia(p+1)=i
    enddo
  enddo


  do i=1,nbc
    ival=maxloc(bcdof)
    do j=1,ne
      if (ja(j).gt.bcdof(ival(1)))  ja(j)=ja(j)-1
    enddo
    do j=bcdof(ival(1)),ni-1
      ia(j)=ia(j+1)
    enddo
    bcdof(ival)=-bcdof(ival)
  enddo
  bcdof=-bcdof

  j=0
  do i=1,ne
    if (ja(i).ne.-1) then
      j=j+1
      ja(j)=ja(i)
    endif
  enddo
  p=p-nbc

  iparm= 0
  iparm(1) = 1 !
  iparm(2) = 2 !
  iparm(3) = 1 !
  iparm(4) = 0 !
  iparm(5) = 0 !
  iparm(6) = 0 !
  iparm(7) = 0 !
  iparm(8) = 9 !
  iparm(9) = 0 !
  iparm(10) = 13
  iparm(11) = 1 
  iparm(12) = 0 
  iparm(13) = 0 
  iparm(14) = 0 
  iparm(15) = 0 
  iparm(16) = 0 
  iparm(17) = 0 
  iparm(18) = -1
  iparm(19) = -1
  iparm(20) = 0 
  pt = 0
  error=0
  msglvl=0
  mtype=11
  maxfct=1
  mnum=1
  ddum=0d0

!msglvl = 1
!iparm(27) = 1

  phase=13
  call pardiso(pt, maxfct, mnum, mtype, phase, n, K%a(1:j), ia(1:p), ja(1:j), &
               idumn, nrhs, iparm, msglvl, REDUCED_FORCE, x, error)

!write(*,*) error

  F=dp
  F(fdof)=x

  phase = -1 ! release internal memory
  call pardiso(pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
               idumn, nrhs, iparm, msglvl, ddum, ddum, error)


  deallocate(dp)
  deallocate(res)
  deallocate(INX)

  deallocate(ia)
  deallocate(ja)
  deallocate(REDUCED_FORCE)
  deallocate(x)
  deallocate(fdof)
  deallocate(idumn)

  return
end subroutine eq_solve_sparse_bc_old

subroutine eq_solve_sparse_bc(K,F,BCdof,BCval)
  implicit none
  type(sparse)                    :: K
  double precision                :: F(:), BCval(:)
  integer                         :: bcdof(:)
  integer                         :: dofnod, NegativeE

  integer                         :: j, ndof

  double precision, allocatable   :: dp(:)
  integer, allocatable            :: INX(:)
  integer                         :: i, status, p, node, nent

! Note that this routine does not allow that 
! bc contains redundant constraints

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bc'
  ndof=size(f)

! would be nice to get rid of the allocation part
  allocate(INX(ndof),stat=status); if (status.ne.0) write(*,*)' allp INX'
  allocate(dp(ndof),stat=status); if (status.ne.0) write(*,*)' allp dp'
  INX=0
  INX(bcdof)=1
  dp=0d0
  dp(bcdof)=bcval

  nent=size(K%ia)-1
  do i=1,nent
    if (INX(i).eq.0) then
      do p=K%ia(i),K%ia(i+1)-1
        j=K%ja(p)
        if (INX(j).eq.1) then
          f(i)=f(i)-K%a(p)*dp(j)
          K%a(p)=0d0
        end if
      end do
    else
      do p=K%ia(i),K%ia(i+1)-1
        j=K%ja(p)
        if (i.eq.j) then
          f(i)=K%a(p)*dp(j)
        else
          K%a(p)=0d0
        end if
      end do
    end if
  end do

  deallocate(dp)
  deallocate(INX)

  call eq_solve_sparse(K,F,NegativeE)

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 9'

  return
end subroutine eq_solve_sparse_bc



subroutine eq_solve_sparse_bcnx(K,F,BCnod,BCval,dofnod, NegativeE)
  implicit none
  type(sparse)                    :: K
  double precision                :: F(:), BCval(:)
  integer                         :: bcnod(:,:), BCdof(size(bcnod,1))
  integer                         :: dofnod

  integer                         :: j, ndof, NegativeE

  double precision, allocatable   :: dp(:)
  integer, allocatable            :: INX(:)
  integer                         :: i, status, p, node, nent

! Note that this routine does not allow that 
! bc contains redundant constraints

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn'

! find bcdof
  do i=1,size(bcnod,1)
    node=bcnod(i,1)
    bcdof(i)=dofnod*(node-1)+bcnod(i,2)
  enddo

  ndof=size(f)

! would be nice to get rid of the allocation part
  allocate(INX(ndof),stat=status); if (status.ne.0) write(*,*)' allp INX'
  allocate(dp(ndof),stat=status); if (status.ne.0) write(*,*)' allp dp'
  INX=0
  INX(bcdof)=1
  dp=0d0
  dp(bcdof)=bcval

  nent=size(K%ia)-1
  do i=1,nent
    if (INX(i).eq.0) then
      do p=K%ia(i),K%ia(i+1)-1
        j=K%ja(p)
        if (INX(j).eq.1) then
          f(i)=f(i)-K%a(p)*dp(j)
          K%a(p)=0d0
        end if
      end do
    else
      do p=K%ia(i),K%ia(i+1)-1
        j=K%ja(p)
        if (i.eq.j) then
          f(i)=K%a(p)*dp(j)
        else
          K%a(p)=0d0
        end if
      end do
    end if
  end do

  deallocate(dp)
  deallocate(INX)

  call eq_solve_sparse(K,F,NegativeE)

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 9'

  return
end subroutine eq_solve_sparse_bcnx



! Solves using symmetric sparse pattern provided in KSYM and cellSym (Mathias Wallin)
subroutine eq_solve_sparse_sym_bcnx(K,F,bcdof,BCval,dofnod,NegativeE,KSYM,cellSym,nbcNew)
  implicit none
  type(sparse)                    :: K, Ksym

  integer                         :: dofnod, cellSym(:,:)
  double precision                :: F(:), BCval(:)
  integer                         :: idumn(size(F,1)), invcellSym(size(cellSym,1),2)
  integer                         :: bcdof(:), nbcNew
  
  integer                         :: j, ndof, NegativeE, sym

  double precision, allocatable   :: dp(:), valsym(:), smallBCval(:)
  integer, allocatable            :: INX(:), smallBCdof(:)
  integer                         :: i, status, p, node, nent 

  allocate(valsym(size(Ksym%a)),stat=status)
 ! write(*,*) 'inne i modelllen, size valsym', size(valsym,1)

! Note that this routine does not allow that 
! bc contains redundant constraints

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn'

  allocate(smallBCdof(nbcNew), stat=status)
  allocate(smallBCval(nbcNew), stat=status)
   
  smallBCval = 0d0

! find bcdof
  do i=1,nbcNew
    smallBCdof(i) = bcdof(i)
    smallBCval(i) = bcval(i)
  enddo

  ndof=size(f)

! would be nice to get rid of the allocation part
  allocate(INX(ndof),stat=status); if (status.ne.0) write(*,*)' allp INX'
  allocate(dp(ndof),stat=status); if (status.ne.0) write(*,*)' allp dp'
  INX=0
  INX(smallBCdof)=1
  dp=0d0
  dp(smallBCdof)=smallBCval

  nent=size(K%ia)-1
  do i=1,nent
    if (INX(i).eq.0) then
      do p=K%ia(i),K%ia(i+1)-1
        j=K%ja(p)
        if (INX(j).eq.1) then
          f(i)=f(i)-K%a(p)*dp(j)
          K%a(p)=0d0
        end if
      end do
    else
      do p=K%ia(i),K%ia(i+1)-1
        j=K%ja(p)
        if (i.eq.j) then
          f(i)=K%a(p)*dp(j)
        else
          K%a(p)=0d0
        end if
      end do
    end if
  end do

  deallocate(dp)
  deallocate(INX)

  deallocate(smallBCdof)
  deallocate(smallBCval)

  !write(*,*) 'Inne i module'
  sym=1
  !write(*,*) 'size cell', size(Cellsym(:,1))

  !write(*,*) 'size KSYM', size(KSYM%a)
 ! write(*,*) 'size K', size(K%a)

  !write(*,*) 'size valsym', size(valSym)

  do i=1,size(KSYM%a)
     valSym(i)=spagetval(K,cellSym(i,1),cellSym(i,2))
  enddo       

 !write(*,*) 'Inne i module 999'
  call spaputval(KSYM,cellSym,valSym)
  !write(*,*) 'size KSYM', KSYM%a
	!write(*,*) 'size K', K%a
 ! write(*,*) 'Inne i module 2' 
 
  call eq_solve_sparse_SYM(KSYM,F,NegativeE,sym,BCdof,nbcNew)

 ! Obtain full sparse pattern again as ouput in K
  do i=1,size(KSYM%a)
     valSym(i)=spagetval(KSYM,cellSym(i,1),cellSym(i,2))
  enddo  

  call spaputval(K,cellSym,valSym)

  invcellSym(:,1) = cellSym(:,2)
  invcellSym(:,2) = cellSym(:,1)

  call spaputval(K,invcellSym,valSym)


!  call eq_solve_sparse(K,F,NegativeE)

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 9'


  return
end subroutine eq_solve_sparse_sym_bcnx

subroutine eq_solve_sparse_bcn(K,F,BCnod,BCval,dofnod)
  implicit none
  type(sparse)                    :: K
  double precision                :: F(:), BCval(:)
  integer                         :: bcnod(:,:), BCdof(size(bcnod,1))
  integer                         :: dofnod

  integer                         :: iparm(64), pt(64)
  integer                         :: maxfct, mnum, mtype, phase, error, msglvl, idum
  integer                         :: nrhs, n, nbc, na, j, j2, ndof, ne
  double precision                :: ddum

  double precision, allocatable   :: REDUCED_FORCE(:), x(:), dp(:), res(:)
  integer, allocatable            :: fdof(:), idumn(:), ia(:), ja(:), INX(:)
  integer                         :: i, status, p, ni, ival(1), nj, node

! Note that this routine does not allow that 
! bc contains redundant constraints

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn'


 

! find bcdof
  do i=1,size(bcnod,1)
    node=bcnod(i,1)
    bcdof(i)=dofnod*(node-1)+bcnod(i,2)
  enddo

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 1'

  nrhs=1
  n=size(F,1)
  nbc=size(bcval)
  ndof=size(f)
  ne=size(K%a)
  n=ndof-nbc

  allocate(dp(ndof),stat=status); if (status.ne.0) write(*,*)' allp dp'
  allocate(res(ndof),stat=status); if (status.ne.0) write(*,*)' allp res'
  allocate(INX(ndof),stat=status); if (status.ne.0) write(*,*)' allp INX'

  allocate(REDUCED_FORCE(ndof-nbc),stat=status); if (status.ne.0) write(*,*)' allp res'
  allocate(x(ndof-nbc),stat=status); if (status.ne.0) write(*,*)' allp x'
  allocate(fdof(ndof-nbc),stat=status); if (status.ne.0) write(*,*)' allp fdof'
  allocate(idumn(ndof-nbc),stat=status); if (status.ne.0) write(*,*)' allp idumn'

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 2'

  fdof=0
  x=0d0
  REDUCED_FORCE=0d0

  dp=0D0
  INX=0
  INX(bcdof)=1

  na=0
  DO i=1,ndof
    if (INX(i)==1) then
      na=na+1
      INX(i)=-1
    ELSE
      INX(i)=na
      fdof(i-na)=i
    END IF    
  END DO

  dp(bcdof)=bcval
  f=f-matmul(K,dp)
  REDUCED_FORCE=F(fdof)

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 3'

! reduce the sparse matrix
! do not change K%ia and K%ja

  nj=size(K%ja)
  allocate(ja(nj),stat=status)
  ja=K%ja

  do i=1,nbc
    do j=1,ne
      if (bcdof(i).eq.ja(j)) then
        ja(j)=-1
      endif
    enddo
    do j=K%ia(bcdof(i)),K%ia(bcdof(i)+1)-1
       ja(j)=-1
    enddo
  enddo

  j=0
  do i=1,ne
    if (ja(i).eq.-1) then
      j=j+1
    else
      K%a(i-j)=K%a(i)
    endif
  enddo
  j2=i-j

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 4'
  
  ni=size(K%ia)
  allocate(ia(ni),stat=status)
  ia=0
  ia(1)=1
  i=1
  do p=1,ni-1
    do j=K%ia(p),K%ia(p+1)-1
      if (ja(j).ne.-1) then
        i=i+1
      else
        j2=j2-1
      endif   
      ia(p+1)=i
    enddo
  enddo

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 5'

  do i=1,nbc
    ival=maxloc(bcdof)
    do j=1,ne
      if (ja(j).gt.bcdof(ival(1)))  ja(j)=ja(j)-1
    enddo
    do j=bcdof(ival(1)),ni-1
      ia(j)=ia(j+1)
    enddo
    bcdof(ival)=-bcdof(ival)
  enddo
  bcdof=-bcdof

  j=0
  do i=1,ne
    if (ja(i).ne.-1) then
      j=j+1
      ja(j)=ja(i)
    endif
  enddo
  p=p-nbc

  iparm= 0
  iparm(1) = 1 !
  iparm(2) = 2 !
  iparm(3) = 1 !
  iparm(4) = 0 !
  iparm(5) = 0 !
  iparm(6) = 0 !
  iparm(7) = 0 !
  iparm(8) = 9 !
  iparm(9) = 0 !
  iparm(10) = 13
  iparm(11) = 1 
  iparm(12) = 0 
  iparm(13) = 0 
  iparm(14) = 0 
  iparm(15) = 0 
  iparm(16) = 0 
  iparm(17) = 0 
  iparm(18) = -1
  iparm(19) = -1
  iparm(20) = 0 
  pt = 0
  error=0
  msglvl=0
  mtype=11 ! mtype=11 <-> real and non-symmetric
  maxfct=1
  mnum=1
  ddum=0d0

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 6'

  phase=13
  call pardiso(pt, maxfct, mnum, mtype, phase, n, K%a(1:j), ia(1:p), ja(1:j), &
               idumn, nrhs, iparm, msglvl, REDUCED_FORCE, x, error)


!  write(*,*) 'Number of neg. eigenvatues', iparm(23)
!  write(*,*) 'Number of pos. eigenvatues', iparm(24)
!  if (iparm(23).gt.0) then
!    write(*,*) 'Critical piont passed !'
!    pause
!  endif




  F=dp
  F(fdof)=x

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 7'

  phase = -1 ! release internal memory
  call pardiso(pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
               idumn, nrhs, iparm, msglvl, ddum, ddum, error)

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 8'

  deallocate(dp)
  deallocate(res)
  deallocate(INX)

  deallocate(ia)
  deallocate(ja)
  deallocate(REDUCED_FORCE)
  deallocate(x)
  deallocate(fdof)
  deallocate(idumn)

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 9'

  return
end subroutine eq_solve_sparse_bcn


subroutine eq_solve_sparse2(K,F)
  implicit none
  type(sparse)                    :: K
  double precision                :: f(:,:)
  integer                         :: iparm(64), pt(64)
  integer                         :: maxfct, mnum, mtype, phase, error, msglvl, idum
  integer                         :: nrhs, n
  double precision                :: x(size(F,1),size(F,2)), ddum
  integer                         :: idumn(size(F,1))
  integer                         :: i
  double precision                :: df(size(F,1)), da(size(F,1)) ! MW Diagonal elements


  nrhs=size(F,2)
  n=size(F,1)

  iparm= 0
  iparm(1) = 1 !
  iparm(2) = 2 !
  iparm(3) = 1 !
  iparm(4) = 0 !
  iparm(5) = 0 !
  iparm(6) = 0 !
  iparm(7) = 0 !
  iparm(8) = 9 !
  iparm(9) = 0 !
  iparm(10) = 13
  iparm(11) = 1 
  iparm(12) = 0 
  iparm(13) = 0 
  iparm(14) = 0  
  iparm(15) = 0 
  iparm(16) = 0 
  iparm(17) = 0 
  iparm(18) = -1
  iparm(19) = -1
  iparm(20) = 0 
  pt = 0
  error=0
  msglvl=0
  mtype=11
  maxfct=1
  mnum=1
  ddum=0d0

  phase=13

  call pardiso(pt, maxfct, mnum, mtype, phase, n, K%a, K%ia, K%ja, &
               idumn, nrhs, iparm, msglvl, f,x, error)



  f=x
  phase = -1 ! release internal memory



!  call pardiso_getdiag(pt, df, da, mnum, error) ! MW Diagonal elements
! Does not exist in old Intel Compiler


! write(*,*) 'Number of negative eigenvalues', iparm(23)
! write(*,*) 'Number of positive eigenvalues', iparm(24)
! if (iparm(23).gt.0) then
!   write(*,*) 'Critical point passed !'
!   pause
! endif


  call pardiso(pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
               idumn, nrhs, iparm, msglvl, ddum, ddum, error)

return
end subroutine eq_solve_sparse2




! Allow for multiple RHS
subroutine eq_solve_sparse2_bcnx(K,F,BCnod,BCval,dofnod)
  implicit none
  type(sparse)                    :: K
  double precision                :: F(:,:), BCval(:)
  integer                         :: bcnod(:,:), BCdof(size(bcnod,1))
  integer                         :: dofnod

  integer                         :: j, ndof

  double precision, allocatable   :: dp(:)
  integer, allocatable            :: INX(:)
  integer                         :: i, status, p, node, nent, nrhs

! Note that this routine does not allow that 
! bc contains redundant constraints

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn'

  nrhs=size(F,2)
  ndof=size(f,1)
  
! find bcdof
  do i=1,size(bcnod,1)
    node=bcnod(i,1)
    bcdof(i)=dofnod*(node-1)+bcnod(i,2)
  enddo
 
! would be nice to get rid of the allocation part
  allocate(INX(ndof),stat=status); if (status.ne.0) write(*,*)' allp INX'
  allocate(dp(ndof),stat=status); if (status.ne.0) write(*,*)' allp dp'
  INX=0
  INX(bcdof)=1
  dp=0d0
  dp(bcdof)=bcval

  nent=size(K%ia)-1
  do i=1,nent
    if (INX(i).eq.0) then
      do p=K%ia(i),K%ia(i+1)-1
        j=K%ja(p)
        if (INX(j).eq.1) then
          f(i,:)=f(i,:)-K%a(p)*dp(j)
          K%a(p)=0d0
        end if
      end do
    else
      do p=K%ia(i),K%ia(i+1)-1
        j=K%ja(p)
        if (i.eq.j) then
          f(i,:)=K%a(p)*dp(j)
        else
          K%a(p)=0d0
        end if
      end do
    end if
  end do

  deallocate(dp)
  deallocate(INX)



  call eq_solve_sparse2(K,F)



  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 9'

  return
end subroutine eq_solve_sparse2_bcnx



subroutine eq_eigenvalue_sparse_bcnx(K,ndof,BCnod,BCval,dofnod,Emin,Emax,L,info)
  implicit none
  type(sparse)                    :: K
  double precision                :: BCval(:), Emin, Emax
  integer                         :: bcnod(:,:), BCdof(size(bcnod,1))
  integer                         :: dofnod

  integer                         :: j, ndof, NegativeE

  double precision, allocatable   :: dp(:)
  integer, allocatable            :: INX(:)
  integer                         :: i, status, p, node, nent
  double precision, allocatable    :: E(:), X(:,:), res(:)
  integer                         :: M0,M, L,info



  M0=L
  M=M0

  allocate(E(ndof),stat=status)
  allocate(X(ndof,M0),stat=status)

! Note that this routine does not allow that 
! bc contains redundant constraints

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn'

! find bcdof
  do i=1,size(bcnod,1)
    node=bcnod(i,1)
    bcdof(i)=dofnod*(node-1)+bcnod(i,2)
  enddo

! would be nice to get rid of the allocation part
  allocate(INX(ndof),stat=status); if (status.ne.0) write(*,*)' allp INX'
  allocate(dp(ndof),stat=status); if (status.ne.0) write(*,*)' allp dp'
  INX=0
  INX(bcdof)=1
  dp=0d0
  dp(bcdof)=bcval

  nent=size(K%ia)-1
  do i=1,nent
    if (INX(i).eq.0) then
      do p=K%ia(i),K%ia(i+1)-1
        j=K%ja(p)
        if (INX(j).eq.1) then
          K%a(p)=0d0
        end if
      end do
    else
      do p=K%ia(i),K%ia(i+1)-1
        j=K%ja(p)
        if (i.eq.j) then
        else
          K%a(p)=0d0
        end if
      end do
    end if
  end do

  deallocate(dp)
  deallocate(INX)

  call eq_eigenvalue_sparse(K,ndof,M,X,Emin,Emax,L,info)

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 9'

  return
end subroutine eq_eigenvalue_sparse_bcnx

subroutine eq_geneigenvalue_sparse_bcnx(K,Mass,ndof,BCnod,BCval,dofnod,Emin,Emax,L,sfactor,phi,info)
  implicit none
  type(sparse)                    :: K, Mass
  double precision                :: BCval(:), Emin, Emax, sfactor(3)
  integer                         :: bcnod(:,:), BCdof(size(bcnod,1))
  integer                         :: dofnod

  integer                         :: j, ndof, NegativeE

  double precision, allocatable   :: dp(:)
  integer, allocatable            :: INX(:)
  integer                         :: i, status, p, node, nent
  double precision, allocatable    :: E(:), X(:,:), res(:),phi(:)
  integer                         :: M0, M, L, info

  M0=L
  M=M0

  allocate(E(ndof),stat=status)
  allocate(X(ndof,M0),stat=status)
  allocate(phi(ndof),stat=status)


! Note that this routine does not allow that 
! bc contains redundant constraints

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn'

! find bcdof
  do i=1,size(bcnod,1)
    node=bcnod(i,1)
    bcdof(i)=dofnod*(node-1)+bcnod(i,2)
  enddo

! would be nice to get rid of the allocation part
  allocate(INX(ndof),stat=status); if (status.ne.0) write(*,*)' allp INX'
  allocate(dp(ndof),stat=status); if (status.ne.0) write(*,*)' allp dp'
  INX=0
  INX(bcdof)=1
  dp=0d0
  dp(bcdof)=bcval

  nent=size(K%ia)-1
  do i=1,nent
    if (INX(i).eq.0) then
      do p=K%ia(i),K%ia(i+1)-1
        j=K%ja(p)
        if (INX(j).eq.1) then
          K%a(p)=0d0
          Mass%a(p)=0d0
        end if
      end do
    else
      do p=K%ia(i),K%ia(i+1)-1
        j=K%ja(p)
        if (i.eq.j) then
        else
          K%a(p)=0d0
          Mass%a(p)=0d0
        end if
      end do
    end if
  end do

  deallocate(dp)
  deallocate(INX)

  call eq_geneigenvalue_sparse(K,Mass,ndof,M,X,Emin,Emax,L,sfactor,phi,info)

  return
end subroutine eq_geneigenvalue_sparse_bcnx

subroutine eq_geneigenvalue3_sparse_bcnx(K,Mass,ndof,BCdof,BCval,dofnod,Emin,Emax,L,sfactor,phi1,phi2,phi3,info)
  implicit none
  type(sparse)                    :: K, Mass
  double precision                :: BCval(:), Emin, Emax, sfactor(3)
  integer                         :: bcdof(:)
  integer                         :: dofnod

  integer                         :: j, ndof, NegativeE

  double precision, allocatable   :: dp(:)
  integer, allocatable            :: INX(:)
  integer                         :: i, status, p, node, nent
  double precision, allocatable    :: E(:), X(:,:), res(:),phi1(:),phi2(:),phi3(:)
  integer                         :: M0, M, L, info

  M0=L
  M=M0

  allocate(E(ndof),stat=status)
  allocate(X(ndof,M0),stat=status)
  allocate(phi1(ndof),stat=status)
  allocate(phi2(ndof),stat=status)
  allocate(phi3(ndof),stat=status)

! Note that this routine does not allow that 
! bc contains redundant constraints

  if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn'


! find bcdof
 ! do i=1,nbcNew
 !   smallBCdof(i) = bcdof(i)
 ! enddo

! would be nice to get rid of the allocation part
  allocate(INX(ndof),stat=status); if (status.ne.0) write(*,*)' allp INX'
  allocate(dp(ndof),stat=status); if (status.ne.0) write(*,*)' allp dp'
  INX=0
  INX(BCdof)=1
  dp=0d0
  dp(BCdof)=bcval

  nent=size(K%ia)-1
  do i=1,nent
    if (INX(i).eq.0) then
      do p=K%ia(i),K%ia(i+1)-1
        j=K%ja(p)
        if (INX(j).eq.1) then
          K%a(p)=0d0
          Mass%a(p)=0d0
        end if
      end do
    else
      do p=K%ia(i),K%ia(i+1)-1
        j=K%ja(p)
        if (i.eq.j) then
        else
          K%a(p)=0d0
          Mass%a(p)=0d0
        end if
      end do
    end if
  end do

  deallocate(dp)
  deallocate(INX)

  call eq_geneigenvalue3_sparse(K,Mass,ndof,M,X,Emin,Emax,L,sfactor,phi1,phi2,phi3,info)

  return
end subroutine eq_geneigenvalue3_sparse_bcnx


subroutine eq_eigenvalue_sparse(K,ndof,M,X,Emin,Emax,L,info)
  implicit none
  type(sparse)                     :: K

  integer                          :: fpm(128)
  character*1 UPLO
  parameter   (UPLO='F')
  double precision                 ::  epsout, Emin,Emax
  integer     loop
  integer                          :: L, info

  integer     			   :: M0,M, ndof, status
  double precision, allocatable    :: E(:), X(:,:), res(:)

  M0=L
  M=M0

  allocate(E(ndof),stat=status)
  allocate(X(ndof,M0),stat=status)
  allocate(RES(ndof),stat=status)

  call feastinit (fpm)
  fpm(1)=1
  fpm(2)=24
  fpm(7)=4
 
  call dfeast_scsrev(UPLO,ndof,K%a, K%ia,K%ja,fpm,epsout,loop,Emin,Emax,M0,E,X,M,res,info)


  print  *,' FEAST OUTPUT INFO ',info
  if(info.ne.0) stop 1

return
end subroutine eq_eigenvalue_sparse



subroutine eq_geneigenvalue_sparse(K,Mass,ndof,M,X,Emin,Emax,L,sfactor,phi,info)
  implicit none
  type(sparse)                     :: K, Mass


  integer                          :: fpm(128), info
  character*1 UPLO
  parameter   (UPLO='F')
  double precision                 ::  epsout, Emin,Emax, sfactor(3)
  integer     loop
  integer     L
  integer     			   :: M0,M, ndof, status
  double precision, allocatable    :: E(:), X(:,:), res(:),phi(:)


  M0=L
  M=M0

  allocate(phi(ndof),stat=status)
  allocate(E(ndof),stat=status)
  allocate(X(ndof,M0),stat=status)
  allocate(RES(ndof),stat=status)

  call feastinit (fpm)
  fpm(1)=1
  fpm(2)=24
  fpm(7)=4
 
  call dfeast_scsrgv(UPLO,ndof,K%a, K%ia,K%ja,Mass%a, Mass%ia,Mass%ja,fpm,epsout,loop,Emin,Emax,M0,E,X,M,res,info)
  sfactor=E(1:3)

  phi=X(:,1)

  print  *,' FEAST OUTPUT INFO ',info
  if(info.ne.0) then
    write(*,*) 'M0:', M0
  endif

!  write(*,*) 'Inne och rotat i eigenvalue gen', E(1:3)
!  write(*,*) 'M became', E(1:3)
!  write(*,*) 'sfactor', sfactor

return
end subroutine eq_geneigenvalue_sparse

subroutine eq_geneigenvalue3_sparse(K,Mass,ndof,M,X,Emin,Emax,L,sfactor,phi1,phi2,phi3,info)
  implicit none
  type(sparse)                     :: K, Mass


  integer                          :: fpm(128), info
  character*1 UPLO
  parameter   (UPLO='F')
  double precision                 ::  epsout, Emin,Emax, sfactor(3)
  integer     loop
  integer     L
  integer     			   :: M0,M, ndof, status
  double precision, allocatable    :: E(:), X(:,:), res(:),phi1(:),phi2(:),phi3(:)


  M0=L
  M=M0

  allocate(phi1(ndof),stat=status)
  allocate(phi2(ndof),stat=status)
  allocate(phi3(ndof),stat=status)
  allocate(E(ndof),stat=status)
  allocate(X(ndof,M0),stat=status)
  allocate(RES(ndof),stat=status)

  call feastinit (fpm)
  fpm(1)=1
  fpm(2)=24
  fpm(4)=10 ! CHANGED NUMBER OF REFINEMENT LOOPS
  fpm(7)=4
 
  call dfeast_scsrgv(UPLO,ndof,K%a, K%ia,K%ja,Mass%a, Mass%ia,Mass%ja,fpm,epsout,loop,Emin,Emax,M0,E,X,M,res,info)
  sfactor=E(1:3)

  phi1=X(:,1)
  phi2=X(:,2)
  phi3=X(:,3)

  print  *,' FEAST OUTPUT INFO ',info
  if(info.ne.0) then
    write(*,*) 'M0:', M0
  endif

!  write(*,*) 'Inne och rotat i eigenvalue gen', E(1:3)
!  write(*,*) 'M became', E(1:3)
!  write(*,*) 'sfactor', sfactor

return
end subroutine eq_geneigenvalue3_sparse


subroutine adjoint(hp,edrho,edphis,edphik,ed,coord,enod,removed,delta0,B,eta,p,Hprojection,dofnod,pdeFilter,rhotilde)
  implicit none

  double precision                 :: edphis(:,:), edphik(:,:), edrho(:,:), rhogp(8), hp(:), rhotilde(:)
  double precision					  :: ed(:,:), coord(:,:)
  double precision, allocatable	  :: dg(:,:,:)
  integer								  :: removed(:), enod(:,:)
  double precision                 :: CRAP98(9,8), Ones(8), GradFATP(6,8), Dgp(6,6,8), Bl(6,8), X3P(6,8), AargP(6,8), TargP(6,8), g, hep(24), tmp24(24), Etmps1(6,8), Etmpk1(6,8), Etmp2(6,8), Etmp3(6,8)
  integer                          :: ie, igp, nelm, ngp, Hprojection, ierr, dofnod
  double precision 					  :: delta0, B, eta, p, gfun
  logical								  :: pdeFilter

  ngp = 8
  nelm = size(ed,2)
  Ones=1D0

  allocate(dg(9,8,nelm), stat=ierr)

  dg=0d0
  dg(1,:,:)=1d0
  dg(5,:,:)=1d0
  dg(9,:,:)=1d0

  ! Calculate the adjoint force vector
  do ie=1,nelm	
		if (pdeFilter) then
         ! Extract gauss point values from nodal quantities
         call fl3d8_gp(rhogp,edrho(:,ie)) 
      else
         rhogp = rhotilde(ie)
      endif

		! Calculate the different E-tilde
		! s = v => T, k = phi => A
		! E(a,v) = T
      call c3dtl8_emu(coord(:,enod(:,ie)),ed(:,ie),edphis(:,ie),Etmps1,crap98,Ones) 
		! E(a,phi) = A
      call c3dtl8_emu(coord(:,enod(:,ie)),ed(:,ie),edphik(:,ie),Etmpk1,crap98,Ones) 

		! E(0,phi)
      call c3dtl8_emu(coord(:,enod(:,ie)),0D0*ed(:,ie),edphik(:,ie),Etmp2,crap98,Ones) 
		! E(v,phi)
      call c3dtl8_emu(coord(:,enod(:,ie)),edphis(:,ie),edphik(:,ie),Etmp3,crap98,Ones)  

      ! The deformation gradient 
		call c3dtl8_d(dg(:,:,ie),coord(:,enod(:,ie)),ed(:,ie))

		! Calculate nabla(f_AT)
		call dfat(GradFATP,dg(:,:,ie),Etmpk1,Etmps1) 

		! Calculate the tangent stiffness 
		call dneohooke('tl',Dgp,dg(:,:,ie))

		! Introduce quantity
		Bl=Etmp3-Etmp2  

		! Thresholding               
      do igp=1,ngp
			if (removed(ie).eq.0) then
            g=gfun(rhogp(igp),delta0,p,eta,B,Hprojection)  
				Dgp(:,:,igp) = g*Dgp(:,:,igp)
				
            X3P(:,igp)=(matmul(Dgp(:,:,igp),Bl(:,igp))+1D0*GradFATP(:,igp)*g) ! Note that 2D0 is replaced by 1D0 since the factor 2D0 is inside sub. dfat
            AargP(:,igp)=matmul(Dgp(:,:,igp),Etmpk1(:,igp))
            TargP(:,igp)=matmul(Dgp(:,:,igp),Etmps1(:,igp))
			else
				Dgp(:,:,igp) = 0d0

            X3P(:,igp)=(matmul(Dgp(:,:,igp),Bl(:,igp))+1D0*GradFATP(:,igp)*0d0)
            AargP(:,igp)=matmul(Dgp(:,:,igp),Etmpk1(:,igp))
            TargP(:,igp)=matmul(Dgp(:,:,igp),Etmps1(:,igp))
			endif 
      enddo

		! First term in h
		call c3dtl8_f(hep,coord(:,enod(:,ie)),ed(:,ie),X3P)

		! The second term is split in two and then added to hep
      call c3dtl8_f(tmp24,coord(:,enod(:,ie)),edphis(:,ie),AargP) 
      hep=hep+tmp24
      call c3dtl8_f(tmp24,coord(:,enod(:,ie)),0D0*edphis(:,ie),AargP) 
      hep=hep-tmp24

		! As well as the third term
		call c3dtl8_f(tmp24,coord(:,enod(:,ie)),edphik(:,ie),TargP) 
      hep=hep+tmp24
      call c3dtl8_f(tmp24,coord(:,enod(:,ie)),0D0*edphik(:,ie),TargP) 
      hep=hep-tmp24

		! Instert in global h
      call insert(hp,hep,enod(:,ie),dofnod)  
  enddo       

  deallocate(dg)   

  return
end subroutine adjoint


subroutine fsk(ffilter,omega,edphis,edphik,edrho,ed,edmu,coord,enod,removed,delta0,B,eta,p,Hprojection,massRho,pdeFilter,rhotilde)
  implicit none

  double precision                 :: edphis(:,:), edphik(:,:), edrho(:,:), edmu(:,:), rhogp(8), ffilter(:), rhotilde(:)
  double precision					  :: ed(:,:), coord(:,:), valm, valf
  double precision, allocatable	  :: dg(:,:,:), es(:,:,:)
  integer								  :: removed(:), enod(:,:)
  double precision                 :: CRAP98(9,8), Ones(8), Dgp(6,6,8), AargP(6,8), gp, esscaled(6,8), tmp6(6), fen(8), fme(8), Etmps1(6,8), Etmpk1(6,8), Etmp2(6,8), Etmp3(6,8), Etmp5(6,8), Scalar(8)
  integer                          :: ie, igp, nelm, ngp, Hprojection, ierr
  double precision 					  :: delta0, B, eta, p, gfun, gfunprim, phiKphi, omega, Hprimfun, massRho
  logical								  :: pdeFilter
  
  ngp = 8
  nelm = size(ed,2) 
  Ones=1D0
  scalar = 0d0

  allocate(dg(9,8,nelm), stat=ierr)
  allocate(es(6,ngp,nelm), stat=ierr)
  
  dg=0d0
  dg(1,:,:)=1d0
  dg(5,:,:)=1d0
  dg(9,:,:)=1d0
  es=0d0

  ! Now calculate fsk
  do ie=1,nelm
		if (pdeFilter) then
         ! Extract gauss point values from nodal quantities
         call fl3d8_gp(rhogp,edrho(:,ie)) 
      else
         rhogp = rhotilde(ie)
      endif  

		! Obtain all Etilde
		! E(a,v) = T
      call c3dtl8_emu(coord(:,enod(:,ie)),ed(:,ie),edphis(:,ie),Etmps1,crap98,Ones) 
		! E(a,phi) = A
      call c3dtl8_emu(coord(:,enod(:,ie)),ed(:,ie),edphik(:,ie),Etmpk1,crap98,Ones) 

		! E(0,phi)
      call c3dtl8_emu(coord(:,enod(:,ie)),0D0*ed(:,ie),edphik(:,ie),Etmp2,crap98,Ones) 
		! E(v,phi)
      call c3dtl8_emu(coord(:,enod(:,ie)),edphis(:,ie),edphik(:,ie),Etmp3,crap98,Ones) 
	   ! E(a,mu)
      call c3dtl8_emu(coord(:,enod(:,ie)),ed(:,ie),edmu(:,ie),Etmp5,crap98,Ones) 

		! The deformation gradient 
		call c3dtl8_d(dg(:,:,ie),coord(:,enod(:,ie)),ed(:,ie))

		! Calculate the tangent stiffness 
		call dneohooke('tl',Dgp,dg(:,:,ie)) 

		! Calculate the 2nd Piola Kirchhoff stress
		call neohooke('2ndPiola',es(:,:,ie),dg(:,:,ie))   

      do igp=1,ngp
			 ! Thresholding, now the derivative of xi is used                 
			 if (removed(ie).eq.0) then
             gp=gfunprim(rhogp(igp),delta0,p,eta,B,Hprojection)   
             esscaled(:,igp)= gp*es(:,igp,ie) 
				 Dgp(:,:,igp) = gp*Dgp(:,:,igp)
			 else
             esscaled(:,igp)= 0d0
				 Dgp(:,:,igp) = 0d0
			 endif 

			 ! Reset
			 Scalar(igp)=0D0
			 phiKphi = 0D0

			 ! Calculate the integral Etilde^T*S
			 Scalar(igp) = Scalar(igp) + dot_product(Etmp5(:,igp),esscaled(:,igp))
			 ! First term in K
			 tmp6=matmul(Dgp(:,:,igp),Etmpk1(:,igp))
			 phiKphi=dot_product(Etmps1(:,igp),tmp6)

		    ! Second term
			 tmp6=Etmp3(:,igp)-Etmp2(:,igp) 
			 phiKphi=phiKphi+dot_product(tmp6,esscaled(:,igp))

			 ! Add to vector
			 Scalar(igp) = Scalar(igp) - phiKphi

      enddo
		! Calculate as an body force
		call fl3d8_bf(fen,coord(:,enod(:,ie)),Scalar)

		! Calculate the derivative of the mass matrix 
		if (removed(ie).eq.0) then
			do igp=1,ngp
         	Scalar(igp)= Hprimfun(rhogp(igp),B,eta)    
      	enddo
		else
			do igp=1,ngp
         	Scalar(igp)=0d0  
      	enddo
		endif
		
		call c3dtl8_m(fme,coord(:,enod(:,ie)),massRho,Scalar,edphis(:,ie),edphik(:,ie))
		
		if (pdeFilter) then
         ! Assemble. ffilter is dfdrhotilde
	 		ffilter(enod(:,ie))=ffilter(enod(:,ie))-(-omega*fme)+fen
      else
         valm = dot_product(fme,ones)
         valf = dot_product(fen,ones)
         ffilter(ie)=-(-omega*valm)+valf
      endif 

  enddo

  deallocate(dg)  
  deallocate(es)

  return
end subroutine fsk

subroutine fskInter(ffilter,omega,edphis,edphik,edrho,ed,edmu,coord,enod,removed,delta0,B,eta,p,Hprojection,massRho,c1,c2)
  implicit none

  double precision                 :: edphis(:,:), edphik(:,:), edrho(:,:), edmu(:,:), rhogp(8), ffilter(:)
  double precision					  :: ed(:,:), coord(:,:)
  double precision, allocatable	  :: dg(:,:,:), es(:,:,:)
  integer								  :: removed(:), enod(:,:)
  double precision                 :: CRAP98(9,8), Ones(8), Dgp(6,6,8), AargP(6,8), gp, esscaled(6,8), tmp6(6), fen(8), fme(8), Etmps1(6,8), Etmpk1(6,8), Etmp2(6,8), Etmp3(6,8), Etmp5(6,8), Scalar(8), dScalar(8)
  integer                          :: ie, igp, nelm, ngp, Hprojection, ierr
  double precision 					  :: delta0, B, eta, p, gfun, gfunprim, phiKphi, omega, Hprimfun, c1, c2, Hfun, massRho

  ngp = 8
  nelm = size(ed,2) 
  Ones=1D0
  scalar = 0d0
  dscalar = 0d0

  allocate(dg(9,8,nelm), stat=ierr)
  allocate(es(6,ngp,nelm), stat=ierr)
  
  dg=0d0
  dg(1,:,:)=1d0
  dg(5,:,:)=1d0
  dg(9,:,:)=1d0
  es=0d0

  ! Now calculate fsk
  do ie=1,nelm
		! Extract gauss point values from nodal quantities, rhotilde
      call fl3d8_gp(rhogp,edrho(:,ie))   

		! Obtain all Etilde
		! E(a,v) = T
      call c3dtl8_emu(coord(:,enod(:,ie)),ed(:,ie),edphis(:,ie),Etmps1,crap98,Ones) 
		! E(a,phi) = A
      call c3dtl8_emu(coord(:,enod(:,ie)),ed(:,ie),edphik(:,ie),Etmpk1,crap98,Ones) 

		! E(0,phi)
      call c3dtl8_emu(coord(:,enod(:,ie)),0D0*ed(:,ie),edphik(:,ie),Etmp2,crap98,Ones) 
		! E(v,phi)
      call c3dtl8_emu(coord(:,enod(:,ie)),edphis(:,ie),edphik(:,ie),Etmp3,crap98,Ones) 
	   ! E(a,mu)
      call c3dtl8_emu(coord(:,enod(:,ie)),ed(:,ie),edmu(:,ie),Etmp5,crap98,Ones) 

		! The deformation gradient 
		call c3dtl8_d(dg(:,:,ie),coord(:,enod(:,ie)),ed(:,ie))

		! Calculate the tangent stiffness 
		call dneohooke('tl',Dgp,dg(:,:,ie)) 

		! Calculate the 2nd Piola Kirchhoff stress
		call neohooke('2ndPiola',es(:,:,ie),dg(:,:,ie))   

      do igp=1,ngp
			 ! Thresholding, now the derivative of xi is used                 
			 if (removed(ie).eq.0) then
             gp=gfunprim(rhogp(igp),delta0,p,eta,B,Hprojection)   
             esscaled(:,igp)= gp*es(:,igp,ie) 
				 Dgp(:,:,igp) = gp*Dgp(:,:,igp)
			 else
             esscaled(:,igp)= 0d0
				 Dgp(:,:,igp) = 0d0
			 endif 

			 ! Reset
			 Scalar(igp)=0D0
			 phiKphi = 0D0

			 ! Calculate the integral Etilde^T*S
			 Scalar(igp) = Scalar(igp) + dot_product(Etmp5(:,igp),esscaled(:,igp))
			 ! First term in K
			 tmp6=matmul(Dgp(:,:,igp),Etmpk1(:,igp))
			 phiKphi=dot_product(Etmps1(:,igp),tmp6)

		    ! Second term
			 tmp6=Etmp3(:,igp)-Etmp2(:,igp) 
			 phiKphi=phiKphi+dot_product(tmp6,esscaled(:,igp))

			 ! Add to vector
			 Scalar(igp) = Scalar(igp) - phiKphi

      enddo
		! Calculate as an body force
		call fl3d8_bf(fen,coord(:,enod(:,ie)),Scalar)

		! Calculate the derivative of the mass matrix 
		do igp=1,ngp
         Scalar(igp)= Hfun(rhogp(igp),B,eta)    
         dScalar(igp)= Hprimfun(rhogp(igp),B,eta) 
      enddo
		call c3dtl8_m(fme,coord(:,enod(:,ie)),massRho,c1,c2,Scalar,dScalar,edphis(:,ie),edphik(:,ie))


		! Assemble. ffilter is dfdrhotilde
	 	ffilter(enod(:,ie))=ffilter(enod(:,ie))-(-omega*fme)+fen
  enddo

  deallocate(dg)  
  deallocate(es)

  return
end subroutine fskInter



end module fem_system
 

