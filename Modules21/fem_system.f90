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
     module procedure eq_solve_sparse_bcnx_multRHS
     module procedure eq_solve_sparse2_bcnx
     module procedure eq_solve_sparse_sym_bcnx_mult
     module procedure eq_solve_sparse_sym_mult
     module procedure eq_solve_sparse_sym_bcnx_single
     module procedure eq_solve_sparse_sym_single
     module procedure eq_solve_sparse2
     module procedure eq_solve_sparse
     module procedure eq_solve_sparse_multRHS
  end interface
  private eq_solve, eq_solve_rs, eq_solve_bc, eq_solve_bcn, eq_solve_sparse_bc
  private eq_solve_sparse_bcnx_multRHS, eq_solve_sparse_bcnx,eq_solve_sparse2_bcnx, eq_solve_sparse2,eq_solve_sparse, eq_solve_sparse_multRHS
  private eq_solve_sparse_sym_bcnx_mult, eq_solve_sparse_sym_bcnx_single, eq_solve_sparse_sym_mult, eq_solve_sparse_sym_single
  private eq_solve_sparse_bcn  ! never used
  
  interface eigenvalue
     module procedure eq_eigenvalue_sparse
     module procedure eq_geneigenvalue_sparse
     module procedure eq_eigenvalue_sparse_bcnx
     module procedure eq_geneigenvalue_sparse_bcnx
  end interface
  private eq_eigenvalue_sparse, eq_eigenvalue_sparse_bcnx, eq_geneigenvalue_sparse, eq_geneigenvalue_sparse_bcnx
  
  !interface sensCal
  !	module procedure dFintdrho
  !	module procedure adjoint
  !	module procedure Sjk
  !end interface
  !private dFintdrho, adjoint, Sjk
  
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
      edof(tmp1:tmp1+dofnod-1)=(/tmp2:tmp2+dofnod-1/) !!!!!!!!!!!!!!!!! Should it be + ????????????
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
    double precision                 ::  da(size(F,1)), df(size(F,1))
  
  
  
    nrhs=1
    n=size(F,1)
  
    iparm= 0
    iparm(1) = 1 !
    iparm(2) = 3 !
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
    iparm(19) = 0
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
  
  
    if (minval(df).lt.0D0) then
      NegativeE=1
    else
      NegativeE=-1
    endif
  
    !write(*,*) 'Nbr of negaivte eigenvalues', NegativeE
  
   ! NegativeE=-1
  
  
    f=x
    phase = -1 ! release internal memory
    call pardiso(pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
                 idumn, nrhs, iparm, msglvl, ddum, ddum, error)
  
  return
  end subroutine eq_solve_sparse
  
  
  
  
  
  subroutine eq_solve_sparse_multRHS(K,F,NegativeE,nrhs)
    implicit none
    type(sparse)                     :: K
    double precision                 :: f(:,:)
    integer                          :: iparm(64), pt(64), NegativeE
    integer                          :: maxfct, mnum, mtype, phase, error, msglvl, idum
    integer                          :: nrhs, n, ki
    double precision                 :: x(size(F,1),size(F,2)), ddum
    integer                          :: idumn(size(F,1))
    double precision                 ::  da(size(F,1)), df(size(F,1))
  
  
  
    !nrhs=size(f,2)
    n=size(F,1)
  
    iparm= 0
    iparm(1) = 1 !
    iparm(2) = 3 !
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
    iparm(19) = 0
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
  
    if (minval(df).lt.0D0) then
      NegativeE=1
    else
      NegativeE=-1
    endif
  
    f=x
  
    phase = -1 ! release internal memory
    call pardiso(pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
                 idumn, nrhs, iparm, msglvl, ddum, ddum, error)
  
  return
  end subroutine eq_solve_sparse_multRHS
  
  
  
  
  
  ! Same as above only with symmetric sparse pattern (Mathias Wallin)
  ! Multiple RHS
  subroutine eq_solve_sparse_SYM_mult(K,F,NegativeE,dum,dum2)
    implicit none
    type(sparse)                     :: K
    double precision                 :: f(:,:)
    integer                          :: iparm(64), pt(64), NegativeE
    integer                          :: maxfct, mnum, mtype, phase, error, msglvl, idum
    integer                          :: nrhs, n, ki, dum, dum2
    double precision                 :: x(size(F,1),size(F,2)), ddum
    integer                          :: idumn(size(F,1))
    double precision                 ::  da(size(F,1)), df(size(F,1))
  
    nrhs=size(F,2)
    n=size(F,1)
  
    iparm= 0
    iparm(1) = 1 !
    iparm(2) = 3 !
    !iparm(3) = 1 !
    iparm(4) = 0 !
    iparm(5) = 0 !
    iparm(6) = 0 !
    iparm(7) = 0 !
    iparm(8) = 9 !
    iparm(9) = 0 !
    iparm(10) = 13!8 ! 13 before
    iparm(11) = 1!0 ! 1 before
    iparm(12) = 0
    iparm(13) = 0 ! 0 before
    iparm(14) = 0
    iparm(15) = 0
    iparm(16) = 0
    iparm(17) = 0
    iparm(18) = -1
    iparm(19) = 0
    iparm(20) = 0
    iparm(56) = 1
    pt = 0
    error=0
    msglvl=0
    maxfct=1
    mnum=1
    ddum=0d0
    mtype=11!-2
  
  
    phase=13
    call pardiso(pt, maxfct, mnum, mtype, phase, n, K%a, K%ia, K%ja, &
                 idumn, nrhs, iparm, msglvl, f,x, error)
  
    call pardiso_getdiag(pt, df, da, mnum, error)
  
  
    !do ki=1,n
    !  if (df(ki).eq.0D0) write(*,*) 'Zero eigenvalue!'!NegativeE=NegativeE+1
    !enddo
    !if (minval(abs(df)).lt.1d0) write(*,*) 'minval(abs(df))', minval(abs(df))
  
    if (minval(df).lt.0D0) then
      write(*,*) 'Negative eigenvalue!'
      NegativeE=1
    !elseif (minval(abs(df)).lt.1D-7) then
    !  write(*,*) 'Zero eigenvalue!'
    !  NegativeE=1
    endif
  
    !write(*,*) 'Nbr of Neg eigenvalues, df'   , NegativeE
  
    if (error.ne.0) then
     write(*,*) 'Problem in pardiso_getdiag, stopping'
     stop
    endif
  
    f=x
    phase = -1 ! release internal memory
    call pardiso(pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
                 idumn, nrhs, iparm, msglvl, ddum, ddum, error)
  
  return
  end subroutine eq_solve_sparse_SYM_mult
  
  
  ! Same as above only with symmetric sparse pattern (Mathias Wallin)
  ! Single RHS
  subroutine eq_solve_sparse_SYM_single(K,F,NegativeE,dum)
    implicit none
    type(sparse)                     :: K
    double precision                 :: f(:)
    integer                          :: iparm(64), pt(64), NegativeE
    integer                          :: maxfct, mnum, mtype, phase, error, msglvl, idum
    integer                          :: nrhs, n, ki, dum
    double precision                 :: x(size(F,1)), ddum
    integer                          :: idumn(size(F,1))
    double precision                 ::  da(size(F,1)), df(size(F,1))
  
    nrhs=1
    n=size(F,1)
  
    iparm= 0
    iparm(1) = 1 !
    iparm(2) = 3 !
    !iparm(3) = 1 !
    iparm(4) = 0 !
    iparm(5) = 0 !
    iparm(6) = 0 !
    iparm(7) = 0 !
    iparm(8) = 9 !
    iparm(9) = 0 !
    iparm(10) = 13!8 ! 13 before
    iparm(11) = 1!0 ! 1 before
    iparm(12) = 0
    iparm(13) = 0 ! 0 before
    iparm(14) = 0
    iparm(15) = 0
    iparm(16) = 0
    iparm(17) = 0
    iparm(18) = -1
    iparm(19) = 0
    iparm(20) = 0
    iparm(56) = 1
    pt = 0
    error=0
    msglvl=0
    maxfct=1
    mnum=1
    ddum=0d0
    mtype=11!-2
  
  
  
    phase=13
    call pardiso(pt, maxfct, mnum, mtype, phase, n, K%a, K%ia, K%ja, &
                 idumn, nrhs, iparm, msglvl, f,x, error)
  
    if (iparm(30).gt.0) then
        write(*,*) 'Nbr of negative or zero pivots, '   , iparm(30)
    endif
  
    call pardiso_getdiag(pt, df, da, mnum, error)
  
    NegativeE=-1
    !do ki=1,n
    !  if (df(ki).eq.0D0) write(*,*) 'Zero eigenvalue!'!NegativeE=NegativeE+1
    !enddo
    !if (minval(abs(df)).lt.1d0) write(*,*) 'minval(abs(df))', minval(abs(df))
  
    if (minval(df).lt.0D0) then
      write(*,*) 'Negative eigenvalue!'
      NegativeE=1
    !elseif (minval(abs(df)).lt.1D-7) then
    !	 write(*,*) 'Zero eigenvalue!'
    !  NegativeE=1
    endif
  
    !write(*,*) 'Nbr of Neg eigenvalues, df'   , NegativeE
  
    if (error.ne.0) then
     write(*,*) 'Problem in pardiso_getdiag, stopping'
     stop
    endif
  
    f=x
    phase = -1 ! release internal memory
    call pardiso(pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
                 idumn, nrhs, iparm, msglvl, ddum, ddum, error)
  
  return
  end subroutine eq_solve_sparse_SYM_single
  
  
  
  
  
  
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
    iparm(2) = 3 !
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
  
  
  
  subroutine eq_solve_sparse_bcnx_multRHS(K,F,BCnod,BCval,dofnod, NegativeE,nrhs)
    implicit none
    type(sparse)                    :: K
    double precision                :: F(:,:), BCval(:)
    integer                         :: bcnod(:,:), BCdof(size(bcnod,1))
    integer                         :: dofnod
  
    integer                         :: j, ndof, NegativeE, nrhs
  
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
  
    ndof=size(f,1)
  
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
  
    call eq_solve_sparse_multRHS(K,F,NegativeE,nrhs)
  
    if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 9'
  
    return
  end subroutine eq_solve_sparse_bcnx_multRHS
  
  
  
  
  ! Solves using symmetric sparse pattern provided in KSYM and cellSym (Mathias Wallin)
  ! Multiple RHS
  subroutine eq_solve_sparse_sym_bcnx_mult(K,F,BCnod,BCval,dofnod, NegativeE,KSYM,cellSym)
    implicit none
    type(sparse)                    :: K, Ksym
  
    integer                         :: dofnod, cellSym(:,:)
    double precision                :: F(:,:), BCval(:)
    integer                         :: idumn(size(F,1))
    integer                         :: bcnod(:,:), BCdof(size(bcnod,1))
  
    integer                         :: j, ndof, NegativeE, sym
  
    double precision, allocatable   :: dp(:), valsym(:)
    integer, allocatable            :: INX(:)
    integer                         :: i, status, p, node, nent
  
    allocate(valsym(size(Ksym%a)),stat=status)
  
   ! write(*,*) 'inne i modelllen, size valsym', size(valsym,1)
  
  ! Note that this routine does not allow that
  ! bc contains redundant constraints
  
    if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn'
  
  ! find bcdof
    do i=1,size(bcnod,1)
      node=bcnod(i,1)
      bcdof(i)=dofnod*(node-1)+bcnod(i,2)
    enddo
  
    ndof=size(f,1)
  
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
  
    sym=1
  
    !do i=1,size(KSYM%a)
    !   valSym(i)=spagetval(K,cellSym(i,1),cellSym(i,2))
    !enddo
  
    !call spaputval(KSYM,cellSym,valSym)
  
    !call eq_solve_sparse_SYM_mult(KSYM,F,NegativeE,sym)
    call eq_solve_sparse_SYM_mult(K,F,NegativeE,sym,sym)
  
    if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 9'
  
    return
  end subroutine eq_solve_sparse_sym_bcnx_mult
  
  
  
  
  
  ! Solves using symmetric sparse pattern provided in KSYM and cellSym (Mathias Wallin)
  ! Single RHS
  subroutine eq_solve_sparse_sym_bcnx_single(K,F,BCnod,BCval,dofnod, NegativeE,KSYM,cellSym)
    implicit none
    type(sparse)                    :: K, Ksym
  
    integer                         :: dofnod, cellSym(:,:)
    double precision                :: F(:), BCval(:)
    integer                         :: idumn(size(F,1))
    integer                         :: bcnod(:,:), BCdof(size(bcnod,1))
  
    integer                         :: j, ndof, NegativeE, sym
  
    double precision, allocatable   :: dp(:), valsym(:)
    integer, allocatable            :: INX(:)
    integer                         :: i, status, p, node, nent
  
    allocate(valsym(size(Ksym%a)),stat=status)
  
   ! write(*,*) 'inne i modelllen, size valsym', size(valsym,1)
  
  ! Note that this routine does not allow that
  ! bc contains redundant constraints
  
    if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn'
  
  ! find bcdof
    do i=1,size(bcnod,1)
      node=bcnod(i,1)
      bcdof(i)=dofnod*(node-1)+bcnod(i,2)
    enddo
  
    ndof=size(f,1)
  
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
  
    sym=1
  
    !do i=1,size(KSYM%a)
    !   valSym(i)=spagetval(K,cellSym(i,1),cellSym(i,2))
    !enddo
  
    !call spaputval(KSYM,cellSym,valSym)
  
    !call eq_solve_sparse_SYM_single(KSYM,F,NegativeE,sym)
    call eq_solve_sparse_SYM_single(K,F,NegativeE,sym)
  
    if (debugg.eq.1) write(*,*)'eq_solve_sparse_bcn 9'
  
    return
  end subroutine eq_solve_sparse_sym_bcnx_single
  
  
  
  
  
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
    iparm(2) = 3 !
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
    iparm(2) = 3 !
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
  
  
  !subroutine eq_geneigenvalue_sparse_bcnx(K,Mass,ndof,BCnod,BCval,dofnod,Emin,Emax,L,sfactor,phi,M,info)
  !subroutine eq_geneigenvalue_sparse_bcnx(K,Mass,ndof,BCnod,BCval,dofnod,Emin,Emax,L,E,X,M,info)
  subroutine eq_geneigenvalue_sparse_bcnx(K,Mass,ndof,BCnod,BCval,dofnod,Emin,Emax,L,Lnew,sfactor,nbrEig,X,Xold,iter,info)
    implicit none
    type(sparse)                    :: K, Mass
    double precision                :: BCval(:), Emin, Emax, Xold(:,:), sfactor(:)
    integer                         :: bcnod(:,:), BCdof(size(bcnod,1))
    integer                         :: dofnod, j, ndof, NegativeE, nbrEig
    double precision, allocatable   :: dp(:)
    integer, allocatable            :: INX(:)
    integer                         :: i, status, p, node, nent
    double precision, allocatable   :: E(:), X(:,:), res(:)
    integer                         :: M0, M, L, info, Lnew, iter
  
    M0=L
    M=M0
  
    allocate(X(ndof,M0),stat=status)
    X = 0d0
  
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
  
    ! K%a - vector of assigned values
    ! K%ja - element i of K%ja is the number of the column of the i:th element of K%a
    ! K%ia - element i of K%ia is the first non-zero element's index on each row in K%a
    ! - last index in K%ia is the number of entries in K%a + 1
  
    ! For each dimension, i.e. row/col in K
    nent=size(K%ia)-1
    do i=1,nent
      ! If no bc is prescribed at this dof
      if (INX(i).eq.0) then
        ! Do for all elements on this row
        do p=K%ia(i),K%ia(i+1)-1
          ! Which column
          j=K%ja(p)
          ! If bc is prescribed at this column dof
          if (INX(j).eq.1) then
            ! Zero this position in K and M
            K%a(p)=0d0
            Mass%a(p)=0d0
          end if
        end do
      ! If bc is described here
      else
      ! Do for all elements on this row
        do p=K%ia(i),K%ia(i+1)-1
          ! Which column
          j=K%ja(p)
          ! If we are on the diagonal, set mass entry to zero and stiffness entry to one such that phi(dof) = 0
          if (i.eq.j) then
            Mass%a(p) = 0d0
            K%a(p) =1d0
          ! Else zero
          else
            K%a(p)=0d0
            Mass%a(p)=0d0
          end if
        end do
      end if
    end do
  
    deallocate(dp)
    deallocate(INX)
  
    !call eq_geneigenvalue_sparse(K,Mass,ndof,M,X,Emin,Emax,L,E,info)
    ! Now using initial guess of basis for eigenvectors
    call eq_geneigenvalue_sparse(K,Mass,ndof,M,X,Xold,Emin,Emax,L,Lnew,sfactor,nbrEig,iter,info)
  
    return
  end subroutine eq_geneigenvalue_sparse_bcnx
  
  
  
  !subroutine eq_geneigenvalue_sparse(K,Mass,ndof,M,X,Emin,Emax,L,sfactor,phi,info)
  !subroutine eq_geneigenvalue_sparse(K,Mass,ndof,M,X,Emin,Emax,L,sfactor,info)
  !subroutine eq_geneigenvalue_sparse(K,Mass,ndof,M,X,Emin,Emax,L,E,info)
  subroutine eq_geneigenvalue_sparse(K,Mass,ndof,M,X,Xold,Emin,Emax,L,Lnew,sfactor,nbrEig,iter,info)
    implicit none
    type(sparse)                     :: K, Mass
    integer                          :: fpm(128), info, nbrEig, ii, iter, Lnew
    character*1 UPLO
    parameter   (UPLO='F')
    double precision                 :: epsout, Emin, Emax, sfactor(:), Xold(:,:), X(:,:)
    integer                          :: loop, L
    integer                          :: M0, M, ndof, status
    double precision, allocatable    :: E(:), res(:)
  
    M0=L
    M=M0
  
    allocate(E(ndof),stat=status)
    allocate(RES(ndof),stat=status)
  
    call feastinit (fpm)
    fpm(1)=1  ! Print run-time output
    fpm(2)=16 !Number of contours 24!8!16
    fpm(3)=9  !Accuracy (error-trace stop criteria)
    fpm(4)=10 !Number of refinement steps
  
    ! If we aren't on the first iteration and if the size of the subspace hasn't changed
    ! Use old eigenspace as input
   ! if ((iter.ne.1).and.(size(Xold,2).eq.L)) then
   !    write(*,*) 'Using the previous span of eigenvectors as basis for initial guess'
   !    fpm(5)=1 ! CHANGED: now inputs old eigenvectors as initial guess of basis
   !    X = Xold
   ! else
   !    fpm(5)=0
   ! endif
  
    call dfeast_scsrgv(UPLO,ndof,K%a,K%ia,K%ja,Mass%a,Mass%ia,Mass%ja,fpm,epsout,loop,Emin,Emax,M0,E,X,M,res,info)
  
    Lnew = M
  
    if (Lnew.gt.nbrEig) then
      sfactor = E(1:nbrEig)
      !write(*,*) 'Eigenvalues', E(1:nbrEig)
    else
      sfactor = E(1:Lnew)
      !write(*,*) 'Eigenvalues', E(1:Lnew)
    endif
  
  
    print  *,' FEAST OUTPUT INFO ',info
    if(info.ne.0) then
      write(*,*) 'M0:', M0
    endif
  
  
  return
  end subroutine eq_geneigenvalue_sparse
  
  
  
  
  
  
  
  end module fem_system