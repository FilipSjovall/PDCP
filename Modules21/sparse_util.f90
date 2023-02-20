include 'mkl_spblas.f90'

module sparse_util

 
! last modified
! M. Ristinmaa 2011-04-13
! M. Ristinmaa 2011-04-15
!    - introduced spaaddval
! M. Ristinmaa 2011-04-19
!    - new code for defining sparse matrix in allocate_sparse_m
!      scales with n, large speed improvments for large systems
!    - added sparemove
! M. Ristinmaa 2011-04-20
!    - sparse_add updated, using mkl to add two sparse matrices with different
!      sparse patterns	
! M. Ristinmaa 2011-04-23
!    - non-square sparse matrix allowed, large part rewritten	
!    - sparse cell implemented
! M. Ristinmaa 2011-04-23
!    - cleaned module
! M. Ristinmaa 2011-04-26
!    - cols updated in *
! M. Ristinmaa 2011-04-27
!    - renamed spadefedof to spatopdef
!    - renamed spadefcell to spacelldef
!    - implemented nod based allocation in spatopdef
! M. Ristinmaa 2011-05-03
!    - change in spatopdef due change of enod and edof
! M. Ristinmaa 2011-01-21
!    - bugg in spaadd :compearing y%ja and z%ja different sizes
! M. Ristinmaa 2012-06-11
!    - bugg in spaadd : if sparse matrix exists deallocate first
! -----------------------------------------------------------------------------

use memory_util
use fem_util
use mkl_spblas

implicit none


type sparse
   double precision, allocatable  :: a(:)
   integer, allocatable           :: ia(:), ja(:)
   integer                        :: cols
end type
! number of rows is found from rows=size(a%ia)-1

interface spaputval
   module procedure spaputval_s
   module procedure spaputval_v
end interface
private spaputval_s, spaputval_v


interface spaaddval
   module procedure spaaddval_s
   module procedure spaaddval_v
end interface
private spaaddval_s, spaaddval_v

interface spagetval
   module procedure spagetval_s
   module procedure spagetval_v
end interface
private spagetval_s, spagetval_v


interface matmul
   module procedure matmul_sparse
   !module procedure matmul_sparse2
end interface
private matmul_sparse!, matmul_sparse2

interface matmulTrans
   module procedure matmul_sparseTrans
end interface
private matmul_sparseTrans

interface matmul_sym
   module procedure matmul_sparse_sym
end interface
private matmul_sparse_sym

interface spatopdef
   module procedure spadefedof1
   module procedure spadefedof2
   module procedure spadefedof3
   module procedure spadefedof4
   module procedure spadefedof5
   module procedure spadefedof6
   module procedure spadefedofs
   module procedure spadefedofm
end interface
private spadefedof1, spadefedof2, spadefedof3, spadefedof4
private spadefedof5, spadefedof6, spadefedofs, spadefedofm

interface spatopdef_sym
   module procedure spadefedof1_sym
   module procedure spadefedof2_sym
   module procedure spadefedof3_sym
   module procedure spadefedof4_sym
   module procedure spadefedof5_sym
   module procedure spadefedof6_sym
   module procedure spadefedofs_sym
   module procedure spadefedofm_sym
end interface
private spadefedof1_sym, spadefedof2_sym, spadefedof3_sym, spadefedof4_sym
private spadefedof5_sym, spadefedof6_sym, spadefedofs_sym, spadefedofm_sym

interface spacelldef
   module procedure spadefcell1
   module procedure spadefcell2
   module procedure spadefcell3
   module procedure spadefcellm
end interface
private spadefcell1, spadefcell2, spadefcell3, spadefcellm


interface operator(*)
   module procedure scalar_mult_sparse
   module procedure sparse_mult_scalar
end interface
private scalar_mult_sparse, sparse_mult_scalar


private bubble

integer, parameter  :: debugg=0
private debugg

!------------------------------------------------------------------------------

contains

subroutine sparemove(K)
  implicit none
  type(sparse)                    :: K

  if (allocated(K%a)) then
     deallocate(K%a)
     deallocate(K%ia)
     deallocate(K%ja)
  endif

  return
end subroutine sparemove


subroutine spa2full(Amat, K)
  implicit none
  type(sparse),     intent(in)    :: K
  double precision, intent(out)   :: Amat(:,:)
  integer                         :: i, j, rows

! dimension check
  if (size(Amat,1).ne.(size(K%ia)-1)) & 
          stop "spa2full matrices are of different sizes"
  if (size(Amat,2).ne.K%cols)         &
          stop "spa2full matrices are of different sizes"

  rows=size(K%ia)-1
  Amat=0.0d0
  do i=1,rows
    do j=K%ia(i),K%ia(i+1)-1
      Amat(i,K%ja(j))=K%a(j)
    enddo
  enddo

  return
end subroutine spa2full


subroutine spa2full_sym(Amat, K)
  implicit none
  type(sparse),     intent(in)    :: K
  double precision, intent(out)   :: Amat(:,:)
  double precision, allocatable   :: AmatTrans(:,:)
  integer                         :: i, j, rows, ierr

  allocate(AmatTrans(size(Amat,1),size(Amat,2)), stat=ierr)

! dimension check
  if (size(Amat,1).ne.(size(K%ia)-1)) & 
          stop "spa2full_sym matrices are of different sizes"
  if (size(Amat,2).ne.K%cols)         &
          stop "spa2full_sym matrices are of different sizes"

  rows=size(K%ia)-1
  Amat=0.0d0
  Amattrans=0.0d0
  do i=1,rows
    do j=K%ia(i),K%ia(i+1)-1
      Amat(i,K%ja(j))=K%a(j)
      if (i.ne.K%ja(j)) then
         AmatTrans(i,K%ja(j))=K%a(j)
      endif
    enddo
  enddo

  Amat = Amat + Transpose(AmatTrans)

  return
end subroutine spa2full_sym



subroutine spaputval_s(K, i, j, val) 
  implicit none
  type(sparse),     intent(inout) :: K
  integer,          intent(in)    :: i, j
  double precision, intent(in)    :: val
  integer                         :: ii
    
! dimension check 
  if (i.gt.(size(K%ia)-1)) stop "spaputval index i outside matrix size"
  if (j.gt.K%cols)         stop "spaputval index j outside matrix size"

  do ii= K%ia(i), K%ia(i+1)-1
    if ( K%ja (ii) == j ) exit
  enddo
  if (ii < K%ia(i+1)) then
    K%a(ii) = val
  else
    stop "spaput Location of index not preallocated"
  endif
  
  return
end subroutine spaputval_s


subroutine spaputval_v(K, iv, val) 
  implicit none
  type(sparse),     intent(inout) :: K
  integer,          intent(in)    :: iv(:,:)
  double precision, intent(in)    :: val(:)
  integer                         :: i
    
  do i=1,size(iv,1)
    call spaputval_s(K,iv(i,1),iv(i,2),val(i))
  enddo
  
  return
end subroutine spaputval_v


subroutine spaaddval_s(K, i, j, val) 
  implicit none
  type(sparse),     intent(inout) :: K
  integer,          intent(in)    :: i, j
  double precision, intent(in)    :: val
  integer                         :: ii
    
! dimension check 
  if (i.gt.(size(K%ia)-1)) stop "spaaddval index i outside matrix size"
  if (j.gt.K%cols)         stop "spaaddval index j outside matrix size"

  do ii= K%ia(i), K%ia(i+1)-1
    if ( K%ja (ii) == j ) exit
  enddo
  if (ii < K%ia(i+1)) then
    K%a(ii) = K%a(ii)+val
  else
    stop "Location of index not preallocated"
  endif
  
  return
end subroutine spaaddval_s


subroutine spaaddval_v(K, iv, val) 
  implicit none
  type(sparse),     intent(inout) :: K
  integer,          intent(in)    :: iv(:,:)
  double precision, intent(in)    :: val(:)
  integer                         :: i
    
  do i=1,size(iv,1)
    call spaaddval_s(K,iv(i,1),iv(i,2),val(i))
  enddo
  
  return
end subroutine spaaddval_v


function spagetval_s(K, i, j) 
  implicit none
  type(sparse),     intent(inout) :: K
  integer,          intent(in)    :: i, j
  double precision                :: spagetval_s
  integer                         :: ii

! dimension check 
  if (i.gt.(size(K%ia)-1)) stop "spagetval index i outside matrix size"
  if (j.gt.K%cols)         stop "spagetval index j outside matrix size"
    
  do ii= K%ia(i), K%ia(i+1)-1
    if ( K%ja (ii) == j ) exit
  enddo
  if (ii < K%ia(i+1)) then
    spagetval_s=K%a(ii)
  else
    spagetval_s=0d0
  endif
  
  return
end function spagetval_s


function spagetval_v(K, iv) 
  implicit none
  type(sparse),     intent(inout) :: K
  integer,          intent(in)    :: iv(:,:)
  double precision                :: spagetval_v(size(iv,1))
  integer                         :: i
    
  do i=1,size(iv,1)
    spagetval_v(i)=spagetval_s(K,iv(i,1),iv(i,2))
  enddo  
  return
end function spagetval_v



! add sparse matrices
subroutine spaadd(X, Y, Z) 
  implicit none
  type(sparse), intent(in)        :: Y
  type(sparse), intent(in)        :: Z
  type(sparse)                    :: X
  integer                         :: n1, n2, status, request, cols, rows, tmp
  integer                         :: last_entry
  double precision                :: dtmp

! dimension check
  if (Y%cols.ne.Z%cols)         stop "spaadd matrices are of different cols"
  if (size(Y%ia).ne.size(Z%ia)) stop "spaadd matrices are of different rows"

  cols=Z%cols       
  rows=size(Z%ia)-1

  n1=size(Z%a)
  n2=size(Z%ia)

! note that mkl_dcsradd allows for the operation Y+beta*op(Y)
! op where op could be transpose
 
!  if ((sum(Z%ia-Y%ia).eq.0).and.(sum(Z%ja-Y%ja).eq.0)) then 
! Bugg: Compearing ja does not work size they can have different sizes
! Z and Y have the same sparse structure
!write(*,*)'3a'
!stop
!     allocate(X%a(n1),stat=status)
!     allocate(X%ja(n1),stat=status)
!     allocate(X%ia(n2),stat=status)
!     X%a = Y%a + Z%a
!     X%ia= Y%ia
!     X%ja= Y%ja
!  else
! Z and Y have different sparse structures
!
! Code should be written such that intel library is notre needed
!
    if (allocated(X%ia)) then
      deallocate(X%ia,X%ja,X%a)
    endif
    allocate(X%ia(rows+1),stat=status)
    request=1
    call mkl_dcsradd('n',request,0,rows,cols,y%a,y%ja,y%ia,1.0d0,& 
                     z%a,z%ja,z%ia,dtmp,tmp,x%ia,,status)
    last_entry=size(x%ia)

    allocate(X%a(x%ia(last_entry)-1),stat=status)
    allocate(X%ja(x%ia(last_entry)-1),stat=status)

    request=2
    call mkl_dcsradd('n',request,0,rows,cols,y%a,y%ja,y%ia,1.0d0, &
                     z%a,z%ja,z%ia,x%a,x%ja,x%ia,,status)
    x%cols=y%cols
  return
end subroutine spaadd


! add 2D arrays of type matrix
function sparse_mult_scalar(Z, val) result(X)
  implicit none
  double precision, intent(in)    :: val
  type(sparse), intent(in)        :: Z
  type(sparse)                    :: X
  integer :: n1, n2, status

  n1=size(Z%a)
  n2=size(Z%ia)

  allocate(X%a(n1),stat=status)
  allocate(X%ja(n1),stat=status)
  allocate(X%ia(n2),stat=status)

  X%a = Z%a*val
  X%ia= Z%ia
  X%ja= Z%ja
  X%cols=Z%cols

end function sparse_mult_scalar

! add 2D arrays of type matrix
function scalar_mult_sparse(val, Z) result(X)
  implicit none
  double precision, intent(in)    :: val
  type(sparse), intent(in)        :: Z
  type(sparse)                    :: X
  integer                         :: n1, n2, status

  n1=size(Z%a)
  n2=size(Z%ia)

  allocate(X%a(n1),stat=status)
  allocate(X%ja(n1),stat=status)
  allocate(X%ia(n2),stat=status)

  X%a = Z%a*val
  X%ia= Z%ia
  X%ja= Z%ja
  X%cols=Z%cols

end function scalar_mult_sparse


function matmul_sparse(K,b) 
  implicit none
  double precision                :: b(:)
  double precision,allocatable    :: matmul_sparse(:)
  type(sparse)                    :: K
  double precision                :: Aij
  integer                         :: i, j, p, nent, ierr
  
! use matrix*vector
! call mkl_dcsrgemv(transa, m, a, ia, ja, x, y)	

  if (size(b).ne.K%cols)   &
       stop "matmul sparse matrix and vector are of different sizes"
  
  nent=size(K%ia)-1

  allocate(matmul_sparse(nent),stat=ierr)
  matmul_sparse=0.0d0
  do i=1,nent
    do p=K%ia(i),K%ia(i+1)-1
      j=K%ja(p)
      Aij=K%a(p)
      matmul_sparse(i)=matmul_sparse(i)+Aij*b(j)
    end do
  end do
end function matmul_sparse


function matmul_sparseTrans(K,b) 
  implicit none
  double precision                :: b(:)
  double precision,allocatable    :: matmul_sparseTrans(:)
  type(sparse)                    :: K
  integer                         :: Krows, request, info, nzmax, ierr, tmp
  character(6)                    :: matdescra
  
  Krows=size(K%ia)-1
  matdescra = 'GXXF'
  
  if (Krows.ne.size(b)) stop "matmulTrans sparse matrix dimensions does not match"
  
  if (allocated(matmul_sparseTrans))  deallocate(matmul_sparseTrans)
  allocate(matmul_sparseTrans(K%cols),stat=ierr)
  matmul_sparseTrans = 0d0
  
  ! use matrix*vector
  call mkl_dcsrmv('t', Krows, K%cols, 1d0, matdescra, K%a, K%ja, K%ia(1), K%ia(2), b, 0d0, matmul_sparseTrans)

end function matmul_sparseTrans


function matmul_sparse_sym(K,b) 
  implicit none
  double precision                :: b(:)
  double precision,allocatable    :: matmul_sparse_sym(:)
  type(sparse)                    :: K
  integer                         :: Krows, request, info, nzmax, ierr, tmp
  character(6)                    :: matdescra
  
  Krows=size(K%ia)-1
  matdescra = 'SUNF'
  
  if (Krows.ne.size(b)) stop "matmulTrans sparse matrix dimensions does not match"
  
  if (allocated(matmul_sparse_sym))  deallocate(matmul_sparse_sym)
  allocate(matmul_sparse_sym(K%cols),stat=ierr)
  matmul_sparse_sym = 0d0
  
  ! use matrix*vector
  call mkl_dcsrmv('n', Krows, K%cols, 1d0, matdescra, K%a, K%ja, K%ia(1), K%ia(2), b, 0d0, matmul_sparse_sym)

end function matmul_sparse_sym




!!! Computes C = A^T*B*A !!!
subroutine spamul_symprod(c,a,b)
  use ISO_C_BINDING
  implicit none
  type(sparse)                   :: a, b, c
  integer                        :: arows, brows, request, info, nzmax, ierr, tmp, i,j
  double precision               :: dtmp
  integer                        :: crows, ccols
  integer(C_INT)                 :: indxing
  TYPE(C_PTR)                    :: crows_start, crows_end, ccol_indx, cvalues
  integer, pointer               :: crows_start_f(:), crows_end_f(:), ccol_indx_f(:)
  double precision, pointer      :: cvalues_f(:)
  type(SPARSE_MATRIX_T)          :: Amkl, Bmkl, Cmkl
  type(MATRIX_DESCR)             :: descrB


!! dimension check
  arows=size(a%ia)-1
  brows=size(b%ia)-1
  if (b%cols.ne.arows) stop "spamul_symprod sparse matrix dimensions do not match"

  descrB%type = SPARSE_MATRIX_TYPE_SYMMETRIC
  descrB%mode = SPARSE_FILL_MODE_UPPER 
  descrB%diag = SPARSE_DIAG_NON_UNIT

  ierr = mkl_sparse_d_create_csr(Amkl,SPARSE_INDEX_BASE_ONE,arows,a%cols,a%ia(1),a%ia(2),a%ja,a%a)
  
  ierr = mkl_sparse_d_create_csr(Bmkl,SPARSE_INDEX_BASE_ONE,brows,b%cols,b%ia(1),b%ia(2),b%ja,b%a)
  
  ierr = mkl_sparse_sypr(SPARSE_OPERATION_TRANSPOSE,Amkl,Bmkl,descrB,Cmkl,SPARSE_STAGE_FULL_MULT)

  if (ierr.ne.SPARSE_STATUS_SUCCESS) then
     write(*,*) 'spamul_symprod: Error'
  endif
  
  ierr = mkl_sparse_order(Cmkl)
  ierr = mkl_sparse_d_export_csr(Cmkl,indxing,crows,ccols,crows_start,crows_end,ccol_indx,cvalues)

  call C_F_POINTER(crows_start, crows_start_f, [crows])
  call C_F_POINTER(crows_end, crows_end_f, [crows])
  call C_F_POINTER(ccol_indx, ccol_indx_f, [crows_end_f(crows)-1])
  call C_F_POINTER(cvalues, cvalues_f, [crows_end_f(crows)-1])
  
  ! Allocate output matrix
  if (allocated(c%a))  deallocate(c%a)
  if (allocated(c%ia)) deallocate(c%ia)
  if (allocated(c%ja)) deallocate(c%ja)
  allocate(c%ia(crows+1),stat=ierr)
  allocate(c%a(crows_end_f(crows)-1),stat=ierr)
  allocate(c%ja(crows_end_f(crows)-1),stat=ierr)

  c%a = cvalues_f
  c%ja = ccol_indx_f
  c%ia = (/crows_start_f(1:crows), crows_end_f(crows)/) 
  c%cols = ccols
 
  ierr = mkl_sparse_destroy(Amkl)
  ierr = mkl_sparse_destroy(Bmkl)
  ierr = mkl_sparse_destroy(Cmkl)

  return
end subroutine spamul_symprod









subroutine spamul(c,a,b)
 ! use ISO_C_BINDING
  implicit none
  type(sparse)                   :: a, b, c
  integer                        :: arows, brows, request, info, nzmax, ierr, tmp
  double precision               :: dtmp
 ! integer                        :: arows, brows, request, info, nzmax, ierr, tmp, i,j
 ! double precision               :: dtmp
 ! 
 ! integer                        :: crows, ccols
 ! integer(C_INT)                 :: indxing
 ! TYPE(C_PTR)                    :: crows_start, crows_end, ccol_indx, cvalues
 ! integer, pointer               :: crows_start_f(:), crows_end_f(:), ccol_indx_f(:)
 ! double precision, pointer      :: cvalues_f(:)
 ! type(SPARSE_MATRIX_T)          :: Amkl, Bmkl, Cmkl

!! dimension check
 ! arows=size(a%ia)-1
 ! brows=size(b%ia)-1
 ! if (a%cols.ne.brows) stop "spamul sparse matrix dimensions do not match"

 ! !ierr = mkl_sparse_d_create_csr(Amkl,SPARSE_INDEX_BASE_ONE,arows,a%cols,a%ia(1:arows),a%ia(2:arows+1),a%ja,a%a)
 ! ierr = mkl_sparse_d_create_csr(Amkl,SPARSE_INDEX_BASE_ONE,arows,a%cols,a%ia,a%ia(2),a%ja,a%a)
 ! 
 ! !ierr = mkl_sparse_d_create_csr(Bmkl,SPARSE_INDEX_BASE_ONE,brows,b%cols,b%ia(1:brows),b%ia(2:brows+1),b%ja,b%a)
 ! ierr = mkl_sparse_d_create_csr(Bmkl,SPARSE_INDEX_BASE_ONE,brows,b%cols,b%ia,b%ia(2),b%ja,b%a)
 ! 
 ! ierr = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE,Amkl,Bmkl,Cmkl)
 ! 
 ! ierr = mkl_sparse_d_export_csr(Cmkl,indxing,crows,ccols,crows_start,crows_end,ccol_indx,cvalues)

 ! call C_F_POINTER(crows_start, crows_start_f, [crows])
 ! call C_F_POINTER(crows_end, crows_end_f, [crows])
 ! call C_F_POINTER(ccol_indx, ccol_indx_f, [crows_end_f(crows)])
 ! call C_F_POINTER(cvalues, cvalues_f, [crows_end_f(crows)])
 ! 
 ! ! Allocate output matrix
 ! if (allocated(c%a))  deallocate(c%a)
 ! if (allocated(c%ia)) deallocate(c%ia)
 ! if (allocated(c%ja)) deallocate(c%ja)
 ! allocate(c%ia(crows+1),stat=ierr)
 ! allocate(c%a(crows_end_f(crows)-1),stat=ierr)
 ! allocate(c%ja(crows_end_f(crows)-1),stat=ierr)

 ! c%a = cvalues_f
 ! c%ja = ccol_indx_f
 ! c%ia = (/crows_start_f(1:crows), crows_end_f(crows)/) 
 ! c%cols = ccols
 ! 
 ! ierr = mkl_sparse_destroy(Amkl)
 ! ierr = mkl_sparse_destroy(Bmkl)
 ! ierr = mkl_sparse_destroy(Cmkl)
 ! 
 ! !do i=1,crows
 ! !  do j=c%ia(i),c%ia(i+1)-1
 ! !     write(*,*) 'row ', i
 ! !     write(*,*) 'col ', j
 ! !     write(*,*) 'value ', c%a(j)
 ! !  enddo
 ! !enddo
  
 
  ! dimension check
  arows=size(a%ia)-1
  brows=size(b%ia)-1
  if (a%cols.ne.brows) stop "spamul sparse matrix dimensions do not match"
  
  if (allocated(c%a))  deallocate(c%a)
  if (allocated(c%ia)) deallocate(c%ia)
  if (allocated(c%ja)) deallocate(c%ja)
  allocate(c%ia(arows+1),stat=ierr)

  request=1
  call mkl_dcsrmultcsr('n', request, 0, arows, a%cols, b%cols, a%a, a%ja, &
                        a%ia, b%a, b%ja, b%ia, dtmp, tmp, c%ia, nzmax, info)
  allocate(c%a(c%ia(arows+1)-1),stat=ierr)
  allocate(c%ja(c%ia(arows+1)-1),stat=ierr)


  request=2
  call mkl_dcsrmultcsr('n', request, 0, arows, a%cols, b%cols, a%a, a%ja, &
                        a%ia, b%a, b%ja, b%ia, c%a, c%ja, c%ia, nzmax, info)
  c%cols=b%cols

  return
end subroutine spamul


subroutine spamulTrans(c,a,b)
  !use ISO_C_BINDING
  implicit none
  type(sparse)                   :: a, b, c
  integer                        :: arows, brows, request, info, nzmax, ierr, tmp
  double precision               :: dtmp
  !integer                        :: arows, brows, request, info, nzmax, ierr, tmp, i,j
  !double precision               :: dtmp
  !
  !integer                        :: crows, ccols
  !integer(C_INT)                 :: indxing
  !TYPE(C_PTR)                    :: crows_start, crows_end, ccol_indx, cvalues
  !integer, pointer               :: crows_start_f(:), crows_end_f(:), ccol_indx_f(:)
  !double precision, pointer      :: cvalues_f(:)
  !type(SPARSE_MATRIX_T)          :: Amkl, Bmkl, Cmkl

! !dimension check
  !arows=size(a%ia)-1
  !brows=size(b%ia)-1
  !if (arows.ne.brows) stop "spamulTrans sparse matrix dimensions do not match"
  !!if (allocated(c%a))  then
  !!write(*,*) 'c%a', c%a
  !!endif
  !
  !!ierr = mkl_sparse_d_create_csr(Amkl,SPARSE_INDEX_BASE_ONE,arows,a%cols,a%ia(1:arows),a%ia(2:arows+1),a%ja,a%a)
  !ierr = mkl_sparse_d_create_csr(Amkl,SPARSE_INDEX_BASE_ONE,arows,a%cols,a%ia,a%ia(2),a%ja,a%a)
  !
  !!ierr = mkl_sparse_d_create_csr(Bmkl,SPARSE_INDEX_BASE_ONE,brows,b%cols,b%ia(1:brows),b%ia(2:brows+1),b%ja,b%a)
  !ierr = mkl_sparse_d_create_csr(Bmkl,SPARSE_INDEX_BASE_ONE,brows,b%cols,b%ia,b%ia(2),b%ja,b%a)
  !
  !ierr = mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE,Amkl,Bmkl,Cmkl)
  !
  !ierr = mkl_sparse_d_export_csr(Cmkl,indxing,crows,ccols,crows_start,crows_end,ccol_indx,cvalues)

  !call C_F_POINTER(crows_start, crows_start_f, [crows])
  !call C_F_POINTER(crows_end, crows_end_f, [crows])
  !call C_F_POINTER(ccol_indx, ccol_indx_f, [crows_end_f(crows)])
  !call C_F_POINTER(cvalues, cvalues_f, [crows_end_f(crows)])
  !
  !! Allocate output matrix
  !if (allocated(c%a))  deallocate(c%a)
  !if (allocated(c%ia)) deallocate(c%ia)
  !if (allocated(c%ja)) deallocate(c%ja)
  !allocate(c%ia(crows+1),stat=ierr)
  !allocate(c%a(crows_end_f(crows)-1),stat=ierr)
  !allocate(c%ja(crows_end_f(crows)-1),stat=ierr)

  !c%a = cvalues_f
  !c%ja = ccol_indx_f
  !c%ia = (/crows_start_f(1:crows), crows_end_f(crows)/) 
  !c%cols = ccols
  !
  !ierr = mkl_sparse_destroy(Amkl)
  !ierr = mkl_sparse_destroy(Bmkl)
  !ierr = mkl_sparse_destroy(Cmkl)
  
  !do i=1,crows
  !  do j=c%ia(i),c%ia(i+1)-1
  !     write(*,*) 'row ', i
  !     write(*,*) 'col ', j
  !     write(*,*) 'value ', c%a(j)
  !  enddo
  !enddo

  !dimension check
  arows=size(a%ia)-1
  brows=size(b%ia)-1
  if (arows.ne.brows) stop "spamulTrans sparse matrix dimensions do not match"
  
  if (allocated(c%a))  deallocate(c%a)
  if (allocated(c%ia)) deallocate(c%ia)
  if (allocated(c%ja)) deallocate(c%ja)

 
  allocate(c%ia(a%cols+1),stat=ierr)
  
  request=1
  call mkl_dcsrmultcsr('t', request, 0, arows, a%cols, b%cols, a%a, a%ja, &
                        a%ia, b%a, b%ja, b%ia, dtmp, tmp, c%ia, nzmax, info)
  
  allocate(c%a(c%ia(a%cols+1)-1),stat=ierr)
  allocate(c%ja(c%ia(a%cols+1)-1),stat=ierr)

  
  request=2
  call mkl_dcsrmultcsr('t', request, 0, arows, a%cols, b%cols, a%a, a%ja, &
                        a%ia, b%a, b%ja, b%ia, c%a, c%ja, c%ia, nzmax, info)
  c%cols=b%cols

  return
end subroutine spamulTrans



subroutine spadefedof1(K,edof)
  implicit none
  type(sparse)                    :: K
  integer                         :: edof(:,:)
  integer                         :: rows,cols

  rows=maxval(edof)
  cols=rows

  call spadefedofm(K,edof,rows,cols)
  
  return
end subroutine spadefedof1


subroutine spadefedof1_sym(K,edof)
  implicit none
  type(sparse)                    :: K
  integer                         :: edof(:,:)
  integer                         :: rows,cols

  rows=maxval(edof)
  cols=rows

  call spadefedofm_sym(K,edof,rows,cols)
  
  return
end subroutine spadefedof1_sym


subroutine spadefedof2(K,edof)
  implicit none
  type(sparse)                    :: K
  integer                         :: edof(:)
  integer, allocatable            :: tmp(:,:)
  integer                         :: rows, cols, ierr

  rows=maxval(edof)
  cols=rows
  allocate(tmp(size(edof),1),stat=ierr)
  tmp(:,1)=edof
  
  call spadefedofm(K,tmp,rows,cols)
  
  return
end subroutine spadefedof2


subroutine spadefedof2_sym(K,edof)
  implicit none
  type(sparse)                    :: K
  integer                         :: edof(:)
  integer, allocatable            :: tmp(:,:)
  integer                         :: rows, cols, ierr

  rows=maxval(edof)
  cols=rows
  allocate(tmp(size(edof),1),stat=ierr)
  tmp(:,1)=edof
  
  call spadefedofm_sym(K,tmp,rows,cols)
  
  return
end subroutine spadefedof2_sym



subroutine spadefedof3(K,edof,rows,cols)
  implicit none
  type(sparse)                    :: K
  integer                         :: edof(:)
  integer, allocatable            :: tmp(:,:)
  integer                         :: rows, cols, ierr

  allocate(tmp(size(edof),1),stat=ierr)
  tmp(:,1)=edof
  
  call spadefedofm(K,tmp,rows,cols)
  
  return
end subroutine spadefedof3


subroutine spadefedof3_sym(K,edof,rows,cols)
  implicit none
  type(sparse)                    :: K
  integer                         :: edof(:)
  integer, allocatable            :: tmp(:,:)
  integer                         :: rows, cols, ierr

  allocate(tmp(size(edof),1),stat=ierr)
  tmp(:,1)=edof
  
  call spadefedofm_sym(K,tmp,rows,cols)
  
  return
end subroutine spadefedof3_sym


subroutine spadefedof4(K,enod,dofnod)
  implicit none
  type(sparse)                    :: K
  integer                         :: enod(:,:), dofnod
  integer, allocatable            :: tmp(:,:)
  integer                         :: ierr, rows, cols, nelm, noel
  integer                         :: ie, ii, tmp1, tmp2

  rows=maxval(enod)*dofnod
  cols=rows

  call spadefedofs(K,enod,dofnod,rows,cols)

  return
end subroutine spadefedof4

subroutine spadefedof4_sym(K,enod,dofnod)
  implicit none
  type(sparse)                    :: K
  integer                         :: enod(:,:), dofnod
  integer, allocatable            :: tmp(:,:)
  integer                         :: ierr, rows, cols, nelm, noel
  integer                         :: ie, ii, tmp1, tmp2

  rows=maxval(enod)*dofnod
  cols=rows

  call spadefedofs_sym(K,enod,dofnod,rows,cols)

  return
end subroutine spadefedof4_sym

subroutine spadefedof5(K,enod,dofnod)
  implicit none
  type(sparse)                    :: K
  integer                         :: enod(:), dofnod
  integer, allocatable            :: tmp(:,:)
  integer                         :: ierr, rows, cols, nelm, noel
  integer                         :: ie, ii, tmp1, tmp2

  rows=maxval(enod)*dofnod
  cols=rows
  allocate(tmp(size(enod),1),stat=ierr)
  tmp(:,1)=enod

  call spadefedofs(K,tmp,dofnod,rows,cols)

  return
end subroutine spadefedof5


subroutine spadefedof5_sym(K,enod,dofnod)
  implicit none
  type(sparse)                    :: K
  integer                         :: enod(:), dofnod
  integer, allocatable            :: tmp(:,:)
  integer                         :: ierr, rows, cols, nelm, noel
  integer                         :: ie, ii, tmp1, tmp2

  rows=maxval(enod)*dofnod
  cols=rows
  allocate(tmp(size(enod),1),stat=ierr)
  tmp(:,1)=enod

  call spadefedofs_sym(K,tmp,dofnod,rows,cols)

  return
end subroutine spadefedof5_sym



subroutine spadefedof6(K,enod,dofnod,rows,cols)
  implicit none
  type(sparse)                    :: K
  integer                         :: enod(:), dofnod
  integer, allocatable            :: tmp(:,:)
  integer                         :: ierr, rows, cols, nelm, noel
  integer                         :: ie, ii, tmp1, tmp2

  allocate(tmp(size(enod),1),stat=ierr)
  tmp(:,1)=enod

  call spadefedofs(K,tmp,dofnod,rows,cols)

  return
end subroutine spadefedof6


subroutine spadefedof6_sym(K,enod,dofnod,rows,cols)
  implicit none
  type(sparse)                    :: K
  integer                         :: enod(:), dofnod
  integer, allocatable            :: tmp(:,:)
  integer                         :: ierr, rows, cols, nelm, noel
  integer                         :: ie, ii, tmp1, tmp2

  allocate(tmp(size(enod),1),stat=ierr)
  tmp(:,1)=enod

  call spadefedofs_sym(K,tmp,dofnod,rows,cols)

  return
end subroutine spadefedof6_sym



subroutine spadefedofs(K,enod,dofnod,rows,cols)
  implicit none
  type(sparse)                    :: K
  integer                         :: enod(:,:), dofnod
  integer, allocatable            :: edof(:,:)
  integer                         :: ierr, rows, cols, nelm, noel
  integer                         :: ie, ii, tmp1, tmp2

  nelm=size(enod,2)
  noel=size(enod,1)
  allocate(edof(noel*dofnod,nelm), stat=ierr)

  do ie=1,nelm
    do ii=1,noel
      tmp1=dofnod*(ii-1)+1
      tmp2=dofnod*(enod(ii,ie)-1)+1
      edof(tmp1:tmp1+dofnod-1,ie)=(/tmp2:(tmp2+dofnod-1)/)
    enddo
    if (debugg.eq.1) then
       write(*,*)'ele ',ie
       write(*,*)'enod ',enod(:,ie)
       write(*,*)'edof ',edof(:,ie)
       if (minval(edof(:,ie)).lt.1) stop
    end if
  enddo
  
  call spadefedofm(K,edof,rows,cols)

  deallocate(edof)

  return
end subroutine spadefedofs



subroutine spadefedofs_sym(K,enod,dofnod,rows,cols)
  implicit none
  type(sparse)                    :: K
  integer                         :: enod(:,:), dofnod
  integer, allocatable            :: edof(:,:)
  integer                         :: ierr, rows, cols, nelm, noel
  integer                         :: ie, ii, tmp1, tmp2

  nelm=size(enod,2)
  noel=size(enod,1)
  allocate(edof(noel*dofnod,nelm), stat=ierr)

  do ie=1,nelm
    do ii=1,noel
      tmp1=dofnod*(ii-1)+1
      tmp2=dofnod*(enod(ii,ie)-1)+1
      edof(tmp1:tmp1+dofnod-1,ie)=(/tmp2:(tmp2+dofnod-1)/)
    enddo
    if (debugg.eq.1) then
       write(*,*)'ele ',ie
       write(*,*)'enod ',enod(:,ie)
       write(*,*)'edof ',edof(:,ie)
       if (minval(edof(:,ie)).lt.1) stop
    end if
  enddo
  
  call spadefedofm_sym(K,edof,rows,cols)

  deallocate(edof)

  return
end subroutine spadefedofs_sym





subroutine spadefedofm(K,edof,rows,cols)
  implicit none
  type(sparse)                    :: K
  integer                         :: edof(:,:), rows, cols

  integer                         :: num(rows), tmp(size(edof,1))
  integer, allocatable            :: ed(:,:)
  integer                         :: nelm, dofe, ndof
  integer                         :: ie, ierr, numax, id, loc, row
  integer                         :: n

  if (debugg.eq.1) write(*,*)'spadefedofm'
  
  if (allocated(K%a))  deallocate(K%a)
  if (allocated(K%ia)) deallocate(K%ia)
  if (allocated(K%ja)) deallocate(K%ja)

  K%cols=cols

  nelm=size(edof,2)
  dofe=size(edof,1)
  ndof=maxval(edof)
  n=maxval(edof)

! check
  if (rows.lt.ndof) stop "row size smaller that largest number in edof"
  if (cols.lt.ndof) stop "column size smaller that largest number in edof"
 
  num=0
  do ie=1,nelm
    num(edof(:,ie))=num(edof(:,ie))+1
  enddo
  ! Numax = maximum number of elements connected to a single node
  numax=maxval(num)
  
! in the first position in ed is the number of elements on the row of ed
! each row corresponds to a dof
! On each row are the dofs associated with the elements that are connected
! to the node of that dof
! I.e. if a corner node, only that elements dofs are on the row
! But otherwise all elements dofs connected to that node are present
  allocate(ed(rows,dofe*numax+1),stat=ierr)
  ed=0
  do ie=1,nelm
    tmp=edof(:,ie)
    do id=1,dofe
      row=tmp(id)
      loc=ed(row,1)+1
      ! Put all dofs on this row of ed
      ed(row,(loc+1):(loc+dofe))=tmp
      ed(row,1)=ed(row,1)+dofe
    enddo
  enddo

  call createsparse(k,ed,rows,numax*dofe) 
 
  deallocate(ed)

  return
end subroutine spadefedofm





subroutine spadefedofm_sym(K,edof,rows,cols)
  implicit none
  type(sparse)                    :: K
  integer                         :: edof(:,:), rows, cols

  integer                         :: num(rows), tmp(size(edof,1))
  integer, allocatable            :: ed(:,:)
  integer                         :: nelm, dofe, ndof
  integer                         :: ie, ierr, numax, id, loc, row
  integer                         :: n

  if (debugg.eq.1) write(*,*)'spadefedofm'
  
  if (allocated(K%a))  deallocate(K%a)
  if (allocated(K%ia)) deallocate(K%ia)
  if (allocated(K%ja)) deallocate(K%ja)

  K%cols=cols

  nelm=size(edof,2)
  dofe=size(edof,1)
  ndof=maxval(edof)
  n=maxval(edof)

! check
  if (rows.lt.ndof) stop "row size smaller that largest number in edof"
  if (cols.lt.ndof) stop "column size smaller that largest number in edof"
 
  num=0
  do ie=1,nelm
    num(edof(:,ie))=num(edof(:,ie))+1
  enddo
  ! Numax = maximum number of elements connected to a single node
  numax=maxval(num)
  
! in the first position in ed is the number of elements on the row of ed
! each row corresponds to a dof
! On each row are the dofs associated with the elements that are connected
! to the node of that dof
! I.e. if a corner node, only that elements dofs are on the row
! But otherwise all elements dofs connected to that node are present
  allocate(ed(rows,dofe*numax+1),stat=ierr)
  ed=0
  do ie=1,nelm
    tmp=edof(:,ie)
    do id=1,dofe
      row=tmp(id)
      loc=ed(row,1)+1
      ! Put all dofs on this row of ed
      ed(row,(loc+1):(loc+dofe))=tmp
      ed(row,1)=ed(row,1)+dofe
    enddo
  enddo

  call createsparse_sym(k,ed,rows,numax*dofe) 
 
  deallocate(ed)

  return
end subroutine spadefedofm_sym







subroutine spadefcell1(K,cell)
  implicit none
  type(sparse)                    :: K
  integer                         :: cell(:,:)
  integer                         :: rows,cols

  rows=maxval(cell(:,1))
  cols=maxval(cell(:,2))

  call spadefcellm(K,cell,rows,cols)
  
  return
end subroutine spadefcell1

subroutine spadefcell2(K,cell)
  implicit none
  type(sparse)                    :: K
  integer                         :: cell(:)
  integer                         :: tmp(1,2)
  integer                         :: rows, cols

  rows=cell(1)
  cols=cell(2)
  tmp(1,:)=cell
  
  call spadefcellm(K,tmp,rows,cols)
  
  return
end subroutine spadefcell2


subroutine spadefcell3(K,cell,rows,cols)
  implicit none
  type(sparse)                    :: K
  integer                         :: cell(:)
  integer                         :: tmp(1,2)
  integer                         :: rows, cols

  tmp(1,:)=cell
  
  call spadefedofm(K,tmp,rows,cols)
  
  return
end subroutine spadefcell3



subroutine spaDefCellm(K,cell,rows,cols)
  implicit none
  type(sparse)                    :: K
  integer                         :: cell(:,:), rows, cols

  integer                         :: num(rows)
  integer, allocatable            :: ed(:,:)
  integer                         :: ie, nelm, ierr, numax, row, loc
  
  if (allocated(K%a))  deallocate(K%a)
  if (allocated(K%ia)) deallocate(K%ia)
  if (allocated(K%ja)) deallocate(K%ja)

  K%cols=cols
  
! check
  if (rows.lt.maxval(cell(:,1)))  &
      stop "row size smaller that largest number in cell"
  if (cols.lt.maxval(cell(:,2)))  &
      stop "column size smaller that largest number in cell"
 
  nelm=size(cell,1)
  num=0
  do ie=1,nelm
    num(cell(ie,1))=num(cell(ie,1))+1
  enddo
  numax=maxval(num)
  
! in the first position in ed is the number of elements on the row
  allocate(ed(rows,numax+1),stat=ierr)
  ed=0
  do ie=1,nelm
    row=cell(ie,1)
    loc=ed(row,1)+1
    ed(row,loc+1)=cell(ie,2)
    ed(row,1)=loc
  enddo
 
  call createsparse(k,ed,rows,numax) 
 
  deallocate(ed)
  
end subroutine spaDefCellm
 

subroutine createsparse(K,ed,rows,numd)
  implicit none
  type(sparse)                    :: K
  integer                         :: ed(:,:), rows, numd
  
  integer                         :: ie, res, inc, id, loc, numax, ierr
  integer, allocatable            :: val(:)
  
  if (debugg.eq.1) write(*,*)'createsparse'

  ! numd = maximum number of elements connected to a single node times the
  ! number of element dofs
  allocate(val(numd),stat=ierr)

  val=0
  !!! This loop is for removing multiple copies of same dof !!!
  do ie=1,rows
    ! if there are elements connected to this dofs node (when are there not??)
    if (ed(ie,1).ne.0) then
      ! Extract all element dofs connected to this dofs node
      val=ed(ie,2:(ed(ie,1)+1))
      ! Sort dofs in ascending order
      call bubble(val,ed(ie,1))
      ! remove redundant/multiple values
      res=val(1)
      inc=1
      do id=2,ed(ie,1)
        if (res.ne.val(id)) then
          inc=inc+1
          val(inc)=val(id)
          res=val(id)
        endif
      enddo
      ed(ie,1)=inc
      ed(ie,2:size(ed,2))=0
      ed(ie,2:(inc+1))=val(1:inc)
    endif
  enddo

  deallocate(val)

  numax=sum(ed(:,1))
  allocate(K%a(numax),stat=ierr)
  allocate(K%ja(numax),stat=ierr)
  allocate(K%ia(rows+1),stat=ierr)

  K%a=0

  loc=1
  K%ia(1)=1

  do ie=1,rows
    inc=ed(ie,1)
    K%ja(loc:(loc+inc-1))=ed(ie,2:inc+1)
    loc=loc+inc
    K%ia(ie+1)=loc
  enddo

end subroutine createsparse  



subroutine createsparse_sym(K,ed,rows,numd)
  implicit none
  type(sparse)                    :: K
  integer                         :: ed(:,:), rows, numd
  
  integer                         :: ie, res, inc, id, loc, numax, ierr
  integer, allocatable            :: val(:)
  
  if (debugg.eq.1) write(*,*)'createsparse'

  ! numd = maximum number of elements connected to a single node times the
  ! number of element dofs
  allocate(val(numd),stat=ierr)

  val=0
  !!! This loop is for removing multiple copies of same dof !!!
  do ie=1,rows
    ! if there are elements connected to this dofs node (when are there not??)
    if (ed(ie,1).ne.0) then
      ! Extract all element dofs connected to this dofs node
      val=ed(ie,2:(ed(ie,1)+1))
      ! Sort dofs in ascending order
      call bubble(val,ed(ie,1))
      ! remove redundant/multiple values
      res=val(1)
      inc=1
      do id=2,ed(ie,1)
        if (res.ne.val(id)) then
          inc=inc+1
          val(inc)=val(id)
          res=val(id)
        endif
      enddo
      ed(ie,1)=inc
      ed(ie,2:size(ed,2))=0
      ed(ie,2:(inc+1))=val(1:inc)
    endif
  enddo

  deallocate(val)
  
  numax = 0
  ! Now only include upper triangular pattern!
  do ie=1,rows
      numax = numax + ed(ie,1)
      id = 1
      ! Remove indices below diagonal
      do while (ed(ie,1+id).lt.ie)
         numax = numax - 1
         id = id + 1
      enddo

      ed(ie,1) = ed(ie,1) - (id-1)
      ed(ie,2:ed(ie,1)+1) = ed(ie,id+1:ed(ie,1)+id)
      ed(ie,(ed(ie,1)+1)+1:size(ed,2)) = 0
  enddo

 ! do ie=1,rows
 !    write(*,*) 'i', ie, 'ed(i,:)', ed(ie,:)
 ! enddo
 ! write(*,*) 'numax', numax
  
  allocate(K%a(numax),stat=ierr)
  allocate(K%ja(numax),stat=ierr)
  allocate(K%ia(rows+1),stat=ierr)

  K%a=0

  loc=1
  K%ia(1)=1
  do ie=1,rows
    inc=ed(ie,1)
    K%ja(loc:(loc+inc-1))=ed(ie,2:inc+1)
    loc=loc+inc
    K%ia(ie+1)=loc
  enddo

end subroutine createsparse_sym  



!return p,q in ascending order
!Subroutine order(p,q)
!  implicit none
!  integer :: p,q,temp
!  if (p.gt.q) then
!    temp=p
!    p=q
!    q=temp
!  end if
!  return
!end subroutine order

!Bubble sorting of integer array Amat
Subroutine bubble(Amat, n)
  implicit none
  integer :: Amat(:), i, j, n, temp
  
  do i=1, n
    do j=n, i+1, -1
        if (Amat(j-1).gt.Amat(j)) then
          temp=Amat(j-1)
          Amat(j-1)=Amat(j)
          Amat(j)=temp
        end if
    end do
  end do
  return
end subroutine bubble





!=======================================================================

!         Create the total (real + imaginary) system matrices
!                          Sparse format

!=======================================================================
subroutine spa_real_imaginary(K_re_im, M_re_im, K, M, ndof)
   implicit none
   type(sparse)      :: K_Re_im, M_re_im, K, M
   integer           :: ndof, ierr

   write(*,*)
   write(*,*) 'Start to create sparse pattern for K and M'
   
   if (allocated(K_Re_Im%a)) deallocate(K_Re_Im%a)
   if (allocated(K_Re_Im%ja)) deallocate(K_Re_Im%ja)
   if (allocated(K_Re_Im%ia)) deallocate(K_Re_Im%ia)

   if (allocated(M_Re_Im%a)) deallocate(M_Re_Im%a)
   if (allocated(M_Re_Im%ja)) deallocate(M_Re_Im%ja)
   if (allocated(M_Re_Im%ia)) deallocate(M_Re_Im%ia)

   allocate(K_Re_Im%a(size(K%a)*2),stat=ierr)
   allocate(M_Re_Im%a(size(M%a)*2),stat=ierr)
   
   K_Re_Im%a = 0d0
   M_Re_Im%a = 0d0

   allocate(K_Re_Im%ja(size(K%ja)*2),stat=ierr)
   allocate(M_Re_Im%ja(size(M%ja)*2),stat=ierr)
   
   K_Re_Im%ja = 0
   M_Re_Im%ja = 0
   
   K_Re_Im%ja(1:size(K%ja)) = K%ja
   K_Re_Im%ja(size(K%ja)+1:size(K_Re_Im%ja)) = ndof + K%ja

   M_Re_Im%ja(1:size(M%ja)) = M%ja
   M_Re_Im%ja(size(M%ja)+1:size(M_Re_Im%ja)) = ndof + M%ja
   
   allocate(K_Re_Im%ia((size(K%ia)-1)*2+1),stat=ierr)
   allocate(M_Re_Im%ia((size(M%ia)-1)*2+1),stat=ierr)

   K_Re_Im%ia = 0
   M_Re_Im%ia = 0
   
   K_Re_Im%ia(1:size(K%ia)-1) = K%ia(1:size(K%ia)-1)
   K_Re_Im%ia(size(K%ia):size(K_Re_Im%ia)-1) = K%ia(size(K%ia))-1+K%ia(1:size(K%ia)-1)

   M_Re_Im%ia(1:size(M%ia)-1) = M%ia(1:size(M%ia)-1)
   M_Re_Im%ia(size(M%ia):size(M_Re_Im%ia)-1) = M%ia(size(M%ia))-1+M%ia(1:size(M%ia)-1)
   
   ! Last element in %ia is the total number of nonzero elements + 1
   K_Re_Im%ia(size(K_Re_Im%ia)) = size(K_Re_Im%a)+1
   M_Re_Im%ia(size(M_Re_Im%ia)) = size(M_Re_Im%a)+1

   K_Re_Im%cols = K%cols*2
   M_Re_Im%cols = M%cols*2
  
   !allocate(K_Re_Imfull(ndof*2,ndof*2),stat=ierr)
   !allocate(M_Re_Imfull(ndof*2,ndof*2),stat=ierr)
   
   !allocate(Kstar_Re_Imfull(dof_m*2,dof_m*2),stat=ierr)
   !allocate(Mstar_Re_Imfull(dof_m*2,dof_m*2),stat=ierr)

 
   write(*,*) 'Sparse pattern created'
   write(*,*)
  

   return
end subroutine spa_real_imaginary











end module sparse_util
