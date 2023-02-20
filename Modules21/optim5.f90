module optim5

! last modified
! M. Wallin

implicit none

!interface mises2
!   module procedure mises1
!   module procedure mises2
!end interface

interface MMAroutine
   module procedure MMAroutine1
end interface

interface MMAsub
   module procedure MMAsub
end interface



!-------------------------------------------------------------------------------

contains

!call MMAroutine1(rho,g0prime,rhoold1,rhoold2,&
!           xlow,xup,alfa,beta,Velm,nelm,Volfrac,geps, f0val, isimp )


!subroutine MMAroutine1(rho,Wp,dFdrho, matptr,rhoold1,rhoold2,&
!           xlow,xupp,alfa,beta,Ve,nelm,volumefrac,geps)

subroutine MMAroutine1(rho,dFdrho,rhoold1,rhoold2,&
           xlow,xupp,alfa,beta,Ve,nelm,volumefrac,geps, f0val, iter)


  use omp_lib
  implicit none
  
!  integer*8           :: matptr
  double precision     :: dFdrho(:), rho(:),xlow(:),xupp(:) 
  
  integer             ::iter,M,N,IYFREE,nelm,ierr
  double precision    ::GEPS,rhoold1(:),rhoold2(:),f0val,xmma(nelm),volumefrac
  double precision    ::Z,ve(:),fmax(1),fval(1),GRADF(1),DSRCH(1),HESSF(1),finish,start
  double precision, allocatable ::xmin(:),xmax(:),xold1(:),xold2(:),xval(:),alfa(:),beta(:)
  double precision, allocatable ::A(:),B(:),C(:),Y(:),ULAM(:),df0dx(:),dfdx(:),P(:)
  double precision, allocatable ::Q(:),P0(:),Q0(:),UU(:)

  M=1 ! Only one constraint
  
!  call matputval(matptr,'dFdphi',dFdrho)
!  call matputval(matptr,'Wptot', Wptot)
!  call matputval(matptr,'rho',rho)
!  call matputval(matptr,'force',force)
!  call matputval(matptr,'am',am)
!  call matcommand(matptr,'mma_routine')

!C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
!C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
!C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
!C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
!C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
!C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
!C    IYFREE(M)


  allocate(xmin(nelm), stat=ierr)
  allocate(xmax(nelm), stat=ierr)
  allocate(xval(nelm), stat=ierr)
  allocate(xold1(nelm), stat=ierr)
  allocate(xold2(nelm), stat=ierr)
  allocate(alfa(nelm), stat=ierr)
  allocate(beta(nelm), stat=ierr)
  allocate(A(M), stat=ierr)
  allocate(B(M), stat=ierr)
  allocate(C(M), stat=ierr)
  allocate(Y(M), stat=ierr)
!  allocate(Z(nelm), stat=ierr)
  allocate(ULAM(M), stat=ierr)
  allocate(df0dx(nelm), stat=ierr)
  allocate(dfdx(nelm), stat=ierr)
  allocate(P(nelm*M), stat=ierr)
  allocate(Q(M*nelm), stat=ierr)
  allocate(P0(nelm), stat=ierr)
  allocate(Q0(nelm), stat=ierr)
  allocate(UU(M), stat=ierr)


  iter=1 ! Varför väljs iter=1 hela tiden ????, Borde uppdateras så att asymptoterna flyttas ?!?!

  dFDX=Ve ! Derivative of volume constraints 
  geps=1d-5 ! Tolerance parameter for the constraints.
 
!  df0dx=-dFdrho ! Derivative of the objective functional

  IYFREE=1     ! Vector of size M: (nbr of constraints)

  xval=rho   ! Design variable
  xmin=1d-3
  xmax=1d0

!  f0val Objective value, now as an input
!  N=nelm

  fval=dot_product(rho,ve)-sum(ve)*volumefrac ! Value of volume constraint
 
  fmax=0d0 ! Right hand side of constraint, i.e. volume constraint

! Plastic work maximizarion  df0dx=-dFdrho 
  df0dx=dFdrho

  xold1=rhoold1
  xold2=rhoold2
   
  A=0d0
  
  C=1000d0 !kolla värde
!  write(*,*) 'rhoold1 innan', rhoold1(1:3)
!  write(*,*) 'rhoold2 innan', rhoold2(1:3)


  call cpu_time(start)
 call MMASUB(ITER,M,Nelm,GEPS,IYFREE,XVAL,XMMA, &
                   XMIN,XMAX,XOLD1,XOLD2,XLOW,XUPP, &
                   ALFA,BETA,A,B,C,Y,Z,ULAM, &
                   F0VAL,FVAL,FMAX,DF0DX,DFDX, &
                   P,Q,P0,Q0,UU,GRADF,DSRCH,HESSF)
 !finish=OMP_GET_WTIME()
 call cpu_time(finish)
! write(*,*) 'MMASUB time',finish-start 
!write(*,*) shape(xmma),shape(xold1),shape(xold2)
rho=xmma
rhoold1=xold1
rhoold2=xold2


!  pause
!  write(*,*) 'rhoold1 efter', rhoold1(1:3)
!  write(*,*) 'rhoold2 efter', rhoold2(1:3)


!write(*,*) 'ieie'
!write(*,*) '------------------------z',z,abs(y)

!pause

  
end subroutine MMAroutine1


!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-!
!********+*********+*********+*********+*********+*********+*********+
!
      SUBROUTINE MMASUB(ITER,M,N,GEPS,IYFREE,XVAL,XMMA, &
                       XMIN,XMAX,XOLD1,XOLD2,XLOW,XUPP, &
                       ALFA,BETA,A,B,C,Y,Z,ULAM, &
                       F0VAL,FVAL,FMAX,DF0DX,DFDX, &
                       P,Q,P0,Q0,UU,GRADF,DSRCH,HESSF)
                       
    
!C       Version "December 2006".
!C    !-----------------------------------------!
!C    !  The author of this subroutine is       !
!C    !  Krister Svanberg <krille@math.kth.se>  !
!C    !-----------------------------------------!
!C
!C    Use of this code is for academic purposes only,
!C    regulated by an agreement with Krister Svanberg.
!C    The code is not to be redistributed.
!C
!C    MMASUB generates and solves the MMA subproblem,
!C    which is of the following form in the variables
!C    x_1,...,x_N, y_1,...,y_M, and z.
!C
!C   minimize h_0(x) + r_0 + z + 0.05*z^2 + sum{c_i*y_i + 0.5*(y_i)^2}
!C
!C subject to h_i(x) - a_i*z - y_i <= b_i ,     i=1,..,M
!C                   alfa_j <= x_j <= beta_j ,  j=1,..,N
!C                             y_i >= 0 ,       i=1,..,M
!C                               z >= 0 .
!C
!C    with h_i(x) = sum{p_ij/(xupp_j-x_j) + q_ij/(x_j-xlow_j)}.
!C
    !  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!C
      INTEGER  :: IYFREE,iter,m,n
     ! DIMENSION XVAL(N),XMMA(N),XMIN(N),XMAX(N), &
     !          XOLD1(N),XOLD2(N),XLOW(N),XUPP(N), &
     !          ALFA(N),BETA(N), &
     !          DF0DX(N),DFDX(M*N), &
     !          P(M*N),Q(M*N),P0(N),Q0(N)!, &
               
      double precision::XVAL(N),XMMA(N),XMIN(N),XMAX(N)
      double precision::XOLD1(N),XOLD2(N),XLOW(N),XUPP(N)
      double precision::ALFA(N),BETA(N)
      double precision::DF0DX(N),DFDX(M*N)
      double precision::P(M*N),Q(M*N),P0(N),Q0(N)!, &
              
               
      double precision::fval(M),fmax(M),UU(M),GRADF(M),DSRCH(M)
      double precision::A(M),B(M),C(M),Y(M),ULAM(M),HESSF(M*(M+1)/2),R0(M)
      double precision::geps,z,f0val
      
!C
!C********+*********+*********+*********+*********+*********+*********+
!C  The sizes of the above areas must be at least as follows:
!C
!C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
!C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
!C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
!C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
!C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
!C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
!C    IYFREE(M)
!C********+*********+*********+*********+*********+*********+*********+
!C
!C***  Input to the subroutine MMASUB:
!C
!C  ITER  = Current iteration number ( =1 the first iteration).
!C     N  = Number of variables x_j in the problem.
!C     M  = Number of constraints in the problem (not including
!C          the simple upper and lower bounds on the variables).
!C  GEPS  = Tolerance parameter for the constraints.
!C          (Used in the termination criteria for the subproblem.)
!C   XVAL(j) = Current value of the variable x_j.
!C   XMIN(j) = Original lower bound for the variable x_j.
!C   XMAX(j) = Original upper bound for the variable x_j.
!C  XOLD1(j) = Value of the variable x_j one iteration ago.
!C  XOLD2(j) = Value of the variable x_j two iterations ago.
!C   XLOW(j) = Current value of the lower asymptot l_j.
!C   XUPP(j) = Current value of the upper asymptot u_j.
!C      A(i) = Coefficient a_i for the minimax variable z.
!C      C(i) = Coefficient c_i for the artificial variable y_i.
!C    F0VAL  = Value of the objective function f_0(x)
!C   FVAL(i) = Value of the i:th constraint function f_i(x).
!C   FMAX(i) = Right hand side of the i:th constraint.
!C  DF0DX(j) = Derivative of f_0(x) with respect to x_j.
!C   DFDX(k) = Derivative of f_i(x) with respect to x_j,
!C             where k = (j-1)*M + i.
!C
!C*** Output from the subroutine MMASUB:
!C
!C   XMMA(j) = Optimal value of x_j in the MMA subproblem.
!C      Y(i) = Optimal value of the "artificial" variable y_i.
!C      Z    = Optimal value of the "minimax" variable z.
!C   ULAM(i) = Optimal value of the dual variable lambda_i.
!C   XLOW(j) = New value on the lower asymptot l_j.
!C   XUPP(j) = New value on the upper asymptot u_j.
!C
!C*** Working areas and their usage in MMASUB:
!C
!C   ALFA(j) = Lower bound for x_j in the MMA subproblem.
!C   BETA(j) = Upper bound for x_j in the MMA subproblem.
!C      P(k) = Coefficient p_ij in the MMA subproblem,
!C             where k = (j-1)*M + i.
!C      Q(k) = Coefficient q_ij in the MMA subproblem,
!C             where k = (j-1)*M + i.
!C     P0(j) = Coefficient p_0j in the MMA subproblem.
!C     Q0(j) = Coefficient q_0j in the MMA subproblem.
!C      B(i) = Right hand side b_i in the MMA subproblem.
!C  GRADF(i) = Gradient component of the dual objective function.
!C  DSRCH(i) = Search direction component in the dual subproblem.
!C  HESSF(k) = Hessian matrix component of the dual function.
!C     UU(i) = Component in a working area.
!C IYFREE(i) = 0 for dual variables which are fixed to zero in
!C               the current subspace of the dual subproblem,
!C           = 1 for dual variables which are "free" in
!C               the current subspace of the dual subproblem.
!C
!C********+*********+*********+*********+*********+*********+*********+
!C
     ! write(*,*) 'in MMAsub'
      
      CALL ASYMPT(ITER,M,N,XVAL,XMIN,XMAX,XOLD1,XOLD2, &
                 XLOW,XUPP,ALFA,BETA)
     ! write(*,*) 'after asympt'
!C
!C****  ASYMPT calculates the asymptotes XLOW(j) and XUPP(j),
!C****  and the bounds ALFA(j) and BETA(j).
!C
      CALL GENSUB(M,N,XVAL,XMIN,XMAX,F0VAL,DF0DX,FMAX,FVAL, &
                 DFDX,P,Q,B,P0,Q0,R0,XLOW,XUPP)
     !   write(*,*) 'after gensub'
!C
!C***** GENSUB generates the MMA subproblem by calculating the
!C***** coefficients P(i,j),Q(i,j),B(i),P0(j),Q0(j) and R0.
!C  
!write(*,*) size(xmma),'size xmma b4 maxim'   
!write(*,*) xmma(1), 'xmma(1)'


      CALL MAXIM(M,N,GEPS,IYFREE,GRADF,DSRCH,HESSF,XMMA,Y,Z, &
                ULAM,UU,XLOW,XUPP,ALFA,BETA,A,B,C,P,Q,P0,Q0)
   ! write(*,*) 'after maxim'
!C
!C***** MAXIM solves the dual problem of the MMA subproblem.
!C***** ULAM = optimal solution of this dual problem.
!C***** XMMA,Y,Z = optimal solution of the MMA subproblem.
!C
      RETURN
      END subroutine MMASUB
!C
!C********+*********+*********+*********+*********+*********+*********+
!C
      SUBROUTINE ASYMPT(ITER,M,N,XVAL,XMIN,XMAX,XOLD1,XOLD2, &
                       XLOW,XUPP,ALFA,BETA)
!C
!C       Version "December 2006".
!C    !-----------------------------------------!
!C    !  The author of this subroutine is       !
!C    !  Krister Svanberg <krille@math.kth.se>  !
!C    !-----------------------------------------!
!C
!C     ASYMPT calculates the asymptotes XLOW and XUPP,
!C     and the bounds ALFA and BETA, for the current subproblem.
!C
!C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
!C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
!C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
!C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
!C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
!C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
!C    IYFREE(M)

      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer ::iter,M,N,J,I
!      DIMENSION XVAL(1),XMIN(1),XMAX(1),XOLD1(1),XOLD2(1), &
!               XLOW(1),XUPP(1),ALFA(1),BETA(1)
       
       double precision :: XVAL(N),XMIN(N),XMAX(N)
       double precision :: XOLD1(N),XOLD2(N),XLOW(N),XUPP(N) 
       double precision :: ALFA(N),BETA(N),ALBEFA,GHINIT,GHDECR,GHINCR,FAK
       double precision :: GMINJ,GMAXJ,HMINJ,HMAXJ,XTEST,XMMJ,XMIJ,XMAJ
               
      !write(*,*) 'ASYMPT 1'
!C
      ALBEFA=0.1
      GHINIT=0.5
!      GHINIT=0.1
! Write(*,*) 'GHINIT=0.1 inne mmasub'
      GHDECR=0.7
      GHINCR=1.2 
      !write(*,*) 'ASYMPT 1.5'
!C      GHINIT=0.1
!C      GHDECR=0.95
!C      GHINCR=1.0
      IF(ITER.GE.3) GOTO 350
!C
      !write(*,*) 'ASYMPT 2'
!C***  Here ITER = 1 or 2 .
      DO 200 J=1,N
      XMMJ=XMAX(J)-XMIN(J)
      IF(XMMJ.LT.0.00001) XMMJ=0.00001
      XLOW(J)=XVAL(J)-GHINIT*XMMJ
      XUPP(J)=XVAL(J)+GHINIT*XMMJ
      !write(*,*) 'ASYMPT 3'
  200 CONTINUE
      GOTO 500
!C
!C***  Here ITER is greater than 2.
  350 CONTINUE
!C
    ! write(*,*) 'ASYMPT 4'
      DO 400 J=1,N
      !write(*,*) 'ASYMPT 5'
      !write(*,*) shape(xval)
      !write(*,*) shape(xold1)
      !write(*,*) shape(xold2)

      XTEST=(XVAL(J)-XOLD1(J))*(XOLD1(J)-XOLD2(J))
     ! write(*,*) 'ASYMPT 5.1'
      FAK=1.0
      IF(XTEST.LT.0.) FAK=GHDECR
      IF(XTEST.GT.0.) FAK=GHINCR
      !write(*,*) 'ASYMPT 5.2'
      XLOW(J)=XVAL(J)-FAK*(XOLD1(J)-XLOW(J))
      XUPP(J)=XVAL(J)+FAK*(XUPP(J)-XOLD1(J))
      !write(*,*) 'ASYMPT 5.3'
      XMMJ=XMAX(J)-XMIN(J)
     ! write(*,*) 'ASYMPT 5.4'
      IF(XMMJ.LT.0.00001) XMMJ=0.00001
      GMINJ = XVAL(J)-10.0*XMMJ
      GMAXJ = XVAL(J)-0.01*XMMJ
     ! write(*,*) 'ASYMPT 5.5'
      HMINJ = XVAL(J)+0.01*XMMJ
      HMAXJ = XVAL(J)+10.0*XMMJ
      IF(XLOW(J).LT.GMINJ) XLOW(J)=GMINJ
      IF(XLOW(J).GT.GMAXJ) XLOW(J)=GMAXJ
      IF(XUPP(J).LT.HMINJ) XUPP(J)=HMINJ
      IF(XUPP(J).GT.HMAXJ) XUPP(J)=HMAXJ
     !write(*,*) 'ASYMPT 6'
  400 CONTINUE
!C
  500 CONTINUE
!C
     ! write(*,*) 5
      DO 600 J=1,N
      XMIJ=XMIN(J)-0.000001
      XMAJ=XMAX(J)+0.000001
      IF(XVAL(J).GE.XMIJ) GOTO 550
      XLOW(J)=XVAL(J)-(XMAJ-XVAL(J))/0.9
      XUPP(J)=XVAL(J)+(XMAJ-XVAL(J))/0.9
      !write(*,*) 6
      GOTO 600
  550 CONTINUE
      IF(XVAL(J).LE.XMAJ) GOTO 600
      XLOW(J)=XVAL(J)-(XVAL(J)-XMIJ)/0.9
      XUPP(J)=XVAL(J)+(XVAL(J)-XMIJ)/0.9
     ! write(*,*) 7
  600 CONTINUE
!C
      DO 700 J=1,N
      ALFA(J)=XLOW(J)+ALBEFA*(XVAL(J)-XLOW(J))
      BETA(J)=XUPP(J)-ALBEFA*(XUPP(J)-XVAL(J))
      IF(ALFA(J).LT.XMIN(J)) ALFA(J)=XMIN(J)
      IF(BETA(J).GT.XMAX(J)) BETA(J)=XMAX(J)
      !write(*,*) 8
  700 CONTINUE
!C
      !write(*,*) 9
      RETURN
      END subroutine ASYMPT
!C
!C********+*********+*********+*********+*********+*********+*********+
!C
      SUBROUTINE GENSUB(M,N,XVAL,XMIN,XMAX,F0VAL,DF0DX,FMAX,FVAL, &
                       DFDX,P,Q,B,P0,Q0,R0,XLOW,XUPP)
!C
!C       Version "December 2006".
!C    !-----------------------------------------!
!C    !  The author of this subroutine is       !
!C    !  Krister Svanberg <krille@math.kth.se>  !
!C    !-----------------------------------------!
!C
!C     GENSUB calculates P( ),Q( ),B( ),P0( ),Q0( ) and R0
!C     for the current subproblem.
!C
        
!C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
!C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
!C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
!C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
!C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
!C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
!C    IYFREE(M)
          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer ::M,N,I,J,MJ,IJ
 !     DIMENSION XVAL(1),XMIN(1),XMAX(1),DF0DX(1),FMAX(1),FVAL(1), &
!               DFDX(1),P(1),Q(1),B(1),P0(1),Q0(1),XLOW(1),XUPP(1)
!C
  
    
      DIMENSION XVAL(N),XMIN(N),XMAX(N), &
               XLOW(N),XUPP(N), &
               DF0DX(N),DFDX(M*N), &
               P(M*N),Q(M*N),P0(N),Q0(N)!, &
               
      double precision::fval(M),fmax(M),UU(M),GRADF(M),DSRCH(M)
      double precision::A(M),B(M),C(M),Y(M),ULAM(M),HESSF(M*(M+1)/2),R0(M)
      
      
      RAA0=0.00001
      R0=F0VAL
      DO 20 I=1,M
      B(I)=FMAX(I)-FVAL(I)
   20 CONTINUE
!C
      DO 50 J=1,N
      MJ=M*(J-1)
      UJLJ=XUPP(J)-XLOW(J)
      UJXJ=XUPP(J)-XVAL(J)
      XJLJ=XVAL(J)-XLOW(J)
      UJXJ2=UJXJ*UJXJ
      XJLJ2=XJLJ*XJLJ
      XMMJ=XMAX(J)-XMIN(J)
      IF(XMMJ.LT.0.00001) XMMJ=0.00001
      P0J=RAA0/XMMJ
      Q0J=RAA0/XMMJ
      IF(DF0DX(J).GT.0.) P0J=P0J+1.001*DF0DX(J)
      IF(DF0DX(J).GT.0.) Q0J=Q0J+0.001*DF0DX(J)
      IF(DF0DX(J).LT.0.) Q0J=Q0J-1.001*DF0DX(J)
      IF(DF0DX(J).LT.0.) P0J=P0J-0.001*DF0DX(J)
      P0J=P0J*UJXJ2
      Q0J=Q0J*XJLJ2
      P0(J)=P0J
      Q0(J)=Q0J
      R0=R0-P0J/UJXJ-Q0J/XJLJ
!C
      DO 40 I=1,M
      IJ=MJ+I
      PIJ=RAA0/XMMJ
      QIJ=RAA0/XMMJ
      DFIJ=DFDX(IJ)
      IF(DFIJ.GT.0.) PIJ=PIJ+1.001*DFIJ
      IF(DFIJ.GT.0.) QIJ=QIJ+0.001*DFIJ
      IF(DFIJ.LT.0.) QIJ=QIJ-1.001*DFIJ
      IF(DFIJ.LT.0.) PIJ=PIJ-0.001*DFIJ
      PIJ=PIJ*UJXJ2
      QIJ=QIJ*XJLJ2
      P(IJ)=PIJ
      Q(IJ)=QIJ
      B(I)=B(I)+PIJ/UJXJ+QIJ/XJLJ
!C
   40 CONTINUE
   50 CONTINUE
!C
      RETURN
      END subroutine GENSUB
!C
!C********+*********+*********+*********+*********+*********+*********+
!C********+*********+*********+*********+*********+*********+*********+
!C
      SUBROUTINE MAXIM(M,N,GEPS,IYFREE,GRADF,DSRCH,HESSF,X,Y,Z,ULAM, &
                      UU,XLOW,XUPP,ALFA,BETA,A,B,C,P,Q,P0,Q0)
!C
!C       Version "December 2006".
!C    !-----------------------------------------!
!C    !  The author of this subroutine is       !
!C    !  Krister Svanberg <krille@math.kth.se>  !
!C    !-----------------------------------------!
!C
!C     MAXIM solves the dual MMA subproblem.
!C     The dual variables are ulam(i), i=1,..,m,
!C     which are required to be non-negative.
!C
!C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
!C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
!C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
!C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
!C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
!C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
!C    IYFREE(M)


     ! IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER ::M,N,ITR,M3,I,IHITY,IGMX
      INTEGER ::IYFREE(M)
 !     DIMENSION GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),X(M*N),Y(M),ULAM(M), &
 !              UU(M),XLOW(N),XUPP(N),ALFA(N),BETA(N), &
 !              A(M),B(M),C(M),P(M*N),Q(M*N),P0(N),Q0(N)
      double precision :: GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),X(M*N),Y(M),ULAM(M)
      double precision :: UU(M),XLOW(N),XUPP(N),ALFA(N),BETA(N)
      double precision :: A(M),B(M),C(M),P(M*N),Q(M*N),P0(N),Q0(N),F,GMX,GEPS,Z
      
!C
     ! write(*,*) 'maxim 1'
      ITR=0
      M3=3*M+30
!C
      DO 10 I=1,M
      ULAM(I)=0.
      IYFREE(I)=1
 10   CONTINUE
!C
   !   write(*,*) 'maxim 2'
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA, &
                A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
!C
   !   write(*,*) 'maxim 3'
      GMX=0.
      DO 20 I=1,M
      IYFREE(I)=0
      IF(GRADF(I).GT.GEPS) IYFREE(I)=1
      IF(GRADF(I).GT.GMX) GMX=GRADF(I)
 20   CONTINUE
!C
      IF(GMX.LE.GEPS) GOTO 100
!C     Vi avbryter optimeringen, ulam=0 ar optimal losning.
!C
 30   CONTINUE
      ITR=ITR+1
      IF(ITR.GT.M3) GOTO 100
!C     Vi avbryter optimeringen pga for manga subspa-anrop.
!C
      CALL SUBSPA(ITR,M,N,GEPS,F,IYFREE,GRADF,DSRCH,HESSF, &
                 X,Y,ULAM,UU,XLOW,XUPP,ALFA,BETA,A,B,C, &
                 P,Q,P0,Q0,IHITY)
!C
      IF(IHITY.EQ.0) GOTO 40
!C     Om ihity = 0 sa ar ulam optimal pa aktuellt underrum.
!C     Om ihity > 0 sa har vi slagit i ett nytt bivillkor.
      IYFREE(IHITY)=0
      ULAM(IHITY)=0.
      GOTO 30
!C
 40   CONTINUE
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA, &
                A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
!C
      GMX=0.
      IGMX=0
      DO 50 I=1,M
      IF(IYFREE(I).EQ.1) GOTO 50
      IF(GRADF(I).LE.GMX) GOTO 50
      GMX=GRADF(I)
      IGMX=I
 50   CONTINUE
!C
      IF(GMX.LE.GEPS) GOTO 100
!C     Om gmx =< geps sa ar ulam optimal losning.
!C     Om gmx > geps sa tar vi bort kravet att ulam(igmx)=0.
      IYFREE(IGMX)=1
      GOTO 30
!C
 100  CONTINUE
!C     Nu ar antingen ulam optimal losning eller itr>m3.
      CALL XYZLAM(M,N,X,Y,Z,ULAM,XLOW,XUPP,ALFA,BETA, &
                 A,B,C,P,Q,P0,Q0,IYFREE)
      IF(ITR.GT.M3) WRITE(*,911)
 911  FORMAT(' ITR GT M3 IN MAXIM')
!C
      RETURN
      END subroutine MAXIM
!C
!C********+*********+*********+*********+*********+*********+*********+
!C
      SUBROUTINE SUBSPA(ITR,M,N,GEPS,F,IYFREE,GRADF,DSRCH,HESSF, &
                       X,Y,ULAM,UU,XLOW,XUPP,ALFA,BETA, &
                       A,B,C,P,Q,P0,Q0,IHITY)
!C
!C       Version "December 2006".
!C    !-----------------------------------------!
!C    !  The author of this subroutine is       !
!C    !  Krister Svanberg <krille@math.kth.se>  !
!C    !-----------------------------------------!
!C
!C    SUBSPA maximizes the dual objective function on the subspace
!C    defined by ulam(i) = 0 for every i such that iyfree(i) = 0.
!C    The first three iterations a steepest ascent method is used,
!C    and after that a Newton method is used.
!C

!C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
!C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
!C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
!C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
!C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
!C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
!C    IYFREE(M)
     ! IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    !  DIMENSION GRADF(M),DSRCH(M),HESSF(M*(M+1)),X(M*N),Y(M), &
    !           ULAM(M),UU(M),XLOW(N),XUPP(N),ALFA(N),BETA(N), &
    !           A(M),B(M),C(M),P(M*N),Q(M*N),P0(M),Q0(N)
      INTEGER ::M,N,I,K,Itemax,itesub,nydim,IHITMX,IOPT,IK,IKRED,IRED,ITR,IHITY
      INTEGER ::IYFREE(M)
      double precision:: GRADF(M),DSRCH(M),HESSF(M*(M+1)),X(M*N),Y(M)
      double precision:: ULAM(M),UU(M),XLOW(N),XUPP(N),ALFA(N),BETA(N)
      double precision:: A(M),B(M),C(M),P(M*N),Q(M*N),P0(M),Q0(N)
      double precision:: geps,F,DSRTOL,TMAX0,TMAX,GTD,T,TOPT,HTRACE,ZZZZ,HESMOD
!C
      IHITY=0
      ITESUB=0
      NYDIM=0
      DSRTOL=-0.0000001*GEPS
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA, &
                A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
!C
      DO 10 I=1,M
      DSRCH(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 10
      NYDIM=NYDIM+1
      DSRCH(I)=GRADF(I)
 10   CONTINUE
!C
      IF(NYDIM.EQ.0) GOTO 100
!C     Vi avbryter med ihity = 0, ty inga variabler ulam(i) ar fria.
      ITEMAX=50+5*NYDIM
!C
 15   ITESUB=ITESUB+1
!C     Har startar en ny iteration.
!C
      TMAX0=1.0D8
      TMAX=TMAX0
      IHITY=0
      GTD=0.
      DO 20 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 20
      GTD=GTD+GRADF(I)*DSRCH(I)
      IF(DSRCH(I).GE.0.) GOTO 20
      IF((-DSRCH(I)).LE.(ULAM(I)/TMAX0)) GOTO 20
      IF(DSRCH(I).GT.DSRTOL) GOTO 20
      T=ULAM(I)/(-DSRCH(I))
      IF(T.GE.TMAX) GOTO 20
      TMAX=T
      IHITY=I
 20   CONTINUE
      IF(TMAX.LT.0.) TMAX=0.
      IF(GTD.GT.0.) GOTO 25
      IHITY=0
      WRITE(*,912)
      GOTO 100
!C     Vi avbryter med ihity = 0, ty dsrch ar ej en ascentriktning.
!C
 25   CONTINUE
      CALL LINSE(M,N,ITESUB,IHITMX,IYFREE,TMAX,TOPT,ULAM,DSRCH, &
                X,Y,UU,XLOW,XUPP,ALFA,BETA,A,B,C,P,Q,P0,Q0)
!C
      DO 30 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 30
      ULAM(I)=ULAM(I)+TOPT*DSRCH(I)
      IF(ULAM(I).LT.0.) ULAM(I)=0.
 30   CONTINUE
!C
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA, &
                A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
!C
      IF(IHITMX.EQ.1.AND.IHITY.GT.0) GOTO 100
!C     Vi avbryter med ihity > 0, ty vi har slagit i det tidigare
!C     inaktiva bivillkoret ulam(ihity) >= 0.
      IHITY=0
      IOPT=1
      DO 40 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 40
      IF(DABS(GRADF(I)).GT.GEPS) IOPT=0
 40   CONTINUE
!C
      IF(IOPT.EQ.1) GOTO 100
!C     Vi avbryter med ihity = 0, ty optimal losning hittad.
      IF(ITESUB.GT.ITEMAX) GOTO 97
!C     Vi avbryter med ihity = 0, ty for manga iterationer.
      IF(ITESUB.GE.3) GOTO 55
!C     Om itesub>=3 sa byter vi fran steepest ascent till Newton.
!C     Om itesub=<2 sa fortsatter vi med steepest ascent.
      DO 50 I=1,M
      DSRCH(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 50
      DSRCH(I)=GRADF(I)
 50   CONTINUE
      GOTO 15
!C
 55   CONTINUE
!C
      CALL HESSI(M,N,ULAM,HESSF,X,Y,ALFA,BETA, &
                A,B,C,P,Q,P0,Q0,XLOW,XUPP,IYFREE)
!C
      IK=0
      IKRED=0
      DO 70 K=1,M
      DO 65 I=K,M
      IK=IK+1
      IF(IYFREE(K).EQ.0) GOTO 65
      IF(IYFREE(I).EQ.0) GOTO 65
      IKRED=IKRED+1
      HESSF(IKRED)=HESSF(IK)
 65   CONTINUE
 70   CONTINUE
!C
      HTRACE=0.
      IKRED=0
      ZZZZ=0.
      DO 73 K=1,NYDIM
      DO 72 I=K,NYDIM
      IKRED=IKRED+1
      IF(I.EQ.K) HTRACE=HTRACE+HESSF(IKRED)
      IF(I.EQ.K) ZZZZ=ZZZZ+1.
 72   CONTINUE
 73   CONTINUE
!C
      HESMOD=0.0001*HTRACE/ZZZZ
      IF(HESMOD.LT.GEPS) HESMOD=GEPS
      IKRED=0
      DO 77 K=1,NYDIM
      DO 76 I=K,NYDIM
      IKRED=IKRED+1
      IF(I.EQ.K) HESSF(IKRED)=HESSF(IKRED)+HESMOD
 76   CONTINUE
 77   CONTINUE
!C
      CALL LDLFAC(NYDIM,GEPS,HESSF,UU)
!C
      IRED=0
      DO 79 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 79
      IRED=IRED+1
      UU(IRED)=GRADF(I)
 79   CONTINUE
!C
      CALL LDLSOL(NYDIM,UU,HESSF,DSRCH)
!C
      DO 80 I=1,M
      UU(I)=DSRCH(I)
 80   CONTINUE
!C
      IRED=0
      DO 85 I=1,M
      DSRCH(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 85
      IRED=IRED+1
      DSRCH(I)=UU(IRED)
 85   CONTINUE
!C
      GOTO 15
!C
 97   CONTINUE
      WRITE(*,911)
!C
 100  CONTINUE
!C
      DO 110 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 110
      IF(ULAM(I).LT.0.) ULAM(I)=0.
 110  CONTINUE
!C
 911  FORMAT(' ITESUB GT ITEMAX IN SUBSPA')
 912  FORMAT(' GTD LE 0 IN SUBSPA')
!C
      RETURN
      END subroutine SUBSPA
!C
!C********+*********+*********+*********+*********+*********+*********+
!C
      SUBROUTINE LINSE(M,N,ITESUB,IHITMX,IYFREE,TMAX,TOPT,ULAM,DSRCH, &
                      X,Y,UU,XLOW,XUPP,ALFA,BETA,A,B,C,P,Q,P0,Q0)
!C
!C       Version "December 2006".
!C    !-----------------------------------------!
!C    !  The author of this subroutine is       !
!C    !  Krister Svanberg <krille@math.kth.se>  !
!C    !-----------------------------------------!
!C
!C     LINSE makes an approximate line search (maximization) in the
!C     direction DSRCH from the point ULAM.
!C     Main input:  ULAM, DSRCH, TMAX.
!C     Main output: TOPT, IHITMX.
!C
     !C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
!C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
!C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
!C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
!C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
!C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
!C    IYFREE(M)
 
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER ::M,N,ITESUB,IHITMX,ITT1,ITT2,ITT3
      INTEGER ::IYFREE(M)
      !DIMENSION ULAM(M),DSRCH(M),X(1),Y(1),UU(M),XLOW(N),XUPP(N), &
      !         ALFA(N),BETA(N),A(M),B(M),C(M),P(M*N),Q(M*N),P0(N),Q0(N)
      double precision ::ULAM(M),DSRCH(M),X(M*N),Y(M),UU(M),XLOW(N),XUPP(N)
      double precision :: ALFA(N),BETA(N),A(M),B(M),C(M),P(M*N),Q(M*N),P0(N),Q0(N)
      double precision :: TMAX,TOPT,DFDTMX,T2,T1,DFDT1,DFDT2,SQT1,SQT2,TM,DFDTM,TKVOT
               

!C
      ITT1=0
      ITT2=0
      ITT3=0
!C
      CALL LINDER(M,N,TMAX,DFDTMX,ULAM,DSRCH,X,Y,UU,XLOW,XUPP, &
                 ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDTMX.GE.0.) GOTO 80
!C     Linjesokningen klar. Optimalt steg ar T=TMAX. IHITMX=1.
      IF(TMAX.GT.1.) GOTO 40
      T2=TMAX
!C
 30   CONTINUE
!C     Nu sker en upprepad minskning av steget.
      ITT1=ITT1+1
      IF(ITT1.GT.13) GOTO 90
      T1=T2/2.
      IF(ITESUB.LE.3) T1=T2/16.
      CALL LINDER(M,N,T1,DFDT1,ULAM,DSRCH,X,Y,UU,XLOW,XUPP, &
                 ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDT1.GT.0.) GOTO 60
      T2=T1
      GOTO 30
!C
 40   CONTINUE
!C     Nu testas enhetssteget, dvs T=1.
      T1=1.
      T2=T1
      CALL LINDER(M,N,T1,DFDT1,ULAM,DSRCH,X,Y,UU,XLOW,XUPP, &
                 ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(ITESUB.GE.6.AND.DFDT1.GE.0.) GOTO 90
!C     Linjesokningen klar. Enhetssteget duger. T=1, IHITMX=0.
      IF(ITESUB.LE.5.AND.DFDT1.GT.0.) GOTO 50
!C     Enhetssteget ar for kort.
      GOTO 30
!C     Enhetssteget ar for langt.
!C
 50   ITT2=ITT2+1
!C     Nu sker en upprepad okning av steget.
      IF(ITT2.GT.10) GOTO 90
      T2=2.*T1
      IF(ITESUB.LE.3) T2=16.*T1
      IF(T2.LT.TMAX) GOTO 55
      T2=TMAX
      GOTO 60
 55   CALL LINDER(M,N,T2,DFDT2,ULAM,DSRCH,X,Y,UU,XLOW,XUPP, &
                 ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDT2.LE.0.) GOTO 60
      T1=T2
      GOTO 50
!C
 60   CONTINUE
!C     Nu sker en upprepad krympning av intervallet T1,T2.
      SQT1=DSQRT(T1)
      SQT2=DSQRT(T2)
 62   ITT3=ITT3+1
      IF(ITT3.GT.10) GOTO 90
      TM=SQT1*SQT2
      CALL LINDER(M,N,TM,DFDTM,ULAM,DSRCH,X,Y,UU,XLOW,XUPP, &
                 ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDTM.GT.0.) GOTO 65
      T2=TM
      TKVOT=T1/T2
      IF(TKVOT.GT.0.97) GOTO 90
!C     Linjesokningen klar. T1 ar approx optimal. IHITMX=0.
      SQT2=DSQRT(T2)
      GOTO 62
 65   T1=TM
      TKVOT=T1/T2
      IF(TKVOT.GT.0.97) GOTO 90
!C     Linjesokningen klar. T1 ar approx optimal. IHITMX=0.
      SQT1=DSQRT(T1)
      GOTO 62
!C
 80   TOPT=TMAX
      IHITMX=1
      GOTO 100
 90   TOPT=T1
      IHITMX=0
      IF(ITT1.GT.13) WRITE(*,911)
      IF(ITT2.GT.10) WRITE(*,912)
      IF(ITT3.GT.10) WRITE(*,913)
 911  FORMAT(' ITT1 GT 13 in LINSE')
 912  FORMAT(' ITT2 GT 10 in LINSE')
 913  FORMAT(' ITT3 GT 10 in LINSE')
 100  CONTINUE
!C
      RETURN
      END subroutine LINSE
!C
!C********+*********+*********+*********+*********+*********+*********+
!C********+*********+*********+*********+*********+*********+*********+
!C
      SUBROUTINE XUPDAT(N,ITER,XMMA,XVAL,XOLD1,XOLD2)
!C
!C       Version "December 2006".
!C    !-----------------------------------------!
!C    !  The author of this subroutine is       !
!C    !  Krister Svanberg <krille@math.kth.se>  !
!C    !-----------------------------------------!
!C
!C     XUPDAT updates XVAL, XOLD1 and if possible XOLD2.
!C
!C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
!C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
!C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
!C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
!C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
!C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
!C    IYFREE(M)
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Integer ::J,N,ITER
      !DIMENSION XMMA(N),XVAL(N),XOLD1(N),XOLD2(N)
      double precision:: XMMA(N),XVAL(N),XOLD1(N),XOLD2(N)
!C
      
      DO 70 J=1,N
      IF(ITER.GE.2) XOLD2(J)=XOLD1(J)
      XOLD1(J)=XVAL(J)
      XVAL(J)=XMMA(J)
 70   CONTINUE
!C
      RETURN
      END subroutine XUPDAT
!C
!C********+*********+*********+*********+*********+*********+*********+
!C
      SUBROUTINE XYZLAM(M,N,X,Y,Z,ULAM,XLOW,XUPP,ALFA,BETA, &
                       A,B,C,P,Q,P0,Q0,IYFREE)
      implicit none
!C
!C       Version "December 2006".
!C    !-----------------------------------------!
!C    !  The author of this subroutine is       !
!C    !  Krister Svanberg <krille@math.kth.se>  !
!C    !-----------------------------------------!
!C
!C     XYZLAM calculates the X,Y,Z that minimize the Lagrange
!C     function, given the vector ULAM of Lagrange multipliers.
!C
!C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
!C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
!C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
!C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
!C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
!C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
!C    IYFREE(M)
      
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER ::IYFREE(1), M,N,I,J,MJ1
   !   DIMENSION X(1),Y(1),ULAM(1),XLOW(1),XUPP(1),ALFA(1),BETA(1), &
   !            A(1),B(1),C(1),P(1),Q(1),P0(1),Q0(1)
      double precision:: Y(m),ULAM(M),XLOW(N),XUPP(N),ALFA(N),BETA(N)
      double precision:: A(M),B(M),C(M),P(M*N),Q(M*N),P0(N),Q0(N)      
      double precision:: PJ,QJ,PIJ,QIJ,SRPJ,SRQJ,XJ,UA,YI,UA1,Z,X(N*m)   
      
     ! write(*,*) 'XYZLAM 1',n,m,m*n
!C
      x(2)=x(1)
      !write(*,*) 'xyzlam 1.1'

    !  write(*,*) x(1)
      DO 10 I=1,M
        IF(IYFREE(I).EQ.0) GOTO 10
        IF(ULAM(I).LT.0.) ULAM(I)=0.
 10   CONTINUE
!C
    !  write(*,*) 'XYZLAM 2'
     ! write(*,*) J
    !  write(*,*) N
      DO 30 J=1,N
       ! write(*,*) shape(P0)
       ! write(*,*) J
       ! write(*,*) PJ

        PJ=P0(J)
      !  write(*,*) 'XYZLAM 2.1'
        QJ=Q0(J)
       ! write(*,*) 'XYZLAM 2.2'
        MJ1=M*(J-1)
      !  write(*,*) 'XYZLAM 2.3'
!C
     !   write(*,*) 'XYZLAM 3'
        DO 20 I=1,M
          IF(IYFREE(I).EQ.0) GOTO 20
          PIJ=P(MJ1+I)
          QIJ=Q(MJ1+I)
          PJ=PJ+ULAM(I)*PIJ
          QJ=QJ+ULAM(I)*QIJ
      !   write(*,*) 'XYZLAM 4'
 20     CONTINUE
!C
      ! write(*,*) 'XYZLAM 4.5'
        SRPJ=DSQRT(PJ)
        SRQJ=DSQRT(QJ)
      !  write(*,*) 'XYZLAM 4.6'
        XJ=(SRPJ*XLOW(J)+SRQJ*XUPP(J))/(SRPJ+SRQJ)
        IF(XJ.LT.ALFA(J)) XJ=ALFA(J)
        IF(XJ.GT.BETA(J)) XJ=BETA(J)
       ! write(*,*) 'XYZLAM 4.7'
       ! write(*,*) J,'J'
       ! write(*,*) shape(X),'shape(X)'
      !  write(*,*) XJ,'xj'
      
     
 

        X(J)=XJ               !!!!!!!
       ! write(*,*) 'XYZLAM 5'
 30   CONTINUE
!C
      UA=0.
      DO 40 I=1,M
        Y(I)=0.
        IF(IYFREE(I).EQ.0) GOTO 40
        UA=UA+ULAM(I)*A(I)
        YI=ULAM(I)-C(I)
        IF(YI.GT.0.) Y(I)=YI
       ! write(*,*) 'XYZLAM 6'
 40   CONTINUE
!C
      Z=0.
      UA1=UA-1.
     ! write(*,*) 'XYZLAM 7'
      IF(UA1.GT.0.) Z=10.*UA1
!C
    !  write(*,*) 'XYZLAM 8'
      RETURN
      END subroutine XYZLAM
!C
!C********+*********+*********+*********+*********+*********+*********+
!C
      SUBROUTINE GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA, &
                      A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
!C
!C       Version "December 2006".
!C    !-----------------------------------------!
!C    !  The author of this subroutine is       !
!C    !  Krister Svanberg <krille@math.kth.se>  !
!C    !-----------------------------------------!
!C
!C     GRADI calculates the gradient GRADF of the dual
!C     objective function, given the vector ULAM of dual variables.
!C
!C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
!C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
!C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
!C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
!C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
!C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
!C    IYFREE(M)
      
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
     ! DIMENSION X(M*N),Y(M),ULAM(M),XLOW(N),XUPP(N),ALFA(N),BETA(N), &
      !         A(M),B(M),C(M),P(M*N),Q(M*N),P0(N),Q0(N),GRADF(M)
      INTEGER ::I,M,J,N,MJ1
      INTEGER ::IYFREE(M)
      double precision ::X(M*N),Y(M),ULAM(M),XLOW(N),XUPP(N),ALFA(N),BETA(N)
      double precision ::A(M),B(M),C(M),P(M*N),Q(M*N),P0(N),Q0(N),GRADF(M)
      double precision ::XJLJ,UJXJ,QIJ,PIJ,Z
      
!C
      !write(*,*) 'gradi 1'
      !write(*,*) shape(x),m*n, 'shape x gradi,m*n'
      !write(*,*) x(1), 'x(1) gradi'
      CALL XYZLAM(M,N,X,Y,Z,ULAM,XLOW,XUPP,ALFA,BETA, &
                 A,B,C,P,Q,P0,Q0,IYFREE)
!C
      !write(*,*) 'gradi 2'
      DO 10 I=1,M
      !write(*,*) 'gradi2.1'
      GRADF(I)=-B(I)-Y(I)-A(I)*Z
 10   CONTINUE
!C
     ! write(*,*) 'gradi2.2'
      DO 30 J=1,N
      MJ1=M*(J-1)
      UJXJ=XUPP(J)-X(J)
      XJLJ=X(J)-XLOW(J)
!C
      !write(*,*) 'gradi2.3'
      DO 20 I=1,M
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      GRADF(I)=GRADF(I)+PIJ/UJXJ+QIJ/XJLJ
     ! write(*,*) 'gradi2.4'
 20   CONTINUE
 30   CONTINUE
!C
    !  write(*,*) 'gradi2.5'
      RETURN
      END subroutine GRADI
!C
!C********+*********+*********+*********+*********+*********+*********+
!C
      SUBROUTINE LINDER(M,N,T,DFDT,ULAM,DSRCH,X,Y,UU,XLOW,XUPP, &
                       ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
!C
!C       Version "December 2006".
!C    !-----------------------------------------!
!C    !  The author of this subroutine is       !
!C    !  Krister Svanberg <krille@math.kth.se>  !
!C    !-----------------------------------------!
!C
!C     LINDER calculates the scalar product DFDT of GRADF and DSRCH
!C     (= the directional derivative) at the point ULAM + T*DSRCH.
!C
!C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
!C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
!C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
!C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
!C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
!C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
!C    IYFREE(M)
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER ::M,N,I,J,MJ1
      INTEGER ::IYFREE(M)
      double precision:: ULAM(M),DSRCH(M),X(M*N),Y(M),UU(M)
      double precision:: XLOW(N),XUPP(N),ALFA(N),BETA(N)
      double precision:: A(M),B(M),C(M),P(M*N),Q(M*N),P0(N),Q0(N)
      double precision:: Z,UJXJ,XJLJ,PIJ,QIJ,T,DFDT
      
!C
      DO 10 I=1,M
      UU(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 10
      UU(I)=ULAM(I)+T*DSRCH(I)
      IF(UU(I).LT.0.) UU(I)=0.
 10   CONTINUE
!C
      CALL XYZLAM(M,N,X,Y,Z,UU,XLOW,XUPP,ALFA,BETA, &
                 A,B,C,P,Q,P0,Q0,IYFREE)
!C
      DO 20 I=1,M
      UU(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 20
      UU(I)=-B(I)-Y(I)-A(I)*Z
   20 CONTINUE
!C
      DO 40 J=1,N
      MJ1=M*(J-1)
      UJXJ=XUPP(J)-X(J)
      XJLJ=X(J)-XLOW(J)
      DO 30 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 30
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      UU(I)=UU(I)+PIJ/UJXJ+QIJ/XJLJ
   30 CONTINUE
   40 CONTINUE
!C
      DFDT=0.
      DO 50 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 50
      DFDT=DFDT+UU(I)*DSRCH(I)
   50 CONTINUE
!C
      RETURN
      END subroutine LINDER
!C
!C********+*********+*********+*********+*********+*********+*********+
!C
      SUBROUTINE HESSI(M,N,ULAM,HESSF,X,Y,ALFA,BETA, &
                      A,B,C,P,Q,P0,Q0,XLOW,XUPP,IYFREE)
!C
!C       Version "December 2006".
!C    !-----------------------------------------!
!C    !  The author of this subroutine is       !
!C    !  Krister Svanberg <krille@math.kth.se>  !
!C    !-----------------------------------------!
!C
!C   HESSI calculates HESSF = minus the reduced Hessian matrix of the
!C   dual objective function, given the vector ULAM of dual variables.
!C
!C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
!C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
!C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
!C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
!C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
!C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
!C    IYFREE(M)
     ! IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER ::M,N,IK,K,I,II,KK,J,MJ1
      INTEGER ::IYFREE(M)
      !DIMENSION ULAM(M),HESSF(M*(M+1)/2),X(M*N),Y(M),ALFA(N),BETA(N),A(M), &
     !          B(M),C(M),P(M*N),Q(M*N),P0(N),Q0(N),XLOW(N),XUPP(N)
      double precision::  ULAM(M),HESSF(M*(M+1)/2),X(M*N),Y(M),ALFA(N),BETA(N),A(M)
      double precision::  B(M),C(M),P(M*N),Q(M*N),P0(N),Q0(N),XLOW(N),XUPP(N)
      double precision:: ulamta,PJ,QJ,PIJ,QIJ,SRPJ,SRQJ,XJ,UJXJ,XJLJ,UJXJ2,XJLJ2,RR
      double precision:: PKJ,QKJ,TTK,TTI 
!C
      IK=0
      DO 12 K=1,M
      DO 11 I=K,M
      IK=IK+1
      HESSF(IK)=0.
 11   CONTINUE
 12   CONTINUE
!C
      ULAMTA=0.
      II=1
      DO 15 I=1,M
      IF(I.GT.1) II=II+M+2-I
      IF(IYFREE(I).EQ.0) GOTO 15
      IF(ULAM(I).LT.0.) ULAM(I)=0.
      IF(ULAM(I).GT.C(I)) HESSF(II)=1.
      ULAMTA=ULAMTA+ULAM(I)*A(I)
 15   CONTINUE
!C
      IF(ULAMTA.LE.1.) GOTO 40
!C
      KK=1
      DO 30 K=1,M
      IF(K.GT.1) KK=KK+M+2-K
      IF(IYFREE(K).EQ.0) GOTO 30
      DO 20 I=K,M
      IF(IYFREE(I).EQ.0) GOTO 20
      IK=KK+I-K
      HESSF(IK)=HESSF(IK)+10.*A(I)*A(K)
 20   CONTINUE
 30   CONTINUE
!C
 40   CONTINUE
!C
      DO 100 J=1,N
      PJ=P0(J)
      QJ=Q0(J)
      MJ1=M*(J-1)
      DO 50 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 50
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      PJ=PJ+ULAM(I)*PIJ
      QJ=QJ+ULAM(I)*QIJ
 50   CONTINUE
!C
      SRPJ=DSQRT(PJ)
      SRQJ=DSQRT(QJ)
      XJ=(SRPJ*XLOW(J)+SRQJ*XUPP(J))/(SRPJ+SRQJ)
      IF(XJ.GE.BETA(J)) GOTO 100
      IF(XJ.LE.ALFA(J)) GOTO 100
!C
      UJXJ=XUPP(J)-XJ
      XJLJ=XJ-XLOW(J)
      UJXJ2=UJXJ**2
      XJLJ2=XJLJ**2
      RR=2.*PJ/UJXJ**3+2.*QJ/XJLJ**3
!C
      KK=1
      DO 80 K=1,M
      IF(K.GT.1) KK=KK+M+2-K
      IF(IYFREE(K).EQ.0) GOTO 80
      PKJ=P(MJ1+K)
      QKJ=Q(MJ1+K)
      TTK=PKJ/UJXJ2-QKJ/XJLJ2
      DO 70 I=K,M
      IF(IYFREE(I).EQ.0) GOTO 70
      IK=KK+I-K
      PIJ=P(MJ1+I)
      QIJ=Q(MJ1+I)
      TTI=PIJ/UJXJ2-QIJ/XJLJ2
      HESSF(IK)=HESSF(IK)+TTI*TTK/RR
 70   CONTINUE
 80   CONTINUE
!C
 100  CONTINUE
!C
      RETURN
      END subroutine HESSI
!C
!C********+*********+*********+*********+*********+*********+*********+
!C
      SUBROUTINE LDLFAC(N,EPS,ADL,E)
!C
!C       Version "December 2006".
!C    !-----------------------------------------!
!C    !  The author of this subroutine is       !
!C    !  Krister Svanberg <krille@math.kth.se>  !
!C    !-----------------------------------------!
!C
!C    LDLFAC makes a factorization of a given symmetric matrix A.
!C    If A is positive definite, then A = L*D*LT.
!C    If A is not positive definite, then A + E = L*D*LT,
!C    where E is a positive semidefinite diagonal matrix such that
!C    A + E is positive definite.
!C    On entry, ADL defines the given matrix A.
!C    On leave, ADL defines the calculated matrices D and L.
!C
!C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
!C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
!C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
!C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
!C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
!C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
!C    IYFREE(M)
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      !DIMENSION ADL(1),E(1)
      Integer::N,JJ,J,KK,JK,L,K,IJ,I,IK
      double precision:: ADL(:),E(:),ADLIJ,EPS
!C
      JJ=1
!C
!write(*,*) 'EDL SIZE E SIZE'
      DO 100 J=1,N
        E(J)=0.
        IF(J.GT.1) JJ=JJ+N+2-J
        IF(J.EQ.1) GOTO 25
        write(*,*) JJ, 'JJ'
        pause
        KK=JJ
        JK=JJ
        DO 20 L=1,J-1
          K=J-L
          JK=JK-N+K
          KK=KK-N+K-1
          ADL(JJ)=ADL(JJ)-ADL(KK)*ADL(JK)*ADL(JK)
 20     CONTINUE
 25     IF(ADL(JJ).GE.EPS) GOTO 35
        E(J)=EPS-ADL(JJ)
        ADL(JJ)=EPS
 35     IF(J.EQ.N) GOTO 100
        IJ=JJ
        DO 50 I=J+1,N
          IJ=IJ+1
          ADLIJ=ADL(IJ)
          IF(J.EQ.1) GOTO 45
          IK=IJ
          JK=JJ
          KK=JJ
          DO 40 L=1,J-1
            K=J-L
            IK=IK-N+K
            JK=JK-N+K
            KK=KK-N+K-1
            ADLIJ=ADLIJ-ADL(KK)*ADL(IK)*ADL(JK)
 40       CONTINUE
 45       ADL(IJ)=ADLIJ/ADL(JJ)
 50     CONTINUE
 100  CONTINUE
!C
 !     write(*,*) JJ, 'JJ I ADL'
 !     pause
      RETURN
      END subroutine LDLFAC
!C
!C********+*********+*********+*********+*********+*********+*********+
!C
      SUBROUTINE LDLSOL(N,B,DL,X)
!C
!C       Version "December 2006".
!C    !-----------------------------------------!
!C    !  The author of this subroutine is       !
!C    !  Krister Svanberg <krille@math.kth.se>  !
!C    !-----------------------------------------!
!C
!C     LDLSOL solves a system of linear equations: A*X = B,
!C     where A has already been factorized as L*D*Ltranspose.
!C
!C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
!C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
!C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
!C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
!C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
!C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
!C    IYFREE(M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(1),DL(1),X(1)
      Integer ::N,JJ,JK,J,L,K,KJ

!C
      JJ=1
!C
      DO 30 J=1,N
      X(J)=B(J)
      IF(J.EQ.1) GOTO 30
      JJ=JJ+N+2-J
      JK=JJ
      DO 20 L=1,J-1
      K=J-L
      JK=JK-N+K
      X(J)=X(J)-DL(JK)*X(K)
 20   CONTINUE
 30   CONTINUE
!C
      JJ=1
      DO 40 J=1,N
      IF(J.GT.1) JJ=JJ+N+2-J
      X(J)=X(J)/DL(JJ)
 40   CONTINUE
!C
      DO 60 L=1,N-1
      J=N-L
      JJ=JJ-N+J-1
      KJ=JJ
      DO 50 K=J+1,N
      KJ=KJ+1
      X(J)=X(J)-DL(KJ)*X(K)
 50   CONTINUE
 60   CONTINUE
!C
      RETURN
      END subroutine LDLSOL
!C
!C********+*********+*********+*********+*********+*********+*********+

	
end module optim5
