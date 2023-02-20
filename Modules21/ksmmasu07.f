C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE MMASUB(ITER,M,N,GEPS,IYFREE,XVAL,XMMA,
C      SUBROUTINE MMASUB(ITER,M,N,GEPS,GMOVE,IYFREE,XVAL,XMMA,
     1                  XMIN,XMAX,XOLD1,XOLD2,XLOW,XUPP,
     2                  ALFA,BETA,A,B,C,Y,Z,ULAM,
     3                  F0VAL,FVAL,FMAX,DF0DX,DFDX,
     4                  P,Q,P0,Q0,UU,GRADF,DSRCH,HESSF)
C
C       Version "December 2006".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C    Use of this code is for academic purposes only,
C    regulated by an agreement with Krister Svanberg.
C    The code is not to be redistributed.
C
C    MMASUB generates and solves the MMA subproblem,
C    which is of the following form in the variables
C    x_1,...,x_N, y_1,...,y_M, and z.
C
C   minimize h_0(x) + r_0 + z + 0.05*z^2 + sum{c_i*y_i + 0.5*(y_i)^2}
C
C subject to h_i(x) - a_i*z - y_i <= b_i ,     i=1,..,M
C                   alfa_j <= x_j <= beta_j ,  j=1,..,N
C                             y_i >= 0 ,       i=1,..,M
C                               z >= 0 .
C
C    with h_i(x) = sum{p_ij/(xupp_j-x_j) + q_ij/(x_j-xlow_j)}.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION XVAL(1),XMMA(1),XMIN(1),XMAX(1),
     1          XOLD1(1),XOLD2(1),XLOW(1),XUPP(1),
     2          ALFA(1),BETA(1),A(1),B(1),C(1),Y(1),ULAM(1),
     3          FVAL(1),FMAX(1),DF0DX(1),DFDX(1),
     4          P(1),Q(1),P0(1),Q0(1),UU(1),
C     5          GRADF(1),DSRCH(1),HESSF(1),GMOVE(1)
     5          GRADF(1),DSRCH(1),HESSF(1)
      INTEGER IYFREE(1)
C
C********+*********+*********+*********+*********+*********+*********+
C  The sizes of the above areas must be at least as follows:
C
C    XVAL(N),XMMA(N),XMIN(N),XMAX(N),
C    XOLD1(N),XOLD2(N),XLOW(N),XUPP(N),
C    ALFA(N),BETA(N),A(M),B(M),C(M),Y(M),ULAM(M),
C    FVAL(M),FMAX(M),DF0DX(N),DFDX(M*N),
C    P(M*N),Q(M*N),P0(N),Q0(N),UU(M),
C    GRADF(M),DSRCH(M),HESSF(M*(M+1)/2),
C    IYFREE(M)
C********+*********+*********+*********+*********+*********+*********+
C
C***  Input to the subroutine MMASUB:
C
C  ITER  = Current iteration number ( =1 the first iteration).
C     N  = Number of variables x_j in the problem.
C     M  = Number of constraints in the problem (not including
C          the simple upper and lower bounds on the variables).
C  GEPS  = Tolerance parameter for the constraints.
C          (Used in the termination criteria for the subproblem.)
C   XVAL(j) = Current value of the variable x_j.
C   XMIN(j) = Original lower bound for the variable x_j.
C   XMAX(j) = Original upper bound for the variable x_j.
C  XOLD1(j) = Value of the variable x_j one iteration ago.
C  XOLD2(j) = Value of the variable x_j two iterations ago.
C   XLOW(j) = Current value of the lower asymptot l_j.
C   XUPP(j) = Current value of the upper asymptot u_j.
C      A(i) = Coefficient a_i for the minimax variable z.
C      C(i) = Coefficient c_i for the artificial variable y_i.
C    F0VAL  = Value of the objective function f_0(x)
C   FVAL(i) = Value of the i:th constraint function f_i(x).
C   FMAX(i) = Right hand side of the i:th constraint.
C  DF0DX(j) = Derivative of f_0(x) with respect to x_j.
C   DFDX(k) = Derivative of f_i(x) with respect to x_j,
C             where k = (j-1)*M + i.
C
C*** Output from the subroutine MMASUB:
C
C   XMMA(j) = Optimal value of x_j in the MMA subproblem.
C      Y(i) = Optimal value of the "artificial" variable y_i.
C      Z    = Optimal value of the "minimax" variable z.
C   ULAM(i) = Optimal value of the dual variable lambda_i.
C   XLOW(j) = New value on the lower asymptot l_j.
C   XUPP(j) = New value on the upper asymptot u_j.
C
C*** Working areas and their usage in MMASUB:
C
C   ALFA(j) = Lower bound for x_j in the MMA subproblem.
C   BETA(j) = Upper bound for x_j in the MMA subproblem.
C      P(k) = Coefficient p_ij in the MMA subproblem,
C             where k = (j-1)*M + i.
C      Q(k) = Coefficient q_ij in the MMA subproblem,
C             where k = (j-1)*M + i.
C     P0(j) = Coefficient p_0j in the MMA subproblem.
C     Q0(j) = Coefficient q_0j in the MMA subproblem.
C      B(i) = Right hand side b_i in the MMA subproblem.
C  GRADF(i) = Gradient component of the dual objective function.
C  DSRCH(i) = Search direction component in the dual subproblem.
C  HESSF(k) = Hessian matrix component of the dual function.
C     UU(i) = Component in a working area.
C IYFREE(i) = 0 for dual variables which are fixed to zero in
C               the current subspace of the dual subproblem,
C           = 1 for dual variables which are "free" in
C               the current subspace of the dual subproblem.
C
C********+*********+*********+*********+*********+*********+*********+
C
      CALL ASYMPT(ITER,M,N,XVAL,XMIN,XMAX,XOLD1,XOLD2,
C      CALL ASYMPT(ITER,M,N,GMOVE,XVAL,XMIN,XMAX,XOLD1,XOLD2,
     1            XLOW,XUPP,ALFA,BETA)
C
C****  ASYMPT calculates the asymptotes XLOW(j) and XUPP(j),
C****  and the bounds ALFA(j) and BETA(j).
C
      CALL GENSUB(M,N,XVAL,XMIN,XMAX,F0VAL,DF0DX,FMAX,FVAL,
     1            DFDX,P,Q,B,P0,Q0,R0,XLOW,XUPP)
C
C***** GENSUB generates the MMA subproblem by calculating the
C***** coefficients P(i,j),Q(i,j),B(i),P0(j),Q0(j) and R0.
C     
      CALL MAXIM(M,N,GEPS,IYFREE,GRADF,DSRCH,HESSF,XMMA,Y,Z,
     1           ULAM,UU,XLOW,XUPP,ALFA,BETA,A,B,C,P,Q,P0,Q0)
C
C***** MAXIM solves the dual problem of the MMA subproblem.
C***** ULAM = optimal solution of this dual problem.
C***** XMMA,Y,Z = optimal solution of the MMA subproblem.
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE ASYMPT(ITER,M,N,XVAL,XMIN,XMAX,XOLD1,XOLD2,
C      SUBROUTINE ASYMPT(ITER,M,N,GMOVE,XVAL,XMIN,XMAX,XOLD1,XOLD2,
     1                  XLOW,XUPP,ALFA,BETA)
C
C       Version "December 2006".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     ASYMPT calculates the asymptotes XLOW and XUPP,
C     and the bounds ALFA and BETA, for the current subproblem.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XVAL(1),XMIN(1),XMAX(1),XOLD1(1),XOLD2(1),
C     1          XLOW(1),XUPP(1),ALFA(1),BETA(1),GMOVE(1)
     1          XLOW(1),XUPP(1),ALFA(1),BETA(1)
C
      ALBEFA=0.1d0
      
      GHDECR=0.7d0
      GHINCR=1.2d0
      GHINIT=0.1d0
      
C     Mathias
C      GHINIT=0.5d0
C      GHDECR=0.95d0
C      GHINCR=1.0d0
      IF(ITER.GE.3) GOTO 350
C
C***  Here ITER = 1 or 2 .
      DO 200 J=1,N
      XMMJ=XMAX(J)-XMIN(J)
C     IF(XMMJ.LT.0.00001) XMMJ=0.00001
      IF(XMMJ.LT.0.00001d0) XMMJ=0.00001d0
      XLOW(J)=XVAL(J)-GHINIT*XMMJ
      XUPP(J)=XVAL(J)+GHINIT*XMMJ
  200 CONTINUE
      GOTO 500
C
C***  Here ITER is greater than 2.
  350 CONTINUE
C
      DO 400 J=1,N
      XTEST=(XVAL(J)-XOLD1(J))*(XOLD1(J)-XOLD2(J))
C     FAK=1.0
      FAK=1d0
C     IF(XTEST.LT.0.) FAK=GHDECR
C     IF(XTEST.GT.0.) FAK=GHINCR
      IF(XTEST.LT.(0d0)) FAK=GHDECR
      IF(XTEST.GT.(0d0)) FAK=GHINCR
      XLOW(J)=XVAL(J)-FAK*(XOLD1(J)-XLOW(J))
      XUPP(J)=XVAL(J)+FAK*(XUPP(J)-XOLD1(J))
      XMMJ=XMAX(J)-XMIN(J)
C     IF(XMMJ.LT.0.00001) XMMJ=0.00001
      IF(XMMJ.LT.(0.00001d0)) XMMJ=0.00001d0
C     GMINJ = XVAL(J)-10.0*XMMJ
C     GMAXJ = XVAL(J)-0.01*XMMJ
C     HMINJ = XVAL(J)+0.01*XMMJ
C     HMAXJ = XVAL(J)+10.0*XMMJ
      GMINJ = XVAL(J)- 10d0 * XMMJ
      GMAXJ = XVAL(J)- 0.01d0*XMMJ
      HMINJ = XVAL(J)+ 0.01d0*XMMJ
      HMAXJ = XVAL(J)+ 10d0 * XMMJ
      IF(XLOW(J).LT.GMINJ) XLOW(J)=GMINJ
      IF(XLOW(J).GT.GMAXJ) XLOW(J)=GMAXJ
      IF(XUPP(J).LT.HMINJ) XUPP(J)=HMINJ
      IF(XUPP(J).GT.HMAXJ) XUPP(J)=HMAXJ
  400 CONTINUE
C
  500 CONTINUE
C
      DO 600 J=1,N
C     XMIJ=XMIN(J)-0.000001
C     XMAJ=XMAX(J)+0.000001
      XMIJ=XMIN(J) - 0.000001d0
      XMAJ=XMAX(J) + 0.000001d0
      IF(XVAL(J).GE.XMIJ) GOTO 550
C     XLOW(J)=XVAL(J)-(XMAJ-XVAL(J))/0.9
C     XUPP(J)=XVAL(J)+(XMAJ-XVAL(J))/0.9
      XLOW(J)=XVAL(J)-(XMAJ-XVAL(J))/0.9d0
      XUPP(J)=XVAL(J)+(XMAJ-XVAL(J))/0.9d0
      GOTO 600
  550 CONTINUE
      IF(XVAL(J).LE.XMAJ) GOTO 600
C     XLOW(J)=XVAL(J)-(XVAL(J)-XMIJ)/0.9
C     XUPP(J)=XVAL(J)+(XVAL(J)-XMIJ)/0.9
      XLOW(J)=XVAL(J)-(XVAL(J)-XMIJ)/0.9d0
      XUPP(J)=XVAL(J)+(XVAL(J)-XMIJ)/0.9d0
  600 CONTINUE
C
      DO 700 J=1,N
C     EXTERNAL MOVE LIMIT
C      IF(XLOW(j).lt.(xval(j)-GMOVE(j))) then
C        xlow(j) = xval(j)-GMOVE(j)
C      ENDIF
C      IF(XUPP(j).gt.(xval(j)+GMOVE(j))) then
C        XUPP(j) = xval(j)+GMOVE(j)
C      ENDIF
      ALFA(J)=XLOW(J)+ALBEFA*(XVAL(J)-XLOW(J))
      BETA(J)=XUPP(J)-ALBEFA*(XUPP(J)-XVAL(J))
      IF(ALFA(J).LT.XMIN(J)) ALFA(J)=XMIN(J)
      IF(BETA(J).GT.XMAX(J)) BETA(J)=XMAX(J)      
  700 CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE GENSUB(M,N,XVAL,XMIN,XMAX,F0VAL,DF0DX,FMAX,FVAL,
     1                  DFDX,P,Q,B,P0,Q0,R0,XLOW,XUPP)
C
C       Version "December 2006".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     GENSUB calculates P( ),Q( ),B( ),P0( ),Q0( ) and R0
C     for the current subproblem.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XVAL(1),XMIN(1),XMAX(1),DF0DX(1),FMAX(1),FVAL(1),
     1          DFDX(1),P(1),Q(1),B(1),P0(1),Q0(1),XLOW(1),XUPP(1)
C
C     RAA0=0.00001
      RAA0=0.00001d0
      R0=F0VAL
      DO 20 I=1,M
      B(I)=FMAX(I)-FVAL(I)
   20 CONTINUE
C
      DO 50 J=1,N
      MJ=M*(J-1)
      UJLJ=XUPP(J)-XLOW(J)
      UJXJ=XUPP(J)-XVAL(J)
      XJLJ=XVAL(J)-XLOW(J)
      UJXJ2=UJXJ*UJXJ
      XJLJ2=XJLJ*XJLJ
      XMMJ=XMAX(J)-XMIN(J)
C     IF(XMMJ.LT.0.00001) XMMJ=0.00001
      IF(XMMJ.LT.(0.00001d0)) XMMJ=0.00001d0
      P0J=RAA0/XMMJ
      Q0J=RAA0/XMMJ
C     IF(DF0DX(J).GT.0.) P0J=P0J+1.001*DF0DX(J)
C     IF(DF0DX(J).GT.0.) Q0J=Q0J+0.001*DF0DX(J)
C     IF(DF0DX(J).LT.0.) Q0J=Q0J-1.001*DF0DX(J)
C     IF(DF0DX(J).LT.0.) P0J=P0J-0.001*DF0DX(J)
      IF(DF0DX(J).GT.(0d0)) P0J=P0J+1.001d0*DF0DX(J)
      IF(DF0DX(J).GT.(0d0)) Q0J=Q0J+0.001d0*DF0DX(J)
      IF(DF0DX(J).LT.(0d0)) Q0J=Q0J-1.001d0*DF0DX(J)
      IF(DF0DX(J).LT.(0d0)) P0J=P0J-0.001d0*DF0DX(J)
      P0J=P0J*UJXJ2
      Q0J=Q0J*XJLJ2
      P0(J)=P0J
      Q0(J)=Q0J
      R0=R0-P0J/UJXJ-Q0J/XJLJ
C
      DO 40 I=1,M
      IJ=MJ+I
      PIJ=RAA0/XMMJ
      QIJ=RAA0/XMMJ
      DFIJ=DFDX(IJ)
C     IF(DFIJ.GT.0.) PIJ=PIJ+1.001*DFIJ
C     IF(DFIJ.GT.0.) QIJ=QIJ+0.001*DFIJ
C     IF(DFIJ.LT.0.) QIJ=QIJ-1.001*DFIJ
C     IF(DFIJ.LT.0.) PIJ=PIJ-0.001*DFIJ
      IF(DFIJ.GT.(0d0)) PIJ=PIJ+1.001d0*DFIJ
      IF(DFIJ.GT.(0d0)) QIJ=QIJ+0.001d0*DFIJ
      IF(DFIJ.LT.(0d0)) QIJ=QIJ-1.001d0*DFIJ
      IF(DFIJ.LT.(0d0)) PIJ=PIJ-0.001d0*DFIJ
      PIJ=PIJ*UJXJ2
      QIJ=QIJ*XJLJ2
      P(IJ)=PIJ
      Q(IJ)=QIJ
      B(I)=B(I)+PIJ/UJXJ+QIJ/XJLJ
C
   40 CONTINUE
   50 CONTINUE
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
