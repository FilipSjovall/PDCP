C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE MAXIM(M,N,GEPS,IYFREE,GRADF,DSRCH,HESSF,X,Y,Z,ULAM,
     1                 UU,XLOW,XUPP,ALFA,BETA,A,B,C,P,Q,P0,Q0)
C
C       Version "December 2006".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     MAXIM solves the dual MMA subproblem.
C     The dual variables are ulam(i), i=1,..,m,
C     which are required to be non-negative.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GRADF(1),DSRCH(1),HESSF(1),X(1),Y(1),ULAM(1),
     1          UU(1),XLOW(1),XUPP(1),ALFA(1),BETA(1),
     2          A(1),B(1),C(1),P(1),Q(1),P0(1),Q0(1)
      INTEGER IYFREE(1)
C
      ITR=0
      M3=3*M+30
C
      DO 10 I=1,M
      ULAM(I)=0.
      IYFREE(I)=1
 10   CONTINUE
C
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA,
     1           A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
C
      GMX=0.
      DO 20 I=1,M
      IYFREE(I)=0
      IF(GRADF(I).GT.GEPS) IYFREE(I)=1
      IF(GRADF(I).GT.GMX) GMX=GRADF(I)
 20   CONTINUE
C
      IF(GMX.LE.GEPS) GOTO 100
C     Vi avbryter optimeringen, ulam=0 ar optimal losning.
C
 30   CONTINUE
      ITR=ITR+1
      IF(ITR.GT.M3) GOTO 100
C     Vi avbryter optimeringen pga for manga subspa-anrop.
C
      CALL SUBSPA(ITR,M,N,GEPS,F,IYFREE,GRADF,DSRCH,HESSF,
     1            X,Y,ULAM,UU,XLOW,XUPP,ALFA,BETA,A,B,C,
     2            P,Q,P0,Q0,IHITY)
C
      IF(IHITY.EQ.0) GOTO 40
C     Om ihity = 0 sa ar ulam optimal pa aktuellt underrum.
C     Om ihity > 0 sa har vi slagit i ett nytt bivillkor.
      IYFREE(IHITY)=0
      ULAM(IHITY)=0.
      GOTO 30
C
 40   CONTINUE
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA,
     1           A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
C
      GMX=0.
      IGMX=0
      DO 50 I=1,M
      IF(IYFREE(I).EQ.1) GOTO 50
      IF(GRADF(I).LE.GMX) GOTO 50
      GMX=GRADF(I)
      IGMX=I
 50   CONTINUE
C
      IF(GMX.LE.GEPS) GOTO 100
C     Om gmx =< geps sa ar ulam optimal losning.
C     Om gmx > geps sa tar vi bort kravet att ulam(igmx)=0.
      IYFREE(IGMX)=1
      GOTO 30
C
 100  CONTINUE
C     Nu ar antingen ulam optimal losning eller itr>m3.
      CALL XYZLAM(M,N,X,Y,Z,ULAM,XLOW,XUPP,ALFA,BETA,
     1            A,B,C,P,Q,P0,Q0,IYFREE)
      IF(ITR.GT.M3) WRITE(*,911)
 911  FORMAT(' ITR GT M3 IN MAXIM')
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE SUBSPA(ITR,M,N,GEPS,F,IYFREE,GRADF,DSRCH,HESSF,
     1                  X,Y,ULAM,UU,XLOW,XUPP,ALFA,BETA,
     2                  A,B,C,P,Q,P0,Q0,IHITY)
C
C       Version "December 2006".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C    SUBSPA maximizes the dual objective function on the subspace
C    defined by ulam(i) = 0 for every i such that iyfree(i) = 0.
C    The first three iterations a steepest ascent method is used,
C    and after that a Newton method is used.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GRADF(1),DSRCH(1),HESSF(1),X(1),Y(1),
     1          ULAM(1),UU(1),XLOW(1),XUPP(1),ALFA(1),BETA(1),
     2          A(1),B(1),C(1),P(1),Q(1),P0(1),Q0(1)
      INTEGER IYFREE(1)
C
      IHITY=0
      ITESUB=0
      NYDIM=0
      DSRTOL=-0.0000001*GEPS
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA,
     1           A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
C
      DO 10 I=1,M
      DSRCH(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 10
      NYDIM=NYDIM+1
      DSRCH(I)=GRADF(I)
 10   CONTINUE
C
      IF(NYDIM.EQ.0) GOTO 100
C     Vi avbryter med ihity = 0, ty inga variabler ulam(i) ar fria.
      ITEMAX=50+5*NYDIM
C
 15   ITESUB=ITESUB+1
C     Har startar en ny iteration.
C
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
C     Vi avbryter med ihity = 0, ty dsrch ar ej en ascentriktning.
C
 25   CONTINUE
      CALL LINSE(M,N,ITESUB,IHITMX,IYFREE,TMAX,TOPT,ULAM,DSRCH,
     1           X,Y,UU,XLOW,XUPP,ALFA,BETA,A,B,C,P,Q,P0,Q0)
C
      DO 30 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 30
      ULAM(I)=ULAM(I)+TOPT*DSRCH(I)
      IF(ULAM(I).LT.0.) ULAM(I)=0.
 30   CONTINUE
C
      CALL GRADI(M,N,X,Y,ULAM,XLOW,XUPP,ALFA,BETA,
     1           A,B,C,P,Q,P0,Q0,GRADF,IYFREE)
C
      IF(IHITMX.EQ.1.AND.IHITY.GT.0) GOTO 100
C     Vi avbryter med ihity > 0, ty vi har slagit i det tidigare
C     inaktiva bivillkoret ulam(ihity) >= 0.
      IHITY=0
      IOPT=1
      DO 40 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 40
      IF(DABS(GRADF(I)).GT.GEPS) IOPT=0
      !IF(DABS(GRADF(I)).GT.GEPS)  then
      !	write(*,*) 'Tolerance, I: ', DABS(GRADF(I)), I
 		!endif
 40   CONTINUE
C
      IF(IOPT.EQ.1) GOTO 100
C     Vi avbryter med ihity = 0, ty optimal losning hittad.
      IF(ITESUB.GT.ITEMAX) GOTO 97
C     Vi avbryter med ihity = 0, ty for manga iterationer.
      IF(ITESUB.GE.3) GOTO 55
C     Om itesub>=3 sa byter vi fran steepest ascent till Newton.
C     Om itesub=<2 sa fortsatter vi med steepest ascent.
      DO 50 I=1,M
      DSRCH(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 50
      DSRCH(I)=GRADF(I)
 50   CONTINUE
      GOTO 15
C
 55   CONTINUE
C
      CALL HESSI(M,N,ULAM,HESSF,X,Y,ALFA,BETA,
     1           A,B,C,P,Q,P0,Q0,XLOW,XUPP,IYFREE)
C
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
C
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
C
      HESMOD=0.0001*HTRACE/ZZZZ
      IF(HESMOD.LT.GEPS) HESMOD=GEPS
      IKRED=0
      DO 77 K=1,NYDIM
      DO 76 I=K,NYDIM
      IKRED=IKRED+1
      IF(I.EQ.K) HESSF(IKRED)=HESSF(IKRED)+HESMOD
 76   CONTINUE
 77   CONTINUE
C
      CALL LDLFAC(NYDIM,GEPS,HESSF,UU)
C
      IRED=0
      DO 79 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 79
      IRED=IRED+1
      UU(IRED)=GRADF(I)
 79   CONTINUE
C
      CALL LDLSOL(NYDIM,UU,HESSF,DSRCH)
C
      DO 80 I=1,M
      UU(I)=DSRCH(I)
 80   CONTINUE
C
      IRED=0
      DO 85 I=1,M
      DSRCH(I)=0.
      IF(IYFREE(I).EQ.0) GOTO 85
      IRED=IRED+1
      DSRCH(I)=UU(IRED)
 85   CONTINUE
C
      GOTO 15
C
 97   CONTINUE
      WRITE(*,911)
C
 100  CONTINUE
C
      DO 110 I=1,M
      IF(IYFREE(I).EQ.0) GOTO 110
      IF(ULAM(I).LT.0.) ULAM(I)=0.
 110  CONTINUE
C
 911  FORMAT(' ITESUB GT ITEMAX IN SUBSPA')
 912  FORMAT(' GTD LE 0 IN SUBSPA')
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
C
      SUBROUTINE LINSE(M,N,ITESUB,IHITMX,IYFREE,TMAX,TOPT,ULAM,DSRCH,
     1                 X,Y,UU,XLOW,XUPP,ALFA,BETA,A,B,C,P,Q,P0,Q0)
C
C       Version "December 2006".
C    !-----------------------------------------!
C    !  The author of this subroutine is       !
C    !  Krister Svanberg <krille@math.kth.se>  !
C    !-----------------------------------------!
C
C     LINSE makes an approximate line search (maximization) in the
C     direction DSRCH from the point ULAM.
C     Main input:  ULAM, DSRCH, TMAX.
C     Main output: TOPT, IHITMX.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ULAM(1),DSRCH(1),X(1),Y(1),UU(1),XLOW(1),XUPP(1),
     1          ALFA(1),BETA(1),A(1),B(1),C(1),P(1),Q(1),P0(1),Q0(1)
      INTEGER IYFREE(1)
C
      ITT1=0
      ITT2=0
      ITT3=0
C
      CALL LINDER(M,N,TMAX,DFDTMX,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDTMX.GE.0.) GOTO 80
C     Linjesokningen klar. Optimalt steg ar T=TMAX. IHITMX=1.
      IF(TMAX.GT.1.) GOTO 40
      T2=TMAX
C
 30   CONTINUE
C     Nu sker en upprepad minskning av steget.
      ITT1=ITT1+1
      IF(ITT1.GT.13) GOTO 90
      T1=T2/2.
      IF(ITESUB.LE.3) T1=T2/16.
      CALL LINDER(M,N,T1,DFDT1,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDT1.GT.0.) GOTO 60
      T2=T1
      GOTO 30
C
 40   CONTINUE
C     Nu testas enhetssteget, dvs T=1.
      T1=1.
      T2=T1
      CALL LINDER(M,N,T1,DFDT1,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(ITESUB.GE.6.AND.DFDT1.GE.0.) GOTO 90
C     Linjesokningen klar. Enhetssteget duger. T=1, IHITMX=0.
      IF(ITESUB.LE.5.AND.DFDT1.GT.0.) GOTO 50
C     Enhetssteget ar for kort.
      GOTO 30
C     Enhetssteget ar for langt.
C
 50   ITT2=ITT2+1
C     Nu sker en upprepad okning av steget.
      IF(ITT2.GT.10) GOTO 90
      T2=2.*T1
      IF(ITESUB.LE.3) T2=16.*T1
      IF(T2.LT.TMAX) GOTO 55
      T2=TMAX
      GOTO 60
 55   CALL LINDER(M,N,T2,DFDT2,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDT2.LE.0.) GOTO 60
      T1=T2
      GOTO 50
C
 60   CONTINUE
C     Nu sker en upprepad krympning av intervallet T1,T2.
      SQT1=DSQRT(T1)
      SQT2=DSQRT(T2)
 62   ITT3=ITT3+1
      IF(ITT3.GT.10) GOTO 90
      TM=SQT1*SQT2
      CALL LINDER(M,N,TM,DFDTM,ULAM,DSRCH,X,Y,UU,XLOW,XUPP,
     1            ALFA,BETA,A,B,C,P,Q,P0,Q0,IYFREE)
      IF(DFDTM.GT.0.) GOTO 65
      T2=TM
      TKVOT=T1/T2
      IF(TKVOT.GT.0.97) GOTO 90
C     Linjesokningen klar. T1 ar approx optimal. IHITMX=0.
      SQT2=DSQRT(T2)
      GOTO 62
 65   T1=TM
      TKVOT=T1/T2
      IF(TKVOT.GT.0.97) GOTO 90
C     Linjesokningen klar. T1 ar approx optimal. IHITMX=0.
      SQT1=DSQRT(T1)
      GOTO 62
C
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
C
      RETURN
      END
C
C********+*********+*********+*********+*********+*********+*********+
