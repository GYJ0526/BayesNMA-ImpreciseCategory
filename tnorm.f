      REAL*8 FUNCTION YTUVN(A,B,LA,LB,iseed)
C  This function generates a N(0,1) random variable
C  subject to the constraint that it be in an interval
C  (a,b), where the endpoints may be finite or infinite.
C  References:  J. Geweke, 1991.  "Efficient Simulation from the Multivariate
C               Normal and Student-t Distributions Subject to LInear
C               Constraints," in E.M. Keramidas (ed.), Computing Science and
C               Statistics: Proceedings of the 23rd Symposium on the
C               Interface, pp. 571-578.  Fairfax, VA: Interface Foundation
C               of North America, Inc.
C  Update history:  Develped in gibbs/tnorm.  Brought into ylib 1/23/92
C                   Modified to allow singleton 2/20/92
C                   Modified to prevent overflow before 10, 3/31/92
C  Inputs:
C  A, B    Endpoints of interval; A < B if LA = LB = .FALSE.
C  LA      .TRUE. if left endpoint is - infinity; in this
C          case A is ignored.
C  LB      .TRUE. if right endpoint is + infinity; in this
C          case B is ignored.
C  Output:
C  GGRNRM  Random variable
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL LA,LB,LFLIP
      DATA EPS,T1,T2,T3,T4/2.0D0,.375D0,2.18D0,.725D0,.45D0/
      F(X)=DEXP(-.5D0*X**2)
      IF(LA.AND.LB)GO TO 160
      LFLIP=.FALSE.
      IF(LA.OR.LB)GO TO 100
      IF(B.LE.A)GO TO 170
C ******* Finite interval
      C1=A
      C2=B
      IF((C1*C2).GT.0.0D0)GO TO 30
C ++++ (A,B) includes 0
      IF((C1.GT.-T1).AND.(C2.LT.T1))GO TO 20
C -- F(A) or F(B) small: full normal with rejection
   10 X=DRNNOF()
      IF(X.LT.C1)GO TO 10
      IF(X.GT.C2)GO TO 10
      GO TO 150
C -- F(A) and F(B) large: uniform importance sampling
   20 CDEL=C2-C1
   25 X=C1+CDEL*DRNUNF()
      IF(DRNUNF().GT.F(X))GO TO 25
      GO TO 150
C ++++ (A,B) excludes 0
C -- Transform to both positive
   30 IF(C1.GT.0.0D0)GO TO 40
      C=C1
      C1=-C2
      C2=-C
      LFLIP=.TRUE.
   40 F1=F(C1)
      F2=F(C2)
      IF(F2.LT.EPS)GO TO 60
      IF((F1/F2).GT.T2)GO TO 60
C  -- F(A)/F(B) not large: uniform importance sampling
   50 CDEL=C2-C1
   55 X=C1+CDEL*DRNUNF()
      IF(DRNUNF().GT.(F(X)/F1))GO TO 55
      GO TO 140
   60 IF(C1.GT.T3)GO TO 80
C -- P(X>A) and F(A)/F(B) large: half-normal with rejection
   70 X=DABS(DRNNOF())
      IF(X.LT.C1)GO TO 70
      IF(X.GT.C2)GO TO 70
      GO TO 140
C -- P(X>A) small, F(A)/F(B) large: exponential importance
C    sampling with rejection
   80 C=C2-C1
   90 Z=-DLOG(DRNUNF())/C1
      IF(Z.GT.C)GO TO 90
      IF(DRNUNF().GT.F(Z))GO TO 90
      X=C1+Z
      GO TO 140
C ****** Half-line interval
  100 C1=A
C -- Transform to bound from below if A = -infinity
      IF(LB)GO TO 110
      C1=-B
      LFLIP=.TRUE.
  110 IF(C1.GT.T4)GO TO 130
C -- A not large: full normal with rejection
  120 X=DRNNOF()
      IF(X.LT.C1)GO TO 120
      GO TO 140
C -- A small: exponential importance sampling
  130 Z=-DLOG(DRNUNF())/C1
      IF(DRNUNF().GT.F(Z))GO TO 130
      X=C1+Z
  140 IF(LFLIP)X=-X
  150 YTUVN=X
      RETURN
C ****** Whole interval
  160 YTUVN=DRNNOF()
      RETURN
C  ***** Singleton
  170 YTUVN=A
      RETURN
      END
