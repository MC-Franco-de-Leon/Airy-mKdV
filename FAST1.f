ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   THIS HAS INTEGRATION ON [0,1]
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      SUBROUTINE FAST(B, N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z),INTEGER(I-N)
      DIMENSION B(2)
      COMMON /CONS/ PII, P7, P7TWO, C22, S22, PI2
      PII = 4.D0 * DATAN(1.D0)
      PI8 = PII / 8.D0
      P7 = 1.D0 / DSQRT(2.D0)
      P7TWO = 2.D0 * P7
      C22 = DCOS(PI8)
      S22 = DSIN(PI8)
      PI2 = 2.D0 * PII
      DO 10 I = 1,15
          M = I
          NT = 2**I
          IF ( N .EQ. NT ) GO TO 20
10    CONTINUE
      WRITE(6,9999)
9999  FORMAT(33H N IS NOT A POWER OF TWO FOR FAST)
      STOP
20    N4POW = M / 2
      IF ( M - N4POW * 2 ) 40, 40, 30
30    NN = 2
      INT = N / NN
      CALL FR2TR(INT, B(1), B(INT+1))
      GO TO 50
40    NN = 1
50    IF ( N4POW .EQ. 0 ) GO TO 70
      DO 60 IT = 1,N4POW
          NN = NN * 4
          INT = N / NN
          CALL FR4TR(INT,NN,B(1),B(INT+1),
     *         B(2*INT+1),B(3*INT+1),
     *         B(1),B(INT+1),B(2*INT+1),B(3*INT+1))
60    CONTINUE
70    CALL FORD1(M, B)
      CALL FORD2(M, B)
      T = B(2)
      B(2) = 0.D0
      B(N+1) = T
      B(N+2) = 0.D0
      DO 80 IT = 4,N,2
          B(IT) = -B(IT)
80    CONTINUE
      RETURN
      END
      SUBROUTINE FSST(B, N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z),INTEGER (I-N)
      DOUBLE PRECISION B(2)
      COMMON /CONST/ PII, P7, P7TWO, C22, S22, PI2
      PII = 4.D0 *  DATAN(1.D0)
      PI8 = PII / 8.D0
      P7 = 1.D0 / DSQRT(2.D0)
      P7TWO = 2.D0 * P7
      C22 = DCOS(PI8)
      S22 = DSIN(PI8)
      PI2 = 2.D0 * PII
      DO 10 I = 1,15
          M = I
          NT = 2**I
          IF ( N .EQ. NT ) GO TO 20
10    CONTINUE
      WRITE(6,9999)
9999  FORMAT(33H N IS NOT A POWER OF TWO FOR FSST)
      STOP
20    B(2) = B(N+1)
      DO 30 I =  4,N,2
          B(I) = -B(I)
30    CONTINUE
      DO 40 I = 1,N
          B(I) = B(I) / DFLOAT(N)
40    CONTINUE
      N4POW = M / 2
      CALL FORD2(M, B)
      CALL FORD1(M, B)
C
      IF ( N4POW .EQ. 0 ) GO TO 60
      NN = 4 * N
      DO 50 IT = 1,N4POW
          NN = NN / 4
          INT = N / NN
          CALL FR4SYN(INT, NN, B(1), B(INT+1), B(2*INT+1), B(3*INT+1),
     *                B(1), B(INT+1), B(2*INT+1), B(3*INT+1))
50    CONTINUE
60    IF ( M - N4POW * 2 ) 80, 80, 70
70    INT = N / 2
      CALL FR2TR(INT, B(1), B(INT+1))
80    RETURN
      END
      SUBROUTINE FR2TR(INT, B0, B1)
      DOUBLE PRECISION B0(2), T, B1(2)
      DO 10 K = 1,INT
          T = B0(K) + B1(K)
          B1(K) = B0(K) - B1(K)
          B0(K) = T
10    CONTINUE
      RETURN
      END
      SUBROUTINE FR4TR(INT, NN, B0, B1, B2, B3, B4, B5, B6, B7)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      INTEGER L(15)
      DOUBLE PRECISION B0(2),B1(2),B2(2),B3(2),B4(2),B5(2),B6(2),B7(2)
      COMMON /CONS/ PII, P7, P7TWO, C22, S22, PI2
      EQUIVALENCE (L15,L(1)),(L14,L(2)),(L13,L(3)),(L12,L(4)),
     *  (L11,L(5)),(L10,L(6)),(L9,L(7)),(L8,L(8)),(L7,L(9)),
     *  (L6,L(10)),(L5,L(11)),(L4,L(12)),(L3,L(13)),(L2,L(14)),
     *  (L1,L(15))
      L(1) = NN / 4
      DO 40 K = 2,15
          IF ( L(K-1) - 2 ) 10, 20, 30
10        L(K-1) = 2
20        L(K) =  2
          GO TO 40
30        L(K) = L(K-1) / 2
40    CONTINUE
C
      PIOVN = PII / DFLOAT(NN)
      JI = 3
      JL = 2
      JR = 2
C
      DO  120 J1 = 2,L1,2
      DO  120 J2 = J1,L2,L1
      DO  120 J3 = J2,L3,L2
      DO  120 J4 = J3,L4,L3
      DO  120 J5 = J4,L5,L4
      DO  120 J6 = J5,L6,L5
      DO  120 J7 = J6,L7,L6
      DO  120 J8 = J7,L8,L7
      DO  120 J9 = J8,L9,L8
      DO  120 J10 = J9,L10,L9
      DO  120 J11 = J10,L11,L10
      DO  120 J12 = J11,L12,L11
      DO  120 J13 = J12,L13,L12
      DO  120 J14 = J13,L14,L13
      DO  120 JTHET = J14,L15,L14
           TH2 = JTHET - 2
           IF ( TH2 ) 50, 50, 90
50         DO 60 K = 1,INT
               T0 = B0(K) + B2(K)
               T1 = B1(K) + B3(K)
               B2(K) = B0(K) - B2(K)
               B3(K) = B1(K) - B3(K)
               B0(K) = T0 + T1
               B1(K) = T0 - T1
60         CONTINUE
C
           IF ( NN - 4 ) 120, 120, 70
70         K0 = INT * 4 + 1
           KL = K0 + INT - 1
           DO 80 K = K0,KL
               PR = P7 * (B1(K) - B3(K))
               PI = P7 * (B1(K) + B3(K))
               B3(K) = B2(K) + PI
               B1(K) = PI - B2(K)
               B2(K) = B0(K) - PR
               B0(K) = B0(K) + PR
80         CONTINUE
           GO TO 120
C
90         ARG = TH2 * PIOVN
           C1 = DCOS(ARG)
           S1 = DSIN(ARG)
           C2 = C1**2 - S1**2
           S2 = C1*S1 + C1*S1
           C3 = C1*C2 - S1*S2
           S3 = C2*S1 + S2*C1
C
           INT4 = INT * 4
           J0 = JR * INT4 + 1
           K0 = JI * INT4 + 1
           JLAST = J0 + INT - 1
           DO 100 J = J0,JLAST
              K = K0 + J - J0
              R1 = B1(J)*C1 - B5(K)*S1
              R5 = B1(J)*S1 + B5(K)*C1
              T2 = B2(J)*C2 - B6(K)*S2
              T6 = B2(J)*S2 + B6(K)*C2
              T3 = B3(J)*C3 - B7(K)*S3
              T7 = B3(J)*S3 + B7(K)*C3
              T0 = B0(J) + T2
              T4 = B4(K) + T6
              T2 = B0(J) - T2
              T6 = B4(K) - T6
              T1 = R1 + T3
              T5 = R5 + T7
              T3 = R1 - T3
              T7 = R5 - T7
              B0(J) = T0 + T1
              B7(K) = T4 + T5
              B6(K) = T0 - T1
              B1(J) = T5 - T4
              B2(J) = T2 - T7
              B5(K) = T6 + T3
              B4(K) = T2 + T7
              B3(J) = T3 - T6
100       CONTINUE
C
          JR = JR + 2
          JI = JI - 2
          IF ( JI - JL ) 110, 110, 120
110       JI = 2 * JR - 1
          JL = JR
120   CONTINUE
      RETURN
      END
      SUBROUTINE FORD1(M, B)
      DOUBLE PRECISION B(2)
      DOUBLE PRECISION T
C
      K = 4
      KL = 2
      N = 2 ** M
      DO 40 J = 4,N,2
          IF ( K - J ) 20, 20, 10
10        T = B(J)
          B(J) = B(K)
          B(K) = T
20        K = K - 2
          IF ( K - KL ) 30, 30, 40
30        K = 2 * J
          KL = J
40    CONTINUE
      RETURN
      END
      SUBROUTINE FORD2(M, B)
      INTEGER L(15)
      DOUBLE PRECISION  B(2), T
      EQUIVALENCE (L15,L(1)),(L14,L(2)),(L13,L(3)),(L12,L(4)),
     *  (L11,L(5)),(L10,L(6)),(L9,L(7)),(L8,L(8)),(L7,L(9)),
     *  (L6,L(10)),(L5,L(11)),(L4,L(12)),(L3,L(13)),(L2,L(14)),
     *  (L1,L(15))
      N = 2 ** M
      L(1) = N
      DO 10 K = 2,M
          L(K) = L(K-1) / 2
10    CONTINUE
      DO 20 K = M,14
          L(K+1) = 2
20    CONTINUE
      IJ = 2
      DO  40 J1 = 2,L1,2
      DO  40 J2 = J1,L2,L1
      DO  40 J3 = J2,L3,L2
      DO  40 J4 = J3,L4,L3
      DO  40 J5 = J4,L5,L4
      DO  40 J6 = J5,L6,L5
      DO  40 J7 = J6,L7,L6
      DO  40 J8 = J7,L8,L7
      DO  40 J9 = J8,L9,L8
      DO  40 J10 = J9,L10,L9
      DO  40 J11 = J10,L11,L10
      DO  40 J12 = J11,L12,L11
      DO  40 J13 = J12,L13,L12
      DO  40 J14 = J13,L14,L13
      DO  40 JI = J14,L15,L14
           IF ( IJ - JI ) 30, 40, 40
30         T = B(IJ-1)
           B(IJ-1) = B(JI-1)
           B(JI-1) = T
           T = B(IJ)
           B(IJ) = B(JI)
           B(JI) = T
40         IJ = IJ + 2
       RETURN
       END
      SUBROUTINE FR4SYN(INT, NN, B0, B1, B2, B3, B4, B5, B6, B7)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z),INTEGER (I-N)
      INTEGER L(15)
      DOUBLE PRECISION B0(2), B1(2), B2(2), B3(2), B4(2), B5(2), B6(2),
     *    B7(2)
      COMMON /CONST/ PII, P7, P7TWO, C22, S22, PI2
      EQUIVALENCE (L15,L(1)),(L14,L(2)),(L13,L(3)),(L12,L(4)),
     *  (L11,L(5)),(L10,L(6)),(L9,L(7)),(L8,L(8)),(L7,L(9)),
     *  (L6,L(10)),(L5,L(11)),(L4,L(12)),(L3,L(13)),(L2,L(14)),
     *  (L1,L(15))
C
      L(1) = NN / 4
      DO 40 K = 2,15
          IF ( L(K-1) - 2 ) 10, 20, 30
10        L(K-1) = 2
20        L(K) = 2
          GO TO 40
30        L(K) = L(K-1)/2
40    CONTINUE
C
      PIOVN = PII / DFLOAT(NN)
      JI = 3
      JL = 2
      JR = 2
C
      DO  120 J1 = 2,L1,2
      DO  120 J2 = J1,L2,L1
      DO  120 J3 = J2,L3,L2
      DO  120 J4 = J3,L4,L3
      DO  120 J5 = J4,L5,L4
      DO  120 J6 = J5,L6,L5
      DO  120 J7 = J6,L7,L6
      DO  120 J8 = J7,L8,L7
      DO  120 J9 = J8,L9,L8
      DO  120 J10 = J9,L10,L9
      DO  120 J11 = J10,L11,L10
      DO  120 J12 = J11,L12,L11
      DO  120 J13 = J12,L13,L12
      DO  120 J14 = J13,L14,L13
      DO  120 JTHET = J14,L15,L14
         TH2 = JTHET - 2
         IF ( TH2 ) 50, 50, 90
50       DO 60 K = 1,INT
             T0 = B0(K) + B1(K)
             T1 = B0(K) - B1(K)
             T2 = B2(K) * 2.D0
             T3 = B3(K) * 2.D0
             B0(K) = T0 + T2
             B2(K) = T0 - T2
             B1(K) = T1 + T3
             B3(K) = T1 - T3
60       CONTINUE
C
         IF ( NN - 4 ) 120, 120, 70
70       K0 = INT * 4 + 1
         KL = K0 + INT - 1
         DO 80 K = K0,KL
             T2 = B0(K) - B2(K)
             T3 = B1(K) + B3(K)
             B0(K) = (B0(K) + B2(K)) * 2.D0
             B2(K) = (B3(K) - B1(K)) * 2.D0
             B1(K) = (T2 + T3) * P7TWO
             B3(K) = (T3 - T2) * P7TWO
80       CONTINUE
         GO TO 120
90       ARG = TH2 * PIOVN
         C1 = DCOS(ARG)
         S1 = -DSIN(ARG)
         C2 = C1**2 - S1**2
         S2 = C1*S1 + C1*S1
         C3 = C1*C2 - S1*S2
         S3 = C2*S1 + S2*C1
C
         INT4 = INT * 4
         J0 = JR * INT4 + 1
         K0 = JI * INT4 + 1
         JLAST = J0 + INT - 1
         DO 100 J = J0,JLAST
             K = K0 + J - J0
             T0 = B0(J) + B6(K)
             T1 = B7(K) - B1(J)
             T2 = B0(J) - B6(K)
             T3 = B7(K) + B1(J)
             T4 = B2(J) + B4(K)
             T5 = B5(K) - B3(J)
             T6 = B5(K) + B3(J)
             T7 = B4(K) - B2(J)
             B0(J) = T0 + T4
             B4(K) = T1 + T5
             B1(J) = (T2+T6)*C1 - (T3+T7)*S1
             B5(K) = (T2+T6)*S1 + (T3+T7)*C1
             B2(J) = (T0-T4)*C2 - (T1-T5)*S2
             B6(K) = (T0-T4)*S2 + (T1-T5)*C2
             B3(J) = (T2-T6)*C3 - (T3-T7)*S3
             B7(K) = (T2-T6)*S3 + (T3-T7)*C3
100      CONTINUE
         JR = JR + 2
         JI = JI - 2
         IF ( JI - JL ) 110, 110, 120
110      JI = 2 * JR - 1
         JL = JR
120   CONTINUE
      RETURN
      END

C
C------------------------------------------------------------------------------.C SUBROUTINE: FFT842
C FAST FOURIER TRANSFORM FOR N=2**M
C COMPLEX INPUT
C-------------------------------------------------------------------------------C
      SUBROUTINE FFT842(IN,N,X,Y)
C
C THIS ROGRAM REPLACES THE VECTOR Z=X+IY BY ITS FINITE
C DISCRETE, COMPLEX FOURIER TRANSFORM IF IN=0.  THE INVERSE TRANSFORM
C IS CALCULATED FOR IN=1.  IT PERFORMS AS MANY BASE 
C 8 ITERATIONS AS POSSIBLE AND THEN FINISHES WITH A BASE 4 ITERATION
C OR A BASE 2 ITERATION IF NEEDED.
C
C THE SUBROUTINE IS CALLED AS SUBROUTINE FFT842(IN,N,X,Y).
C THE INTEGER N (A POWER OF 2), THE N REAL LOCATION ARRAY X,
C AND THE N REAL LOCATION ARRAY Y MUST BE SUPPLIES TO THE SUBROUTINE.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z),INTEGER(I-N)
      INTEGER L(15)
      DOUBLE PRECISION X(2), Y(2)
      COMMON /CON2/ PI2, P7
      EQUIVALENCE (L15,L(1)), (L14,L(2)), (L13,L(3)), (L12,L(4)),
     *    (L11,L(5)), (L10,L(6)), (L9,L(7)), (L8,L(8)), (L7,L(9)),
     *    (L6,L(10)), (L5,L(11)), (L4,L(12)), (L3,L(13)), (L2,L(14)),
     *    (L1,L(15))
      PI2 = 8.D0*DATAN(1.D0)
      P7 = 1.D0/DSQRT(2.D0)
      DO 10 I=1,15
        M = I
        NT = 2**I
        IF (N.EQ.NT) GO TO 20
 10   CONTINUE
      WRITE(6,9999)
 9999 FORMAT (35H N IS NOT A POWER OF TWO FOR FFT842)
      STOP
 20   N2POW = M
      NTHPO = N
      FN = NTHPO
      IF (IN.EQ.1) GO TO 40
      DO 30 I=1,NTHPO
        Y(I) = -Y(I)
 30   CONTINUE
 40   N8POW = N2POW/3
      IF (N8POW.EQ.0) GO TO 60
C
C RADIX 8 PASSES, IF ANY.
C
      DO 50 IPASS=1,N8POW
        NXTLT = 2**(N2POW-3*IPASS)
        LENGT = 8*NXTLT
        CALL R8TX(NXTLT,NTHPO,LENGT,X(1),X(NXTLT+1),X(2*NXTLT+1),
     *      X(3*NXTLT+1),X(4*NXTLT+1),X(5*NXTLT+1),X(6*NXTLT+1),
     *      X(7*NXTLT+1),Y(1),Y(NXTLT+1),Y(2*NXTLT+1),Y(3*NXTLT+1),
     *      Y(4*NXTLT+1),Y(5*NXTLT+1),Y(6*NXTLT+1),Y(7*NXTLT+1))
 50   CONTINUE
C
C IS THERE A FOUR FACTOR LEFT
C
 60   IF (N2POW-3*N8POW-1) 90,70,80
C
C GO THROUGH THE BASE 2 ITERATION
C
C
 70   CALL R2TX(NTHPO,X(1),X(2),Y(1),Y(2))
      GO TO 90
C
C GO THROUGH THE BASE 4 ITERATION
C
 80   CALL R4TX(NTHPO,X(1),X(2),X(3),
     &      X(4),Y(1),Y(2),Y(3),Y(4))
C
 90   DO 110 J=1,15
        L(J) = 1
        IF (J-N2POW) 100,100,110
 100    L(J) = 2**(N2POW+1-J)
 110  CONTINUE
      IJ = 1
      DO 130 J1=1,L1
      DO 130 J2=J1,L2,L1
      DO 130 J3=J2,L3,L2
      DO 130 J4=J3,L4,L3
      DO 130 J5=J4,L5,L4
      DO 130 J6=J5,L6,L5
      DO 130 J7=J6,L7,L6
      DO 130 J8=J7,L8,L7
      DO 130 J9=J8,L9,L8
      DO 130 J10=J9,L10,L9
      DO 130 J11=J10,L11,L10
      DO 130 J12=J11,L12,L11
      DO 130 J13=J12,L13,L12
      DO 130 J14=J13,L14,L13
      DO 130 JI=J14,L15,L14
         IF (IJ-JI) 120, 130, 130
 120     R = X(IJ)
         X(IJ) = X(JI)
         X(JI) = R
         FI = Y(IJ)
         Y(IJ) = Y(JI)
         Y(JI) = FI
 130     IJ = IJ +1
       IF (IN.EQ.1) GO TO 150
       DO 140 I=1,NTHPO
         Y(I) = -Y(I)
 140   CONTINUE
       GO TO 170
 150   DO 160 I=1,NTHPO
         X(I) = X(I)/FN
         Y(I) = Y(I)/FN
 160   CONTINUE
 170   RETURN
       END
C
C------------------------------------------------------------------------
C SUBROUTINE: R2TX
C RADIX 2 ITERATION SUBROUTINE
C-------------------------------------------------------------------------------C
      SUBROUTINE R2TX(NTHPO,CR0,CR1,CI0,CI1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z),INTEGER (I-N)
      DOUBLE PRECISION CR0(2),CR1(2),CI0(2),CI1(2)
      DO 10 K=1,NTHPO,2
         R1 = CR0(K) + CR1(K)
         CR1(K) = CR0(K) -CR1(K)
         CR0(K) = R1
         FI1 = CI0(K) + CI1(K)
         CI1(K) = CI0(K) - CI1(K)
         CI0(K) = FI1
 10   CONTINUE
      RETURN
      END
C
C-------------------------------------------------------------------
C SUBROUTINE: R4TX
C RADIX 4 ITERATION SUBROUTINE
C-------------------------------------------------------------------------------C
      SUBROUTINE R4TX(NTHPO,CR0,CR1,CR2,CR3,CI0,CI1,CI2,CI3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z),INTEGER (I-N)
      DOUBLE PRECISION CR0(2),CR1(2),CR2(2),CR3(2),CI0(2),CI1(2),CI2(2),
     *    CI3(2)
      DO 10 K=1,NTHPO,4
        R1 = CR0(K) + CR2(K)
        R2 = CR0(K) - CR2(K)
        R3 = CR1(K) + CR3(K)
        R4 = CR1(K) - CR3(K)
        FI1 = CI0(K) + CI2(K)
        FI2 = CI0(K) - CI2(K)
        FI3 = CI1(K) + CI3(K)
        FI4 = CI1(K) - CI3(K)
        CR0(K) = R1 + R3
        CI0(K) = FI1 + FI3
        CR1(K) = R1 - R3
        CI1(K) = FI1 - FI3
        CR2(K) = R2 - FI4
        CI2(K) = FI2 + R4
        CR3(K) = R2 + FI4
        CI3(K) = FI2 - R4
 10   CONTINUE
      RETURN
      END
C
C--------------------------------------------------------------------
C SUBROUTINE: R8TX
C RADIX 8 ITERATION SUBROUTINE
C-------------------------------------------------------------------------------C
      SUBROUTINE R8TX(NXTLT,NTHPO,LENGT,CR0,CR1,CR2,CR3,CR4,CR5,CR6,CR7,
     *    CI0,CI1,CI2,CI3,CI4,CI5,CI6,CI7)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z),INTEGER (I-N)
      DOUBLE PRECISION CR0(2),CR1(2),CR2(2),CR3(2),CR4(2),CR5(2),
     *    CR6(2),CR7(2),
     *    CI0(2),CI1(2),CI2(2),CI3(2),CI4(2),CI5(2),CI6(2),CI7(2)
      COMMON /CON2/ PI2, P7
C
      SCALE = PI2/DFLOAT(LENGT)
      DO 30 J=1,NXTLT
        ARG = DFLOAT(J-1)*SCALE
        C1 = DCOS(ARG)
        S1 = DSIN(ARG)
        C2 = C1**2 - S1**2
        S2 = C1*S1 + C1*S1
        C3 = C1*C2 - S1*S2
        S3 = C2*S1 + S2*C1
        C4 = C2**2 - S2**2
        S4 = C2*S2 + C2*S2
        C5 = C2*C3 - S2*S3
        S5 = C3*S2 + S3*C2
        C6 = C3**2 - S3**2
        S6 = C3*S3 + C3*S3
        C7 = C3*C4 - S3*S4
        S7 = C4*S3 + S4*C3
        DO 20 K=J,NTHPO,LENGT
          AR0 = CR0(K) + CR4(K)
          AR1 = CR1(K) + CR5(K)
          AR2 = CR2(K) + CR6(K)
          AR3 = CR3(K) + CR7(K)
          AR4 = CR0(K) - CR4(K)
          AR5 = CR1(K) - CR5(K)
          AR6 = CR2(K) - CR6(K)
          AR7 = CR3(K) - CR7(K)
          AI0 = CI0(K) + CI4(K)
          AI1 = CI1(K) + CI5(K)
          AI2 = CI2(K) + CI6(K)
          AI3 = CI3(K) + CI7(K)
          AI4 = CI0(K) - CI4(K)
          AI5 = CI1(K) - CI5(K)
          AI6 = CI2(K) - CI6(K)
          AI7 = CI3(K) - CI7(K)
          BR0 = AR0 + AR2
          BR1 = AR1 + AR3
          BR2 = AR0 - AR2
          BR3 = AR1 - AR3
          BR4 = AR4 - AI6
          BR5 = AR5 - AI7
          BR6 = AR4 + AI6
          BR7 = AR5 + AI7
          BI0 = AI0 + AI2
          BI1 = AI1 + AI3
          BI2 = AI0 - AI2
          BI3 = AI1 - AI3
          BI4 = AI4 + AR6
          BI5 = AI5 + AR7
          BI6 = AI4 - AR6
          BI7 = AI5 - AR7
          CR0(K) = BR0 + BR1
          CI0(K) = BI0 + BI1
          IF (J.LE.1) GO TO 10
          CR1(K) = C4*(BR0-BR1) - S4*(BI0-BI1)
          CI1(K) = C4*(BI0-BI1) + S4*(BR0-BR1)
          CR2(K) = C2*(BR2-BI3) - S2*(BI2+BR3)
          CI2(K) = C2*(BI2+BR3) + S2*(BR2-BI3)
          CR3(K) = C6*(BR2+BI3) - S6*(BI2-BR3)
          CI3(K) = C6*(BI2-BR3) + S6*(BR2+BI3)
          TR = P7*(BR5-BI5)
          TI = P7*(BR5+BI5)
          CR4(K) = C1*(BR4+TR) - S1*(BI4+TI)
          CI4(K) = C1*(BI4+TI) + S1*(BR4+TR)
          CR5(K) = C5*(BR4-TR) - S5*(BI4-TI)
          CI5(K) = C5*(BI4-TI) + S5*(BR4-TR)
          TR = -P7*(BR7+BI7)
          TI = P7*(BR7-BI7)
          CR6(K) = C3*(BR6+TR) - S3*(BI6+TI)
          CI6(K) = C3*(BI6+TI) + S3*(BR6+TR)
          CR7(K) = C7*(BR6-TR) - S7*(BI6-TI)
          CI7(K) = C7*(BI6-TI) + S7*(BR6-TR)
          GO TO 20
 10       CR1(K) = BR0 - BR1
          CI1(K) = BI0 - BI1
          CR2(K) = BR2 - BI3
          CI2(K) = BI2 + BR3
          CR3(K) = BR2 + BI3
          CI3(K) = BI2 - BR3
          TR = P7*(BR5-BI5)
          TI = P7*(BR5+BI5)
          CR4(K) = BR4 + TR
          CI4(K) = BI4 + TI
          CR5(K) = BR4 - TR
          CI5(K) = BI4 - TI
          TR = -P7*(BR7+BI7)
          TI = P7*(BR7-BI7)
          CR6(K) = BR6 + TR
          CI6(K) = BI6 + TI
          CR7(K) = BR6 - TR
          CI7(K) = BI6 - TI
 20     CONTINUE
 30   CONTINUE
      RETURN
      END
        subroutine fdiff(x,xi,n,h)
        implicit double precision (a-h,o-z)
        parameter(lda=33000)
        double precision x(1),xi(1)
        double precision a(lda),h
        pi=4.d0*datan(1.d0)
        p2=2.d0*pi
        n1=n-1
        n2=n1/2
        do 10 j=1,n1
        a(j)=x(j)
 10     continue
        call fast(a,n1)
        xi(1)=0.d0
        xi(2)=0.d0
        do 20 j=2,n2+1
        k=2*j-1
        xi(k)=-a(k+1)*dfloat(j-1)*dble(p2)
        xi(k+1)=a(k)*dfloat(j-1)*dble(p2)
 20     continue
        call fsst(xi,n1)
        xi(n)=xi(1)
        return
        end
                
        subroutine fdiff2(x,xi,n,h)
        implicit double precision (a-h,o-z)
        parameter(lda=33000)
        double precision x(0),xi(0)
        double precision a(lda),h
        pi=4.d0*datan(1.d0)
        p2=2.d0*pi
        n1=n-1
        n2=n1/2
        do 10 j=0,n1-1
        a(j)=x(j)
 10     continue
        a(n1)=a(0)
        call fast(a,n1)
        xi(0)=0.d0
        xi(1)=0.d0
        do 20 j=2,n2+1
        k=2*(j-1)
        xi(k)=-a(k+1)*dfloat(j-1)*dble(p2)
        xi(k+1)=a(k)*dfloat(j-1)*dble(p2)
 20     continue
        call fsst(xi,n1)

c        xi(n1)=xi(0)
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine filterfdiff(x,xi,n,h)
        implicit double precision (a-h,o-z)
        parameter(lda=33000)
        double precision x(1),xi(1)
        double precision a(lda),h,almbd
        pi=4.d0*datan(1.d0)
        p2=2.d0*pi
        n1=n-1
        n2=n1/2
        almbd=.5
        do 10 j=1,n1
        a(j)=x(j)
 10     continue
        call fast(a,n1)
        xi(1)=0.d0
        xi(2)=0.d0
        do 20 j=2,n2+1
        k=2*j-1 
          auxrho=dble(0)
          avar=dble(j-1)/n2
          if (avar .le. almbd) then
             auxrho=1
          else if ((avar .gt. almbd) .and. (avar .le. 1)) then
             aux3=1-((avar-almbd)/(1-almbd))**2
             aquotient=exp(-dfloat(1)/aux3)
             auxrho=exp(dfloat(1))*aquotient          
             
c              ax=1+((-avar+almbd)/(1-almbd))
c              auxrho=(35-84*ax+70*ax**2-20*ax**3)*ax**4
c              write(6,*)'normalizado ', auxrho
          end if
          
c          auxcoef1=(a(k+1)*dfloat(j-1)*dble(p2))**2
c          auxcoef2=(a(k)*dfloat(j-1)*dble(p2))**2
c          auxs1=dsqrt(auxcoef1+auxcoef2)
c          write(6,*)'spect norm before', auxs1

         xi(k)=-a(k+1)*dfloat(j-1)*dble(p2)*auxrho
         xi(k+1)=a(k)*dfloat(j-1)*dble(p2)*auxrho
        
c          auxs2=sqrt( (xi(k))**2+(xi(k+1))**2)
c          write(6,*)'spect norm after', auxs2
        
 20     continue
         xi(n1)=0
         xi(n1+1)=0
        call fsst(xi,n1)
        xi(n)=xi(1)
        return
        end

        SUBROUTINE SYNTHCV(Z,ZI,E,N,L)
        PARAMETER(LDA=33000)
        COMPLEX*16 Z(N),ZI(L)
        DOUBLE PRECISION E(L)
        COMPLEX*16 CXP1(LDA),CXP2(LDA),U1(LDA),U2(LDA)
        COMPLEX*16 ZJ1,ZJ2
        DOUBLE PRECISION A(LDA),B(LDA)
        N1=N-1
        N2=N1/2
        N21=N2-1
        N22=N2-2
        DO 1 K=1,N1
        A(K)=DREAL(Z(K))/DFLOAT(N1)
        B(K)=DIMAG(Z(K))/DFLOAT(N1)
1       CONTINUE
        CALL FFT842(0,N1,A,B)
        DO 5 K=1,L
        CXP1(K)=DCMPLX(DCOS(E(K)),DSIN(E(K)))
        CXP2(K)=DCONJG(CXP1(K))
        U1(K)=DCMPLX(A(N2),B(N2))
        U2(K)=DCMPLX(A(N2+2),B(N2+2))
 5      CONTINUE
        DO 10 J=1,N22
        ZJ1=DCMPLX(A(N2-J),B(N2-J))
        ZJ2=DCMPLX(A(N2+J+2),B(N2+J+2))
        DO 10 K=1,L
        U1(K)=ZJ1+CXP1(K)*U1(K)
        U2(K)=ZJ2+CXP2(K)*U2(K)
 10     CONTINUE
        DO 20 K=1,L
        ZI(K)=DCMPLX(A(1),B(1))+DCMPLX(A(N2+1),B(N2+1))*DCOS(N2*E(K))+
     *        CXP1(K)*U1(K)+CXP2(K)*U2(K)
 20     CONTINUE
        RETURN
        END
        subroutine fintg(x,xi,n,h)
        implicit double precision (a-h,o-z)
        parameter(lda=33000)
        double precision x(1),xi(1)
        double precision a(lda),h
        pi=4.d0*datan(1.d0)
        p2=2.d0*pi
        n1=n-1
        n2=n1/2
        do 10 j=1,n1
        a(j)=x(j)
 10     continue
        call fast(a,n1)
        xi(1)=0.d0
        xi(2)=0.d0
        do 20 j=2,n2+1
        k=2*j-1
        xi(k)=a(k+1)/(p2*dfloat(j-1))
        xi(k+1)=-a(k)/(p2*dfloat(j-1))
 20     continue
        call fsst(xi,n1)
        xi(n)=xi(1)
        do 30 j=1,n 
        xi(j)=(a(1)/dfloat(n1))*(j-1)*h+xi(j)-xi(n)
 30     continue
        return
        end
        SUBROUTINE SYNTH(X,XI,E,N,A,B)
        implicit double precision (a-h,o-z)
        double precision X(1),XI
        DOUBLE PRECISION E
        COMPLEX*16 CXP1,CXP2,U1,U2
        COMPLEX*16 ZJ1,ZJ2
        DOUBLE PRECISION A(1),B(1)
        N1=N-1
        N2=N1/2
        N21=N2-1
        N22=N2-2
C        DO 1 K=1,N1
C        A(K)=X(K)/DFLOAT(N1)
C        B(K)=0.D0
C1       CONTINUE
C        CALL FFT842(0,N1,A,B)
        CXP1=DCMPLX(DCOS(E),DSIN(E))
        CXP2=DCONJG(CXP1)
        U1=DCMPLX(A(N2),B(N2))
        U2=DCMPLX(A(N2+2),B(N2+2))
        DO 10 J=1,N22
        ZJ1=DCMPLX(A(N2-J),B(N2-J))
        ZJ2=DCMPLX(A(N2+J+2),B(N2+J+2))
        U1=ZJ1+CXP1*U1
        U2=ZJ2+CXP2*U2
 10     CONTINUE
        XI=A(1)+A(N2+1)*DCOS(N2*E)+DREAL(CXP1*U1+CXP2*U2)
        RETURN
        END


