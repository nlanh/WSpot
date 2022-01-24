      PROGRAM WSpot
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NMAX = 10000)
      CHARACTER NOTE 
      DIMENSION DRI(NMAX),DV0(NMAX),DEXPF(NMAX),DERDEXPF(NMAX),
     &          DVSO(NMAX),DVCOU(NMAX),DVNU(NMAX),DVCE(NMAX),DPOT(NMAX),
     &          WF(NMAX)
      DATA FSC,AMN,HC,RCOMP/137.036D0,931.494D0,197.327D0,1.4138/
      DATA PI/3.1415926535D0/ ,EPS/1D-3/
      E2 = HC/FSC
      HM = HC*HC/2/AMN
      OPEN(1, FILE='BOUND.INP', FORM='FORMATTED', STATUS='OLD')
      OPEN(2, FILE='BOUND.OUT', FORM='FORMATTED')
      OPEN(3, FILE='BOUND.POT', FORM='FORMATTED')
      OPEN(4, FILE='BOUND.WF', FORM='FORMATTED')
      
      PRINT *, 'THE CALCULATION FOR BOUND STATE'
      PRINT *, 'WRITTEN BY N.L.ANH AND B.M.LOC'
      WRITE(2,*) 'THE CALCULATION FOR BOUND STATE'
      WRITE(2,*) 'WRITTEN BY N.L.ANH AND B.M.LOC'
      READ(1,*) NOTE
      READ(1,*) N, H
      READ(1,*) NOTE
      READ(1,*) A1, A2, KZ1, KZ2
      N1 = INT(A1 + 0.5D0)
      N2 = INT(A2 + 0.5D0)
      LZ = 2 - MOD(N1+N2,2)
      LZ1 = 2 - MOD(N1,2)
      LZ2 = 2 - MOD(N2,2)
      AMA = A1 + A2
      RMU = A1*A2/AMA
      RM = RMU/HM
      ZE = KZ1*KZ2*E2            
      READ(1,*) NOTE
      READ(1,*) N0, LF, JF
      JF = JF*LZ2
      LS = 1*LZ2
      PRINT 2001, N1+N2, KZ1+KZ2, N1+N2-KZ1-KZ2
      WRITE(2,2001) N1+N2, KZ1+KZ2, N1+N2-KZ1-KZ2
      
      IF(KZ2.EQ.1) PRINT 2002, 'PROTON', N1, KZ1, N1-KZ1
      IF(KZ2.EQ.1) WRITE(2,2002) 'PROTON', N1, KZ1, N1-KZ1
      IF(KZ2.EQ.0) PRINT 2002, 'NEUTRON', N1, KZ1, N1-KZ1
      IF(KZ2.EQ.0) WRITE(2,2002) 'NEUTRON', N1, KZ1, N1-KZ1
 2001 FORMAT('THE NUCLEUS WITH A =', I3, '; Z =', I2, '; N =', I2) 
 2002 FORMAT(A7,' + CORE NUCLEUS WITH A =' I3, '; Z =',I2, '; N =', I2) 

! READ THE POTENTIAL
      PRINT*, 'THE POTENTIAL'
      WRITE(2,*)'THE POTENTIAL'
      READ(1,*) NOTE
      READ(1,*) KPOT
      IF (KPOT.EQ.0) THEN
       READ(1,*) NOTE
       READ(1,*) V0, V1, A0, R0, RC
       PRINT*, 'CENTRAL POTENTIAL IS OF WS FORM' 
       WRITE(2,*)'CENTRAL POTENTIAL IS OF WS FORM'
       PRINT 2003, V0, V1, R0, A0    
       WRITE(2,2003) V0, V1, R0, A0    
      END IF
      
      READ(1,*) NOTE
      READ(1,*) KSO
      IF (KSO.EQ.0) THEN
       READ(1,*) NOTE
       READ(1,*) VSO, ASO, RSO, QSO
       PRINT*, 'SPIN-ORBIT POTENTIAL IS OF DERIVATIVE OF WS'
       WRITE(2,*)'SPIN-ORBIT POTENTIAL IS OF DERIVATIVE OF WS'
       PRINT 2004, VSO, RSO, ASO 
       WRITE(2,2004) VSO, RSO, ASO  
      END IF
      
      IF (RC.GE.0.0000001) THEN
       PRINT*, 
     &      'COULOMB POTENTIAL IS USED FOR UNIFORM CHARGED DISTRIBUTION'
       WRITE(2,*)
     &      'COULOMB POTENTIAL IS USED FOR UNIFORM CHARGED DISTRIBUTION'
       PRINT 2005, RC
       WRITE(2,2005) RC
      END IF
! END READ THE POTENTIAL

 2003 FORMAT('WITH V0 =', F8.2, ' MeV;     V1 =', F8.2, ' MeV;   R0 = ', 
     &        F5.3, ' fm;   a0 =', F5.2, ' fm')
 2004 FORMAT('WITH VSO =', F6.2,' MeV;     RSO = ', F5.3,' fm;   aSO =', 
     &        F5.2, ' fm')
 2005 FORMAT('WITH THE COULOMB RADIUS RC =', F6.3, ' fm')
 
! COMPUTING THE FINAL POTENTIAL (POTF)
      UMAX = 0 
      ZF = (JF*(JF + 2) - LS*(LS + 2) - 4*LF*(LF + 1))/4.0 
      DO I = 1, N
       DRI(I) = I*H
      END DO
      
      IF (KPOT.EQ.0) THEN
       DO I = 1, N
        IF (KZ2.EQ.1) THEN
         DV0(I) = (V0+V1*(AMA-2*(KZ1+KZ2))/AMA)/(1+EXP((DRI(I)-R0)/A0)) 
        ELSE IF (KZ2.EQ.0) THEN
         DV0(I) = (V0-V1*(AMA-2*(KZ1+KZ2))/AMA)/(1+EXP((DRI(I)-R0)/A0))
        END IF
       END DO
      END IF
      
      IF (KSO.EQ.0) THEN
       DO I = 1, N
        DEXPF(I) = 1/(1+QSO*EXP((DRI(I) - RSO)/ASO))
       END DO
       CALL DER1(DEXPF, DERDEXPF, N, H)
       DO I = 1, N
        DVSO(I) = RCOMP**2*VSO*DERDEXPF(I)/DRI(I)
       END DO
      END IF
      
      IF (RC.GE.0.0000001) THEN
       DO I = 1, N
        IF (DRI(I).LE.RC) THEN
         DVCOU(I) = ZE*(3 - (DRI(I)/RC)**2)/2/RC
        ELSE
         DVCOU(I) = ZE/DRI(I)
        END IF
       END DO
      END IF
      
      DO I = 1, N
       DVNU(I) = DV0(I) + ZF*DVSO(I)
       DVCE(I) = LF*(LF+1)/DRI(I)**2
       DPOT(I) = (DVNU(I)+DVCOU(I))*RM + DVCE(I)
      END DO
      
C  PRINT POTENTIALS
      WRITE(3,2006)'R (fm)', 'V0+V1 (MeV)', 'VSO*LSI (MeV)', 
     &              'VNUCL (MeV)', 'VCOUL (MeV)', 'VTOTAL (MeV)', 
     &              'VCENTRI (MeV)', 'VEFF (MeV)'
      DO I = 1, N
       WRITE(3,2007) DRI(I), DV0(I), ZF*DVSO(I), DVNU(I), DVCOU(I), 
     &               DVNU(I)+DVCOU(I), DVCE(I)/RM, DPOT(I)/RM
      END DO
      

C COMPUTING THE FINAL ENERGY (EBOUND) AND WAVE FUNCTION (WFF)
      S2 = 0
      UMAX = 0.d0
      CALL NUM1L(N, H, E, S2, DPOT, WF, N0, EPS)
      EBOUND = (E + UMAX)/RM 
      IF(E.GE.0) THEN
       PRINT *, 'THERE IS NO BOUND STATE'
       WRITE(2,*) 'THERE IS NO BOUND STATE'
      ELSE
       PRINT 2008, 'BOUND-STATE ENERGY EB =', EBOUND, ' MEV'
       WRITE(2,2008) 'BOUND-STATE ENERGY EB =', EBOUND, ' MEV'
      END IF
      



! PRINT WAVE FUNCTION      
      WRITE(4,*) 'BOUND WAVE FUNCTION'
      WRITE(4,2009) 'R (fm)', 'PHI', 'PHI**2'
      DO I = 1, N
       WRITE(4,2010) DRI(I), WF(I), WF(I)**2
      END DO

 2006 FORMAT(A10, 7(A15))
 2007 FORMAT(F10.3, 7(F15.6))
 2008 FORMAT(A24, F8.3, A4)
 2009 FORMAT(A10, A15, A15)
 2010 FORMAT(F10.3, F15.6, F15.6)
 
      END PROGRAM
      
      
*COUFRA 
      SUBROUTINE COUFRA(RHO,ETA,MINL,MAXL,FC,FCP,GC,GCP)
C*** FONCTIONS COULOMBIENNES CALCULEES EN R = RHO PAR LA METHODE DES FRA
C*** CONTINUES DE STEED. MINL ET MAXL CORRESPONDENT AUX VRAIES VALEURS D
C*** VOIR BARNETT, FENG, STEED ET GOLDFARB, CPC 1974 *******************
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 K,K1,K2,K3,K4,M1,M2,M3,M4
      DIMENSION FC(MAXL),FCP(MAXL),GC(MAXL),GCP(MAXL)
      SAVE  
      DATA ACCUR,STEP/1.D-7,100.0D0/
      PACE = STEP 
      ACC = ACCUR 
      R = RHO 
      KTR = 1 
      LMAX = MAXL 
      LMIN1 = MINL+1  
      XLL1 = MINL*LMIN1 
      ETA2 = ETA**2 
      TURN = ETA+SQRT(ETA2+XLL1)
      IF(R.LT.TURN.AND.ABS(ETA).GE.1.D-6) KTR = -1  
      KTRP = KTR
      GO TO 2 
    1 R = TURN
      TF = F
      TFP = FP
      LMAX = MINL 
      KTRP = 1
    2 ETAR = ETA*R
   21 RHO2=R*R
      PL = LMAX+1 
      PMX = PL+0.5D0  
C** FRACTION CONTINUE POUR FP(MAXL)/F(MAXL) ; XL=F ; XLPRIME=FP ********
      FP = ETA/PL+PL/R
      DK = ETAR+ETAR  
      DEL = 0
      D = 0
      F = 1 
      K = (PL*PL-PL+ETAR)*(PL+PL-1) 
      IF(PL*PL+PL+ETAR.NE.0.) GO TO 3 
      R = R*1.0000001D0 
      GO TO 2 
    3 H = (PL*PL+ETA2)*(1-PL*PL)*RHO2 
      K = K+DK+PL*PL*6
      D = 1/(D*H+K) 
      DEL = DEL*(D*K-1) 
      IF(PL.LT.PMX) DEL = -R*(PL*PL+ETA2)*(PL+1)*D/PL 
      PL = PL+1 
      FP = FP+DEL 
      IF(D.LT.0) F = -F
      IF(PL.GT.20000.0D0) GO TO 11
      IF(ABS(DEL/FP).GE.ACC) GO TO 3
      FP = F*FP 
      IF(LMAX.EQ.MINL) GO TO 5  
      FC(LMAX+1) = F  
      FCP(LMAX+1) = FP
C*** RECURRENCE ARRIERE POUR F ET FP ; GC,GCP UTILISES POUR STOCKAGE ***
      L = LMAX
      DO 4 LP=LMIN1,LMAX
      PL = L
      GC(L+1) = ETA/PL+PL/R 
      GCP(L+1) = SQRT(ETA2+PL*PL)/PL
      FC(L) =(GC(L+1)*FC(L+1)+FCP(L+1))/GCP(L+1)
      FCP(L) = GC(L+1)*FC(L)-GCP(L+1)*FC(L+1) 
    4 L = L-1 
      F = FC(LMIN1) 
      FP = FCP(LMIN1) 
    5 IF(KTRP.EQ.-1) GO TO 1
C*** MEME CALCUL POUR R = TURN SI RHO.LT.TURN 
C*** P + I.Q CALCULE EN MINL , EQUATION (32)
      P = 0
      Q = R-ETA 
      PL = 0 
      AR = -(ETA2+XLL1) 
      AI = ETA
      BR = Q+Q
      BI = 2
      WI = ETA+ETA
      DR = BR/(BR*BR+BI*BI) 
      DI = -BI/(BR*BR+BI*BI)
      DP = -(AR*DI+AI*DR) 
      DQ = AR*DR-AI*DI
    6 P = P+DP
      Q = Q+DQ
      PL = PL+2 
      AR = AR+PL
      AI = AI+WI
      BI = BI+2 
      D = AR*DR-AI*DI+BR
      DI = AI*DR+AR*DI+BI 
      T = 1/(D*D+DI*DI) 
      DR = T*D
      DI = -T*DI
      H = BR*DR-BI*DI-1
      K = BI*DR+BR*DI 
      T = DP*H-DQ*K 
      DQ = DP*K+DQ*H  
      DP = T
      IF(PL.GT.46000.0D0) GO TO 11
      IF(ABS(DP)+ABS(DQ).GE.(ABS(P)+ABS(Q))*ACC) GO TO 6
      P = P/R 
      Q = Q/R 
C*** CALCUL DE FP,G,GP, ET NORMALISATION DE F EN L = MINL **************
      G = (FP-P*F)/Q  
      GP = P*G-Q*F
      W = 1/SQRT(FP*G-F*GP) 
      G = W*G 
      GP = W*GP 
      IF(KTR.EQ.1) GO TO 8
      F = TF
      FP = TFP
      LMAX = MAXL 
C*** CALCUL DE G(MINL) ET GP(MINL) PAR INTEGRATION RUNGE-KUTTA A PARTIR 
C***         VOIR FOX ET MAYERS(1968) PG 202
      IF(RHO.LT.0.2D0*TURN) PACE = 999.0D0
      R3=1.0D0/3.0D0  
      H = (RHO-TURN)/(PACE+1) 
      H2 = H/2
      I2 = INT(PACE+0.001D0)
      ETAH = ETA*H
      H2LL = H2*XLL1  
      S = (ETAH+H2LL/R)/R-H2
    7 RH2 = R+H2
      T = (ETAH+H2LL/RH2)/RH2-H2
      K1 = H2*GP
      M1 = S*G
      K2 = H2*(GP+M1) 
      M2 = T*(G+K1) 
      K3 = H*(GP+M2)  
      M3 = T*(G+K2) 
      M3 = M3+M3
      K4 = H2*(GP+M3) 
      RH = R+H
      S = (ETAH+H2LL/RH)/RH-H2  
      M4 = S*(G+K3) 
      G = G+(K1+K2+K2+K3+K4)*R3 
      GP = GP+(M1+M2+M2+M3+M4)*R3 
      R = RH
      I2 = I2-1 
      IF(ABS(GP).GT.1.D300) GO TO 11
      IF(I2.GE.0) GO TO 7 
      W = 1/(FP*G-F*GP) 
C*** RECURRENCE AVANT A PARTIR DE GC(MINL) ET GCP(MINL) 
C*** RENORMALISATION DE FC ET FCP POUR CHAQUE VALEUR DE L **************
    8 GC(LMIN1) = G 
      GCP(LMIN1) = GP 
      IF(LMAX.EQ.MINL) GO TO 10 
      DO 9 L=LMIN1,LMAX 
      T = GC(L+1) 
      GC(L+1) = (GC(L)*GC(L+1)-GCP(L))/GCP(L+1) 
      GCP(L+1) = GC(L)*GCP(L+1)-GC(L+1)*T 
      FC(L+1) = W*FC(L+1) 
    9 FCP(L+1) = W*FCP(L+1) 
      FC(LMIN1) = W*FC(LMIN1) 
      FCP(LMIN1) = W*FCP(LMIN1) 
      RETURN
   10 FC(LMIN1) = W*F 
      FCP(LMIN1) = W*FP 
      RETURN
   11 W = 0
      G = 0
      GP = 0 
      GO TO 8 
      END 
      
*DEPHASE
      SUBROUTINE DEPHASE(MAX,H,W,Y,EPS,DELTA,L,ETA,QK,RFIN) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION FG(2),DFC(101),GC(101),DGC(101),FC(101),W(1),Y(1) 
      DD=0 
      PAS=H*EPS 
      QQ=QK*QK
      LL=L+1
      L1=L*LL 
      R=2*ETA*QK
      RMAX=MAX*H
      DO20 K=2,MAX
      RMAX=RMAX-H 
      IF(ABS(W(MAX-K+1)-(L1/RMAX+R)/RMAX).LT.EPS)GOTO20 
      IF(K.GT.3)GOTO21
      PRINT2001 
      STOP  
   20 CONTINUE
      RMAX=H*(MAX-2)  
   21 H2 = H**2 
      H212 = H2/12 
      AA=H212*(QQ-W(1)) 
      Y(1)=H**LL
      B0=0  
      B1=Y(1)*(1+AA)  
      DO22 K=2,MAX
      AA=H212*(QQ-W(K)) 
      B2=12*Y(K-1)-10*B1-B0 
      Y(K)=B2/(1+AA)  
      B0=B1 
      B1=B2 
   22 CONTINUE
      R=RMAX
      N1=INT((RMAX+1D-5)/H)+1 
      CALL COUFRA(QK*R,ETA,L,L,FC,DFC,GC,DGC) 
      FG(1)=FC(LL)
      FG(2)=GC(LL)
      DO23 N=N1,MAX 
      R=R+H 
      CALL COUFRA(QK*R,ETA,L,L,FC,DFC,GC,DGC) 
      C=FG(2)*Y(N)-GC(LL)*Y(N-1)
      D=FC(LL)*Y(N-1)-FG(1)*Y(N)
      Z=FG(1)*GC(LL)-FG(2)*FC(LL) 
      DEL=ABS(Z)/SQRT(C*C+D*D)  
      IF(ABS(1-DD/DEL).LT.PAS) GO TO 24 
      DD = DEL
      FG(1)=FC(LL)
      FG(2)=GC(LL)
   23 CONTINUE
      PRINT2000,L,QK  
   24 RFIN=R
      D=D/C 
      IF(Y(N)*(FC(LL)+GC(LL)*D).LT.0.)DEL=-DEL
      DELTA=ATAN(D) 
      DO25 K=1,MAX
   25 Y(K)=Y(K)*DEL 
      RETURN
 2000 FORMAT(/27H PAS DE CONVERGENCE POUR L=,I3,5H   K=,F10.6/) 
 2001 FORMAT(/46H LE POTENTIEL N'EST JAMAIS PUREMENT COULOMBIEN)
      END
      
*SIXJ 
      SUBROUTINE SIXJ(J1,J2,J3,L1,L2,L3,Q)
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION M(7),M1(4),M2(4),M3(4),FT(0:1000)
      LOGICAL ICE 
      COMMON/TEUFJ/ICE
      SAVE  
      DATA LMEM,(FT(I),I=0,10)/9,2*1.0D0,2.0D0,6.0D0,24.0D0,
     1 120.0D0,720.0D0,5040.0D0,40320.0D0,362880.0D0,3628800.0D0/  
      ENTRYSIXJI(J1,J2,J3,L1,L2,L3,Q) 
      I1=J1 
      I2=J2 
      I3=J3 
      K1=L1 
      K2=L2 
      K3=L3 
      IS=0  
   24 Q=0  
      M(1)=I1+I2+I3 
      M(2)=I1+K2+K3 
      M(3)=K1+I2+K3 
      M(4)=K1+K2+I3 
      DO 17 I=1,4 
      IF(MOD(M(I),2).EQ.1) GO TO 8
   17 CONTINUE
      L=MAX(I1+I2+K1+K2,I1+I3+K1+K3,I2+I3+K2+K3) 
      L=L/2 
      IF(L.LE.LMEM) GO TO 6 
      DO 10 I=LMEM,L  
   10 FT(I+1)=FT(I)*(I+1) 
      LMEM=L
    6 IF(I1.LT.ABS(I2-I3).OR.I1.GT.I2+I3) RETURN 
      IF(I1.LT.ABS(K2-K3).OR.I1.GT.K2+K3) RETURN 
      IF(K1.LT.ABS(I2-K3).OR.K1.GT.I2+K3) RETURN 
      IF(K1.LT.ABS(K2-I3).OR.K1.GT.K2+I3) RETURN 
      IF(I1) 8,2,1
    2 IF(I2.LT.0) GO TO 8 
    9 IF(I3.LT.0) GO TO 8 
   14 IF(K1.LT.0) GO TO 8 
   19 IF(K2.LT.0) GO TO 8 
   23 IF(K3.LT.0) GO TO 8 
   27 Q=SQRT(1.0D0/(I2+1)/(K2+1))  
      IS=(I2+K2+K1)/2+IS
      IF(MOD(IS,2).EQ.1) Q=-Q 
      RETURN
    1 IF(I1.GT.1) GO TO 3 
      IF(I2.LT.0) RETURN
   12 IF(I3.LT.0) RETURN
   16 IF(K1.LT.0) RETURN
   21 IF(K2.LT.0) RETURN
   25 IF(K3.LT.0) RETURN
   28 IF(I2.LT.I3) GO TO 4
      IC=I2 
      I2=I3 
      I3=IC 
      IC=K2 
      K2=K3 
      K3=IC 
    4 IF(K2.GT.K3) GO TO 5
      I11=I1+K1+I2-K2 
      I11=I11/2 
      I12=I11-I2+K2 
      Q=SQRT(I11*I12*1./I3/(I3+1)/K3/(K3+1))
      IS =I11+K2+IS 
      IF(MOD(IS,2).EQ.1) Q=-Q 
      RETURN
    5 I11=K3-K1+I2
      I11=I11/2+1 
      I12=I11+K1+1
      Q=SQRT(I11*I12*1.0D0/I3/(I3+1)/K2/(K2+1))
      IS =I12-1+IS
      IF(MOD(IS ,2).EQ.1) Q=-Q  
      RETURN
    3 IF(I2.GE.I1) GO TO 7
      IF(I2.LT.0) GO TO 8 
      IC=I2 
      I2=I1 
      I1=IC 
      IC=K1 
      K1=K2 
      K2=IC 
      IF(I1.EQ.0) GO TO 9 
      IF(I1.EQ.1) GO TO 12
    7 IF(I3.GE.I1) GO TO 13 
      IF(I3.LT.0) GO TO 8 
      IC=I3 
      I3=I1 
      I1=IC 
      IC=K3 
      K3=K1 
      K1=IC 
      IF(I1.EQ.0) GO TO 14
      IF(I1.EQ.1) GO TO 16
   13 IF(K1.GE.I1) GO TO 18 
      IF(K1.LT.0) GO TO 8 
      IC=K1 
      K1=I1 
      I1=IC 
      IC=K2 
      K2=I2 
      I2=IC 
      IF(I1.EQ.0) GO TO 19
      IF(I1.EQ.1) GO TO 21
   18 IF(K2.GE.I1) GO TO 22 
      IF(K2.LT.0) GO TO 8 
      IC=K2 
      K2=I1 
      I1=IC 
      IC=K1 
      K1=I2 
      I2=IC 
      IF (I1.EQ.0) GO TO 23 
      IF(I1.EQ.1) GO TO 25
   22 IF(K3.GE.I1) GO TO 26 
      IF(K3.LT.0) GO TO 8 
      IC=K3 
      K3=I1 
      I1=IC 
      IC=K1 
      K1=I3 
      I3=IC 
      IF(I1.EQ.0) GO TO 27
      IF(I1.EQ.1) GO TO 28
   26 M1(4)=I3
      M1(1)=I3
      M1(3)=K3
      M1(2)=K3
      M2(2)=I1
      M2(1)=I1
      M2(4)=K1
      M2(3)=K1
      M3(3)=I2
      M3(1)=I2
      M3(4)=K2
      M3(2)=K2
      M(1)=I1+I2+I3 
      M(2)=I1+K2+K3 
      M(3)=K1+I2+K3 
      M(4)=K1+K2+I3 
      Q1=1 
      DO 11 I=1,4 
      M(I)=M(I)/2 
   11 Q1=FT(M(I)-M1(I))*FT(M(I)-M2(I))*FT(M(I)-M3(I))*Q1/FT(M(I)+1) 
      Q1=SQRT(Q1) 
      M1(1)=I1+K1 
      M1(2)=I2+K2 
      M1(3)=I3+K3 
      IC=M1(1)+M1(2)  
      M(5)=IC/2 
      IC=M1(2)+M1(3)  
      M(6)=IC/2 
      IC=M1(1)+M1(3)  
      M(7)=IC/2 
      MAXZ=MIN(M(5),M(6),M(7)) 
      MINZ=MAX(M(1),M(2),M(3),M(4))
      X=0  
      DO 15 I=MINZ,MAXZ 
      Q2=1 
      DO 20 J=1,7 
      IJ=I-M(J) 
      IF(J.GT.4) IJ=-IJ 
   20 Q2=Q2*FT(IJ)
      Q2=FT(I+1)/Q2 
   15 X=-X+Q2 
      Q=X*Q1
      IS=MAXZ+IS
      IF(MOD(IS,2).EQ.1) Q=-Q 
      RETURN
    8 PRINT 1010,J1,J2,J3,L1,L2,L3
      ICE=.FALSE. 
 1010 FORMAT(10H ERREUR 6J,2(3X,3I3)) 
      RETURN
      ENTRY RACAH(J1,J2,J3,L1,L2,L3,Q)
      IS=(J1+J2+J3+L1)/2
      I1=J1 
      I2=J2 
      I3=L2 
      K1=L1 
      K2=J3 
      K3=L3 
      GO TO 24
      END 
      
*NUM1L
      SUBROUTINE NUM1L(N,H,E,S2,U,S,NO,EPS) 
CC  VERSION CORRIGEE LE 21 NOV 72 
C*****INTEGRATION THE SCHROEDINGER EQUATION BY THE NUMEROV METHOD
C*****POUR E NEGATIF  
C*****RECHERCHE DE L"ENERGIE PROPRE PAR LA METHODE DE RAPHSON-NEWTON
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION U(1),S(1) 
      DATA RAP1,RAP2/0,0/ 
      H12=H*H/12 
C*****CONTROLE DES CONDITIONS ASYMPTOTIQUES 
      IF(E.GT.0) E=0
      DEI=0
      EPSS=.1D-10 
      IF(U(N-1).GT.EPSS) GO TO 10 
      DEI=U(N-1)-EPSS 
      DO 8 K=1,N
    8 U(K)=U(K)-DEI 
   10 U(N)=U(N-1) 
C*****CALCUL DU NOMBRE D"ETATS LIES PAR INTEGRATION A ENERGIE NULLE 
      S(1)=1.D-10 
      B0=0 
      AA=H12*U(1) 
      IF (S2) 16,18,16
   16 B0=-S(1)*AA 
   18 B1=S(1)*(1-AA) 
      DO 38 K=2,N 
      B2=12*S(K-1)-10*B1-B0 
      IF (ABS(B2).LT.1.D+10) GO TO 22 
      B2=B2*1.D-20
      B1=B1*1.D-20
   22 AA=H12*U(K) 
      S(K)=B2/(1-AA) 
      B0=B1 
   38 B1=B2 
      DO 42 K=5,N 
      N0=K  
      IF(U(K).LT.0) GO TO 44 
   42 CONTINUE
   44 NEL=0 
      DO 52 K=N0,N
      IF (S(K-1)*S(K)) 46,50,52 
   46 NEL=NEL+2 
      GO TO 52
   50 NEL=NEL+1 
   52 CONTINUE
      NEL=NEL/2 
      IF(NEL.GT.NO) GO TO 64
      IF(NEL.EQ.NO) GO TO 60
   62 NO=-1 
      RETURN
   60 RAP1=S(N-1)/S(N)
      RAP2=EXP(H*SQRT(U(N-1)-E))
      IF(RAP1.LT.RAP2) GO TO 62 
C*****CALCUL DE EMIN ET EMAX ENTRE LESQUELLES SE TROUVE L"ENERGIE PROPRE
   64 UMIN=U(1) 
      DO 70 K=2,N 
      IF(U(K).LT.UMIN) UMIN=U(K)
   70 CONTINUE
      EMIN=UMIN 
      EMAX=0 
C*****DEBUT DE LA RECHERCHE DE L"ENERGIE PROPRE DANS L"INTERVALLE MAXIMU
      TE=EMAX-EMIN
C*****REJET DE L"ENERGIE D"ESSAI E PROPOSEE SI ELLE EST A L"EXTERIEUR DE
C*****BORNES (EMIN,EMAX)
      IF((E.LT.EMIN).OR.(E.GT.EMAX)) E=EMIN+TE/2 
      E1=EMIN 
      E2=EMAX 
      J=2 
      I=1 
      GO TO 102 
C*****REDUCTION DES BORNES EMIN ET EMAX 
   90 EMIN=E1 
      EMAX=E2 
      TE=EMAX-EMIN
      J=2 
   98 I=1 
  100 E=EMIN+TE*I/J 
  102 DE=0 
  104 E=E+DE
      IF(E.GT.0) GO TO 204 
      S(N)=1.D-10 
      N1=N-1
      RAP2=EXP(H*SQRT((U(N-1)+U(N))/2-E))
      S(N1)=S(N)*RAP2 
      AA=H12*(U(N1)-E)
      B0=S(N)*(1-AA) 
      B1=S(N1)*(1-AA)
      N1=N-2
      DO 138 KAUX=1,N1
      K=N1-KAUX+1 
      B2=12*S(K+1)-10*B1-B0 
      AA=H12*(U(K)-E) 
      S(K)=B2/(1-AA) 
      B0=B1 
      B1=B2 
      IF(U(K).LT.E) GO TO 140 
  138 CONTINUE
  140 N1=K  
C*****NORMALISATION DE LA FONCTION D"ONDE A S(N1) 
      DO 146 KAUX=N1,N
      K=N-KAUX+N1 
  146 S(K)=S(K)/S(N1) 
C*****DEBUT DE L"INTEGRATION VERS L"EXTERIEUR JUSQU"A N1
      S(1)=1.D-10 
      B0=0 
      AA=H12*(U(1)-E) 
      IF(S2) 156,158,156
  156 B0=-S(1)*AA 
  158 B1=S(1)*(1-AA) 
      DO 170 K=2,N1 
      B2=12*S(K-1)-10*B1-B0 
      AA=H12*(U(K)-E) 
      S(K)=B2/(1-AA) 
      B0=B1 
  170 B1=B2 
C*****NORMALISATION DE LA FONCTION A S(N1)
      DO 174 K=1,N1 
  174 S(K)=S(K)/S(N1) 
C*****CALCUL DE LA CORECTION D"ENERGIE
      SOM=0
      DO 180 K=1,N
  180 SOM=SOM+S(K)*S(K) 
      DE=((-S(N1-1)+2-S(N1+1))/(H*H)+U(N1)-E)/SOM  
      IF(ABS(DE).GT.EPS) GO TO 104
C*****CALCUL DU NOMBRE DE NOEUDS DE L"ETAT PROPRE TROUVE
      DO 182 K=5,N
      IF(U(K).LT.E) GO TO 184 
  182 CONTINUE
  184 N0=K  
      NEL=0 
      DO 192 K=N0,N1  
      IF(S(K-1)*S(K)) 186,190,192 
  186 NEL=NEL+2 
      GO TO 192 
  190 NEL=NEL+1 
  192 CONTINUE
      NEL=NEL/2 
C*****L"ETAT PROPRE TROUVE EST-IL LE BON  
      IF(NEL-NO) 198,214,202
  198 IF(E.GT.E1) E1=E
      GO TO 204 
  202 IF(E.LT.E2) E2=E
  204 I=I+2 
      IF (I.LE.J)  GO TO 100
      J=2*J 
      IF(ABS(E1-EMIN).GT.EPS.OR.ABS(EMAX-E2).GT.EPS) GO TO 90 
      GO TO 98
C*****NORMALISATION DE LA FONCTION PROPRE 
  214 SOM=1/SQRT(SOM*H)
      DO 218 K=1,N
  218 S(K)=S(K)*SOM 
      E=E+DEI 
      RETURN
C*****DEBUT FORMATS 
 2000 FORMAT(/,35X,56HL"ETAT DEMANDE N"EST PAS LIE. RETOUR DE NUM1L AVEC
     1 NO=-1,/) 
C*****FIN FORMATS 
      END 
      
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      INTEGER N,NMAX
      REAL*8 YP1,YPN,X(N),Y(N),Y2(N)
      PARAMETER (NMAX=12000)
      INTEGER I,K
      REAL*8 P,QN,SIG,UN,U(NMAX)
      IF (YP1.GT.0.99D40) THEN
        Y2(1)=0.0D0
        U(1)=0.0D0
      ELSE
        Y2(1)=-0.50D0
C        Y2(1)=-0.0D0
        U(1)=(3.0D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      END IF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.0D0
        Y2(I)=(SIG-1.0D0)/P
        U(I)=(6.0D0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-
     & (Y(I)-Y(I-1))/(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT.0.99D40) THEN
        QN=0.0D0
        UN=0.0D0
      ELSE
        QN=0.5D0
        UN=(3.0D0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      END IF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.0D0)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END

      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      INTEGER N
      REAL*8 X,Y,XA(N),Y2A(N),YA(N)
      INTEGER K,KHI,KLO
      REAL*8 A,B,H
      KLO=1
      KHI=N
2     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        END IF
C      GOTO 1
      GOTO 2
      END IF
      H=XA(KHI)-XA(KLO)
C      IF(H.EQ.0.0D0) PAUSE 'BAD XA INPUT IN SPLINT'
      IF(H.EQ.0.0D0) THEN
       A=0
       B=0
      ELSE
       A=(XA(KHI)-X)/H
       B=(X-XA(KLO))/H
       Y=A*YA(KLO)+B*YA(KHI)+((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*
     &(H**2)/6.0D0
      END IF
      RETURN
      END
                                           
      SUBROUTINE DER1(Y,X,N,H)                                            
C NUMERICAL DERIVATION OF THE FUNCTION Y KNOWN AT N POINTS.             
C THIS SUBROUTINE REQUIRES AT LEAST 7 POINTS.                           
C AND RETURNS THE VALUE OF H*(D(Y)/DR) IN X.                            
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(N),Y(N)
      N3=N-3
      X(1)=(-147.D0*Y(1)+360.D0*Y(2)-450.D0*Y(3)+400.D0*Y(4)-225.D0*Y(5)
     1+72.D0*Y(6)-10.D0*Y(7))/60.D0                                     
      X(2)=(-10.D0*Y(1)-77.D0*Y(2)+150.D0*Y(3)-100.D0*Y(4)+50.D0*Y(5)-15
     1.D0*Y(6)+2.D0*Y(7))/60.D0                                         
      X(3)=(2.D0*Y(1)-24.D0*Y(2)-35.D0*Y(3)+80.D0*Y(4)-30.D0*Y(5)+8.D0*Y
     1(6)-Y(7))/60.D0                                                   
      DO I=4,N3                                                       
      X(I)=(45.D0*(Y(I+1)-Y(I-1))-9.D0*(Y(I+2)-Y(I-2))+Y(I+3)-Y(I-3))/60
     1.D0
      END DO                                                               
      X(N-2)=(Y(N-6)-8.D0*Y(N-5)+30.D0*Y(N-4)-80.D0*Y(N3)+35.D0*Y(N-2)+2
     14.D0*Y(N-1)-2.D0*Y(N))/60.D0                                      
      X(N-1)=(-2.D0*Y(N-6)+15.D0*Y(N-5)-50.D0*Y(N-4)+100.D0*Y(N3)-150.D0
     1*Y(N-2)+77.D0*Y(N-1)+10.D0*Y(N))/60.D0                            
      X(N)=(10.D0*Y(N-6)-72.D0*Y(N-5)+225.D0*Y(N-4)-400.D0*Y(N3)+450.D0*
     1Y(N-2)-360.D0*Y(N-1)+147.D0*Y(N))/60.D0
      DO I = 1, N
       X(I) = X(I)/H
      END DO
      RETURN
      END
