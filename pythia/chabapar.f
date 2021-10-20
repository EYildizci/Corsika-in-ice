      SUBROUTINE CHABAPAR(EC,NY)
c  charmed baryon partonic process
      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER PYCOMP
      CHARACTER CHAF*16

      DOUBLE PRECISION m1,m2,m12
      DOUBLE PRECISION MP,M0,DM,MC,MLAMBDAC,MCHARM,MN,MM1,MM2,C,D
      INTEGER COL
      DOUBLE PRECISION LH(25)
      CHARACTER CHANNEL(0:25)*30
      CHARACTER CHAU*16
      CHARACTER LTAR*5
      INTEGER IKC(200)
      INTEGER KFC(-5000:5000)

      DOUBLE PRECISION A
      DOUBLE PRECISION XGLIST(100),ECMLIST(100),QMASS(100)
      INTEGER KFLIST(100)
      INTEGER KLIST(10000,5)
      CHARACTER TLIST(100)*5
      DOUBLE PRECISION PLIST(10000,5),VLIST(10000,5)
      

      COMMON/CONVERTER/KFC

      
      COMMON/ENER/MP,EMIN,EP,ELAB,PXLAB,PYLAB,PZLAB
      COMMON/EMAXLH/EMAXB,EMAXM,IMAXB,IMAXM,IKMAXB,IKMAXM


      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT4/CHAF(500,2)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDATR/MRPY(6),RRPY(100)

C     MAKING TAUS STABLE
C     TAU-
      KC=PYCOMP(15)
      MDCY(KC,1)=0
C     TAU+
      KC=PYCOMP(-15)
      MDCY(KC,1)=0
C     MAKING MUONS STABLE
C     MU-
      KC=PYCOMP(13)
      MDCY(KC,1)=0
C     MU+
      KC=PYCOMP(-13)
      MDCY(KC,1)=0

C     MAKING NUCLEONS STABLE
C     P  
      KC=PYCOMP(2212)
      MDCY(KC,1)=0
C     N     
      KC=PYCOMP(2112)
      MDCY(KC,1)=0

C     MAKING PIONS STABLE
C     PI0
      KC=PYCOMP(111)
      MDCY(KC,1)=0
C     PI-
      KC=PYCOMP(-211)
      MDCY(KC,1)=0
C     PI+
      KC=PYCOMP(211)
      MDCY(KC,1)=0
C     MAKING KAONS STABLE 
C     K0     
C     KC=PYCOMP(311)
C     MDCY(KC,1)=0
C     K0B     
C     KC=PYCOMP(-311)
C     MDCY(KC,1)=0
C     KL     
      KC=PYCOMP(130)
      MDCY(KC,1)=0
C     KS     
      KC=PYCOMP(310)
      MDCY(KC,1)=0
C     K+     
      KC=PYCOMP(321)
      MDCY(KC,1)=0
C     K-     
      KC=PYCOMP(-321)
      MDCY(KC,1)=0

C-------------------------
      KC=PYCOMP(411)
      MDCY(KC,1)=0
      KC=PYCOMP(421)
      MDCY(KC,1)=0
      KC=PYCOMP(431)
      MDCY(KC,1)=0
      KC=PYCOMP(4122)
      MDCY(KC,1)=0
      KC=PYCOMP(4232)
      MDCY(KC,1)=1

      KC=PYCOMP(4132)
      MDCY(KC,1)=1

      KC=PYCOMP(4332)
      MDCY(KC,1)=1

C     D- BAR
      KC=PYCOMP(-411)
      MDCY(KC,1)=0
C     D0 BAR
      KC=PYCOMP(-421)
      MDCY(KC,1)=0
C     Ds BAR
      KC=PYCOMP(-431)
      MDCY(KC,1)=0
C     LAMBDA_C BAR
      KC=PYCOMP(-4122)
      MDCY(KC,1)=0
C     XI_C_- BAR
      KC=PYCOMP(-4232)
      MDCY(KC,1)=1
C     XI_C_0 BAR
      KC=PYCOMP(-4132)
      MDCY(KC,1)=1
C     OMEGA_C_0 BAR
      KC=PYCOMP(-4332)
      MDCY(KC,1)=1

C     FORCED DECAY
C                SIGMAC++
      KC = PYCOMP(4222)
      MDCY(KC,1) = 1
C                SIGMAC+
      KC = PYCOMP(4212)
      MDCY(KC,1) = 1
C                SIGMAC0
      KC = PYCOMP(4112)
      MDCY(KC,1) = 1
C                A-SIGMAC++
      KC = PYCOMP(-4222)
      MDCY(KC,1) = 1
C                A-SIGMAC+
      KC = PYCOMP(-4212)
      MDCY(KC,1) = 1
C                A-SIGMAC0
      KC = PYCOMP(-4112)
      MDCY(KC,1) = 1
C                ETAC
      KC = PYCOMP(441)
      MDCY(KC,1) = 1
C-------------------------
C                K*(892)0
      KC = PYCOMP(313)
      MDCY(KC,1) = 1
C                K*(892)0
      KC = PYCOMP(-313)
      MDCY(KC,1) = 1
C-------------------------

c     DELTAS
      KC = PYCOMP(2114)
      MDCY(KC,1) = 1
      KC = PYCOMP(-2114)
      MDCY(KC,1) = 1
      KC = PYCOMP(1114)
      MDCY(KC,1) = 1
      KC = PYCOMP(-1114)
      MDCY(KC,1) = 1
      KC = PYCOMP(-2214)
      MDCY(KC,1) = 1
      KC = PYCOMP(2214)
      MDCY(KC,1) = 1
C-------------------------


      KC=PYCOMP(2212)
      MP=PMAS(KC,1)
      KC=PYCOMP(4122)
      MLAMBDAC=PMAS(KC,1)
      MCHARM=1.27D0


      MSTP(111)=0
      MSTP(125)=1
      MSEL=0
      DO IE=1,500
         MSUB(IE)=0
      ENDDO

*     Hard QCD processes:
*--------------------
*     11: f + f' -> f + f' (QCD) 
      MSUB(11)=1
*     12: f + fbar -> f' + fbar'
      MSUB(12)=1
*     13: f + fbar -> g + g 
      MSUB(13)=1
*     28: f + g -> f + g 
      MSUB(28)=1
*     53: g + g -> f + fbar
      MSUB(53)=1
*     68: g + g -> g + g
      MSUB(68)=1
*     95: Low-pT scattering	(keep always active)
      MSUB(95)=1

	

      A = 1.D0
      ZB = 0.D0

      EP = EC

C     ASSIGN SEA PAIRS ENERGIES

      NP = 7
      NN = 7
      FAIL = 0.D0
      IF (NY.GT.1) THEN
 123     CONTINUE
         FAIL = 0.D0
         ZB = 0.D0
         DO 124 INW=1,NY-1
            XGLIST(INW) = GXCBP(1.0)
            ZB = ZB+XGLIST(INW)
C     ASSIGN IDS
            IF (PYR(0).LE.2.D0/2.3D0)THEN
               KFLIST(INW) = 111
               QMASS(INW) = 2.D0*PYMASS(1)
               IF (PYR(0).LE.0.5D0) THEN
                  KFLIST(INW) = 221
                  QMASS(INW) = 2.D0*PYMASS(2)
               ENDIF
            ELSE
               KFLIST(INW) = 331
               QMASS(INW) = 2.D0*PYMASS(3)
            ENDIF
            
C     ECMS
            ECMLIST(INW) = 2.D0*(QMASS(INW)**2+(MP*EP*XGLIST(INW)))/
     *           DSQRT(2.D0*MP*EP*XGLIST(INW)+MP**2+QMASS(INW)**2)

            
            IF ( ECMLIST(INW) .LE. PARP(2) ) THEN
               ZB = ZB-XGLIST(INW)
               TLIST(INW) = 'F'
               GOTO 124
            ENDIF

            IF (PYR(0) .LT. 1.D0*NP/(1.D0*NP+1.D0*NN)) THEN
               TLIST(INW) = 'P+'
               NP = NP-1
            ELSE
               TLIST(INW) = 'n0'
               NP = NP-1
            ENDIF


 124     CONTINUE
         IF (ZB.GE.1.D0) GOTO 123
         IF (DSQRT(2.D0*MP*(MP+((MLAMBDAC-MCHARM)/
     *       MLAMBDAC*(EP*(1.D0-ZB))))).LE. PARP(2)) GOTO 123
        

      ENDIF


C     REMOVE QQBAR PAIRS ENERGY: ENERGY FOR QCD COLLISION
      EPR = EP*(1.D0-ZB)

C..   To stop hadronization:

C-----------------------
      EAUX = 0.D0
      M = 1
      IF (NY .GT. 1 ) THEN
         DO INW=1,NY-1

            IF (TLIST(INW).EQ.'F') GOTO 500
            
            CALL PYINIT('CMS','pi0',TLIST(INW),ECMLIST(INW))            

            MINT(11) = KFLIST(INW)
            VINT(3) = QMASS(INW)

            CALL PYEVNT
            CALL PYEXEC
C            CALL PYLIST(1)
            DO IPART=9,N
               IF (K(IPART,1).EQ.1) THEN
                  CALL BOOSTCBP(IPART,EP*XGLIST(INW),MP)
                  DO COL=1,5
                     KLIST(M,COL) = K(IPART,COL)
                     PLIST(M,COL) = P(IPART,COL)
                     VLIST(M,COL) = V(IPART,COL)
                  ENDDO

                  
                  EAUX = EAUX + PLIST(M,4)
                  M = M+1
               ENDIF
            ENDDO

 500        CONTINUE
         ENDDO
      ENDIF

      IF (PYR(0) .LT. 1.D0*NP/(1.D0*NP+1.D0*NN)) THEN
         LTAR = 'P+'
      ELSE
         LTAR = 'n0'
      ENDIF

C MAIN INTERACTION      
*=======================================================================

      KC=PYCOMP(2212)
      MP=PMAS(KC,1)
      KC=PYCOMP(4122)
      MLAMBDAC=PMAS(KC,1)
      MCHARM=1.27D0

      EP=(MLAMBDAC-MCHARM)/MLAMBDAC*EPR
      ECM = DSQRT(2.D0*MP*(MP+EP))
      EPCOR=ECM/2.D0*MLAMBDAC/(MLAMBDAC-MCHARM)
      DELTAE=EPCOR-ECM/2D0
      DELTAPZ=EPCOR*DSQRT(1.D0-MLAMBDAC**2/EPCOR**2)
     .     -ECM/2D0*DSQRT(1.D0-(2D0*MP)**2/ECM**2)


      CALL PYINIT('CMS','P+',LTAR,ECM)
C..   To init u to c conversion:
      CALL KFCINITCBP
 777  CONTINUE
      IC=0
      ig=0

      CALL PYEVNT
c      CALL PYLIST(1)
*=======================================================================
C     U-C SUBSTITUTION ROUTINE


c..   for hard:
      IK1=9
      PMAX=0D0
      DO IK=9,N
         IF (K(IK,2).NE.22) THEN
            IF (K(IK,2).NE.21.and.K(IK,2).GT.0) THEN
               IF(K(IK,2).NE.KFC(K(IK,2)).AND.P(IK,3).GT.0.D0) THEN
                  IF(K(IK,2).NE.4) THEN
                     IC=IC+1
                     IKC(IC)=IK
                     IF (P(IK,3).GT.PMAX) THEN
                        PMAX=P(IK,3)
                        IKMAX=IK
                        ICMAX=IC
                     END IF
                  END IF

               END IF
            END IF
         END IF

      END DO


      IF ( IC .EQ. 0 ) GOTO 777

C..   To modify an entry to replace charm:
c..   random:
      ICR=INT(IC*PYR(0))+1
      IK=IKC(ICR)
      IF (IC.GE.3) THEN
         ICR=ICMAX
         IK=IKMAX
      END IF

      DO I=IK,9,-1
         IF (K(I-1,1).EQ.1) THEN
	    IK1=I
	    GOTO 33
         END IF
      END DO
 33   DO I=IK,N
         IF (K(I,1).EQ.1) THEN
	    IK2=I
	    GOTO 34
         END IF
      END DO
 34   CONTINUE
      

c..   to include both cu_0 and cu_1:
      if (K(IK,2).eq.4203) then
         ii=INT(4*PYR(0))+1
         if (ii.le.3) K(IK,2)=4201
      end if

      IF (IK.EQ.IK1) THEN
         IKO=IK2
      ELSE
         IKO=IK1
      END IF

c..........................................................
      IF (IK.EQ.9.AND.K(IK,1).EQ.1) THEN
         K(IK,2)=KFC(K(IK,2))
         KC=PYCOMP(K(IK,2))
         M0=PMAS(KC,1)
         P(IK,5)=M0
         P(IK,3)=P(IK,3)+DELTAPZ
         P(IK,4)=DSQRT(M0**2+P(IK,1)**2+P(IK,2)**2+P(IK,3)**2)
         GOTO 77
      END IF
c..........................................................

      D0=0.D0
      D1=0.D0
      D2=0.D0
      D3=0.D0
      DO I=IK1,IK2
         D0=D0+P(I,4)
         D1=D1+P(I,1)
         D2=D2+P(I,2)
         D3=D3+P(I,3)
      END DO

****************
*     (a1)
      DM=DSQRT(D0**2-D1**2-D2**2-D3**2)
****************


      K(IK,2)=KFC(K(IK,2))
      KC=PYCOMP(K(IK,2))
      M0=PMAS(KC,1)

****************
*     (b)
      d0=d0+DELTAE
      d3=d3+DELTAPZ
****************

      m1=M0
      x1=P(IK,1)
      y1=P(IK,2)

      m2=P(IKO,5)
      x2=P(IKO,1)
      y2=P(IKO,2)
      
      z3=0.d0
      e3=0.d0
      DO I=IK1+1,IK2-1
         P(I,1)=P(I,1)*MCHARM/MLAMBDAC
         P(I,2)=P(I,2)*MCHARM/MLAMBDAC
         P(I,3)=P(I,3)*MCHARM/MLAMBDAC
         P(I,4)=P(I,4)*MCHARM/MLAMBDAC
         z3=z3+P(I,3)
         e3=e3+P(I,4)
      END DO


      z1a=(4*d0**2*d3-4*d3**3-8*d0*d3*e3+4*d3*e3**2+4*d3*m1**2-
     -     4*d3*m2**2+4*d3*x1**2-4*d3*x2**2+4*d3*y1**2-4*d3*y2**2- 
     -     4*d0**2*z3+12*d3**2*z3+8*d0*e3*z3-4*e3**2*z3-4*m1**2*z3+
     -     4*m2**2*z3-4*x1**2*z3+4*x2**2*z3-4*y1**2*z3+4*y2**2*z3- 
     -     12*d3*z3**2+4*z3**3-dsqrt((-4*d0**2*d3+4*d3**3+8*d0*d3*e3-
     -     4*d3*e3**2-4*d3*m1**2+4*d3*m2**2-4*d3*x1**2+4*d3*x2**2- 
     -     4*d3*y1**2+4*d3*y2**2+4*d0**2*z3-12*d3**2*z3-8*d0*e3*z3+
     -     4*e3**2*z3+4*m1**2*z3-4*m2**2*z3+4*x1**2*z3- 
     -     4*x2**2*z3+4*y1**2*z3-4*y2**2*z3+12*d3*z3**2-4*z3**3)**2- 
     -     4*(4*d0**2-4*d3**2-8*d0*e3+4*e3**2+8*d3*z3-4*z3**2)*
     -     (-d0**4+2*d0**2*d3**2-d3**4+4*d0**3*e3-4*d0*d3**2*e3-
     -     6*d0**2*e3**2+2*d3**2*e3**2+4*d0*e3**3-e3**4+ 
     -     2*d0**2*m1**2+2*d3**2*m1**2-4*d0*e3*m1**2+2*e3**2*m1**2-   
     -     m1**4+2*d0**2*m2**2-2*d3**2*m2**2-4*d0*e3*m2**2+ 
     -     2*e3**2*m2**2+2*m1**2*m2**2-m2**4+2*d0**2*x1**2+ 
     -     2*d3**2*x1**2-4*d0*e3*x1**2+2*e3**2*x1**2-2*m1**2*x1**2+ 
     -     2*m2**2*x1**2-x1**4+2*d0**2*x2**2-2*d3**2*x2**2- 
     -     4*d0*e3*x2**2+2*e3**2*x2**2+2*m1**2*x2**2-2*m2**2*x2**2+ 
     -     2*x1**2*x2**2-x2**4+2*d0**2*y1**2+2*d3**2*y1**2- 
     -     4*d0*e3*y1**2+2*e3**2*y1**2-2*m1**2*y1**2+2*m2**2*y1**2- 
     -     2*x1**2*y1**2+2*x2**2*y1**2-y1**4+2*d0**2*y2**2-
     -     2*d3**2*y2**2-4*d0*e3*y2**2+2*e3**2*y2**2+2*m1**2*y2**2- 
     -     2*m2**2*y2**2+2*x1**2*y2**2-2*x2**2*y2**2+2*y1**2*y2**2-
     -     y2**4-4*d0**2*d3*z3+4*d3**3*z3+8*d0*d3*e3*z3- 
     -     4*d3*e3**2*z3-4*d3*m1**2*z3+4*d3*m2**2*z3-4*d3*x1**2*z3+
     -     4*d3*x2**2*z3-4*d3*y1**2*z3+4*d3*y2**2*z3+2*d0**2*z3**2- 
     -     6*d3**2*z3**2-4*d0*e3*z3**2+2*e3**2*z3**2+2*m1**2*z3**2-
     -     2*m2**2*z3**2+2*x1**2*z3**2-2*x2**2*z3**2+2*y1**2*z3**2- 
     -     2*y2**2*z3**2+4*d3*z3**3-z3**4)))/(2.*(4*d0**2-4*d3**2-
     -     8*d0*e3+4*e3**2+8*d3*z3-4*z3**2))

      z1b=(4*d0**2*d3-4*d3**3-8*d0*d3*e3+4*d3*e3**2+4*d3*m1**2-
     -     4*d3*m2**2+4*d3*x1**2-4*d3*x2**2+4*d3*y1**2-4*d3*y2**2- 
     -     4*d0**2*z3+12*d3**2*z3+8*d0*e3*z3-4*e3**2*z3-4*m1**2*z3+
     -     4*m2**2*z3-4*x1**2*z3+4*x2**2*z3-4*y1**2*z3+4*y2**2*z3- 
     -     12*d3*z3**2+4*z3**3+dsqrt((-4*d0**2*d3+4*d3**3+8*d0*d3*e3-
     -     4*d3*e3**2-4*d3*m1**2+4*d3*m2**2-4*d3*x1**2+4*d3*x2**2- 
     -     4*d3*y1**2+4*d3*y2**2+4*d0**2*z3-12*d3**2*z3-8*d0*e3*z3+
     -     4*e3**2*z3+4*m1**2*z3-4*m2**2*z3+4*x1**2*z3- 
     -     4*x2**2*z3+4*y1**2*z3-4*y2**2*z3+12*d3*z3**2-4*z3**3)**2- 
     -     4*(4*d0**2-4*d3**2-8*d0*e3+4*e3**2+8*d3*z3-4*z3**2)*
     -     (-d0**4+2*d0**2*d3**2-d3**4+4*d0**3*e3-4*d0*d3**2*e3-
     -     6*d0**2*e3**2+2*d3**2*e3**2+4*d0*e3**3-e3**4+ 
     -     2*d0**2*m1**2+2*d3**2*m1**2-4*d0*e3*m1**2+2*e3**2*m1**2-   
     -     m1**4+2*d0**2*m2**2-2*d3**2*m2**2-4*d0*e3*m2**2+ 
     -     2*e3**2*m2**2+2*m1**2*m2**2-m2**4+2*d0**2*x1**2+ 
     -     2*d3**2*x1**2-4*d0*e3*x1**2+2*e3**2*x1**2-2*m1**2*x1**2+ 
     -     2*m2**2*x1**2-x1**4+2*d0**2*x2**2-2*d3**2*x2**2- 
     -     4*d0*e3*x2**2+2*e3**2*x2**2+2*m1**2*x2**2-2*m2**2*x2**2+ 
     -     2*x1**2*x2**2-x2**4+2*d0**2*y1**2+2*d3**2*y1**2- 
     -     4*d0*e3*y1**2+2*e3**2*y1**2-2*m1**2*y1**2+2*m2**2*y1**2- 
     -     2*x1**2*y1**2+2*x2**2*y1**2-y1**4+2*d0**2*y2**2-
     -     2*d3**2*y2**2-4*d0*e3*y2**2+2*e3**2*y2**2+2*m1**2*y2**2- 
     -     2*m2**2*y2**2+2*x1**2*y2**2-2*x2**2*y2**2+2*y1**2*y2**2-
     -     y2**4-4*d0**2*d3*z3+4*d3**3*z3+8*d0*d3*e3*z3- 
     -     4*d3*e3**2*z3-4*d3*m1**2*z3+4*d3*m2**2*z3-4*d3*x1**2*z3+
     -     4*d3*x2**2*z3-4*d3*y1**2*z3+4*d3*y2**2*z3+2*d0**2*z3**2- 
     -     6*d3**2*z3**2-4*d0*e3*z3**2+2*e3**2*z3**2+2*m1**2*z3**2-
     -     2*m2**2*z3**2+2*x1**2*z3**2-2*x2**2*z3**2+2*y1**2*z3**2- 
     -     2*y2**2*z3**2+4*d3*z3**3-z3**4)))/(2.*(4*d0**2-4*d3**2-
     -     8*d0*e3+4*e3**2+8*d3*z3-4*z3**2))

      z2a=d3-z3-z1a
      z2b=d3-z3-z1b

      e1a=dsqrt(m1**2+x1**2+y1**2+z1a**2)
      e2a=dsqrt(m2**2+x2**2+y2**2+z2a**2)
      e1b=dsqrt(m1**2+x1**2+y1**2+z1b**2)
      e2b=dsqrt(m2**2+x2**2+y2**2+z2b**2)


      z1=z1a
      z2=z2a
      e1=e1a
      e2=e2a
      
      if (e1a.lt.e1b) then
         z1=z1b
         z2=z2b
         e1=e1b
         e2=e2b
      end if


      P(IK,3)=z1
      P(IK,4)=e1
      P(IK,5)=m0

      P(IKO,3)=z2
      P(IKO,4)=e2

 77   P(1,3)=EPCOR*DSQRT(1.D0-MLAMBDAC**2/EPCOR**2)
      P(1,4)=EPCOR
      P(1,5)=MLAMBDAC

*=======================================================================



C     --------------------------
      CALL PYEXEC
C      CALL PYLIST(1)
C      EAUX=0
      DO IPART=9,N
         IF (K(IPART,1).EQ.1) THEN
            CALL BOOSTCBP(IPART,EP,MP)

            DO COL=1,5
               KLIST(M,COL) = K(IPART,COL)
               PLIST(M,COL) = P(IPART,COL)
               VLIST(M,COL) = V(IPART,COL)
            ENDDO

            EAUX = EAUX + PLIST(M,4)
            M = M+1
         ENDIF
      ENDDO
      
C     CORSIKA READS PYJETS
      DO III=1,M-1
         DO COL=1,5
            K(III,COL) = KLIST(III,COL)
            P(III,COL) = PLIST(III,COL)
            V(III,COL) = VLIST(III,COL)
         ENDDO
      ENDDO
      
      N=M-1
      
      IF (EAUX.GT.EC) THEN
c     PRINT*,'REESCALED ',EAUX,'  TO  ',EC
c     SHOULD NOT HAPPEN BUT SOMETIMES EXCEEDS BY A SMALL AMOUNT
         DO III=1,N
            P(III,4) = P(III,4)*EC/EAUX
            P(III,3) = DSQRT(P(III,4)**2-P(III,1)**2-
     *           P(III,2)**2-P(III,5)**2)
         ENDDO
      ENDIF


      END
     
      
*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

	SUBROUTINE BOOSTCBP(IK,DE0,DM0)

      	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 	DOUBLE PRECISION MP

	COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
	COMMON/ENER/MP,EMIN,EP,ELAB,PXLAB,PYLAB,PZLAB

C	GAMMA=DSQRT(EP/(2.D0*MP))
	GAMMA=DSQRT(DE0/(2.D0*DM0))
	BETA=DSQRT(1.D0-GAMMA**(-2))
	ELAB =GAMMA*(P(IK,4)+BETA*P(IK,3))
	PXLAB=P(IK,1)
	PYLAB=P(IK,2)
	PZLAB=GAMMA*(BETA*P(IK,4)+P(IK,3))

        P(IK,3) = PZLAB
        P(IK,4) = ELAB


	RETURN
	END

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

	SUBROUTINE KFCINITCBP

        IMPLICIT DOUBLE PRECISION(A-H,O-Z)

	integer i,j,k,s2

	INTEGER KFC(-5000:5000)
	COMMON/CONVERTER/KFC

* quarks:

c	print*,'quarks:'
	do i=1,6
	  kf=i
	  kfc(kf)=kf
	  if (i.eq.2) kfc(kf)=4
c	  call showCBP(kf,kfc(kf))
	end do

* diquarks:

c	print*,'diquarks:'
	do s2=0,2,2
	  do j=1,3
	    do i=j,3
	      kf=1000*i+100*j+s2+1
	      kfc(kf)=kf
	      if (i.eq.2) kfc(kf)=1000*4+100*j+s2+1
	      if (j.eq.2.and.i.ne.2) kfc(kf)=1000*i+100*4+s2+1
c	  call showCBP(kf,kfc(kf))
	    end do
	  end do
	end do

* mesons:

c	print*,'mesons:'
	do s2=0,2,2
	  do j=1,3
	    do i=1,3
	      kf=(100*maxCBP(i,j)+10*minCBP(i,j)+s2+1)*
     .           sign(1,i-j)*(-1)**maxCBP(i,j)
	      if (i.eq.j) kf=(-1)**i*kf
	      kfc(kf)=kf
	      if (i.eq.2) 
     .          kfc(kf)=(100*maxCBP(4,j)+10*minCBP(4,j)+s2+1)*
     .          sign(1,4-j)*(-1)**maxCBP(4,j)
	      if (j.eq.2.and.i.ne.2) 
     .           kfc(kf)=(100*maxCBP(i,4)+10*minCBP(i,4)+s2+1)*
     .           sign(1,i-4)*(-1)**maxCBP(i,4)
c	  call showCBP(kf,kfc(kf))
	    end do
	  end do
	end do

* baryons:

c	print*,'baryons:'
	do s2=1,3,2
	  do k=1,3
	    do j=k,3
	       do i=j,3
		 if (i.eq.3.and.j.eq.2.and.k.eq.1.and.s2.eq.1) then
		   kf=1000*3+100*2+10*1+s2+1
	           kfc(kf)=1000*4+100*3+10*1+s2+1
c		   call showCBP(kf,kfc(kf))
		   kf=1000*3+100*1+10*2+s2+1
	           kfc(kf)=1000*4+100*1+10*3+s2+1
		   goto 1
		 end if
		 if (i.eq.3.and.j.eq.2.and.k.eq.1.and.s2.eq.2) then
		   kf=1000*3+100*2+10*1+s2+1
	           kfc(kf)=1000*4+100*3+10*1+s2+1
c		   call showCBP(kf,kfc(kf))
		   kf=1000*3+100*1+10*2+s2+1
	           kfc(kf)=1000*4+100*1+10*3+s2+1
		   goto 1
		 end if
		 if (i.eq.3.and.j.eq.2.and.k.eq.2.and.s2.eq.1) then
		   kf=1000*3+100*2+10*2+s2+1
	           kfc(kf)=1000*4+100*3+10*2+s2+1
		   goto 1
		 end if
	         kf=1000*i+100*j+10*k+s2+1
	         kfc(kf)=kf
	         if(i.eq.2.and.j.eq.2.and.s2.eq.1) then
		   kfc(kf)=1000*4+100*k+10*2+s2+1
		   goto 1
		 end if 
		 if(i.eq.2) then
		   kfc(kf)=1000*4+100*j+10*k+s2+1
		   goto 1
		 end if 
		 if(j.eq.2) then
		   kfc(kf)=1000*4+100*i+10*k+s2+1
		   goto 1
		 end if 
		 if(k.eq.2) then
		   kfc(kf)=1000*4+100*i+10*j+s2+1
		   goto 1
		 end if 
1	  	 continue
c	  CALL showCBP(kf,kfc(kf))
	        end do
	    end do
	  end do
	end do
C.. To replace pip and rho0 as well:
	kfc(111)=421
	kfc(113)=423

	end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

	integer function maxCBP(i,j)

	  maxCBP=i
	  if (j.gt.i) maxCBP=j
	
	end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

	integer function minCBP(i,j)

	  minCBP=i
	  if (j.lt.i) minCBP=j
	
	end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

	subroutine showCBP(kf1,kf2)
	character cha1*16,cha2*16

	CALL PYNAME(kf1,cha1)
	CALL PYNAME(kf2,cha2)

c	if (kf1.ne.kf2) print*,kf1,cha1,'->',kf2,cha2

	end

*=======================================================================

      DOUBLE PRECISION FUNCTION GXCBP(A)
c  function to generate the energy from projectile parton of a 
c  sea qq_bar pair (including normalization)
      DOUBLE PRECISION Y,YR,PYR
      EXTERNAL PYR
      
 10   CONTINUE
      GXCBP = PYR(0)
      YR = PYR(0)*1000.D0
      Y = 4.77497509534447750D-02*(1.D0-GXCBP)**4/(GXCBP+1.0D-10)
      IF ( YR.GT.Y ) GOTO 10

      RETURN  
      END


