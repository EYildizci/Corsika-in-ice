      SUBROUTINE BOBAPAR(EC,NY)
c  bottom baryon partonic process

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER PYCOMP
      CHARACTER CHAF*16

cdh   real*8 m1,m2,m12
      DOUBLE PRECISION m1,m2,m12
cdh   REAL*8 MP,M0,DM,MC,MSIGMAB,MBOTTOM
      DOUBLE PRECISION MP,M0,DM,MC,MSIGMAB,MBOTTOM
      INTEGER COL
      CHARACTER CHANNEL(0:23)*30
      CHARACTER CHAU*16

      CHARACTER LTAR*5
      INTEGER IKB(200)
      INTEGER IKM(6)

      INTEGER KFB(-5000:5000)
      INTEGER KCODES(0:21)
      INTEGER KCODES_SUB(0:21)
      INTEGER KSAMECODES(0:17)


      DOUBLE PRECISION A
      DOUBLE PRECISION XGLIST(100),ECMLIST(100),QMASS(100)
      INTEGER KFLIST(100)
      INTEGER KLIST(10000,5)
      CHARACTER TLIST(100)*5
      DOUBLE PRECISION PLIST(10000,5),VLIST(10000,5)

      COMMON/CONVERTERBBP/KFB,KCODES,KCODES_SUB,KSAMECODES
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


C     CHARMED HADRONS
C-------------------------
C     D+
      KC=PYCOMP(411)
      MDCY(KC,1)=0
C     D0
      KC=PYCOMP(421)
      MDCY(KC,1)=0
C     DS+
      KC=PYCOMP(431)
      MDCY(KC,1)=0
C     LC
      KC=PYCOMP(4122)
      MDCY(KC,1)=0
C     XIC-
      KC=PYCOMP(4232)
      MDCY(KC,1)=0
C     XIC0
      KC=PYCOMP(4132)
      MDCY(KC,1)=0
C     OMEGAC
      KC=PYCOMP(4332)
      MDCY(KC,1)=0

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
      MDCY(KC,1)=0
C     XI_C_0 BAR
      KC=PYCOMP(-4132)
      MDCY(KC,1)=0
C     OMEGA_C_0 BAR
      KC=PYCOMP(-4332)
      MDCY(KC,1)=0

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
c$$$C                ETAC
c$$$      KC = PYCOMP(-441)
c$$$      MDCY(KC,1) = 1

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
C     K*(892)0
      KC=PYCOMP(313)
      MDCY(KC,1)=1
C     K*(892)0 bar
      KC=PYCOMP(-313)
      MDCY(KC,1)=1

C-------------------------

C     BOTTOM BARYONS
C-------------------------
C     B+
      KC=PYCOMP(521)
      MDCY(KC,1)=0
C     B0
      KC=PYCOMP(511)
      MDCY(KC,1)=0
C     BS+
      KC=PYCOMP(531)
      MDCY(KC,1)=0
C     LB0
      KC=PYCOMP(5122)
      MDCY(KC,1)=0
C     XIB-
      KC=PYCOMP(5132)
      MDCY(KC,1)=1
C     XIB0
      KC=PYCOMP(5232)
      MDCY(KC,1)=1
C     OMEGAB
      KC=PYCOMP(5332)
      MDCY(KC,1)=1

C     B- BAR
      KC=PYCOMP(-521)
      MDCY(KC,1)=0
C     B0 BAR
      KC=PYCOMP(-511)
      MDCY(KC,1)=0
C     BS- BAR
      KC=PYCOMP(-531)
      MDCY(KC,1)=0
C     LB0 BAR
      KC=PYCOMP(-5122)
      MDCY(KC,1)=0
C     XIB+ BAR
      KC=PYCOMP(-5132)
      MDCY(KC,1)=1
C     XIB0 BAR
      KC=PYCOMP(-5232)
      MDCY(KC,1)=1
C     OMEGAB BAR
      KC=PYCOMP(-5332)
      MDCY(KC,1)=1

C     FORCED DECAY
C     SIGMAB+
      KC = PYCOMP(5222)
      MDCY(KC,1) = 1
C     SIGMAB0
      KC = PYCOMP(5212)
      MDCY(KC,1) = 1
c     SigmaB-
      KC = PYCOMP(5112)
      MDCY(KC,1) = 1
C     A-SIGMAB0
      KC = PYCOMP(-5212)
      MDCY(KC,1) = 1
C     A-SIGMAB-
      KC = PYCOMP(-5222)
      MDCY(KC,1) = 1
c     A-SigmaB+
      KC = PYCOMP(-5112)
      MDCY(KC,1) = 1
c$$$C                ETAB
c$$$      KC = PYCOMP(551)
c$$$      MDCY(KC,1) = 1
c$$$C                ANTI-ETAB
c$$$      KC = PYCOMP(-551)
c$$$      MDCY(KC,1) = 1

C-------------------------

      
      MSEL=0
      DO IE=1,500
         MSUB(IE)=0
      ENDDO

* HARD QCD PROCESSES:
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
*     95: Low-pT scattering	(KEEP ALWAYS ACTIVE)
      MSUB(95)=1


      A = 1.D0
      ZB = 0.D0

      KC=PYCOMP(2212)
      MP=PMAS(KC,1)

      KC=PYCOMP(5222)
      MSIGMAB=PMAS(KC,1)

      MBOTTOM=4.20D0

      EP = EC

C     ASSIGN SEA PAIRS ENERGIES

      NP = 7
      NN = 7
      FAIL = 0
      IF (NY.GT.1) THEN
 123     CONTINUE
         FAIL = 0.D0
         ZB = 0.D0
         DO 124 INW=1,NY-1
            XGLIST(INW) = GXBBP(1.0)
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

            IF (PYR(0). LT. 1.D0*NP/(1.D0*NP+1.D0*NN)) THEN
               TLIST(INW) = 'P+'
               NP = NP-1
            ELSE
               TLIST(INW) = 'n0'
               NP = NP-1
            ENDIF


 124     CONTINUE
         IF (ZB.GE.1.D0) GOTO 123
         IF (DSQRT(2*MP*(MP+((MSIGMAB-MBOTTOM)/MSIGMAB*(EP*(1.D0-ZB)))))
     *       .LE. PARP(2)) GOTO 123

      ENDIF


      EPR = EP*(1.D0-ZB)
C     REMOVE QQBAR PAIRS ENERGY: ENERGY FOR QCD COLLISION

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

            DO IPART=9,N
               IF (K(IPART,1).EQ.1) THEN
                  CALL BOOSTBBP(IPART,EP*XGLIST(INW),MP)
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

      EP=(MSIGMAB-MBOTTOM)/MSIGMAB*EPR
      ECM=DSQRT(2.D0*MP*(MP+EP))
      EPCOR=ECM/2.D0*MSIGMAB/(MSIGMAB-MBOTTOM)

      DELTAE=EPCOR-ECM/2.D0
      DELTAPZ=EPCOR*DSQRT(1.D0-MSIGMAB**2/EPCOR**2)
     .     -ECM/2.D0*DSQRT(1.D0-(2D0*MP)**2/ECM**2)

C..   TO STOP HADRONIZATION:
      MSTP(111)=0
C-----------------------
      CALL PYINIT('CMS','P+',LTAR,ECM)

C..   To init u to c conversion:
      CALL KFBINITBBP

C     EVENT LOOP

 777  CONTINUE

      IB=0
      ig=0
      CALL PYEVNT
c      CALL PYLIST(1)
      

      EMAX=0.D0
      EMAXM=0.D0
      EMAXB=0.D0
      IH=0


      IK1=9
      PMAX=0D0
      DO IK=9,N
c     si no es un foton
         IF (K(IK,2).NE.22) THEN
c     si no es un gluon y es particula
            IF (K(IK,2).NE.21.and.K(IK,2).GT.0) THEN
c     si es sustituida y pz>0
               IF(K(IK,2).NE.KFB(K(IK,2)).AND.P(IK,3).GT.0d0) THEN
                  IB=IB+1
                  IKB(IB)=IK
                  IF (P(IK,3).GT.PMAX) THEN
                     PMAX=P(IK,3)
                     IKMAX=IK
                     IBMAX=IB
                  END IF
c     PRINT*,IK,P(IK,3)
               END IF
            END IF
         END IF
         
      END DO
      
      IF ( IB .EQ. 0) GOTO 777
c     if no candidates, skipt to next event


      
C..   To modify an entry to replace charm:
c..   random:
      IBR=INT(IB*PYR(0))+1
      IK=IKB(IBR)
      
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

c     bd_0 [5101]:(1/4) and bd_1[5103]:(3/4)
     
      if (K(IK,2).eq.5103) then
         ii=INT(4*PYR(0))+1
         if (ii.eq.1) K(IK,2)=5101
      end if

      IF (IK.EQ.IK1) THEN
         IKO=IK2
      ELSE
         IKO=IK1
      END IF
      
c..........................................................
      IF (IK.EQ.9.AND.K(IK,1).EQ.1) THEN
         K(IK,2)=KFB(K(IK,2))
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

      K(IK,2)=KFB(K(IK,2))
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
      
      z3=0.D0
      e3=0.D0
      DO I=IK1+1,IK2-1
         P(I,1)=P(I,1)*MBOTTOM/MSIGMAB
         P(I,2)=P(I,2)*MBOTTOM/MSIGMAB
         P(I,3)=P(I,3)*MBOTTOM/MSIGMAB
         P(I,4)=P(I,4)*MBOTTOM/MSIGMAB
         z3=z3+P(I,3)
         e3=e3+P(I,4)
      END DO

c     math	Solve[Sqrt[x1^2 + y1^2 + z1^2 + m1^2] + 
c     math         Sqrt[x2^2 + y2^2 + (d3 - z1 - z3)^2 + m2^2] + e3 == d0, {z1}]

      z1a=(4*d0**2*d3-4*d3**3-8*d0*d3*e3+4*d3*e3**2+4*d3*m1**2-
     -     4*d3*m2**2+4*d3*x1**2-4*d3*x2**2+4*d3*y1**2-4*d3*y2**2- 
     -     4*d0**2*z3+12*d3**2*z3+8*d0*e3*z3-4*e3**2*z3-4*m1**2*z3+
     -     4*m2**2*z3-4*x1**2*z3+4*x2**2*z3-4*y1**2*z3+4*y2**2*z3- 
     -     12*d3*z3**2+4*z3**3-
     -     dsqrt((-4*d0**2*d3+4*d3**3+8*d0*d3*e3-
     -     4*d3*e3**2-4*d3*m1**2+4*d3*m2**2-4*d3*x1**2+4*d3*x2**2- 
     -     4*d3*y1**2+4*d3*y2**2+4*d0**2*z3-12*d3**2*z3-8*d0*e3*z3+
     -     4*e3**2*z3+4*m1**2*z3-4*m2**2*z3+4*x1**2*z3- 
     -     4*x2**2*z3+4*y1**2*z3-4*y2**2*z3+12*d3*z3**2-4*z3**3)**2- 
     -     4*(4*d0**2-4*d3**2-8*d0*e3+4*e3**2+8*d3*z3-4*z3**2)*
     -     (-d0**4+2*d0**2*d3**2-d3**4+4*d0**3*e3-4*d0*d3**2*e3-
     -     6*d0**2*e3**2+2*d3**2*e3**2+4*d0*e3**3-e3**4+ 
     -     2*d0**2*m1**2+2*d3**2*m1**2-4*d0*e3*m1**2+2*e3**2*m1**2-
     -     m1**4+   
     -     2*d0**2*m2**2-2*d3**2*m2**2-4*d0*e3*m2**2+ 
     -     2*e3**2*m2**2+2*m1**2*m2**2-m2**4+2*d0**2*x1**2+
     -     2*d3**2*x1**2- 
     -     4*d0*e3*x1**2+2*e3**2*x1**2-2*m1**2*x1**2+ 
     -     2*m2**2*x1**2-x1**4+2*d0**2*x2**2-2*d3**2*x2**2-
     -     4*d0*e3*x2**2+ 
     -     2*e3**2*x2**2+2*m1**2*x2**2-2*m2**2*x2**2+ 
     -     2*x1**2*x2**2-x2**4+2*d0**2*y1**2+2*d3**2*y1**2-
     -     4*d0*e3*y1**2+ 
     -     2*e3**2*y1**2-2*m1**2*y1**2+2*m2**2*y1**2- 
     -     2*x1**2*y1**2+2*x2**2*y1**2-y1**4+2*d0**2*y2**2-
     -     2*d3**2*y2**2-
     -     4*d0*e3*y2**2+2*e3**2*y2**2+2*m1**2*y2**2- 
     -     2*m2**2*y2**2+2*x1**2*y2**2-2*x2**2*y2**2+2*y1**2*y2**2-
     -     y2**4-
     -     4*d0**2*d3*z3+4*d3**3*z3+8*d0*d3*e3*z3- 
     -     4*d3*e3**2*z3-4*d3*m1**2*z3+4*d3*m2**2*z3-4*d3*x1**2*z3+
     -     4*d3*x2**2*z3-4*d3*y1**2*z3+4*d3*y2**2*z3+2*d0**2*z3**2- 
     -     6*d3**2*z3**2-4*d0*e3*z3**2+2*e3**2*z3**2+2*m1**2*z3**2-
     -     2*m2**2*z3**2+2*x1**2*z3**2-2*x2**2*z3**2+2*y1**2*z3**2- 
     -     2*y2**2*z3**2+4*d3*z3**3-z3**4)))/(2.*(4*d0**2-
     -     4*d3**2-8*d0*e3+4*e3**2+8*d3*z3-4*z3**2))

      z1b=(4*d0**2*d3-4*d3**3-8*d0*d3*e3+4*d3*e3**2+4*d3*m1**2-
     -     4*d3*m2**2+4*d3*x1**2-4*d3*x2**2+4*d3*y1**2-4*d3*y2**2- 
     -     4*d0**2*z3+12*d3**2*z3+8*d0*e3*z3-4*e3**2*z3-4*m1**2*z3+
     -     4*m2**2*z3-4*x1**2*z3+4*x2**2*z3-4*y1**2*z3+4*y2**2*z3- 
     -     12*d3*z3**2+4*z3**3+
     -     dsqrt((-4*d0**2*d3+4*d3**3+8*d0*d3*e3-
     -     4*d3*e3**2-4*d3*m1**2+4*d3*m2**2-4*d3*x1**2+4*d3*x2**2- 
     -     4*d3*y1**2+4*d3*y2**2+4*d0**2*z3-12*d3**2*z3-8*d0*e3*z3+
     -     4*e3**2*z3+4*m1**2*z3-4*m2**2*z3+4*x1**2*z3- 
     -     4*x2**2*z3+4*y1**2*z3-4*y2**2*z3+12*d3*z3**2-4*z3**3)**2- 
     -     4*(4*d0**2-4*d3**2-8*d0*e3+4*e3**2+8*d3*z3-4*z3**2)*
     -     (-d0**4+2*d0**2*d3**2-d3**4+4*d0**3*e3-4*d0*d3**2*e3-
     -     6*d0**2*e3**2+2*d3**2*e3**2+4*d0*e3**3-e3**4+ 
     -     2*d0**2*m1**2+2*d3**2*m1**2-4*d0*e3*m1**2+2*e3**2*m1**2-
     -     m1**4+   
     -     2*d0**2*m2**2-2*d3**2*m2**2-4*d0*e3*m2**2+ 
     -     2*e3**2*m2**2+2*m1**2*m2**2-m2**4+2*d0**2*x1**2+
     -     2*d3**2*x1**2- 
     -     4*d0*e3*x1**2+2*e3**2*x1**2-2*m1**2*x1**2+ 
     -     2*m2**2*x1**2-x1**4+2*d0**2*x2**2-2*d3**2*x2**2-
     -     4*d0*e3*x2**2+ 
     -     2*e3**2*x2**2+2*m1**2*x2**2-2*m2**2*x2**2+ 
     -     2*x1**2*x2**2-x2**4+2*d0**2*y1**2+2*d3**2*y1**2-
     -     4*d0*e3*y1**2+ 
     -     2*e3**2*y1**2-2*m1**2*y1**2+2*m2**2*y1**2- 
     -     2*x1**2*y1**2+2*x2**2*y1**2-y1**4+2*d0**2*y2**2-
     -     2*d3**2*y2**2-
     -     4*d0*e3*y2**2+2*e3**2*y2**2+2*m1**2*y2**2- 
     -     2*m2**2*y2**2+2*x1**2*y2**2-2*x2**2*y2**2+2*y1**2*y2**2-
     -     y2**4-
     -     4*d0**2*d3*z3+4*d3**3*z3+8*d0*d3*e3*z3- 
     -     4*d3*e3**2*z3-4*d3*m1**2*z3+4*d3*m2**2*z3-4*d3*x1**2*z3+
     -     4*d3*x2**2*z3-4*d3*y1**2*z3+4*d3*y2**2*z3+2*d0**2*z3**2- 
     -     6*d3**2*z3**2-4*d0*e3*z3**2+2*e3**2*z3**2+2*m1**2*z3**2-
     -     2*m2**2*z3**2+2*x1**2*z3**2-2*x2**2*z3**2+2*y1**2*z3**2- 
     -     2*y2**2*z3**2+4*d3*z3**3-z3**4)))/(2.*(4*d0**2-
     -     4*d3**2-8*d0*e3+4*e3**2+8*d3*z3-4*z3**2))

      z2a=d3-z3-z1a
      z2b=d3-z3-z1b
      
      
      e1a=dsqrt(m1**2+x1**2+y1**2+z1a**2)
      e2a=dsqrt(m2**2+x2**2+y2**2+z2a**2)
      e1b=dsqrt(m1**2+x1**2+y1**2+z1b**2)
      e2b=dsqrt(m2**2+x2**2+y2**2+z2b**2)

c     print*,'a) ',d0,'= ',e1a+e2a+e3,'e1 =',e1a,'e2 =',e2a
c     print*,'b) ',d0,'= ',e1b+e2b+e3,'e1 =',e1b,'e2 =',e2b

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

c     print	print*,'check a) ',dsqrt(x1**2+y1**2+z1a**2+m1**2)+ 
c     print     .         dsqrt(x2**2+y2**2+(d3-z1a-z3)**2+m2**2)+e3-d0
c     print	print*,'check b) ',dsqrt(x1**2+y1**2+z1b**2+m1**2)+ 
c     print     .         dsqrt(x2**2+y2**2+(d3-z1b-z3)**2+m2**2)+e3-d0

      P(IK,3)=z1
      P(IK,4)=e1
      P(IK,5)=m0

      P(IKO,3)=z2
      P(IKO,4)=e2

 77   P(1,3)=EPCOR*DSQRT(1.D0-MSIGMAB**2/EPCOR**2)
      P(1,4)=EPCOR
      P(1,5)=MSIGMAB

*=======================================================================
 11   CONTINUE

C..   To restart event with hadronization:
      CALL PYEXEC
C------------------------
      DO IPART=9,N
         IF (K(IPART,1).EQ.1) THEN
            CALL BOOSTBBP(IPART,EP,MP)

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
c         PRINT*,'REESCALED ',EAUX,'  TO  ',EC
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

	SUBROUTINE BOOSTBBP(IK,DE0,DM0)

      	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
cdh     REAL*8 MP
 	DOUBLE PRECISION MP

	COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
	COMMON/ENER/MP,EMIN,EP,ELAB,PXLAB,PYLAB,PZLAB

	GAMMA=DSQRT(DE0/(2.D0*DM0))
	BETA=DSQRT(1.D0-GAMMA**(-2))
	ELAB =GAMMA*(P(IK,4)+BETA*P(IK,3))
	PXLAB=P(IK,1)
	PYLAB=P(IK,2)
	PZLAB=GAMMA*(BETA*P(IK,4)+P(IK,3))

        P(IK,3) = PZLAB
        P(IK,4) = ELAB

c	print*,ELAB,PXLAB,PYLAB,PZLAB
c	print*,ELAB**2-PXLAB**2-PYLAB**2-PZLAB**2,
c     .       P(IK,4)**2-P(IK,1)**2-P(IK,2)**2-P(IK,3)**2,P(IK,5)**2

	RETURN
	END

*=======================================================================
*     2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      SUBROUTINE KFBINITBBP

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      INTEGER I,J,K,S2

      INTEGER KFB(-5000:5000)
      INTEGER KCODES(0:21)
      INTEGER KCODES_SUB(0:21)
      INTEGER KSAMECODES(0:17)
      COMMON/CONVERTERBBP/KFB,KCODES,KCODES_SUB,KSAMECODES

*     QUARKS:

C     PRINT*,'QUARKS:'
C     PRINT*,' '
      DO I=1,6
         KF=I
         KFB(KF)=KF
         IF (I.EQ.1) KFB(KF)=5
C     CALL SHOW(KF,KFB(KF))
      END DO

*     DIQUARKS:

C     PRINT*,' '
C     PRINT*,'DIQUARKS:'
C     PRINT*,' '
      DO S2=0,2,2
         DO J=1,4
	    DO I=J,4
               KF=1000*I+100*J+S2+1
               KFB(KF)=KF
               IF (I.EQ.1) KFB(KF)=1000*5+100*J+S2+1
               IF (J.EQ.1.AND.I.NE.1) KFB(KF)=1000*5+100*I+S2+1
C     CALL SHOW(KF,KFB(KF))
	    END DO
         END DO
      END DO

*     MESONS:

C     PRINT*,' '
C     PRINT*,'MESONS:'
C     PRINT*,' '
      DO S2=0,2,2
         DO J=1,4
	    DO I=1,4
               KF=(100*maxBBP(I,J)+10*minBBP(I,J)+S2+1)*
     .              SIGN(1,I-J)*(-1)**maxBBP(I,J)
               IF (I.EQ.J) KF=(-1)**I*KF
               KFB(KF)=KF
               IF (I.EQ.1) 
     .              KFB(KF)=-1*(100*5+10*J+S2+1)
               IF (J.EQ.1.AND.I.NE.1) 
     .              KFB(KF)=(100*5+10*I+S2+1)
               
C     CALL SHOW(KF,KFB(KF))
	    END DO
         END DO
      END DO

*     BARYONS:

C     PRINT*,' '
C     PRINT*,'BARYONS:'
C     PRINT*,' '


*     BARYONS:
      
      DATA KCODES/2212,2112,3122,3212,3112,3312,4122,4212,4112
     .     ,4132,4312,4412,2214,2114,1114,3214,3114,3314
     .     ,4214,4114,4314,4414/

      DATA KCODES_SUB/5222,5122,5232,5322,5312,5332,5242,5422,5412
     .     ,5342,5432,5442,5224,5214,5114,5324,5314,5334
     .     ,5424,5414,5434,5444/
      
      DATA KSAMECODES/2224,3222,3224,3322,3324,3334,4222,4224,4232,4322
     .     ,4324,4332,4334,4422,4424,4432,4434,4444/

      DO I=0,21
         KFB(KCODES(I))=KCODES_SUB(I)
         KFB(-KCODES(I))=-KCODES_SUB(I)
C     CALL SHOW(KCODES(I),KFB(KCODES(I)))
      END DO

      DO I=0,17
C     PRINT*,I," ",KSAMECODES(I)
         
         KFB(KSAMECODES(I))=KSAMECODES(I)
         KFB(-KSAMECODES(I))=-KSAMECODES(I)
C     CALL SHOW(KSAMECODES(I),KFB(KSAMECODES(I)))
      END DO
      
      KFB(111)=511
      KFB(113)=513
      
      END

*=======================================================================
*     2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      integer function maxBBP(i,j)

      maxBBP=i
      if (j.gt.i) maxBBP=j
      
      end

*=======================================================================
*     2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      integer function minBBP(i,j)

      minBBP=i
      if (j.lt.i) minBBP=j
      
      end

*=======================================================================
*     2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      subroutine showBBP(kf1,kf2)
      character cha1*16,cha2*16

      CALL PYNAME(kf1,cha1)
      CALL PYNAME(kf2,cha2)

c      if (kf1.ne.kf2) print*,kf1,cha1,'->',kf2,cha2

      end

*=======================================================================

      DOUBLE PRECISION FUNCTION GXBBP(A)
c  function to generate the energy from projectile parton of a 
c  sea qq_bar pair (including normalization)
      DOUBLE PRECISION Y,YR,PYR
      EXTERNAL PYR

 10   CONTINUE
      GXBBP = PYR(0)
      YR = PYR(0)*1000.D0
      Y = 4.77497509534447750D-02*(1.D0-GXBBP)**4/(GXBBP+1.0D-10)
      IF ( YR.GT.Y ) GOTO 10

      RETURN  
      END


