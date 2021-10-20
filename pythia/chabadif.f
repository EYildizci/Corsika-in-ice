      SUBROUTINE CHABADIF(EC)
c  charmed baryon diffractive process

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER PYCOMP
      CHARACTER CHAF*16

      DOUBLE PRECISION m1,m2,m12
      DOUBLE PRECISION MP,M0,DM,MC,MLAMBDAC,MCHARM
      INTEGER BINES, LH(25)
      CHARACTER CHANNEL(0:25)*30
      CHARACTER CHAU*16

      INTEGER IKC(2)
      INTEGER KFC(-5000:5000)
      COMMON/CONVERTER/KFC

      COMMON/ENER/MP,EMIN,EP,ELAB,PXLAB,PYLAB,PZLAB
      COMMON/BINES/BIN(0:25,0:140),NBIN
      COMMON/EMAXLH/EMAXB,EMAXM,IMAXB,IMAXM,IKMAXB,IKMAXM

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT4/CHAF(500,2)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDATR/MRPY(6),RRPY(100)

      
      NBIN=140
      EMIN=0.1D0
      DO I=1,25
         LH(I)=0
      END DO
      DO I=1,NBIN
         BIN(19,I)=0.D0
      END DO


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
C     SIGMAC++
      KC = PYCOMP(4222)
      MDCY(KC,1) = 1
C     SIGMAC+
      KC = PYCOMP(4212)
      MDCY(KC,1) = 1
C     SIGMAC0
      KC = PYCOMP(4112)
      MDCY(KC,1) = 1
C     A-SIGMAC++
      KC = PYCOMP(-4222)
      MDCY(KC,1) = 1
C     A-SIGMAC+
      KC = PYCOMP(-4212)
      MDCY(KC,1) = 1
C     A-SIGMAC0
      KC = PYCOMP(-4112)
      MDCY(KC,1) = 1
C     ETAC
      KC = PYCOMP(441)
      MDCY(KC,1) = 1
c$$$  C                ETAC
c$$$  KC = PYCOMP(-441)
c$$$  MDCY(KC,1) = 1

C-------------------------
      NEV=1
      
C     SELECT TYPE OF EVENTS TO BE GENERATED
      
      MSEL=0


      DO IE=1,500
         MSUB(IE)=0
      ENDDO

*     Soft QCD processes (excluding elastic):
*----------------------------------------
*     92: Single diffractive (XB)
      MSUB(92)=1
*     93: Single diffractive (AX)
      MSUB(93)=1
*     94: Double  diffractive
      MSUB(94)=1
      
      
      KC=PYCOMP(2212)
      MP=PMAS(KC,1)
      KC=PYCOMP(4122)
      MLAMBDAC=PMAS(KC,1)
      MCHARM=1.27D0

      EP=(MLAMBDAC-MCHARM)/MLAMBDAC*EC
      ECM=DSQRT(2.D0*MP*EP)

C..   To stop hadronization:
      MSTP(111)=0
C-----------------------
      CALL PYINIT('CMS','P+','P+',ECM)

C..   To init u to c conversion:
      CALL KFCINITCBD

C     EVENT LOOP


      iuu1=0
      iud0=0
      iud1=0
      icu0=0
      icu1=0
      icd0=0
      icd1=0



 777  CONTINUE
      icase=0
      iC=0
      ig=0
      CALL PYEVNT

c     CALL PYLIST(1)


      EMAX=0.D0
      EMAXM=0.D0
      EMAXB=0.D0
      IH=0


c..   for diffractive:
      DO IK=5,N
         if(k(ik,2).eq.21) ig=ik
         IF (K(IK,3).EQ.3.AND.K(IK,2).NE.21) THEN
            IF(K(IK,2).NE.KFC(K(IK,2))) THEN
               IC=IC+1
               IKC(IC)=IK
            END IF
            if (K(IK,2).eq.2203) iuu1=iuu1+1
            if (K(IK,2).eq.2101) iud0=iud0+1
            if (K(IK,2).eq.2103) iud1=iud1+1

         END IF

      END DO
      
      IF ( IC .EQ. 0) GOTO 777
      
C..   To modify an entry to replace charm:
c..   random:
      ICR=INT(IC*PYR(0))+1
c     PRINT*,'ICR: ',ICR
      IK=IKC(ICR)

      if (P(3,5).lt.1.938D0.and.ic.eq.2) then
         IK=IKC(1)
         if (K(IKC(2),2).gt.1000) IK=IKC(2)
      end if

      K(IK,2)=KFC(K(IK,2))
      if (K(IK,2).eq.4203) then
         ii=INT(4*PYR(0))+1
         if (ii.le.3) K(IK,2)=4201
      end if
      KC=PYCOMP(K(IK,2))
      M0=PMAS(KC,1)

c     PRINT*,'ENTRY ',IK,' WILL BE REPLACED'

      if (K(IK,2).eq.4201) icu0=icu0+1
      if (K(IK,2).eq.4203) icu1=icu1+1
      if (K(IK,2).eq.4101) icd0=icd0+1
      if (K(IK,2).eq.4103) icd1=icd1+1

      EPCOR=ECM/2d0*MLAMBDAC/(MLAMBDAC-MCHARM)
      
*=======================================================================

      if (P(3,5).gt.mp.and.P(3,5).lt.1.938D0) then

         icase=2

         IK1=IK

         IK2=6
         if(IK.eq.6) IK2=5
         m2=P(IK2,5)

*     conversion exceptions:
         IF (K(IK1,2).EQ.4122.AND.K(IK2,2).EQ.111) THEN
            K(IK1,2)=4212
         END IF
         IF (K(IK1,2).EQ.4114) K(IK1,2)=4112
         IF (K(IK1,2).EQ.4214) K(IK1,2)=4212
         IF (K(IK1,2).EQ.4224) K(IK1,2)=4222

         KC=PYCOMP(K(IK1,2))
         M0=PMAS(KC,1)
         m1=m0

         d0=EPCOR+P(3,4)-P(1,4)	
         d1=P(3,1)-P(1,1)
         d2=P(3,2)-P(1,2)
         d3=EPCOR*DSQRT(1.D0-MLAMBDAC**2/EPCOR**2)+P(3,3)-P(1,3)
         dm=dsqrt(d0**2-d1**2-d2**2-d3**2)

         q0=P(3,4)-P(1,4)
         q1=P(3,1)-P(1,1)
         q2=P(3,2)-P(1,2)
         q3=P(3,3)-P(1,3)

         m12=m1+m2
         if (m12.gt.dm) then
            dm=(m1+m2)*1.00001D0
         end if

         bx=d1/d0
         by=d2/d0
         bz=d3/d0
         bt=dsqrt(d1**2+d2**2+d3**2)/d0
         gm=d0/dm

         e1=(dm**2+m1**2-m2**2)/2.d0/dm
         e2=(dm**2+m2**2-m1**2)/2.d0/dm
         pp=dsqrt((dm**2-(m1+m2)**2)*(dm**2-(m1-m2)**2))/2.d0/dm

         P(IK1,1)=gm*bx*e1+(gm-1.d0)*bx*bz/bt**2*pp
         P(IK1,2)=gm*by*e1+(gm-1.d0)*by*bz/bt**2*pp
         P(IK1,3)=gm*bz*e1+(1.d0+(gm-1.d0)*bz**2/bt**2)*pp
         P(IK1,4)=gm*e1   +gm*bz*pp
         P(IK1,5)=m1

         P(IK2,1)=gm*bx*e2-(gm-1.d0)*bx*bz/bt**2*pp
         P(IK2,2)=gm*by*e2-(gm-1.d0)*by*bz/bt**2*pp
         P(IK2,3)=gm*bz*e2-(1.d0+(gm-1.d0)*bz**2/bt**2)*pp
         P(IK2,4)=gm*e2   -gm*bz*pp
         P(IK2,5)=m2

         P(3,1)=d1
         P(3,2)=d2
         P(3,3)=d3
         P(3,4)=d0
         P(3,5)=dm

         P(1,3)=EPCOR*DSQRT(1.D0-MLAMBDAC**2/EPCOR**2)
         P(1,4)=EPCOR
         P(1,5)=MLAMBDAC

         goto 11

      end if

*=======================================================================

      if (K(3,2).eq.2212) then

         icase=1

         IK=5
         K(IK,2)=4122
         P(IK,3)=EPCOR+P(2,3)-P(4,3)
         P(IK,4)=DSQRT(P(IK,1)**2+P(IK,2)**2+P(IK,3)**2+MLAMBDAC**2)
         P(IK,5)=MLAMBDAC

         CALL PYNAME(K(IK,2),CHAU)
         CALL PYEXEC
         goto 22

      end if

*=======================================================================

      IF ( IC .EQ. 1) THEN
         GOTO 777
      END IF


      ICO=3-ICR
c     PRINT*,'ICO: ',ICO
      d0=EPCOR+P(3,4)-P(1,4)	
      d1=      P(3,1)-P(1,1)
      d2=      P(3,2)-P(1,2)
      d3=EPCOR*DSQRT(1.D0-MLAMBDAC**2/EPCOR**2)+P(3,3)-P(1,3)
      dm=dsqrt(d0**2-d1**2-d2**2-d3**2)

      q0=P(3,4)-P(1,4)
      q1=P(3,1)-P(1,1)
      q2=P(3,2)-P(1,2)
      q3=P(3,3)-P(1,3)
      
      IK=IKC(ICR)
      m1=M0
      x1=P(IK,1)
      y1=P(IK,2)

      IK=IKC(ICO)
      m2=P(IK,5)
      x2=P(IK,1)
      y2=P(IK,2)
      
      z3=0.d0
      e3=0.d0
      if(ig.ne.0) then
         z3=P(ig,3)*MCHARM/MLAMBDAC
         e3=P(ig,4)*MCHARM/MLAMBDAC
         P(ig,3)=z3
         P(ig,4)=e3
      end if


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

      P(3,3)=d3
      P(3,4)=d0
      P(3,5)=dsqrt(P(3,4)**2-P(3,1)**2-P(3,2)**2-P(3,3)**2)

      IK=IKC(ICR)
      P(IK,3)=z1
      P(IK,4)=e1
      P(IK,5)=m0

      IK=IKC(ICO)
      P(IK,3)=z2
      P(IK,4)=e2

      P(1,3)=EPCOR*DSQRT(1.D0-MLAMBDAC**2/EPCOR**2)
      P(1,4)=EPCOR
      P(1,5)=MLAMBDAC

*=======================================================================
 11   continue

C..   To restart event with hadronization:
      CALL PYEXEC
c     CALL PYLIST(1)
C------------------------


 22   continue


      EMAX=0.D0
      EAUX = 0.D0
      EAUX2 = 0.D0
      DO IK=3,N
         IF ( K(IK,1) .EQ. 1) THEN
            CALL BOOSTCBD(IK)
            EAUX = EAUX + P(IK,4)
         ENDIF
      ENDDO

C     IF THE SUM OF THE OUTCOMING PARTS IS .GT. THE INITIAIL E (SOMEHOW)
c     SHOULD NOT HAPPEN BUT SOMETIMES EXCEEDS BY A SMALL AMOUNT
      IF ( EAUX .GT. EC) THEN

         DO IK=3,N
            IF ( K(IK,1) .EQ. 1) THEN
C     RESCALE ENERGY AND MOMENTUM / MASS IS THE SAME
               KC=PYCOMP(K(IK,2))
               
               P(IK,4) = P(IK,4)*EC/EAUX
               P(IK,3) = DSQRT(P(IK,4)**2-P(IK,1)**2-
     .              P(IK,2)**2-PMAS(KC,1)**2)
               EAUX2 = EAUX2 + P(IK,4)
            ENDIF
         ENDDO


      ENDIF


      END

*=======================================================================
*     2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      SUBROUTINE BOOSTCBD(IK)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 MP

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/ENER/MP,EMIN,EP,ELAB,PXLAB,PYLAB,PZLAB

      GAMMA=DSQRT(EP/(2.D0*MP))
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
*     2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      SUBROUTINE KFCINITCBD

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      integer i,j,k,s2

      INTEGER KFC(-5000:5000)
      COMMON/CONVERTER/KFC

*     quarks:

c     print*,'quarks:'
      do i=1,3
         kf=i
         kfc(kf)=kf
         if (i.eq.2) kfc(kf)=4
      end do

*     diquarks:

c     print*,'diquarks:'
      do s2=0,2,2
         do j=1,3
	    do i=j,3
               kf=1000*i+100*j+s2+1
               kfc(kf)=kf
               if (i.eq.2) kfc(kf)=1000*4+100*j+s2+1
               if (j.eq.2.and.i.ne.2) kfc(kf)=1000*i+100*4+s2+1
	    end do
         end do
      end do
      

*     mesons:

c     print*,'mesons:'
      do s2=0,2,2
         do j=1,3
	    do i=1,3
               kf=(100*maxCBD(i,j)+10*minCBD(i,j)+s2+1)*
     .              sign(1,i-j)*(-1)**maxCBD(i,j)
               if (i.eq.j) kf=(-1)**i*kf
               kfc(kf)=kf
               if (i.eq.2) 
     .              kfc(kf)=(100*maxCBD(4,j)+10*minCBD(4,j)+s2+1)*
     .              sign(1,4-j)*(-1)**maxCBD(4,j)
               if (j.eq.2.and.i.ne.2) 
     .              kfc(kf)=(100*maxCBD(i,4)+10*minCBD(i,4)+s2+1)*
     .              sign(1,i-4)*(-1)**maxCBD(i,4)
	    end do
         end do
      end do

*     baryons:

c     print*,'baryons:'
      do s2=1,3,2
         do k=1,3
	    do j=k,3
	       do i=j,3
                  if (i.eq.3.and.j.eq.2.and.k.eq.1.and.s2.eq.1) then
                     kf=1000*3+100*2+10*1+s2+1
                     kfc(kf)=1000*4+100*3+10*1+s2+1
                     kf=1000*3+100*1+10*2+s2+1
                     kfc(kf)=1000*4+100*1+10*3+s2+1
                     goto 1
                  end if
                  if (i.eq.3.and.j.eq.2.and.k.eq.1.and.s2.eq.2) then
                     kf=1000*3+100*2+10*1+s2+1
                     kfc(kf)=1000*4+100*3+10*1+s2+1
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

	integer function maxCBD(i,j)

	  maxCBD=i
	  if (j.gt.i) maxCBD=j
	
	end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

	integer function minCBD(i,j)

	  minCBD=i
	  if (j.lt.i) minCBD=j
	
	end

