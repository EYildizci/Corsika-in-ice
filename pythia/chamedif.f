      SUBROUTINE CHAMEDIF(EC)
c  charmed meson diffractive process

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER PYCOMP
      CHARACTER CHAF*16

      DOUBLE PRECISION m1,m2,m12
      DOUBLE PRECISION MP,M0,DM,MC,MD,MCHARM,MPI
      CHARACTER CHANNEL(0:25)*30
      CHARACTER CHAU*16

      INTEGER IKC(2)
      INTEGER KFC(-5000:5000)
      COMMON/CONVERTER/KFC

      COMMON/MYVAR/ EPCOR,PXC,PYC,PZC,PXP,PYP,PZP
      COMMON/ENER/MP,EMIN,EP,ELAB,PXLAB,PYLAB,PZLAB
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

      EMIN=0.1D0

      PRINT*, 'N ',N

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
      
      DO IE=1,500
         IF(MSUB(IE).EQ.1) PRINT*,'MSUB(',IE,')=',MSUB(IE)
      ENDDO

C     INITIALIZE PYTHIA 

      KC=PYCOMP(2212)
      MP=PMAS(KC,1)
      KC=PYCOMP(411)
      MD=PMAS(KC,1)
      KC=PYCOMP(211)
      MPI=PMAS(KC,1)
      MCHARM=1.27D0

      EP=(MD-MCHARM)/MD*EC
      ECM=DSQRT(2.D0*MP*EP)

C..   To stop hadronization:
      MSTP(111)=0
C-----------------------
      CALL PYINIT('CMS','PI+','P+',ECM)
      PRINT*,'MSUB AFTER PYINIT'
      DO IE=1,500
         IF(MSUB(IE).EQ.1) PRINT*,'MSUB(',IE,')=',MSUB(IE)
      ENDDO

C..   To init u to c conversion:
      CALL KFCINITCMD

C     EVENT LOOP

      NEVC=0

 777  CONTINUE
      icase=0
      IC=0
      ig=0
      CALL PYEVNT
c      CALL PYLIST(1)

      PRINT*,'MSUB AFTER PYEVNT'
      DO IE=1,500
         IF(MSUB(IE).EQ.1) PRINT*,'MSUB(',IE,')=',MSUB(IE)
      ENDDO


      EMAX=0.D0
      EMAXM=0.D0
      EMAXB=0.D0
      IH=0



c..for diffractive:
      DO IK=5,N
         if(k(ik,2).eq.21) ig=ik
         IF (K(IK,3).EQ.3.AND.K(IK,2).NE.21) THEN
            IF(K(IK,2).NE.KFC(K(IK,2))) THEN
c..   pi0 and rho0 replaced only half of the times:
               IF (K(IK,2).EQ.111.OR.K(IK,2).EQ.113) THEN
                  IF (PYR(0).GT.0.5D0) GOTO 33
               END IF
               IC=IC+1
               IKC(IC)=IK
 33            CONTINUE
            END IF
            if (K(IK,2).EQ.-1) IKC(2)=IK
            CALL PYNAME(K(IK,2),CHAU)

         END IF
         
      END DO

      IF ( IC .EQ. 0 ) GOTO 777

C..   To modify an entry to replace charm:
c..   random:
      ICR=INT(IC*PYR(0))+1
      IK=IKC(ICR)
      K(IK,2)=KFC(K(IK,2))
      KC=PYCOMP(K(IK,2))
      M0=PMAS(KC,1)
      PRINT*,'ENTRY ',IK,' WILL BE REPLACED'
      EPCOR=ECM/2.D0*MD/(MD-MCHARM)
      
*=======================================================================

      if (P(3,5).gt.mpi.and.P(3,5).lt.1.13957D0) then

         icase=2

         PRINT*,'icase ',icase
         IK1=IK

         IK2=6
         if(IK.eq.6) IK2=5
         m2=P(IK2,5)

c     print*,K(IK1,2),K(IK2,2)


         KC=PYCOMP(K(IK1,2))
         M0=PMAS(KC,1)
         m1=m0


         d0=EPCOR+P(3,4)-P(1,4)	
         d1=P(3,1)-P(1,1)
         d2=P(3,2)-P(1,2)
         d3=EPCOR*DSQRT(1.D0-MD**2/EPCOR**2)+P(3,3)-P(1,3)
         dm=dsqrt(d0**2-d1**2-d2**2-d3**2)

         PRINT*,'D0,D1,D2,D3,Dm: ',d0,' ',d1,' ',d2,' ',d3,' ',dm
         PRINT*,'EPCOR: ',EPCOR
         PRINT*,'M1,M2: ',m1,' ',m2
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

         e1=(dm**2+m1**2-m2**2)/2.D0/dm
         e2=(dm**2+m2**2-m1**2)/2.D0/dm
         pp=dsqrt((dm**2-(m1+m2)**2)*(dm**2-(m1-m2)**2))/2.D0/dm

         P(IK1,1)=gm*bx*e1+(gm-1.D0)*bx*bz/bt**2*pp
         P(IK1,2)=gm*by*e1+(gm-1.D0)*by*bz/bt**2*pp
         P(IK1,3)=gm*bz*e1+(1.D0+(gm-1.D0)*bz**2/bt**2)*pp
         P(IK1,4)=gm*e1   +gm*bz*pp
         P(IK1,5)=m1

         P(IK2,1)=gm*bx*e2-(gm-1.D0)*bx*bz/bt**2*pp
         P(IK2,2)=gm*by*e2-(gm-1.D0)*by*bz/bt**2*pp
         P(IK2,3)=gm*bz*e2-(1.D0+(gm-1.D0)*bz**2/bt**2)*pp
         P(IK2,4)=gm*e2   -gm*bz*pp
         P(IK2,5)=m2

         P(3,1)=d1
         P(3,2)=d2
         P(3,3)=d3
         P(3,4)=d0
         P(3,5)=dm

         P(1,3)=EPCOR*DSQRT(1.D0-MD**2/EPCOR**2)
         P(1,4)=EPCOR
         P(1,5)=MD

         goto 11

      end if

*=======================================================================

      if (K(3,2).eq.211) then

         icase=1


         IK=5
         K(IK,2)=411
         P(IK,3)=EPCOR+P(2,3)-P(4,3)
         P(IK,4)=DSQRT(P(IK,1)**2+P(IK,2)**2+P(IK,3)**2+MD**2)
         P(IK,5)=MD

         CALL PYNAME(K(IK,2),CHAU)
         CALL PYEXEC
         goto 22

      end if

*=======================================================================



      IF ( IC .EQ. 1 ) THEN
         GOTO 777
      END IF

      icase=3

      ICO=3-ICR

      d0=EPCOR+P(3,4)-P(1,4)	
      d1=      P(3,1)-P(1,1)
      d2=      P(3,2)-P(1,2)
      d3=EPCOR*DSQRT(1.D0-MD**2/EPCOR**2)+P(3,3)-P(1,3)
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
         z3=P(ig,3)*MCHARM/MD
         e3=P(ig,4)*MCHARM/MD
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
     -     2*d0**2*m1**2+2*d3**2*m1**2-4*d0*e3*m1**2+2*e3**2*m1**2   
     -     -m1**4+2*d0**2*m2**2-2*d3**2*m2**2-4*d0*e3*m2**2+ 
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

      P(1,3)=EPCOR*DSQRT(1.D0-MD**2/EPCOR**2)
      P(1,4)=EPCOR
      P(1,5)=MD

*=======================================================================
 11   continue

C..   To restart event with hadronization:
      CALL PYEXEC
c      CALL PYLIST(1)
C------------------------


 22   continue

      EAUX = 0.D0
      EAUX2 = 0.D0
      DO IK=3,N
         IF ( K(IK,1) .EQ. 1) THEN
            CALL BOOSTCMD(IK)
            EAUX = EAUX + P(IK,4)
         ENDIF
      ENDDO

C     IF THE SUM OF THE OUTCOMING PARTS IS .GT. THE INITIAIL E (SOMEHOW)
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
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

	SUBROUTINE BOOSTCMD(IK)

      	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 	DOUBLE PRECISION MP

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

c	print*,ELAB,PXLAB,PYLAB,PZLAB
c	print*,ELAB**2-PXLAB**2-PYLAB**2-PZLAB**2,
c     .       P(IK,4)**2-P(IK,1)**2-P(IK,2)**2-P(IK,3)**2,P(IK,5)**2

	RETURN
	END

*=======================================================================
*     2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      SUBROUTINE KFCINITCMD

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      integer i,j,k,s2

      INTEGER KFC(-5000:5000)
      COMMON/CONVERTER/KFC

*     quarks:

c     print*,'quarks:'
      do i=-3,3
         kf=i
         kfc(kf)=kf
         if (i.eq.2) kfc(kf)=4
         if (i.eq.-2) kfc(kf)=-4
c     call showCMD(kf,kfc(kf))
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
c     call showCMD(kf,kfc(kf))
	    end do
         end do
      end do

*     mesons:

c     print*,'mesons:'
      do s2=0,2,2
         do j=1,3
	    do i=1,3
               kf=(100*maxCMD(i,j)+10*minCMD(i,j)+s2+1)*
     .              sign(1,i-j)*(-1)**maxCMD(i,j)
               if (i.eq.j) kf=(-1)**i*kf
               kfc(kf)=kf
               if (i.eq.2) 
     .              kfc(kf)=(100*maxCMD(4,j)+10*minCMD(4,j)+s2+1)*
     .              sign(1,4-j)*(-1)**maxCMD(4,j)
               if (j.eq.2.and.i.ne.2) 
     .              kfc(kf)=(100*maxCMD(i,4)+10*minCMD(i,4)+s2+1)*
     .              sign(1,i-4)*(-1)**maxCMD(i,4)
c     call showCMD(kf,kfc(kf))
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
c     call showCMD(kf,kfc(kf))
                     kf=1000*3+100*1+10*2+s2+1
                     kfc(kf)=1000*4+100*1+10*3+s2+1
                     goto 1
                  end if
                  if (i.eq.3.and.j.eq.2.and.k.eq.1.and.s2.eq.2) then
                     kf=1000*3+100*2+10*1+s2+1
                     kfc(kf)=1000*4+100*3+10*1+s2+1
c     call showCMD(kf,kfc(kf))
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
c	  CALL showCMD(kf,kfc(kf))
	        end do
	    end do
	  end do
	end do
C.. To replace pip and rho0 as well:
	kfc(111)=421
	kfc(113)=423

	end

*=======================================================================
*     2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      integer function maxCMD(i,j)

      maxCMD=i
      if (j.gt.i) maxCMD=j
      
      end

*=======================================================================
*     2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      integer function minCMD(i,j)

      minCMD=i
      if (j.lt.i) minCMD=j
      
      end

*=======================================================================
*     2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      subroutine showCMD(kf1,kf2)
      character cha1*16,cha2*16

      CALL PYNAME(kf1,cha1)
      CALL PYNAME(kf2,cha2)

      if (kf1.ne.kf2) print*,kf1,cha1,'->',kf2,cha2

      end

