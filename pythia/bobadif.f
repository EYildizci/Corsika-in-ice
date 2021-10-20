      SUBROUTINE BOBADIF(EC)
c  bottom baryon diffractive process
      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER PYCOMP
      CHARACTER CHAF*16

cdh   real*8 m1,m2,m12
      DOUBLE PRECISION m1,m2,m12
cdh   REAL*8 MP,M0,DM,MC,msigmab,mbottom
      DOUBLE PRECISION MP,M0,DM,MC,msigmab,mbottom
      INTEGER BINES, LH(25)
      CHARACTER CHANNEL(0:27)*30
      CHARACTER CHAU*16

      INTEGER IKB(2)

      INTEGER KFB(-5000:5000)
      INTEGER KCODES(0:21)
      INTEGER KCODES_SUB(0:21)
      INTEGER KSAMECODES(0:17)

      COMMON/CONVERTERBBD/KFB,KCODES,KCODES_SUB,KSAMECODES
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
      MDCY(KC,1)=1
C     XIC0
      KC=PYCOMP(4132)
      MDCY(KC,1)=1
C     OMEGAC
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
      KC = PYCOMP(-2214)
      MDCY(KC,1) = 1
C-------------------------
C     K*(892)0
      KC=PYCOMP(313)
      MDCY(KC,1)=1
C     K*(892)0 bar
      KC=PYCOMP(-313)
      MDCY(KC,1)=1

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
C     A-SIGMAB0
      KC = PYCOMP(-5212)
      MDCY(KC,1) = 1
C     A-SIGMAB-
      KC = PYCOMP(-5222)
      MDCY(KC,1) = 1
c$$$C                ETAB
c$$$      KC = PYCOMP(551)
c$$$      MDCY(KC,1) = 1
c$$$C                ANTI-ETAB
c$$$      KC = PYCOMP(-551)
c$$$      MDCY(KC,1) = 1


C-------------------------

C   SELECT TYPE OF EVENTS TO BE GENERATED

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

C     INITIALIZE PYTHIA 

      
      KC=PYCOMP(2212)
      MP=PMAS(KC,1)

      KC=PYCOMP(5222)
      msigmab=PMAS(KC,1)

      mbottom=4.20D0

      EP=(msigmab-mbottom)/msigmab*EC
      ECM=DSQRT(2.D0*MP*EP)


C..   To stop hadronization:
      MSTP(111)=0
C-----------------------
      CALL PYINIT('CMS','P+','P+',ECM)

C..   To init d to b conversion:
      CALL KFBINITBBD

C     EVENT LOOP

      NEVB=0

 777  CONTINUE
      icase=0
      IB=0
      ig=0
      CALL PYEVNT

******************************************
c     IF (IEV.EQ.5) GOTO 11
******************************************

      EMAX=0.D0
      EMAXM=0.D0
      EMAXB=0.D0
      IH=0

c     recorre las particulas entre el sistema post-colision
c     y antes de la hadronizacion
c     si encuentra un gluon, mete en ig el numero de la particula
c     que es el gluon
c     si la particula procede de 3 y no es un gluon
c     y si es sustituida por KFBINIT (porque tiene algun q=d)
c     aumenta en 1 IB y mete en IKB(1), IKB(2) el numero de particula
c     que cumple la condicion

c..   for diffractive:
      DO IK=5,N
c     si la particula es un gluon de p diff, ig = #de part del gluon
         if(k(ik,2).eq.21.and.K(IK,3).EQ.3) ig=ik
c     print*,"gluon en ",ig
c     si su madre es 3 y no es un gluon
         IF (K(IK,3).EQ.3.AND.K(IK,2).NE.21) THEN
c     si sustituida por otra en el algoritmo
            IF(K(IK,2).NE.KFB(K(IK,2))) THEN
               IB=IB+1
               IKB(IB)=IK
            END IF

         END IF

      END DO

      IF ( IB .EQ. 0 ) GOTO 777

C..   To modify an entry to replace bottom:
c..   random:
      IBR=INT(IB*PYR(0))+1
      IK=IKB(IBR)
      
c     si la masa difractiva less than 1.938 gev y hay dos particulas
c     que pueden ser sustituidas
      if (P(3,5).lt.1.938d0.and.ib.eq.2) then 
         IK=IKB(1)
         if (K(IKB(2),2).gt.1000) IK=IKB(2)
      end if


      K(IK,2)=KFB(K(IK,2))
      KC=PYCOMP(K(IK,2))
      M0=PMAS(KC,1)

c     energia del primario corregida al cambiar la particula
      EPCOR=ECM/2.d0*msigmab/(msigmab-mbottom)
      
*=======================================================================
c     decay a dos cuerpos
      if (P(3,5).gt.mp.and.P(3,5).lt.1.938d0) then

         icase=2

         IK1=IK

         IK2=6
         if(IK.eq.6) IK2=5
         m2=P(IK2,5)

         
c     sustituye estados de J=3/2 por los corr. de J=1/2
         IF (K(IK1,2).EQ.5224) K(IK1,2)=5222
         IF (K(IK1,2).EQ.5214) K(IK1,2)=5212
         IF (K(IK1,2).EQ.5114) K(IK1,2)=5112
         
         KC=PYCOMP(K(IK1,2))
         M0=PMAS(KC,1)
c     print*,'M0  of ',K(IK1,2),' = ',M0
         m1=m0

         d0=EPCOR+P(3,4)-P(1,4)	
         d1=P(3,1)-P(1,1)
         d2=P(3,2)-P(1,2)
         d3=EPCOR*DSQRT(1.D0-msigmab**2/EPCOR**2)+P(3,3)-P(1,3)
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

         P(1,3)=EPCOR*DSQRT(1.D0-msigmab**2/EPCOR**2)
         P(1,4)=EPCOR
         P(1,5)=msigmab

         goto 11

      end if

*=======================================================================
      
      if (K(3,2).eq.2212) then

         icase=1

         IK=5
         K(IK,2)=5222
         P(IK,3)=EPCOR+P(2,3)-P(4,3)
         P(IK,4)=DSQRT(P(IK,1)**2+P(IK,2)**2+P(IK,3)**2+msigmab**2)
         P(IK,5)=msigmab

         CALL PYEXEC


         goto 22

      end if

*=======================================================================

      if(ig.eq.0) then
         IK1 = IKB(IBR)
         IK2 = -IK1+11
      end if

      if(ig.eq.6) then
         IK1 = IKB(IBR)
         IK2 = -IK1+12
      end if

      d0=EPCOR+P(3,4)-P(1,4)	
      d1=      P(3,1)-P(1,1)
      d2=      P(3,2)-P(1,2)
      d3=EPCOR*DSQRT(1.D0-msigmab**2/EPCOR**2)+P(3,3)-P(1,3)
      dm=dsqrt(d0**2-d1**2-d2**2-d3**2)

      q0=P(3,4)-P(1,4)
      q1=P(3,1)-P(1,1)
      q2=P(3,2)-P(1,2)
      q3=P(3,3)-P(1,3)
      
      m1=M0
c     print*,"m1= ",m1
      x1=P(IK1,1)
      y1=P(IK1,2)

      m2=P(IK2,5)
c     print*,"m2= ",m2
      x2=P(IK2,1)
      y2=P(IK2,2)
      
      z3=0d0
      e3=0d0
      if(ig.ne.0) then
         z3=P(ig,3)*mbottom/msigmab
         x3=P(ig,1)
         y3=P(ig,2)
         e3=dsqrt(x3**2+y3**2+z3**2)

         P(ig,1)=x3
         P(ig,2)=y3
         P(ig,3)=z3
         P(ig,4)=e3
      end if

c     math	Solve[Sqrt[x1^2 + y1^2 + z1^2 + m1^2] + 
c     math         Sqrt[x2^2 + y2^2 + (d3 - z1 - z3)^2 + m2^2] + e3 == d0, {z1}]

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

c     print*,'a) ',d0,' = ',e1a+e2a+e3,' e1 = ',e1a,' e2 = ',e2a
c     print*,'b) ',d0,' = ',e1b+e2b+e3,' e1 = ',e1b,' e2 = ',e2b
c     print*, 'z1a = ',z1a,' z1b = ',z1b
c     print*, 'z2a = ',z2a,' z2b = ',z2b
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

c     print*,'check a) ',dsqrt(x1**2+y1**2+z1a**2+m1**2)+ 
c     .         dsqrt(x2**2+y2**2+(d3-z1a-z3)**2+m2**2)+e3-d0
c     print*,'check b) ',dsqrt(x1**2+y1**2+z1b**2+m1**2)+ 
c     .         dsqrt(x2**2+y2**2+(d3-z1b-z3)**2+m2**2)+e3-d0

c     print*, 'z1 = ',z1,' e1 = ',e1
c     print*, 'z2 = ',z2,' e2 = ',e2
      if(ig.ne.0) then
c     print*,"gx= ",x3," gy= ",y3," gz= ",z3," gE= ",e3
      end if
      P(3,3)=d3
      P(3,4)=d0
      P(3,5)=dsqrt(P(3,4)**2-P(3,1)**2-P(3,2)**2-P(3,3)**2)

c     IK=IKB(IBR)
      P(IK1,3)=z1
      P(IK1,4)=e1
      P(IK1,5)=m0
c     print*,K(IK1,2),"E= ",P(IK1,4)
c     IK=IKB(ICO)
      P(IK2,3)=z2
      P(IK2,4)=e2
c     print*,K(IK2,2),"E= ",P(IK2,4)

      P(1,3)=EPCOR*DSQRT(1.D0-msigmab**2/EPCOR**2)
      P(1,4)=EPCOR
      P(1,5)=msigmab

*=======================================================================
 11   continue

C..   To restart event with hadronization:
      CALL PYEXEC
      if (MSTU(23).gt.0.OR.MSTU(27).gt.0) THEN
c     print*,'a) ',d0,'= ',e1a+e2a+e3,'e1 =',e1a,'e2 =',e2a
c     print*,'b) ',d0,'= ',e1b+e2b+e3,'e1 =',e1b,'e2 =',e2b
         MSTU(27)=0
         MSTU(23)=0
      END IF
      CALL PYLIST(1)
      
C------------------------

 22   continue

      EMAX=0.D0
      EAUX = 0.D0
      EAUX2 = 0.D0
      DO IK=3,N
         IF ( K(IK,1) .EQ. 1) THEN
            CALL BOOSTBBD(IK)
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

	SUBROUTINE BOOSTBBD(IK)

      	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
cdh     REAL*8 MP
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

c$$$*=======================================================================
c$$$*2345678901--------2---------3---------4---------5---------6---------7--
c$$$*=======================================================================
      SUBROUTINE KFBINITBBD

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      integer i,j,k,s2

      INTEGER KFB(-5000:5000)
      INTEGER KCODES(0:21)
      INTEGER KCODES_SUB(0:21)
      INTEGER KSAMECODES(0:17)
      COMMON/CONVERTERBBD/KFB,KCODES,KCODES_SUB,KSAMECODES

* quarks:

c      print*,'quarks:'
c      print*,' '
      do i=1,4
         kf=i
         kfb(kf)=kf
         if (i.eq.1) kfb(kf)=5
c         call show(kf,kfb(kf))
      end do

*     diquarks:

c      print*,' '
c      print*,'diquarks:'
c      print*,' '
      do s2=0,2,2
         do j=1,4
	    do i=j,4
               kf=1000*i+100*j+s2+1
               kfb(kf)=kf
               if (i.eq.1) kfb(kf)=1000*5+100*j+s2+1
               if (j.eq.1.and.i.ne.1) kfb(kf)=1000*5+100*i+s2+1
c               call show(kf,kfb(kf))
	    end do
         end do
      end do

*     mesons:

c      print*,' '
c      print*,'mesons:'
c      print*,' '
      do s2=0,2,2
         do j=1,4
	    do i=1,4
               kf=(100*maxBBD(i,j)+10*minBBD(i,j)+s2+1)*
     .              sign(1,i-j)*(-1)**maxBBD(i,j)
               if (i.eq.j) kf=(-1)**i*kf
               kfb(kf)=kf
               if (i.eq.1) 
     .              kfb(kf)=-1*(100*5+10*j+s2+1)
               if (j.eq.1.and.i.ne.1) 
     .              kfb(kf)=(100*5+10*i+s2+1)
               
c               call show(kf,kfb(kf))
	    end do
         end do
      end do

*     baryons:

c      print*,' '
c      print*,'baryons:'
c      print*,' '


* baryons:
      
      DATA KCODES/2212,2112,3122,3212,3112,3312,4122,4212,4112
     .      ,4132,4312,4412,2214,2114,1114,3214,3114,3314
     .      ,4214,4114,4314,4414/

      DATA KCODES_SUB/5222,5122,5232,5322,5312,5332,5242,5422,5412
     .     ,5342,5432,5442,5224,5214,5114,5324,5314,5334
     .     ,5424,5414,5434,5444/
      
      DATA KSAMECODES/2224,3222,3224,3322,3324,3334,4222,4224,4232,4322
     .     ,4324,4332,4334,4422,4424,4432,4434,4444/

      do i=0,21
         kfb(KCODES(i))=KCODES_SUB(i)
         kfb(-KCODES(i))=-KCODES_SUB(i)
c         call show(KCODES(i),kfb(KCODES(i)))
      end do

      do i=0,17
c         print*,i," ",KSAMECODES(i)
         
         kfb(KSAMECODES(i))=KSAMECODES(i)
         kfb(-KSAMECODES(i))=-KSAMECODES(i)
c         call show(KSAMECODES(i),kfb(KSAMECODES(i)))
      end do

      kfb(111)=511
      kfb(113)=513

      end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

	integer function maxBBD(i,j)

	  maxBBD=i
	  if (j.gt.i) maxBBD=j
	
	end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

	integer function minBBD(i,j)

	  minBBD=i
	  if (j.lt.i) minBBD=j
	
	end


