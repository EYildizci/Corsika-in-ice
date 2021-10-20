      SUBROUTINE BOBADIF(EC)
c  bottom baryon diffractive process
      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER PYCOMP
      CHARACTER CHAF*16

      DOUBLE PRECISION m1,m2,m12
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
cdh   COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
      COMMON/PYDAT4/CHAF(500,2)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDATR/MRPY(6),RRPY(100)


C  MDCY MUST BE SET FOR ALL PARTICLES BEFORE CALLING THIS SUBROUTINE
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

cc    recorre las particulas entre el sistema post-colision
cc    y antes de la hadronizacion
c     runs over the particles in the system after collision
c     and before the hadronization
cc    si encuentra un gluon, mete en ig el numero de la particula
cc    que es el gluon
c     if it finds a gluon, it sets in ig the number of the particle
c     which is the gluon.
cc    si la particula procede de 3 y no es un gluon
cc    y si es sustituida por KFBINIT (porque tiene algun q=d)
cc    aumenta en 1 IB y mete en IKB(1), IKB(2) el numero de particula
cc    que cumple la condicion
c     if the particle procedes from 3 and is not a gluon
c     and if it is substituted by KFBINIT (as it has sometimes q=d)
c     increase IB  by 1 and put in IKB(1), IKB(2) the number of the
c     particle which fulfils the condition.

c..   for diffractive:
      DO IK=5,N
cc    si la particula es un gluon de p diff, ig = #de part del gluon
c     if the particle is a gluon with p diff, ig = # of part of gluon
         if(k(ik,2).eq.21.and.K(IK,3).EQ.3) ig=ik
c     print*,"gluon en ",ig
cc    si su madre es 3 y no es un gluon
c     if the mother particle is 3 and is not a gluon
         IF (K(IK,3).EQ.3.AND.K(IK,2).NE.21) THEN
cc    si sustituida por otra en el algoritmo
c     if it is substututed by other in the algorithm
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
      
cc    si la masa difractiva less than 1.938 gev y hay dos particulas
cc    que pueden ser sustituidas
c     if the diffractive mass is less than 1.938 GeV and has two particles
c     which can be replaced
      if (P(3,5).lt.1.938d0.and.ib.eq.2) then 
         IK=IKB(1)
         if (K(IKB(2),2).gt.1000) IK=IKB(2)
      end if


      K(IK,2)=KFB(K(IK,2))
      KC=PYCOMP(K(IK,2))
      M0=PMAS(KC,1)

cc    energia del primario corregida al cambiar la particula
c     energy of primary corrected after change of particle
      EPCOR=ECM/2.d0*msigmab/(msigmab-mbottom)
      
*=======================================================================
c     decay a dos cuerpos  (two body decay)
      if (P(3,5).gt.mp.and.P(3,5).lt.1.938d0) then

         icase=2

         IK1=IK

         IK2=6
         if(IK.eq.6) IK2=5
         m2=P(IK2,5)

         
cc    sustituye estados de J=3/2 por los corr. de J=1/2
c     replace states with J=3/2 by the corrected with J=1/2
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
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================
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


*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================
      SUBROUTINE BOBAPAR(EC,NY)
c  bottom baryon partonic process

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER PYCOMP
      CHARACTER CHAF*16

      DOUBLE PRECISION m1,m2,m12
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
cdh   COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
      COMMON/PYDAT4/CHAF(500,2)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDATR/MRPY(6),RRPY(100)

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
*2345678901--------2---------3---------4---------5---------6---------7--
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
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      integer function maxBBP(i,j)

      maxBBP=i
      if (j.gt.i) maxBBP=j
      
      end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      integer function minBBP(i,j)

      minBBP=i
      if (j.lt.i) minBBP=j
      
      end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      subroutine showBBP(kf1,kf2)
      character cha1*16,cha2*16

      CALL PYNAME(kf1,cha1)
      CALL PYNAME(kf2,cha2)

c      if (kf1.ne.kf2) print*,kf1,cha1,'->',kf2,cha2

      end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
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


*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================
      SUBROUTINE BOMEDIF(EC)
c  bottom meson diffractive process
      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER PYCOMP
      CHARACTER CHAF*16

      DOUBLE PRECISION m1,m2,m12
      DOUBLE PRECISION MP,M0,DM,MC,MSIGMAB,MBOTTOM,MBP,MPI
      INTEGER BINES, LH(25)
      CHARACTER CHANNEL(0:26)*30
      CHARACTER CHAU*16

      INTEGER IKB(2)

      INTEGER KFB(-5000:5000)
      INTEGER KCODES(0:21)
      INTEGER KCODES_SUB(0:21)
      INTEGER KSAMECODES(0:17)
      COMMON/CONVERTERBMF/KFB,KCODES,KCODES_SUB,KSAMECODES

      COMMON/ENER/MP,EMIN,EP,ELAB,PXLAB,PYLAB,PZLAB
      COMMON/EMAXLH/EMAXB,EMAXM,IMAXB,IMAXM,IKMAXB,IKMAXM

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
cdh   COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
      COMMON/PYDAT4/CHAF(500,2)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDATR/MRPY(6),RRPY(100)


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

C     INITIALIZE PYTHIA 


      
      KC=PYCOMP(2212)
      MP=PMAS(KC,1)

      KC=PYCOMP(521)
      MBP=PMAS(KC,1)

      KC=PYCOMP(211)
      MPI=PMAS(KC,1)
      mbottom=4.20D0

      EP=(mBp-mbottom)/mBp*EC
      ECM=DSQRT(2.D0*MP*EP)

C..   To stop hadronization:
      MSTP(111)=0
C-----------------------
      CALL PYINIT('CMS','PI+','P+',ECM)

C..   To init d to b conversion:
      CALL KFBINITBMF

C     EVENT LOOP

 777  CONTINUE
c     print*,'EVENT #',IEV
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

cc    recorre las particulas entre el sistema post-colision
cc    y antes de la hadronizacion
c     runs over the particles in the system after collision
c     and before the hadronization
cc    si encuentra un gluon, mete en ig el numero de la particula
cc    que es el gluon
c     if it finds a gluon, it sets in ig the number of the particle
c     which is the gluon.
cc    si la particula procede de 3 y no es un gluon
cc    y si es sustituida por KFBINIT (porque tiene algun q=d)
cc    aumenta en 1 IB y mete en IKB(1), IKB(2) el numero de particula
cc    que cumple la condicion
c     if the particle procedes from 3 and is not a gluon
c     and if it is substituted by KFBINIT (as it has sometimes q=d)
c     increase IB  by 1 and put in IKB(1), IKB(2) the number of the
c     particle which fulfils the condition.

c..   for diffractive:
      DO IK=5,N
cc    si la particula es un gluon del sistema diff,ig = #de part del gluon
c     if the particle is a gluon with p diff, ig = # of part of gluon
         if(k(ik,2).eq.21.and.K(IK,3).EQ.3) ig=ik
c     print*,"gluon en ",ig
cc    si su madre es 3 y no es un gluon
c    if the mother particle is 3 and is not a gluon
         IF (K(IK,3).EQ.3.AND.K(IK,2).NE.21) THEN
cc    si sustituida por otra en el algoritmo
c     if it is substututed by other in the algorithm
            IF(K(IK,2).NE.KFB(K(IK,2))) THEN

               if(k(ik,2).eq.111.or.k(ik,2).eq.113) then
cdh               if(rnbmf(1).gt.0.5) then 
                  if(pyr(1).gt.0.5d0) then 
c     print*, "a 33!"
                     goto 33
                  end if
               end if
               
               IB=IB+1
               IKB(IB)=IK
c     CALL PYNAME(K(IK,2),CHAU)
c     PRINT*,IK,' ',K(IK,2),' ',CHAU

 33            continue

c     print*,"IB ",IB," IKB(",IB,")=",IKB(IB)                
            END IF


         END IF

      END DO

      IF ( IB .EQ. 0) GOTO 777


      DO i=1,IB
c     print*,"IKB(",IB,")= ",IKB(IB)
      END DO
C..   To modify an entry to replace bottom:
c..   random:
      IBR=INT(IB*PYR(0))+1
      IK=IKB(IBR)
c     print*,'IBR: ',IBR        

cc    si la masa difractiva less than 1.13957 gev y hay dos particulas
cc    que pueden ser sustituidas
c     if the diffractive mass is less than 1.938 GeV and has two particles
c     which can be replaced
      if (P(3,5).lt.1.13957d0.and.ib.eq.2) then 
         IK=IKB(1)
         if (K(IKB(2),2).gt.1000) IK=IKB(2)
      end if

*********************************************
c	if (iev.eq.5) then
c	  IBR=2
c	  IK=IKB(IBR)
c	end if
*********************************************

c     print*,IK,' ','will be replaced'
      K(IK,2)=KFB(K(IK,2))
c..   to include both cu_0 and cu_1:
c$$$  if (K(IK,2).eq.4203) then
c$$$  ii=INT(4*RN(1))+1
c$$$  if (ii.le.3) K(IK,2)=4201
c$$$  end if
      KC=PYCOMP(K(IK,2))
      M0=PMAS(KC,1)

cc    energia del primario corregida al cambiar la particulaA
c     energy of primary corrected after change of particle
      EPCOR=ECM/2.D0*mBp/(mBp-mbottom)
	
*=======================================================================
c     decay a dos cuerpos (two body decay)
      if (P(3,5).gt.mpi.and.P(3,5).lt.1.13957D0) then

	icase=2

	IK1=IK

	IK2=6
	if(IK.eq.6) IK2=5
	m2=P(IK2,5)

c     print*,' ',K(IK1,2),' ',K(IK2,2)
        
*     conversion exceptions:
cc    no existen conversiones que den lugar a violacion de isospin
cc    en bottom no hacen falta el siguiente if
c     no conversions exist which will cause a violation of isospin
c     in bottom are not necessary in the following if
        
c$$$  IF (K(IK1,2).EQ.4122.AND.K(IK2,2).EQ.111) THEN
c$$$  K(IK1,2)=4212
c$$$  END IF
        
cc    sustituye estados de J=3/2 por los corr. de J=1/2
c     replace states with J=3/2 by the corrected with J=1/2
	IF (K(IK1,2).EQ.5224) K(IK1,2)=5222
	IF (K(IK1,2).EQ.5214) K(IK1,2)=5212
	IF (K(IK1,2).EQ.5114) K(IK1,2)=5112
        
	KC=PYCOMP(K(IK1,2))
	M0=PMAS(KC,1)
	m1=m0

	d0=EPCOR+P(3,4)-P(1,4)	
	d1=P(3,1)-P(1,1)
	d2=P(3,2)-P(1,2)
	d3=EPCOR*DSQRT(1.D0-mBp**2/EPCOR**2)+P(3,3)-P(1,3)
	dm=dsqrt(d0**2-d1**2-d2**2-d3**2)

	q0=P(3,4)-P(1,4)
	q1=P(3,1)-P(1,1)
	q2=P(3,2)-P(1,2)
	q3=P(3,3)-P(1,3)
c	write(6,999) q1,q2,q3,q0,dsqrt(-q0**2+q1**2+q2**2+q3**2)

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

	P(1,3)=EPCOR*DSQRT(1.D0-mBp**2/EPCOR**2)
	P(1,4)=EPCOR
	P(1,5)=mBp

	goto 11

      end if

*=======================================================================

      if (K(3,2).eq.211) then
         
         icase=1

         IK=5
         K(IK,2)=521
         P(IK,3)=EPCOR+P(2,3)-P(4,3)
         P(IK,4)=DSQRT(P(IK,1)**2+P(IK,2)**2+P(IK,3)**2+mBp**2)
         P(IK,5)=mBp

         CALL PYEXEC

         goto 22

      end if

*=======================================================================

      
      if(ig.eq.0) then
         IK1 = IKB(IBR)
c     print*, "0 gluones ",IKB(IBR)
         IK2 = -IK1+11
      end if

      if(ig.eq.6) then
         IK1 = IKB(IBR)
c     print*, "un gluon en 6 ",IKB(IBR)
         IK2 = -IK1+12
      end if

c     print*,'(IK1,IK2):(',IK1,',',IK2,')'
      
      d0=EPCOR+P(3,4)-P(1,4)	
      d1=      P(3,1)-P(1,1)
      d2=      P(3,2)-P(1,2)
      d3=EPCOR*DSQRT(1.D0-mBp**2/EPCOR**2)+P(3,3)-P(1,3)
      dm=dsqrt(d0**2-d1**2-d2**2-d3**2)

      q0=P(3,4)-P(1,4)
      q1=P(3,1)-P(1,1)
      q2=P(3,2)-P(1,2)
      q3=P(3,3)-P(1,3)
c      write(6,999) q1,q2,q3,q0,dsqrt(-q0**2+q1**2+q2**2+q3**2)
      
c     IK=IKB(IBR)
      m1=M0
c     print*,"m1= ",m1
      x1=P(IK1,1)
      y1=P(IK1,2)

c     IK=IKB(ICO)
      m2=P(IK2,5)
c     print*,"m2= ",m2
      x2=P(IK2,1)
      y2=P(IK2,2)
      
      z3=0.d0
      e3=0.d0
      if(ig.ne.0) then
         z3=P(ig,3)*mbottom/mBp
         x3=P(ig,1)
         y3=P(ig,2)
         e3=dsqrt(x3**2+y3**2+z3**2)

c     e3=P(ig,4)*mbottom/mBp
c     z3=P(ig,3)
c     e3=P(ig,4)
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


c$$$	print*,'a) ',d0,' = ',e1a+e2a+e3,' e1 = ',e1a,' e2 = ',e2a
c$$$	print*,'b) ',d0,' = ',e1b+e2b+e3,' e1 = ',e1b,' e2 = ',e2b
c$$$        print*, 'z1a = ',z1a,' z1b = ',z1b
c$$$        print*, 'z2a = ',z2a,' z2b = ',z2b

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

c$$$	print*,'check a) ',dsqrt(x1**2+y1**2+z1a**2+m1**2)+ 
c$$$     .         dsqrt(x2**2+y2**2+(d3-z1a-z3)**2+m2**2)+e3-d0
c$$$	print*,'check b) ',dsqrt(x1**2+y1**2+z1b**2+m1**2)+ 
c$$$     .         dsqrt(x2**2+y2**2+(d3-z1b-z3)**2+m2**2)+e3-d0

c$$$        print*, 'z1 = ',z1,' e1 = ',e1
c$$$        print*, 'z2 = ',z2,' e2 = ',e2

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

      P(1,3)=EPCOR*DSQRT(1.D0-mBp**2/EPCOR**2)
      P(1,4)=EPCOR
      P(1,5)=mBp

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
c      CALL PYLIST(1)

C------------------------ 

 22     continue

        EAUX = 0D0
        EAUX2 = 0D0
        DO IK=3,N
           IF ( K(IK,1) .EQ. 1) THEN
              CALL BOOSTBMD(IK)
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
     .                P(IK,2)**2-PMAS(KC,1)**2)
                 EAUX2 = EAUX2 + P(IK,4)
              ENDIF
           ENDDO

        ENDIF
        
	END

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================
	SUBROUTINE BOOSTBMD(IK)

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
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================
      SUBROUTINE KFBINITBMF

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      integer i,j,k,s2

      INTEGER KFB(-5000:5000)
      INTEGER KCODES(0:21)
      INTEGER KCODES_SUB(0:21)
      INTEGER KSAMECODES(0:17)
      COMMON/CONVERTERBMF/KFB,KCODES,KCODES_SUB,KSAMECODES

*     quarks:

c     print*,'quarks:'
c     print*,' '
      do i=-4,4
         kf=i
         kfb(kf)=kf
         if (i.eq.1) kfb(kf)=5
         if (i.eq.-1) kfb(kf)=-5
c     call show(kf,kfb(kf))
      end do

*     diquarks:

c     print*,' '
c     print*,'diquarks:'
c     print*,' '
      do s2=0,2,2
         do j=1,4
	    do i=j,4
               kf=1000*i+100*j+s2+1
               kfb(kf)=kf
               if (i.eq.1) kfb(kf)=1000*5+100*j+s2+1
               if (j.eq.1.and.i.ne.1) kfb(kf)=1000*5+100*i+s2+1
c     call show(kf,kfb(kf))
	    end do
         end do
      end do

*     mesons:

c     print*,' '
c     print*,'mesons:'
c     print*,' '
      do s2=0,2,2
         do j=1,4
	    do i=1,4
               kf=(100*maxBMF(i,j)+10*minBMF(i,j)+s2+1)*
     .              sign(1,i-j)*(-1)**maxBMF(i,j)
               if (i.eq.j) kf=(-1)**i*kf
               kfb(kf)=kf
               if (i.eq.1) 
     .              kfb(kf)=-1*(100*5+10*j+s2+1)
               if (j.eq.1.and.i.ne.1) 
     .              kfb(kf)=(100*5+10*i+s2+1)
               
c     call show(kf,kfb(kf))
	    end do
         end do
      end do

*     baryons:

c     print*,' '
c     print*,'baryons:'
c     print*,' '

*     baryons:
      
      DATA KCODES/2212,2112,3122,3212,3112,3312,4122,4212,4112
     .     ,4132,4312,4412,2214,2114,1114,3214,3114,3314
     .     ,4214,4114,4314,4414/

      DATA KCODES_SUB/5222,5122,5232,5322,5312,5332,5242,5422,5412
     .     ,5342,5432,5442,5224,5214,5114,5324,5314,5334
     .     ,5424,5414,5434,5444/
      
      DATA KSAMECODES/2224,3222,3224,3322,3324,3334,4222,4224,4232,4322
     .     ,4324,4332,4334,4422,4424,4432,4434,4444/

      do i=0,21
         kfb(KCODES(i))=KCODES_SUB(i)
         kfb(-KCODES(i))=-KCODES_SUB(i)
c     call show(KCODES(i),kfb(KCODES(i)))
      end do

      do i=0,17
c     print*,i," ",KSAMECODES(i)
         
         kfb(KSAMECODES(i))=KSAMECODES(i)
         kfb(-KSAMECODES(i))=-KSAMECODES(i)
c     call show(KSAMECODES(i),kfb(KSAMECODES(i)))
      end do

      kfb(111)=511
      kfb(113)=513

      end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

	integer function maxBMF(i,j)

	  maxBMF=i
	  if (j.gt.i) maxBMF=j
	
	end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

	integer function minBMF(i,j)

	  minBMF=i
	  if (j.lt.i) minBMF=j
	
	end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================
      SUBROUTINE BOMEPAR(EC,NY)
c  bottom meson partonic process

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER PYCOMP
      INTEGER COL
      CHARACTER CHAF*16

      DOUBLE PRECISION m1,m2,m12
      DOUBLE PRECISION MP,M0,DM,MC,MB,MBOTTOM,MPI
      INTEGER BINES, LH(25)
      CHARACTER CHANNEL(0:25)*30
      CHARACTER CHAU*16
      CHARACTER LTAR*5

      INTEGER IKB(200)
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

      COMMON/CONVERTERBMP/KFB,KCODES,KCODES_SUB,KSAMECODES

      COMMON/ENER/MP,EMIN,EP,ELAB,PXLAB,PYLAB,PZLAB
      COMMON/EMAXLH/EMAXB,EMAXM,IMAXB,IMAXM,IKMAXB,IKMAXM

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
cdh   COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
      COMMON/PYDAT4/CHAF(500,2)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)

      COMMON/PYDATR/MRPY(6),RRPY(100)


C   SELECT TYPE OF EVENTS TO BE GENERATED

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
*     95: Low-pT scattering	(keep always active)
      MSUB(95)=1

	
      A = 1.D0
      ZB = 0.D0
      EP = EC
      KC=PYCOMP(2212)
      MP=PMAS(KC,1)

      KC=PYCOMP(521)
      MB=PMAS(KC,1)

      KC=PYCOMP(211)
      MPI=PMAS(KC,1)

      MBOTTOM=4.20D0

C     ASSIGN SEA PAIRS ENERGIES

      NP = 7
      NN = 7
      FAIL = 0.D0
      IF (NY.GT.1) THEN
 123     CONTINUE
         FAIL = 0.D0
         ZB = 0.D0
         DO 124 INW=1,NY-1
            XGLIST(INW) = GXBMP(1.0)
            ZB = ZB+XGLIST(INW)
C     ASSIGN IDS
            IF (PYR(0).LE.2.D0/2.3)THEN
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
         IF (DSQRT(2.D0*MP*(MP+((MB-MBOTTOM)/MB*(EP*(1.D0-ZB))))).
     *        LE. PARP(2)) GOTO 123

         

      ENDIF


      EPR = EP*(1.D0-ZB)
     
C     REMOVE QQBAR PAIRS ENERGY: ENERGY FOR QCD COLLISION

C..   To stop hadronization:
      MSTP(111)=0
      MSTP(125)=1
C-----------------------
      EAUX = 0.D0
      M = 1
      IF (NY .GT. 1) THEN
         DO INW=1,NY-1

            IF (TLIST(INW).EQ.'F') GOTO 500

            CALL PYINIT('CMS','pi0',TLIST(INW),ECMLIST(INW))            

            MINT(11) = KFLIST(INW)
            VINT(3) = QMASS(INW)

            CALL PYEVNT
            CALL PYEXEC

            DO IPART=9,N
               IF (K(IPART,1).EQ.1) THEN
                  CALL BOOSTBMP(IPART,EP*XGLIST(INW),MP)

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

C-----MAIN INTERACTION------------------

      EP=(MB-MBOTTOM)/MB*EPR
      ECM=DSQRT(2.D0*MP*(MP+EP))

      EPCOR=ECM/2.D0*MB/(MB-MBOTTOM)
      DELTAE=EPCOR-ECM/2.D0
      DELTAPZ=EPCOR*DSQRT(1.D0-MB**2/EPCOR**2)
     .     -ECM/2D0*DSQRT(1.D0-(2.D0*MP)**2/ECM**2)

C..   To stop hadronization:
      MSTP(111)=0
C-----------------------
      CALL PYINIT('CMS','Pi+',LTAR,ECM)

C..   To init u to c conversion:
      CALL KFBINITBMP

C     EVENT LOOP

 777  CONTINUE
      IB=0
      ig=0
      CALL PYEVNT


      EMAX=0.D0
      EMAXM=0.D0
      EMAXB=0.D0
      IH=0
      
c..   for hard:
      IK1=9
c!!!! PMAX=0D0
      PMAX=-10000.D0
      DO IK=9,N
c     si no es un foton (if it is not a photon)
         IF (K(IK,2).NE.22) THEN
c     si no es un gluon (if it is not a gluon)
            IF (K(IK,2).NE.21) THEN  
cc    si es cambiada en el algoritmo y pz>0
c     if it is changed in the algorithm and pz > 0
               IF(K(IK,2).NE.KFB(K(IK,2)).AND.P(IK,3).GT.0.D0) THEN
c     si es particula (if it is a particle)
                  IF(K(IK,2).GT.0) THEN
                     
                     IB=IB+1
                     IKB(IB)=IK
c     CALL PYNAME(K(IK,2),CHAU)
c     PRINT*,IK,' ',K(IK,2),' ',CHAU

                     IF (P(IK,3).GT.PMAX) THEN
                        PMAX=P(IK,3)
                        IKMAX=IK
                        IBMAX=IB
                     END IF
c     PRINT*,IK,P(IK,3)
c     END IF
c     si es antiparticule (if it is an antiparticle)
                  ELSE IF(K(IK,2).LT.0) THEN
c     si su madre es 1,3,5 o 7 (if its mother is 1,3,5 or 7)
                     IF(K(IK,3).eq.1.OR.K(IK,3).eq.3.OR.K(IK,3).eq.5
     .                    .OR.K(IK,3).eq.7) THEN                      
                        IB=IB+1
                        IKB(IB)=IK
c     CALL PYNAME(K(IK,2),CHAU)
c     PRINT*,IK,' ',K(IK,2),' ',CHAU
                        IF (P(IK,3).GT.PMAX) THEN
                           PMAX=P(IK,3)
                           IKMAX=IK
                           IBMAX=IB
                        END IF

                     END IF
                  END IF
               END IF
            END IF
         END IF
         
      END DO
      
      IF (IB .EQ. 0) GOTO 777
      
C..   To modify an entry to replace bottom:
c..   random:
      IBR=INT(IB*PYR(0))+1
      IK=IKB(IBR)
c$$$  IF (IB.GE.3) THEN
c$$$  IBR=IBMAX
c$$$  IK=IKMAX
c$$$  END IF

cc    busca los extremos de la string
c     look for the ends of the string
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

*********************************************
c     if (iev.eq.20) then
c     IBR=1
c     IK=IKB(IBR)
c     end if
*********************************************

c     print*,IK,'will be replaced', IK1,IK2

c     to include both bd_0 and bd_1:
c     bd_0 [5101]:(1/4) and bd_1[5103]:(3/4)     
      if (K(IK,2).eq.5103) then
         ii=INT(4.D0*PYR(0))+1
c     si ii =1
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
c     print	PRINT*,'STRING MASS =',DM

****************
*     (a2)
c     DM=DM*DSQRT(1.5d0/0.33d0)
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
         P(I,1)=P(I,1)*MBOTTOM/MB
         P(I,2)=P(I,2)*MBOTTOM/MB
         P(I,3)=P(I,3)*MBOTTOM/MB
         P(I,4)=P(I,4)*MBOTTOM/MB
         z3=z3+P(I,3)
         e3=e3+P(I,4)
      END DO

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

c     print	print*,'a) ',d0,'= ',e1a+e2a+e3,'e1 =',e1a,'e2 =',e2a
c     print	print*,'b) ',d0,'= ',e1b+e2b+e3,'e1 =',e1b,'e2 =',e2b

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

 77   P(1,3)=EPCOR*DSQRT(1D0-MB**2/EPCOR**2)
      P(1,4)=EPCOR
      P(1,5)=MB

*=======================================================================

C..   To restart event with hadronization:
      CALL PYEXEC
C      CALL PYLIST(1)
      DO IPART=9,N
         IF (K(IPART,1).EQ.1) THEN
            CALL BOOSTBMP(IPART,EP,MP)
            
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

	SUBROUTINE BOOSTBMP(IK,DE0,DM0)

      	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
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
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================
      SUBROUTINE KFBINITBMP

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      integer i,j,k,s2

      INTEGER KFB(-5000:5000)
      INTEGER KCODES(0:21)
      INTEGER KCODES_SUB(0:21)
      INTEGER KSAMECODES(0:17)
      COMMON/CONVERTERBMP/KFB,KCODES,KCODES_SUB,KSAMECODES

* quarks:

c      print*,'quarks:'
c      print*,' '
      do i=-5,5
         kf=i
         kfb(kf)=kf
         if (i.eq.1) kfb(kf)=5
         if (i.eq.-1) kfb(kf)=-5
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
               kf=(100*maxBMP(i,j)+10*minBMP(i,j)+s2+1)*
     .              sign(1,i-j)*(-1)**maxBMP(i,j)
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

	integer function maxBMP(i,j)

	  maxBMP=i
	  if (j.gt.i) maxBMP=j
	
	end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

	integer function minBMP(i,j)

	  minBMP=i
	  if (j.lt.i) minBMP=j
	
	end


*=======================================================================

      DOUBLE PRECISION FUNCTION GXBMP(A)
c  function to generate the energy from projectile parton of a 
c  sea qq_bar pair (including normalization)
      DOUBLE PRECISION Y,YR,PYR
      EXTERNAL PYR

 10   CONTINUE
      GXBMP = PYR(0)
      YR = PYR(0)*1000.D0
      Y = 4.77497509534447750D-02*(1.D0-GXBMP)**4/(GXBMP+1.0D-10)
      IF ( YR.GT.Y ) GOTO 10

      RETURN  
      END


*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================
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
cdh   COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
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
*2345678901--------2---------3---------4---------5---------6---------7--
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
*2345678901--------2---------3---------4---------5---------6---------7--
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


*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================
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
cdh   COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
      COMMON/PYDAT4/CHAF(500,2)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDATR/MRPY(6),RRPY(100)

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
*2345678901--------2---------3---------4---------5---------6---------7--
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

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================
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
cdh   COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
      COMMON/PYDAT4/CHAF(500,2)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDATR/MRPY(6),RRPY(100)

      EMIN=0.1D0

cdh   PRINT*, 'N ',N
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
*2345678901--------2---------3---------4---------5---------6---------7--
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
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      integer function maxCMD(i,j)

      maxCMD=i
      if (j.gt.i) maxCMD=j
      
      end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      integer function minCMD(i,j)

      minCMD=i
      if (j.lt.i) minCMD=j
      
      end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      subroutine showCMD(kf1,kf2)
      character cha1*16,cha2*16

      CALL PYNAME(kf1,cha1)
      CALL PYNAME(kf2,cha2)

      if (kf1.ne.kf2) print*,kf1,cha1,'->',kf2,cha2

      end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================
      SUBROUTINE CHAMEPAR(EC,NY)
c  charmed meson partonic process
      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER PYCOMP
      CHARACTER CHAF*16
      
      DOUBLE PRECISION m1,m2,m12
      DOUBLE PRECISION MPI,M0,DM,MC,MD,MCHARM,MP
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
cdh   COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
      COMMON/PYDAT4/CHAF(500,2)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDATR/MRPY(6),RRPY(100)


C     SELECT TYPE OF EVENTS TO BE GENERATED

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

      KC=PYCOMP(211)
      MPI=PMAS(KC,1)

      KC=PYCOMP(2212)
      MP=PMAS(KC,1)

      KC=PYCOMP(411)
      MD=PMAS(KC,1)
      MCHARM=1.27D0

      EP = EC

C     ASSIGN SEA PAIRS ENERGIES

      NP = 7
      NN = 7
      FAIL = 0.D0
      IF (NY.GT.1) THEN
 123     CONTINUE
         FAIL = 0.D0
         ZB = 0.
         DO 124 INW=1,NY-1
            XGLIST(INW) = GXCMP(1.0)
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
         IF (DSQRT(2.D0*MP*(MP+((MD-MCHARM)/MD*(EP*(1.D0-ZB))))).
     *        LE. PARP(2)) GOTO 123

      ENDIF


      EPR = EP*(1.D0-ZB)
     
C     REMOVE QQBAR PAIRS ENERGY: ENERGY FOR QCD COLLISION

C..   To stop hadronization:
      MSTP(111)=0
      MSTP(125)=1
C-----------------------
      EAUX = 0.D0
      M = 1
      IF (NY .GT. 1) THEN
         DO INW=1,NY-1

            IF (TLIST(INW).EQ.'F') GOTO 500

            CALL PYINIT('CMS','pi0',TLIST(INW),ECMLIST(INW))            

            MINT(11) = KFLIST(INW)
            VINT(3) = QMASS(INW)

            CALL PYEVNT
            CALL PYEXEC
c            CALL PYLIST(1)
            DO IPART=9,N
               IF (K(IPART,1).EQ.1) THEN
                  CALL BOOSTCMP(IPART,EP*XGLIST(INW),MP)

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

C-----------------------

      EP=(MD-MCHARM)/MD*EPR
      ECM=DSQRT(2.D0*MP*(MP+EP))

      EPCOR=ECM/2.D0*MD/(MD-MCHARM)
      DELTAE=EPCOR-ECM/2.D0
      DELTAPZ=EPCOR*DSQRT(1.D0-MD**2/EPCOR**2)
     .     -ECM/2.D0*DSQRT(1.D0-(2D0*MP)**2/ECM**2)

C..   To stop hadronization:
      MSTP(111)=0
      MSTP(125)=1
C-----------------------
      CALL PYINIT('CMS','Pi+',LTAR,ECM)
C..   To init u to c conversion:
      CALL KFCINITCMP
 777  CONTINUE

      IC=0
      ig=0
      
      CALL PYEVNT
c      CALL PYLIST(1)


      EMAXM=0.D0
      EMAXB=0.D0
      IH=0

c..   for hard:
      IK1=9
c!!!! PMAX=0D0
      PMAX=-10000.D0
c      PRINT*,'N before loop ',N
      DO IK=9,N
c     si no es un foton (if it is not a photon)
         IF (K(IK,2).NE.22) THEN
cc    si no es un gluon y no es antiparticula 
c     if it is not a gluon and not an antiparticle
	    IF (K(IK,2).NE.21.and.K(IK,2).GT.0) THEN
c!!!! IF(K(IK,2).NE.KFC(K(IK,2)).AND.P(IK,3).GT.0.D0) THEN
               IF(K(IK,2).NE.KFC(K(IK,2))) THEN
c!!!! IF(K(IK,2).NE.4) THEN
                  IF(K(IK,2).NE.4.AND.K(IK,2).NE.5) THEN
                     IC=IC+1
                     IKC(IC)=IK
c                     PRINT*,'IK inside loop ',IK
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

      IF ( IC .EQ. 0) GOTO 777

C..   To modify an entry to replace charm:
c..   random:
c      PRINT*,'IK  1',IK
      ICR=INT(IC*PYR(0))+1
      IK=IKC(ICR)
c      PRINT*,'IK 2 ',IK
      IF (IC.GE.3) THEN
         ICR=ICMAX
         IK=IKMAX
      END IF
      
c      PRINT*,'ENTRY ',IK,' WILL BE REPLACED'
      
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

*********************************************
c     if (iev.eq.20) then
c     ICR=1
c     IK=IKC(ICR)
c     end if
*********************************************

c..   to include both cu_0 and cu_1:
      if (K(IK,2).eq.4203) then
         ii=INT(4.D0*PYR(0))+1
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

****************
*     (a2)
c     DM=DM*DSQRT(1.5d0/0.33d0)
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
         P(I,1)=P(I,1)*MCHARM/MD
         P(I,2)=P(I,2)*MCHARM/MD
         P(I,3)=P(I,3)*MCHARM/MD
         P(I,4)=P(I,4)*MCHARM/MD
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

 77   P(1,3)=EPCOR*DSQRT(1.D0-MD**2/EPCOR**2)
      P(1,4)=EPCOR
      P(1,5)=MD

*=======================================================================

      CALL PYEXEC
c      CALL PYLIST(1)

C------------------------
C      EAUX=0
      DO IPART=9,N
         IF (K(IPART,1).EQ.1) THEN
            CALL BOOSTCMP(IPART,EP,MP)
            
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

	SUBROUTINE BOOSTCMP(IK,DE0,DM0)

      	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 	REAL*8 MP

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

	RETURN
	END

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================
      SUBROUTINE KFCINITCMP

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
c     call showCMP(kf,kfc(kf))
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
c     call showCMP(kf,kfc(kf))
	    end do
         end do
      end do

*     mesons:

c     print*,'mesons:'
      do s2=0,2,2
         do j=1,3
	    do i=1,3
               kf=(100*maxCMP(i,j)+10*minCMP(i,j)+s2+1)*
     .              sign(1,i-j)*(-1)**maxCMP(i,j)
               if (i.eq.j) kf=(-1)**i*kf
               kfc(kf)=kf
               if (i.eq.2) 
     .              kfc(kf)=(100*maxCMP(4,j)+10*minCMP(4,j)+s2+1)*
     .              sign(1,4-j)*(-1)**maxCMP(4,j)
               if (j.eq.2.and.i.ne.2) 
     .              kfc(kf)=(100*maxCMP(i,4)+10*minCMP(i,4)+s2+1)*
     .              sign(1,i-4)*(-1)**maxCMP(i,4)
c     call showCMP(kf,kfc(kf))
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
c     call showCMP(kf,kfc(kf))
                     kf=1000*3+100*1+10*2+s2+1
                     kfc(kf)=1000*4+100*1+10*3+s2+1
                     goto 1
                  end if
                  if (i.eq.3.and.j.eq.2.and.k.eq.1.and.s2.eq.2) then
                     kf=1000*3+100*2+10*1+s2+1
                     kfc(kf)=1000*4+100*3+10*1+s2+1
c     call showCMP(kf,kfc(kf))
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
c	  CALL showCMP(kf,kfc(kf))
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

	integer function maxCMP(i,j)

	  maxCMP=i
	  if (j.gt.i) maxCMP=j
	
	end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

	integer function minCMP(i,j)

	  minCMP=i
	  if (j.lt.i) minCMP=j
	
	end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

	subroutine showCMP(kf1,kf2)
	character cha1*16,cha2*16

	CALL PYNAME(kf1,cha1)
	CALL PYNAME(kf2,cha2)

	if (kf1.ne.kf2) print*,kf1,cha1,'->',kf2,cha2

	end

*=======================================================================
*2345678901--------2---------3---------4---------5---------6---------7--
*=======================================================================

      DOUBLE PRECISION FUNCTION GXCMP(A)
c  function to generate the energy from projectile parton of a 
c  sea qq_bar pair (including normalization)
      DOUBLE PRECISION Y,YR,PYR
      EXTERNAL PYR

 10   CONTINUE
      GXCMP = PYR(0)
      YR = PYR(0)*1000.D0
      Y = 4.77497509534447750D-02*(1.D0-GXCMP)**4/(GXCMP+1.0D-10)
      IF ( YR.GT.Y ) GOTO 10

      RETURN  
      END

