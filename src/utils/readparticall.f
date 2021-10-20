c=======================================================================
c
c  r e a d p a r t i c a l l . f
c  ----------------------------- 
c      read all particle data files of a parallel corsika simulation 
c      and print tabular of detected particles. 
c-----------------------------------------------------------------------
c  compile+link:
c      gfortran -fbounds-check readparticall.f -o readparticall 
c      f77 -fbounds-check readparticall.f -o readparticall
c      ifort -C -check bounds readparticall.f -o readparticall
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c  running:
c     ls -1 DAT000182-* | grep t -v | grep n -v > readparticall.i000182
c     ./readparticall < readparticall.i000182 > readparticall.out000182   
c-----------------------------------------------------------------------
c     input-files:
c          unit=*: names of particle data files.
c          unit=3: current corsika particle data file.
c     output-file: 
c          unit=*: protocol output.
c-----------------------------------------------------------------------
c                                      juergen.oehlschlaeger@kit.edu
c=======================================================================
  
      program readparticall
 
      implicit double precision (a-h,o-z), integer (i-n) 

      parameter (lenthin=6552,lenstnd=5733,nmax=1000) 

      character cdata(1000)*120,cdat*120,cblk*120
      character cout*120,crunh*4,cevte*4,qpatext(1:200)*19
      double precision aatm(5),batm(5),catm(5)
      double precision parmas(0:101),phead(30)
      dimension lpdat(0:21),qpartic(200),jpartic(0:201,32)
      real pdata(lenthin),qdata(936),prunh,pevte
      equivalence (crunh,prunh),(cevte,pevte)
      common /atmos/aatm,batm,catm
      common /integ/lobs,lsho,ishu,irec,isho,isub,lpdat,jpartic
      common /utabl/pdata,parmas,phead,cpi180,c180pi,qpartic
      data crunh/'RUNH'/,cevte/'EVTE'/

c--initialize names of particles----------------------------------------
      data (qpatext(i),i=1,60)/
     c' gamma             ',' positron          ',' electron          ',
     c' _stck_in_         ',' muon+             ',' muon-             ',
     c' pi0               ',' pi+               ',' pi-               ',
     c' K0long            ',' K+                ',' K-                ',
     c' neutron           ',' proton            ',' anti proton       ',
     c' K0short           ','                   ',' Lambda            ',
     c' Sigma+            ',' Sigma0            ',' Sigma-            ',
     c' Xi0               ',' Xi-               ',' Omega-            ',
     c' anti neutron      ',' anti Lambda       ',' anti Sigma-       ',
     c' anti Sigma0       ',' anti Sigma+       ',' anti Xi0          ',
     c' anti Xi+          ',' anti Omega+       ',16*'                ',
     c'                   ',' omega             ',' rho0              ',
     c' rho+              ',' rho-              ',' Delta++           ',
     c' Delta+            ',' Delta0            ',' Delta-            ',
     c' anti Delta--      ',' anti Delta-       ',' anti Delta0       '/
      data (qpatext(i),i=61,129)/
     c' anti Delta+       ',' K*0               ',' K*+               ',
     c' K*-               ',' anti K*0          ',' electron neutrino ',
     c' anti elec neutrino',' muon neutrino     ',' anti muon neutrino',
     c'                   ',' eta=>2*gamma      ',' eta=>3*pi0        ',
     c' eta=>pi+pi-pi0    ',' eta=>pi+pi-gamma  ',' addi muon+        ',
     c' addi muon-        ','                   ',37*'                ',
     c'                   ',' D0                ',' D+                ',
     c' anti D-           ',' anti D0           ',' Ds+               ',
     c' anti Ds-          ',' eta c             ',' D*0               ',
     c' D*+               ',' anti D*-          ',' anti D*0          ',
     c' D*s+              ',' anti D*s-         ','                   '/
      data (qpatext(i),i=130,200)/
     c' J/psi             ',' tau +             ',' tau -             ',
     c' tau neutrino      ',' anti tau neutrino ','                   ',
     c'                   ',' Lambda c +        ',' Xi c +            ',
     c' Xi c 0            ',' Sigma c ++        ',' Sigma c +         ',
     c' Sigma c 0         ',' Xi c prime +      ',' Xi c prime 0      ',
     c' Omega c 0         ','                   ','                   ',
     c'                   ',' anti Lambda c -   ',' anti Xi c -       ',
     c' anti Xi c 0       ',' anti Sigma c --   ',' anti Sigma c -    ',
     c' anti Sigma c 0    ',' anti Xi c prime - ',' anti Xi c prime 0 ',
     c' anti Omega c 0    ','                   ','                   ',
     c'                   ',' Sigma c * ++      ',' Sigma c * +       ',
     c' Sigma c * 0       ',7*'                 ',' anti Sigma c * -- ',
     c' anti Sigma c * -  ',' anti Sigma c * 0  ',25*'                ',
     c' Cherenkov photon  ','                   '/
      qpatext(176) = ' Stnd              '
      qpatext(196) = ' Thin              '

c--initialize some more quantities--------------------------------------
      cblk='                                                  '
      cdat=cblk
      cpi = 4.d0 * atan(1.d0)
      cpi180 = cpi/180.d0
      c180pi = 180.d0/cpi
      isho = 0
      iend = 0
      irec = 0
 
c--read all names of particle data files--------------------------------
      do  ifil=1,nmax
         read(*,'(a)',end=428,err=427) cdat
         cdata(ifil) = cdat(1:16)
         if ( ifil .eq. 1 ) read(cdata(ifil)(4:9),'(i6)') irun
      enddo
      write(*,'(6x,''maximum number of particle data files reached.'')')
      goto 429
  427 continue
      write(*,'(6x,''ERROR in reading names of particle data files.'')')
      goto 429
  428 continue
      ! write(*,'(6x,''END condition reading particle data files.'')')
  429 continue
      nfil = ifil - 1
      write(*,'(2x)')
 
c----------many big particle data files for one shower------------------
      do  444  ifil=1,nfil
      if (ifil.gt.1) close(unit=3)
      irec = 0
      idat = index(cdata(ifil),'DAT')
      ilen = index(cdata(ifil),' ') - 1
      write(*,'(6x,a)') cdata(ifil)(idat:ilen)
* - - - - - - read data record to get kind of simulation:
      open(unit=3,file=cdata(ifil)(idat:ilen),status='old',
     +     form='unformatted',access='sequential')
      itype = 2
      read(unit=3,err=496,end=432) (pdata(i),i=1,lenstnd)
      close(unit=3)
* - - - - - - detect `standard` simulation instead of `thinning`:
      if ( 217433.0 .lt. pdata(273+1) .and.
     +                   pdata(273+1) .lt. 217433.2 ) then
         itype = 0
         cout = cdata(ifil)(1:ilen)//'.ascistnd'
         lenrec = lenstnd
      elseif ( 217433.0 .lt. pdata(312+1) .and.
     +                       pdata(312+1) .lt. 217433.2 ) then
         itype = 1
         cout = cdata(ifil)(1:ilen)//'.ascithin'
         lenrec = lenthin
      else
         write(*,*) '    ______________________________________________'
         write(*,*) '    ERROR: this corsika syntax should not occur.'
         goto 498
      endif
      lenblk = lenrec / 21
      lenpar = lenblk / 39 
      do  i=0,21
         lpdat(i) = i * lenblk
      enddo
      do  i=1,936 ! keep first two subblocks and some particles.
         qdata(i) = pdata(i)
      enddo
* - - - - - - read data record with lenrec words - - - -
      open(unit=3,file=cdata(ifil)(idat:ilen),status='old',
     +     form='unformatted',access='sequential')
  431 continue
      read(unit=3,err=496,end=432) (pdata(i),i=1,lenrec)
      irec = irec + 1
      ! if (mod(irec,200) .eq. 0) write(*,*) '       irec =',irec
      call blwork(iend,lenrec)
      if ( iend .le. 0 ) goto 431
  432 continue
* - - - - - - end of corsika particle data file reached: 
      write(*,'(11x,''isub ='',i3,'' (rune)      irec ='',i8,
     +        ''      irwc ='',i12)') isub,irec,819*irec
  443 continue
      if ( qdata(274)+qdata(547) .lt. 1. ) then
      !  write(*,'(8x,14(''_''),/,8x,''`thinning run`'',/)')
      else
      !  write(*,'(8x,14(''_''),/,8x,''`standard run`'',/)')
      endif
  444 continue
  445 continue
      close(unit=3)
      goto 499
 
c--end of data----------------------------------------------------------
  496 continue
      inull = 0
      do  l=1,5733,7
         if ( pdata(l) .eq. 0.d0 ) inull = inull + 1
      enddo
      do  i=1,5733,273
         write(*,'(f12.0,1p,3e14.6,''   i='',i4)') (pdata(l),l=i,i+3),i
      enddo
      write(*,*)
     +'    ____________________________________________________________'
      write(*,*) '    ERROR: simulation type of corsika is `standard`.'
      write(*,'(9x,''estimated number of particles'',f15.0,
     +   '' (as 1 sh.)'')') 741.d0+(819.d0*(irec-1))-inull  
      write(*,*) '        subblocks EVTE and RUNE missing!'
      goto 499 
  497 continue
      inull = 0
      do  l=1,6552,8
         if ( pdata(l) .eq. 0.d0 ) inull = inull + 1
      enddo
      do  i=1,6552,312
         write(*,'(f12.0,1p,3e14.6,''   i='',i4)') (pdata(l),l=i,i+3),i
      enddo
      write(*,*) '    _________________________________________________'
      write(*,*) '    ERROR: simulation type of corsika is `thinning`.'
      write(*,'(9x,''estimated number of particles'',f15.0,
     +   '' (as 1 sh.)'')') 741.d0+(819.d0*(irec-1))-inull  
      write(*,*) '        subblocks EVTE and RUNE missing!'
      goto 499 
  498 continue
      write(*,*) '    _________________________________________________'
      write(*,*) '    ERROR: missing some input parameters.'
  499 continue
c--print out tabulars---------------------------------------------------
      npot = 16 ! logarithmic bins per decade.
      fpot = 10.d0**(1.d0/npot) ! factor per bin.
      if ( fpot .lt. 0. ) then
         write(*,'(10x,''photons'')')
         write(*,'(10x,''-------'')')
         do  jp=-4,4
            ja = npot * (jp+4)
            write(*,'(1p,e7.0,1x,0p,15f8.5,:,/)')
     +        10.d0**jp,(fpot**i,i=1,npot-1)
           write(*,'(16i8)') (jpartic(i,1),i=ja,ja+npot-1)
         enddo
         write(*,'(10x,''positrons'')')
         write(*,'(10x,''---------'')')
         do  jp=-4,4
            ja = npot * (jp+4)
            write(*,'(1p,e7.0,1x,0p,15f8.5,:,/)')
     +        10.d0**jp,(fpot**i,i=1,npot-1)
            write(*,'(16i8)') (jpartic(i,2),i=ja,ja+npot-1)
         enddo
         write(*,'(10x,''electrons'')')
         write(*,'(10x,''---------'')')
         do  jp=-4,4
            ja = npot * (jp+4)
            write(*,'(1p,e7.0,1x,0p,15f8.5,:,/)')
     +        10.d0**jp,(fpot**i,i=1,npot-1)
            write(*,'(16i8)') (jpartic(i,3),i=ja,ja+npot-1)
         enddo
      endif
      ! print total counts of particles:
      write(*,'(/,6x,''code name'',18x,''particles '',a)')
     +   cdata(nfil)(idat:ilen)
      write(*,'(6x,''---- ----'',18x,''--------- '',16(''-''))')
      do  l=1,200
         if ( qpartic(l) .gt. 0. ) write(*,'(i10,a19,f13.0)')
     +      l,qpatext(l),qpartic(l)
      enddo
      write(*,'(2x)')
      stop
      end
c=======================================================================
c
c     block data initialization
c
c=======================================================================
      block data blinit

      implicit double precision (a-h,o-z), integer (i-n)

      parameter (lenthin=6552)

      double precision aatm(5),batm(5),catm(5)
      double precision parmas(0:101),phead(30)
      dimension lpdat(0:21),qpartic(200),jpartic(0:201,32)
      real pdata(lenthin)
      common /atmos/aatm,batm,catm
      common /integ/lobs,lsho,ishu,irec,isho,isub,lpdat,jpartic
      common /utabl/pdata,parmas,phead,cpi180,c180pi,qpartic
      data aatm /-186.5562d0,  -94.919d0, 0.61289d0, 0.d0, .01128292d0/
      data batm /1222.6562d0,1144.9069d0, 1305.5948d0, 540.1778d0,0.d0/
      data catm /994186.38d0,878153.55d0,636143.04d0,772170.16d0,1.d-9/
      data parmas/
     + 2*0.0000d0 , 0.000511d0 , 0.000511d0 , 0.000000d0 , 0.105658d0,
     + 0.105658d0 , 0.134973d0 , 0.139568d0 , 0.139568d0 , 0.497671d0,
     + 0.493646d0 , 0.493646d0 , 0.939566d0 , 0.938272d0 , 0.938272d0,
     + 0.497671d0 , 0.5488d0   , 1.11563d0  , 1.18937d0  , 1.192550d0,
     + 1.19743d0  , 1.3149d0   , 1.32132d0  , 1.67243d0  , 0.939566d0,
     + 1.11563d0  , 1.18937d0  , 1.19255d0  , 1.19743d0  , 1.3149d0,
     + 1.32132d0  , 1.67243d0  , 1.7841d0   , 1.7841d0   , 1.8693d0,
     + 1.8693d0   , 1.8645d0   , 1.8645d0   , 1.9693d0   , 1.9693d0,
     + 2.2852d0   , 80.6d0     , 80.6d0     , 91.161d0   , 1.8770d0,
     + 2.817d0    , 3.755d0    , 0.0d0      , 0.0d0      , 0.0000d0,
     + 0.7669d0   , 0.7681d0   , 0.7681d0   , 1.2309d0   , 1.2323d0,
     + 1.2336d0   , 1.2349d0   , 1.2309d0   , 1.2323d0   , 1.2336d0,
     + 1.2349d0   , 0.0d0      , 0.0d0      , 0.0d0      , 0.0000d0,
     + 0.0d0      , 0.0d0      , 0.0d0      , 0.0d0      , 0.0000d0,
     + 0.5488d0   , 0.5488d0   , 0.5488d0   , 0.5488d0   , 0.105658d0,
     + 0.105658d0 , 0.0d0      , 0.0d0      , 0.0d0      , 22*0.0d0/
      end
c=======================================================================
c
c     analyze the contents of all 21 subblocks in a record
c
c=======================================================================
c
c        primary    energy    runnumb   simnrsh     #hadr    #muons
c        #gammas     #elec   #elecnkg    obslv    theta    phi
c         h1wkm    h1wgcmq  ( emymax    ehasum    ehamax    nha )
c
c=======================================================================
 
      subroutine blwork(iend,lenrec)

      implicit double precision (a-h,o-z), integer (i-n)

      parameter (lenthin=6552)
 
      double precision aatm(5),batm(5),catm(5)
      double precision parmas(0:101),phead(30)
      dimension lpdat(0:21),qpartic(200),jpartic(0:201,32)
      real pdata(lenthin)
      common /integ/lobs,lsho,ishu,irec,isho,isub,lpdat,jpartic
      common /utabl/pdata,parmas,phead,cpi180,c180pi,qpartic
      isub = 0
      lenblk = lenrec / 21
      lenpar = lenblk / 39

c-----------loop over subblocks-----------------------------------------
      do  948  lia=1,lenrec,lenblk
      isub = isub + 1
      if ( 211285.2 .le. pdata(lia).and.pdata(lia) .le. 211285.4 ) then
c----------------subblock run header------------------------------------
         lobs = int(pdata(lia+4))
         if (isho.le.0) then
            phead(6) = lobs
            do  908  i=1,lobs
               phead(6+i) = 1.e-2 * pdata(lia+4+i)
  908       continue
            if (lobs.ge.1.and.lobs.le.10) then
               if (lobs.gt.1) write(*,
     +            '(6x,''observation level:'',f10.2,'' meter'')')
     +            phead(6+lobs)
            else
               write(*,'(6x,''nr of observation level undefined.'')')
               ! iend = iend + 9876543
               goto 948 ! 949
            endif
         endif
c - - - - - - - - - - - observation level (meter) - - - - - - - - - - -
         phead(17) = 1.e-2 * pdata(lia+lobs)
c - - - - - - - - - - - run number  - - - - - - - - - - - - - - - - - -
         phead(18) = pdata(lia+1)
      elseif ( 217433.0 .le. pdata(lia).and.pdata(lia) .le. 217433.2 )
     +   then
c----------------subblock event header----------------------------------
         if (iend.lt.0) goto 948
         isho = isho + 1
c- - - - - - - - - - simulated shower number - - - - - - - - - - - - -
         phead(1) = pdata(lia+1)
c - - - - - - - - - - height of first interaction (km)  - - - - - - - -
         phead(2) = 1.e-5 * pdata(lia+6)
c - - - - - - - - - - height of first interaction (grams/cm2)- - - - -
         pdata(lia+4) = thickgr(dble(pdata(lia+6)))
         phead(3) = pdata(lia+4)
c - - - - - - - - - - primary particle code - - - - - - - - - - - - - -
         phead(4) = pdata(lia+2)
c - - - - - - - - - - primary particle energy in gev  - - - - - - - - -
         phead(5) = pdata(lia+3)
c - - - - - - - - - - - phi angle this shower in radian - - - - - - - -
         phead(19) = pdata(lia+11)
c - - - - - - - - - - theta angle this shower in radian - - - - - - - -
         phead(20) = pdata(lia+10)
         timev0 = ( phead(2) * 1.e+5 - phead(6+lobs)*100. ) / 
     /      cos(phead(20)) / 29.9792458e0
         pdata(lia+37) = timev0
         do  912  l=lia,lia+lenblk-1,lenpar
**          write(9,'(1p,8e14.6)') (pdata(i+l),i=0,lenpar-1)
  912    continue
c----------------subblock longi information-----------------------------
      elseif (52815.2 .le. pdata(lia).and.pdata(lia) .le. 52815.4) then
         do  915  l=lia,lia+lenblk-1,lenpar
**          write(9,'(1p,8e14.6)') (pdata(i+l),i=0,lenpar-1)
  915    continue
c----------------subblock event end-------------------------------------
      elseif ( 3397.3 .le. pdata(lia).and.pdata(lia) .le. 3397.5 ) then 
cc       write(*,*) '       isub =',isub,' (evte)    irec =',irec,
cc  +      '    isho =',isho
c----------------subblock run end---------------------------------------
      elseif ( 3301.3 .le. pdata(lia).and.pdata(lia) .le. 3301.5 ) then
         goto 949
      else
c-----------subblock with particle data---------------------------------
         do  917  l=lia,lia+lenblk-1,lenpar
            icode = int(1.e-3*(pdata(l)*1.0000001))
            if ( icode.ge.1 .and. icode.le.200 ) then               
               qpartic(icode) = qpartic(icode) + 1
               if ( icode .le. 32 ) then
                  equad = parmas(icode)*parmas(icode)
     +                  + pdata(l+1)*pdata(l+1)
     +                  + pdata(l+2)*pdata(l+2) + pdata(l+3)*pdata(l+3)
                  if ( equad .gt. 0 ) then 
                     epart = sqrt(equad) 
                     eplog = log10(epart) + 4. ! multiply by 1.e4
                     jlog = eplog * 16. ! 16 log bins per decade
                     if ( jlog .lt. 0 ) jlog = 0
                     if ( jlog .gt. 201 ) jlog = 201 
                     jpartic(jlog,icode) = jpartic(jlog,icode) + 1 
                     ! write(*,*) sngl(epart), sngl(eplog), icode, jlog  
                  endif  
               endif
            endif
  917    continue
      endif
      if (iend.gt.0) goto 949
  948 continue
 
c-----------all records analyzed or end of showers reached--------------
  949 continue
      return
      end
c=======================================================================
 
      double precision function thickgr( arg )
 
c-----------------------------------------------------------------------
c  calculates thickness (g/cm**2) of atmosphere depending on height (cm)
c  argument:    arg    = height in cm
c-----------------------------------------------------------------------
 
      double precision aatm(5),batm(5),catm(5),arg
      common /atmos/aatm,batm,catm
 
      if     ( arg .lt. 4.d5 ) then
         thickgr = aatm(1) + batm(1) * exp( -arg/catm(1) )
      elseif ( arg .lt. 1.d6 ) then
         thickgr = aatm(2) + batm(2) * exp( -arg/catm(2) )
      elseif ( arg .lt. 4.d6 ) then
         thickgr = aatm(3) + batm(3) * exp( -arg/catm(3) )
      elseif ( arg .lt. 1.d7 ) then
         thickgr = aatm(4) + batm(4) * exp( -arg/catm(4) )
      else
         thickgr = aatm(5) - arg * catm(5)
      endif
 
      return
      end
