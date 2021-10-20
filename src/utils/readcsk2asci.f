c=======================================================================
c
c  r e a d c s k 2 a s c i . f
c  --------------------------- 
c      write corsika shower data as readable ascii file, needs 
c      factor 3.5 of disk space as of binary particle data file.
c      This is impossible for 64bit simulations on all machines.
c      End of ascii print out may be incomplete because of
c      differencies in reading particle data of 32bit simulations
c      on 64bit machines.
c-----------------------------------------------------------------------
c      gfortran -fbounds-check readcsk2asci.f -o readcsk2asci 
c      f77 -fbounds-check readcsk2asci.f -o readcsk2asci
c      ifort -C -check bounds readcsk2asci.f -o readcsk2asci
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c      input-files:
c           unit=3: current corsika particle data file.
c           unit=*: number of showers and file name(s):
c           ------------------------------------------------
c                       1          'total_number_of_showers'
c                       1          'total_number_of_files'
c            '/lxdata/d2lx14/joe/DAT045216'
c                       1          
c           ------------------------------------------------
c     output-files: 
c           unit=*: protocol output.
c           unit=9: ascii file named DATnnnnnn.ascithin or ....ascistnd.
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c=======================================================================
  
      program readcsk2asci
 
      implicit double precision (a-h,o-z), integer (i-n) 

      parameter (lenthin=6552,lenstnd=5733) 

      character cout*200,crunh*4,cevte*4
      character cdata(50)*200,cdat*200,cblk*200
      double precision aatm(5),batm(5),catm(5)
      double precision parmas(0:101),phead(30)
      dimension lpdat(0:21)
      real pdata(lenthin),qdata(936),prunh,pevte
      equivalence (crunh,prunh),(cevte,pevte)
      common /atmos/aatm,batm,catm
      common /integ/lobs,lsho,ishu,irec,isho,isub,lpdat
      common /utabl/pdata,parmas,phead,cpi180,c180pi
      data crunh/'RUNH'/,cevte/'EVTE'/

c--initialize some quantities-------------------------------------------
      cblk='                                                  '
      cdat=cblk
      cpi = 4.d0 * atan(1.d0)
      cpi180 = cpi/180.d0
      c180pi = 180.d0/cpi
      isho = 0
      iend = 0
      irec = 0
 
c--read run parameters including file names-----------------------------
      read(*,*,end=498,err=498) lsho
      read(*,*,end=498,err=498) nfil
      ! write(*,'(22x,''total number of showers:'',i7)') lsho
      ! write(*,'(22x,''total number of files  :'',i7)') nfil
      if ( nfil .gt. 50 ) goto 498
      do  ifil=1,nfil
         read(*,*,end=498,err=498) cdat
         read(*,*,end=498,err=498) nfsh
         cdata(ifil) = cdat
      enddo
 
c----------one or more showers in big disk file-------------------------
      do  444  ifil=1,nfil
      if (ifil.gt.1) close(unit=3)
      if (ifil.gt.1) close(unit=9)
      irec = 0
      idat = index(cdata(ifil),'DAT')
      ilen = index(cdata(ifil),' ') - 1
      write(*,'(/,8x,''shower file '',a)') cdata(ifil)(1:ilen)
* - - - - - - read data record with lenstnd words to test run - - - -
      open(unit=3,file=cdata(ifil)(1:ilen),
     +     status='old',form='unformatted')
      ishift = 2
      read(unit=3,err=496,end=432) (pdata(i),i=1,lenstnd)
      close(unit=3)
      write(*,*) '              pdata(1:3)',(pdata(i),i=1,3)
      write(*,*) '              pdata(4:6)',(pdata(i),i=4,6)
* - - - - - - check on reading 32bit simulation on 64bit machine:
      if ( 211285.2 .lt. pdata(1) .and.
     +     pdata(1) .lt. 211285.4 ) then
         ishift = 0
      elseif ( 211285.2 .lt. pdata(2) .and.
     +     pdata(2) .lt. 211285.4 ) then
         ishift = -1
         ! shift contents one element to the left: 
         do  i=1,936-1
            pdata(i) = pdata(i+1)
         enddo
      elseif ( 2.0202 .lt. pdata(3) .and. pdata(3) .lt. 9.9999 ) then 
         ! check version number instead of testing (273+1) and (312+1).
         ishift = 1
         ! shift contents one element to the right:
         do  i=936,2,-1
            pdata(i) = pdata(i-1)
         enddo
         pdata(1) = prunh ! 211285.2812500000000000;
      endif
* - - - - - - detect `standard` simulation instead of `thinning`:
      if ( 217433.0 .lt. pdata(273+1) .and.
     +                   pdata(273+1) .lt. 217433.2 ) then
         cout = cdata(ifil)(1:ilen)//'.ascistnd'
         lenrec = lenstnd
      elseif ( 217433.0 .lt. pdata(312+1) .and.
     +                       pdata(312+1) .lt. 217433.2 ) then
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
      write(*,'(20x,a)') cout(1:ilen+9)
      if ( ishift .eq. -1 ) goto 443 
* - - - - - - check on `joe` path - - - -
      if ( index(cout(1:ilen+9),'joe') .le. 0 ) then
         idat = index(cout(1:ilen+9),'DAT')
         cout(1:ilen-idat+1+9) = cout(idat:ilen+9)
         ilen = ilen-idat+1 
         do  i=ilen+10,200
            cout(i:i) = ' '
         enddo
      endif
* - - - - - - read data record with lenrec words - - - -
      open(unit=9,file=cout(1:ilen+9),status='unknown',form='formatted')
      open(unit=3,file=cdata(ifil),status='old',form='unformatted')
  431 continue
      read(unit=3,err=496,end=432) (pdata(i),i=1,lenrec)
      irec = irec + 1
* - - - - - - check on original 64bit simulation:
      pdatone = pdata(1)
      if ( irec .gt. 1 .and. ( 3.21d-41 .lt. pdatone
     +                 .and. pdatone .lt. 3.23d-41 ) ) goto 443
      if (mod(irec,200) .eq. 0) write(*,*) '       irec =',irec
* - - - - - - - - - - -
      if ( ishift .eq. 1 ) then
* - - - - - - shift data corresp. to a 32bit simul. on a 64bit machine:
         do  i=lenrec,2,-1
            pdata(i) = pdata(i-1)
         enddo
         if ( irec .gt. 1 ) then
            pdata(1) = 177177.1
            ipart = int(pdata(lenpar+1) * 1.000001e-3)
            if ( ipart .ge. 1 .and. ipart .le. 3 )
     +         pdata(1) = 1000. + mod(pdata(lenpar+1),1000.)
            if ( ipart .eq. 0 ) pdata(1) = 1091. ! fixed generation.
            if ( ipart .eq. 5 ) pdata(1) = 70000. + pdata(lenpar+1)
            if ( ipart .eq. 6 ) pdata(1) = 70000. + pdata(lenpar+1)
            if ( ipart .eq. 75 ) pdata(1) = 177177.1
            if ( ipart .eq. 76 ) pdata(1) = 177177.1
            if ( 3301.2 .lt. pdata(lenblk+1) .and.
     +         pdata(lenblk+1) .lt. 3301.4 ) pdata(1) = pevte ! 3397.39185
         else ! irec = 1.
            pdata(1) = prunh ! 211285.281
         endif
      elseif ( ishift .eq. -1 ) then
* - - - - - - shift data corresp. to a 64bit simul. on a 32bit machine:
         do  i=1,lenrec-1
            pdata(i) = pdata(i+1)
         enddo
         pdata(lenrec) = -pdata(lenrec-lenpar) 
         ! check negative elements carefully.
      endif
* - - - - - - - - - - -
      call blwork(iend,lenrec)
      if ( iend .le. 0 ) goto 431
  432 continue
* - - - - - - end of corsika particle data file reached: 
      write(*,*) '       isub =',isub,' (rune)    irec =',irec,
     +      '    irwc =',819*irec
      if ( 0 .lt. isub .and. isub .lt. 21 ) then
         do  i=1+isub*39,819
            write(9,'(1p,8e14.6)') (0.,l=0,lenpar-1) 
         enddo
      endif
  443 continue
* - - - - - - distinguish contents of pdata vector:
      write(*,*)
     +'    ____________________________________________________________'
      if ( ishift .eq. 1 ) then
        write(*,*) '    Fortran: impossible to read more than one rec'
        write(*,*) '        of the simulation, better use C/C++ code.'
      elseif ( ishift .eq. -1 ) then
        write(*,*) '    Fortran: is impossible to read more than one'
        write(*,*) '        record of a 64bit simulation; better use a'
        write(*,*) '        32bit run to get the full ascii output or'
        write(*,*) '        use a C/C++-program to read particle data.'
        write(*,*) '     rec=1.01',(pdata(i),i=1,6)
        write(*,*) '         1.02',(pdata(i),i=lenblk+1,lenblk+6)
        write(*,*) '         1.03',(pdata(i),i=lenblk*2+1,lenblk*2+6)
        write(*,*) '         1.04',(pdata(i),i=lenblk*2+lenpar+1,
     +                                         lenblk*2+lenpar+6)
        open(unit=9,file=cout(1:ilen+9),
     +       status='unknown',form='formatted')
        do  l=1,lenblk*3,lenpar
           write(9,'(1p,8e14.6)') (pdata(i+l),i=0,lenpar-1)
        enddo
        close(unit=9)
      endif
* - - - - - - closing print of kind of simulation:
      if ( ishift .ne. 0 ) write(*,'(9x,''ishift ='',i3)') ishift 
      if ( qdata(274)+qdata(547) .lt. 1. ) then
         write(*,'(/,8x,''`thinning run`'',/)')
      else
         write(*,'(/,8x,''`standard run`'',/)')
      endif
  444 continue
      close(unit=3)
      goto 499
 
c--end of data----------------------------------------------------------
  496 continue
      write(*,*) '       irec =',irec
      write(*,*)
     +'    ____________________________________________________________'
      write(*,'(5x,''ERROR: the estimated number of particles is'',
     +   i11)') 819*irec-80
      write(*,'(9x,''and subblocks EVTE and RUNE not found.'')')
      if ( qdata(274)+qdata(547) .lt. 1. ) then
         write(*,'(/,8x,''`thinning run`'',/)')
      else
         write(*,'(/,8x,''`standard run`'',/)')
      endif
      goto 499 
  498 continue
      write(*,*)
     +'    ____________________________________________________________'
      write(*,'(5x,''ERROR: missing some input parameters.'',/)')
  499 continue
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
      dimension lpdat(0:21)
      real pdata(lenthin)
      common /atmos/aatm,batm,catm
      common /integ/lobs,lsho,ishu,irec,isho,isub,lpdat
      common /utabl/pdata,parmas,phead,cpi180,c180pi
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
 
      double precision parmas(0:101),phead(30)
      dimension lpdat(0:21)
      real pdata(lenthin)
      common /integ/lobs,lsho,ishu,irec,isho,isub,lpdat
      common /utabl/pdata,parmas,phead,cpi180,c180pi
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
               iend = iend + 9876543
               goto 949
            endif
         endif
c - - - - - - - - - - - observation level (meter) - - - - - - - - - - -
         phead(17) = 1.e-2 * pdata(lia+lobs)
c - - - - - - - - - - - run number  - - - - - - - - - - - - - - - - - -
         phead(18) = pdata(lia+1)
         do  909  l=lia,lia+lenblk-1,lenpar
            write(9,'(1p,8e14.6)') (pdata(i+l),i=0,lenpar-1)
  909    continue
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
         pdata(lia+4) = sngl(thickgr(dble(pdata(lia+6))))
         phead(3) = pdata(lia+4)
c - - - - - - - - - - primary particle code - - - - - - - - - - - - - -
         phead(4) = pdata(lia+2)
c - - - - - - - - - - primary particle energy in gev  - - - - - - - - -
         phead(5) = pdata(lia+3)
c - - - - - - - - - - - phi angle this shower in radian - - - - - - - -
         phead(19) = pdata(lia+11)
c - - - - - - - - - - theta angle this shower in radian - - - - - - - -
         phead(20) = pdata(lia+10)
         timev0 = ( phead(2) * 1.d+5 - phead(6+lobs)*100. ) / 
     /      cos(phead(20)) / 29.9792458d0
         pdata(lia+37) = sngl(timev0)
         do  912  l=lia,lia+lenblk-1,lenpar
            write(9,'(1p,8e14.6)') (pdata(i+l),i=0,lenpar-1)
  912    continue
c----------------subblock longi information-----------------------------
      elseif (52815.2 .le. pdata(lia).and.pdata(lia) .le. 52815.4) then
         do  915  l=lia,lia+lenblk-1,lenpar
            write(9,'(1p,8e14.6)') (pdata(i+l),i=0,lenpar-1)
  915    continue
c----------------subblock event end-------------------------------------
      elseif ( 3397.3 .le. pdata(lia).and.pdata(lia) .le. 3397.5 ) then 
         if (iend.lt.0) then ! ignore first showers.
            iend = iend + 1
            isho = iend
            goto 948
         endif
         do  913  l=lia,lia+lenblk-1,lenpar
            write(9,'(1p,8e14.6)') (pdata(i+l),i=0,lenpar-1)
  913    continue
         write(*,*) '       isub =',isub,' (evte)    irec =',irec,
     +      '    isho =',isho
c----------------subblock run end---------------------------------------
      elseif ( 3301.3 .le. pdata(lia).and.pdata(lia) .le. 3301.5 ) then
         iend = iend + iend + isho
         do  914  l=lia,lia+lenblk-1,lenpar
            write(9,'(1p,8e14.6)') (pdata(i+l),i=0,lenpar-1)
  914    continue
         goto 949
      else
c-----------subblock with particle data---------------------------------
         do  917  l=lia,lia+lenblk-1,lenpar
            if ( irec .le. 2 )
     +         write(9,'(1p,8e14.6)') (pdata(i+l),i=0,lenpar-1)
            ! icode = int(1.e-3*(pdata(l)*1.0000001))
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
