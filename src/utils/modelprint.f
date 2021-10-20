c=======================================================================
c
c  m o d e l p r i n t . f
c  ----------------------- 
c     display high and low energy models as texts of the given
c     corsika simulation.
c-----------------------------------------------------------------------
c compilation:
c      f77 -fbounds-check -m32 modelprint.f -o modelprint
c      gfortran -fbounds-check modelprint.f -o modelprint 
c      ifort -C modelprint.f -o modelprint
c execution:
c      ./modelprint 
c      < enter particle data file name >
c-----------------------------------------------------------------------
c      input-files:
c           unit=*: name of particle data file.
c           unit=3: current corsika particle data file.
c     output-files: 
c           unit=*: protocol output.
c-----------------------------------------------------------------------
c           runh=211285.281   evth=217433.078
c           long=52815.2969   evte=3397.39185   rune=3301.33252
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c=======================================================================
  
      program modelprint
 
      implicit double precision (a-h,o-z), integer (i-n) 

      parameter (lenthin=6552,lenstnd=5733) 

      character cout*120,crunh*4,cevte*4
      character cdat*120,cblk*120,qpatext(200)*20
      character cmodlow(3)*8,cmodhig(0:7)*8,cmodqgs(4)*5,cmodsib(4)*5
      double precision aatm(5),batm(5),catm(5)
      double precision parmas(0:101),phead(30)
      dimension lpdat(0:21)
      real pdata(lenthin),qdata(936),prunh,pevte
      equivalence (crunh,prunh),(cevte,pevte)
      common /atmos/aatm,batm,catm
      common /integ/lobs,lsho,ishu,irec,isho,isub,lpdat
      common /utabl/pdata,parmas,phead,cpi180,c180pi
      common /chrgd/ichargd(52),qpatext
      data crunh/'RUNH'/,cevte/'EVTE'/
      data cmodlow/'gheisha ','Urqmd   ','fluka   '/
      data cmodhig/' hdpm   ',' venus  ',' sibyll ','  qgsjet',
     +  '  dpmjet',' nexus  ',' epos   ','        '/
      data cmodqgs/'old  ','01c  ','-II-3','-II-4'/
      data cmodsib/'-1.6 ','-2.1 ','-2.3 ','     '/

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
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - - read next name of a corsika data file:
   10 continue
      read(*,*,err=498,end=498) cdat
      i = 120 + 1
   11 continue
      i = i - 1
      if ( cdat(i:i) .eq. ' ' ) goto 11
      lenchf = i
      idat = index(cdat(1:lenchf),'DAT') - 1
      if ( idat .lt. 0 ) then
         if ( cdat(1:1).eq.'/' ) then
            i = lenchf
   12       continue
            i = i - 1
            if ( cdat(i:i) .ne. '/' ) goto 12
            idat = i ! postion of last slash in the file name.
         endif
      endif
      itpa = index(cdat(1:lenchf),'.part')
      iusc = index(cdat(1:lenchf),'_')

c----------one or more showers in big disk file-------------------------
      ilost = 0
      iswit = 0
      irec = 0
      idat = index(cdat,'DAT')
      ilen = index(cdat,' ') - 1
* - - - - - - read data record with lenstnd words - - - -
      open(unit=3,file=cdat,status='old',form='unformatted')
      ishift = 2
      itype = 2
      read(unit=3,err=496,end=496) (pdata(i),i=1,lenstnd)
      close(unit=3)
* - - - - - - check on reading 32bit simulation on 64bit machine:
      if ( 211285.2 .lt. pdata(1) .and.
     +     pdata(1) .lt. 211285.4 ) then
         ishift = 0
      elseif ( 211285.2 .lt. pdata(2) .and.
     +     pdata(2) .lt. 211285.4 ) then
         ishift = -1
         do  i=1,936-1
            pdata(i) = pdata(i+1)
         enddo
      elseif ( 2.0202 .lt. pdata(3) .and. pdata(3) .lt. 9.9999 ) then 
         ! check version number instead of testing (273+1) and (312+1).
         ishift = 1
         do  i=936,2,-1
            pdata(i) = pdata(i-1)
         enddo
         pdata(1) = prunh ! 211285.281 
      endif
* - - - - - - detect `standard` simulation instead of `thinning`:
      if ( 217433.0 .lt. pdata(273+1) .and.
     +                   pdata(273+1) .lt. 217433.2 ) then
         itype = 0
         lenrec = lenstnd
      elseif ( 217433.0 .lt. pdata(312+1) .and.
     +                       pdata(312+1) .lt. 217433.2 ) then
         itype = 1
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
      if ( ishift .eq. -1 ) goto 443 
* - - - - - - check on original 64bit simulation:
      pdatone = pdata(1)
      if ( irec .gt. 1 .and. ( 3.21d-41 .lt. pdatone
     +                 .and. pdatone .lt. 3.23d-41 ) ) goto 443
* - - - - - - - - - - -
      if ( ishift .eq. 1 ) then
* - - - - - - shift data corresp. to a 32bit simul. on a 64bit machine:
         do  i=lenrec,2,-1
            pdata(i) = pdata(i-1)
         enddo
         if ( irec .gt. 1 ) then
            pdata(1) = 177177.1
            ipart = pdata(lenpar+1) * 1.000001e-3
            if ( ipart .ge. 1 .and. ipart .le. 3 )
     +         pdata(1) = 1000. + mod(pdata(lenpar+1),1000.)
            if ( ipart .eq. 0 ) pdata(1) = 1091. ! fixed generation.
            if ( ipart .eq. 5 ) pdata(1) = 70000. + pdata(lenpar+1)
            if ( ipart .eq. 6 ) pdata(1) = 70000. + pdata(lenpar+1)
            if ( ipart .eq. 75 ) pdata(1) = 177177.1
            if ( ipart .eq. 76 ) pdata(1) = 177177.1
            if ( 3301.2 .lt. pdata(lenblk+1) .and.
     +         pdata(lenblk+1) .lt. 3301.4 ) pdata(1) = pevte ! 3397.39185
            if ( pdata(1) .eq. 177177.1 ) then
               ilost = ilost + 1
            else
               iswit = iswit + 1
            endif
         else ! irec = 1.
            pdata(1) = prunh ! 211285.281
         endif
      elseif ( ishift .eq. -1 ) then
* - - - - - - shift data corresp. to a 64bit simul. on a 32bit machine:
         do  i=1,lenrec-1
            pdata(i) = pdata(i+1)
         enddo
         pdata(lenrec) = pdata(lenrec-lenpar)
         if ( 217433.0 .lt. pdata(lenrec-lenblk+1) .and.
     +         pdata(lenrec-lenblk+1) .lt. 217433.2 )
     +         pdata(lenrec) = -3.333333
         if ( 3397.3 .lt. pdata(lenrec-lenblk+1) .and.
     +         pdata(lenrec-lenblk+1) .lt. 3397.5 )
     +         pdata(lenrec) = -7.777777
         if ( 3301.2 .lt. pdata(lenrec-lenblk+1) .and.
     +         pdata(lenrec-lenblk+1) .lt. 3301.4 )
     +         pdata(lenrec) = -9.999999
      endif
  443 continue
* - - - - - - end of corsika particle data file reached: 
      modlow = pdata(lenblk+75)
      modhig = pdata(lenblk+76)
      modsib = pdata(lenblk+139)
      modqgs = pdata(lenblk+141)
      write(*,*) '    '
      if ( modhig .eq. 3 ) then
         write(*,*) '    ',cmodlow(modlow)
     +               ,'  ',cmodhig(modhig),cmodqgs(modqgs)
      elseif ( modhig .eq. 2 ) then
         write(*,*) '    ',cmodlow(modlow)
     +               ,'  ',cmodhig(modhig),cmodsib(modsib)
      else
         write(*,*) '    ',cmodlow(modlow),'  ',cmodhig(modhig)
      endif
      write(*,*) '    ______________________________________________'
      if ( qdata(274)+qdata(547) .lt. 1. ) then
         write(*,'(8x,''`thinning run`'',/)')
      else
         write(*,'(8x,''`standard run`'',/)')
      endif
  444 continue
  445 continue
      close(unit=3)
      goto 499
 
c--end of data----------------------------------------------------------
  496 continue
      write(*,*) '        irec = 1  pdata(273:275)',(pdata(i),i=273,275)
      write(*,*) '    ______________________________________________'
      write(*,*) '    ERROR: simulation type of corsika is `standard`.'
      goto 499 
  497 continue
      write(*,*) '        irec = 1  pdata(312:314)',(pdata(i),i=312,314)
      write(*,*) '    ______________________________________________'
      write(*,*) '    ERROR: simulation type of corsika is `thinning`.'
      goto 499 
  498 continue
      write(*,*) '    ______________________________________________'
      write(*,*) '    ERROR: missing some input parameters.'
  499 continue
      stop
      end
