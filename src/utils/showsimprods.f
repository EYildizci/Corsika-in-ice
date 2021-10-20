c=======================================================================
c
c  s h o w s i m p r o d s . f
c  ---------------------------
c     reading corsika particle data files of special user `simprod`
c     with file names like `pre16m316_17` or `sie17m178_09` a.s.o. 
c     a one line per file tabular will be created with the following
c     quantities 
c         Primary, lg(E), theta, phi, nsh, runnr, size,
c             obslvme, h1stme, thilev, thiwmax, lg(thirad),
c             verspgm, models, rundate, Xmagn, Zmagn;
c     on 32bit machines all 64bit simulations will be detected
c     correctly and in the list marked by `_64`; simulation protocol
c     files are assumed to be `DATnnnnnn.lst`,
c     otherwise additional conditions must be implemented.
c-----------------------------------------------------------------------
c compilation:
c     gfortran -fbounds-check showsimprods.f -o showsimprods
c     f77 -fbounds-check showsimprods.f -o showsimprods
c execution:
c     ls -l hee* | ./showsimprods
c     ls -l [c,f,g,h,m,p,s]?[e,E]* | ./showsimprods > showsimprods.table
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
c         Primary   lg(E)  theta    phi  nsh  runnr    size  obslvme  h1stme   .....
c proton _64   14   17.50   37.9 -138.0   1  000044    90.8  -500.00  -24880.  .....
c _stackin_     4   15.09    0.0    0.0   1  169051     2.4   194.00   18765.  .....
c Fluorine   1909   16.00    0.0    0.0   1  199070     3.9  1416.51   22222.  .....
c Iron       5626   16.00    0.0    0.0   1  199079     3.9  1413.82   18350.  .....
c proton       14   16.00   30.0   -3.3   1  199080     3.4  1429.25   22224.  .....
c Manganese  5525   16.00   30.0   -3.3   1  199082    10.1  1428.13   22222.  .....
c-----------------------------------------------------------------------
c        ncsumm(16)/ 28117, 15811, 8891, 5000,   2812, 1581, 889, 500,
c                      281,   158,   89,   50,     28,   16,   9,   5/
c-----------------------------------------------------------------------
c     1.0000E+14....1.7783E+14      28117 sh  =  18 * 1500 sh  + 1117 sh
c     1.7783E+14....3.1623E+14      15811 sh  =  15 * 1000 sh  +  811 sh
c     3.1623E+14....5.6234E+14       8891 sh  =  22 *  400 sh  +   91 sh
c     5.6234E+14....1.0000E+15       5000 sh  =  25 *  200 sh
c     1.0000E+15....1.7783E+15       2812 sh  =  28 *  100 sh  +   12 sh
c     1.7783E+15....3.1623E+15       1581 sh  =  26 *   60 sh  +   21 sh
c     3.1623E+15....5.6234E+15        889 sh  =  29 *   30 sh  +   19 sh
c     5.6234E+15....1.0000E+16        500 sh  =  25 *   20 sh
c     1.0000E+16....1.7783E+16        281 sh  =  25 *   11 sh  +    6 sh
c     1.7783E+16....3.1623E+16        158 sh  =  31 *    5 sh  +    3 sh
c     3.1623E+16....5.6234E+16         89 sh  =  44 *    2 sh  +    1 sh
c     5.6234E+16....1.0000E+17         50 sh  =  50 *    1 sh
c     1.0000E+17....1.7783E+17         28 sh  =  28 *    1 sh
c     1.7783E+17....3.1623E+17         16 sh  =  16 *    1 sh
c     3.1623E+17....5.6234E+17          9 sh  =   9 *    1 sh
c     5.6234E+17....1.0000E+18          5 sh  =   5 *    1 sh
c-----------------------------------------------------------------------

      program showsimprods

      implicit real (a-h,o-z), integer(i-n) 
      parameter (nfmx=40000,nchx=160)

      character cdata(nfmx)*160,cdat*160,czeile*80,crunh*4
      character qpatext(1:200)*19,chemical(100)*12,chfmti*5,chaspec*1
      character chuser*7

      dimension fsize(nfmx),nfsho(nfmx),ntskf(nfmx),nbits(nfmx)
      dimension pdata(5733),qrunh(312),qevth(312)
      dimension nevent(16),ncount(16),nclast(16)

      equivalence (crunh,prunh)
      data crunh/'RUNH'/, chaspec/' '/, chuser/'simprod'/

      data chemical/
     c ' Hydrogen   ',' Helium     ',' Lithium    ',' Beryllium  ',
     c ' Boron      ',' Carbon     ',' Nitrogen   ',' Oxygen     ',
     c ' Fluorine   ',' Neon       ',' Sodium     ',' Magnesium  ',
     c ' Aluminium  ',' Silicon    ',' Phosphorus ',' Sulfur     ',
     c ' Chlorine   ',' Argon      ',' Potassium  ',' Calcium    ',
     c ' Scandium   ',' Titanium   ',' Vanadium   ',' Chromium   ',
     c ' Manganese  ',' Iron       ',' Cobalt     ',' Nickel     ',
     c ' Copper     ',' Zinc       ',' Gallium    ',' Germanium  ',
     c ' Arsenic    ',' Selenium   ',' Bromine    ',' Krypton    ',
     c ' Rubidium   ',' Strontium  ',' Yttrium    ',' Zirconium  ',
     c ' Niobium    ',' Molybdenum ',' Technetium ',' Ruthenium  ',
     c ' Rhodium    ',' Palladium  ',' Silver     ',' Cadmium    ',
     c ' Indium     ',' Tin        ',' Antimony   ',' Tellurium  ',
     c ' Iodine     ',' Xenon      ',' Caesium    ',' Barium     ',
     c ' Lanthanum  ',' Cerium     ',' Praseodym. ',' Neodymium  ',
     c ' Promethium ',' Samarium   ',' Europium   ',' Gadolinium ',
     c ' Terbium    ',' Dysprosium ',' Holmium    ',' Erbium     ',
     c ' Thulium    ',' Ytterbium  ',' Lutetium   ',' Hafnium    ',
     c ' Tantalum   ',' Tungsten   ',' Rhenium    ',' Osmium     ',
     c ' Iridium    ',' Platinum   ',' Gold       ',' Mercury    ',
     c ' Thallium   ',' Lead       ',' Bismuth    ',' Polonium   ',
     c ' Astatine   ',' Radon      ',' Francium   ',' Radium     ',
     c ' Actinium   ',' Thorium    ',' Protactin. ',' Uranium    ',
     c ' Neptunium  ',' Plutonium  ',' Americium  ',' Curium     ',
     c ' Berkelium  ',' Californium',' Einsteinium','            '/

      data (qpatext(i),i=1,60)/
     c' gamma             ',' positron          ',' electron          ',
     c' _stackin_         ',' muon+             ',' muon-             ',
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
     c' eta=>pi+pi-pi0    ',' eta=>pi+pi-gamma  ',40*'                ',
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

      data nevent/
     +   1500, 1000, 400, 200, 100, 60, 30, 20, 11, 5, 2, 1, 1, 1, 1, 1/
      data ncount/
     +    18,   15,  22,  25,  28, 26, 29, 25, 25, 31,44,50,28,16, 9, 5/
      data nclast/
     +   1117,  811,  91,   0, 12, 21, 19,  0,  6,  3, 1, 0, 0, 0, 0, 0/

      qpatext(176) = ' Stnd              '
      qpatext(196) = ' Thin              '
      cpi = 4.d0 * atan(1.d0)
      cpi180 = cpi/180.d0
      c180pi = 180.d0/cpi

c--read first file name and fix position of file size-------------------
      ifil = 1
      read(*,'(a)',end=433,err=433) cdat
      if ( index(cdat,'-rw') .le. 0 ) then
         write(*,'(/,14x,''Use of `ls -1 ...` not supported, but'',
     +      '' `ls -l ...` is ok.'',/)') 
         goto 445
      endif
      ibl1 = index(cdat,' ')
      ibl2 = index(cdat(ibl1+2:nchx),' ') ! after -rwx-`s 2 blanks.
      ibl2 = ibl2 + ibl1 + 1
      ibl3 = index(cdat(ibl2+1:nchx),' ')
      ibl3 = ibl3 + ibl2
      ibl4 = index(cdat(ibl3+1:nchx),' ')
      ibl4 = ibl4 + ibl3
      imon = index(cdat,'Jan') ! test string of date and time.
      if ( imon .eq. 0 ) imon = index(cdat,'Feb')
      if ( imon .eq. 0 ) imon = index(cdat,'Mar')
      if ( imon .eq. 0 ) imon = index(cdat,'Apr')
      if ( imon .eq. 0 ) imon = index(cdat,'May')
      if ( imon .eq. 0 ) imon = index(cdat,'Jun')
      if ( imon .eq. 0 ) imon = index(cdat,'Jul')
      if ( imon .eq. 0 ) imon = index(cdat,'Aug')
      if ( imon .eq. 0 ) imon = index(cdat,'Sep')
      if ( imon .eq. 0 ) imon = index(cdat,'Oct')
      if ( imon .eq. 0 ) imon = index(cdat,'Nov')
      if ( imon .eq. 0 ) imon = index(cdat,'Dec')
      if ( imon .eq. 0 ) imon = ibl4 + index(cdat(ibl4+1:nchx),'-') - 4
      iright= imon - 2
      ileft = ibl4 + 1
      write(chfmti(1:5),'(''(i'',i2,'')'')') iright-ileft+1
      read(cdat(ileft:iright-3),chfmti) isize
      fsize(ifil) = 1.d-3 * isize ! now in Mbytes.
      if ( fsize(ifil) .lt. 0.1 ) fsize(ifil) = 0.1

c--check structure of file names (with / or without):
      iddot = index(cdat,':')
      isla = index(cdat,'/')
      islb = 0
      if ( isla .gt. 0 ) then
         if ( iddot .gt. 0 ) isla = iddot + 4
         islb = nchx
  401    continue
         islb = islb - 1
         if ( cdat(islb:islb) .ne. '/' ) goto 401         
      endif

c--check length of particle data file name:
      ilen = nchx
  404 continue 
      ilen = ilen - 1
      if ( cdat(ilen:ilen) .eq. ' ' ) goto 404
      itask = index(cdat,'DAT')
      if ( itask .le. 0 ) then
      if ( chuser .eq. 'simprod' ) then
         if ( index(cdat,'pre').gt.0 .or. index(cdat,'prE').gt.0 ) then
            ipre = index(cdat,'pre') 
            if ( ipre .le. 0 ) ipre = index(cdat,'prE')
            if ( isla .gt. 0 ) ipre = isla
            cdata(ifil) = cdat(ipre:ilen)
            icsk = ilen - ipre + 1 
         elseif (index(cdat,'fee').gt.0.or.index(cdat,'feE').gt.0) then
            ifee = index(cdat,'fee')
            if ( ifee .le. 0 ) ifee = index(cdat,'feE')
            if ( isla .gt. 0 ) ifee = isla
            cdata(ifil) = cdat(ifee:ilen)
            icsk = ilen - ifee + 1 
         elseif (index(cdat,'coe').gt.0.or.index(cdat,'coE').gt.0) then
            icoe = index(cdat,'coe')
            if ( icoe .le. 0 ) icoe = index(cdat,'coE')
            if ( isla .gt. 0 ) icoe = isla
            cdata(ifil) = cdat(icoe:ilen)
            icsk = ilen - icoe + 1 
         elseif (index(cdat,'hee').gt.0.or.index(cdat,'heE').gt.0) then
            ihee = index(cdat,'hee') 
            if ( ihee .le. 0 ) ihee = index(cdat,'heE')
            if ( isla .gt. 0 ) ihee = isla
            cdata(ifil) = cdat(ihee:ilen)
            icsk = ilen - ihee + 1 
         elseif (index(cdat,'sie').gt.0.or.index(cdat,'siE').gt.0) then
            isie = index(cdat,'sie')
            if ( isie .le. 0 ) isie = index(cdat,'siE')
            if ( isla .gt. 0 ) isie = isla
            cdata(ifil) = cdat(isie:ilen)
            icsk = ilen - isie + 1 
         endif

c--check name of file to get number of simulated showers:
         if ( cdata(ifil)(icsk-2:icsk-2) .eq. '_' ) then
            read(cdata(ifil)(icsk-8:icsk-7),'(i2)') mpot
            read(cdata(ifil)(icsk-5:icsk-3),'(i3)') mant
            read(cdata(ifil)(icsk-1:icsk),'(i2)') msimu
         elseif ( cdata(ifil)(icsk-3:icsk-3) .eq. '_' ) then
            read(cdata(ifil)(icsk-9:icsk-8),'(i2)') mpot
            read(cdata(ifil)(icsk-6:icsk-4),'(i3)') mant
            read(cdata(ifil)(icsk-2:icsk),'(i3)') msimu
         endif
         is = 4. * log10(1.001e-2*mant*10.**(mpot-9)) - 19.
         isho = 1
         if ( msimu .le. ncount(is) ) then
            isho = nevent(is)
         else 
            isho = nclast(is)
         endif
         nfsho(ifil) = isho
         naires = 0
      endif
      endif

c--read run parameters including file names (more than 1 file)----------
      do  ifil=2,nfmx
         read(*,'(a)',end=433,err=433) cdat
         read(cdat(ileft:iright-3),chfmti) isize
         fsize(ifil) = 1.d-3 * isize ! now in Mbytes.
         if ( fsize(ifil) .lt. 0.1 ) fsize(ifil) = 0.1
c--check length of particle data file name:
         ilen = nchx 
  411    continue
         ilen = ilen - 1
         if ( cdat(ilen:ilen) .eq. ' ' ) goto 411
         itask = index(cdat,'DAT')
         if ( itask .le. 0 ) then
         if ( chuser .eq. 'simprod' ) then
          if ( index(cdat,'pre').gt.0 .or. index(cdat,'prE').gt.0 ) then
            ipre = index(cdat,'pre')            
            if ( ipre .le. 0 ) ipre = index(cdat,'prE')
            if ( isla .gt. 0 ) ipre = isla
            cdata(ifil) = cdat(ipre:ilen)
            icsk = ilen - ipre + 1
          elseif (index(cdat,'fee').gt.0.or.index(cdat,'feE').gt.0) then
            ifee = index(cdat,'fee') 
            if ( ifee .le. 0 ) ifee = index(cdat,'feE')
            if ( isla .gt. 0 ) ifee = isla
            cdata(ifil) = cdat(ifee:ilen)
            icsk = ilen - ifee + 1
          elseif (index(cdat,'coe').gt.0.or.index(cdat,'coE').gt.0) then
            icoe = index(cdat,'coe') 
            if ( icoe .le. 0 ) icoe = index(cdat,'coE')
            if ( isla .gt. 0 ) icoe = isla
            cdata(ifil) = cdat(icoe:ilen)
            icsk = ilen - icoe + 1
          elseif (index(cdat,'hee').gt.0.or.index(cdat,'heE').gt.0) then
            ihee = index(cdat,'hee') 
            if ( ihee .le. 0 ) ihee = index(cdat,'heE')
            if ( isla .gt. 0 ) ihee = isla
            cdata(ifil) = cdat(ihee:ilen)
            icsk = ilen - ihee + 1
          elseif (index(cdat,'sie').gt.0.or.index(cdat,'siE').gt.0) then
            isie = index(cdat,'sie')
            if ( isie .le. 0 ) isie = index(cdat,'siE')
            if ( isla .gt. 0 ) isie = isla
            cdata(ifil) = cdat(isie:ilen)
            icsk = ilen - isie + 1
          endif
c--check name of file to get number of simulated showers:
            if ( cdata(ifil)(icsk-2:icsk-2) .eq. '_' ) then
               read(cdata(ifil)(icsk-8:icsk-7),'(i2)') mpot
               read(cdata(ifil)(icsk-5:icsk-3),'(i3)') mant
               read(cdata(ifil)(icsk-1:icsk),'(i2)') msimu
            elseif ( cdata(ifil)(icsk-3:icsk-3) .eq. '_' ) then
               read(cdata(ifil)(icsk-9:icsk-8),'(i2)') mpot
               read(cdata(ifil)(icsk-6:icsk-4),'(i3)') mant
               read(cdata(ifil)(icsk-2:icsk),'(i3)') msimu
            endif
            is = 4. * log10(1.001e-2 * mant * 10.**(mpot-9)) - 19.
            isho = 1
            if ( msimu .le. ncount(is) ) then
               isho = nevent(is)
            else
               isho = nclast(is)
            endif
            nfsho(ifil) = isho
            naires = 0
         endif
         endif
      enddo ! end-of loop ifil=2,nfmx.
  433 continue
      nfil = ifil - 1

c--print title lines----------------------------------------------------
      write(*,'(/,9x,''primary   lg(E)  theta   phi   nsh'',
     +   3x,''runnr    size  obslvme  h1stme  thilev  wmax  thirad'',
     +   2x,''verspgm    models   rundate  Xmagn  Zmagn'',/)')

c--work on all particle data files--------------------------------------
      nstop = 0
      do  444  ifil=1,nfil
      if ( ifil .gt. 1 ) close(unit=3)
      nstop = nstop + 1
      open(unit=3,file=cdata(ifil),status='unknown',form='unformatted')
      read(unit=3,end=444,err=444) pdata ! 5733 elements.

c - - - - - check on 64 or 32 bit and standard or thinning corsika:
      ibit = 32
      if ( abs(pdata(1)).lt.1.e-6 .and.
     +   ( pdata(2).ge.211284. .and. pdata(2).le.211286. ) ) ibit = 64
      if ( pdata(273).ge.217432. .and. pdata(273).le.217434.) ibit = -64
      if ( pdata(312).ge.217432. .and. pdata(312).le.217434.) ibit = -64
      nbits(ifil) = ibit
      ibit = ibit / 64
      isubr = 312
      if ( pdata(ibit+274).ge.217432. .and. pdata(ibit+274).le.217434.)
     +   isubr=273
      qrunh(1) = prunh
      do  i=2,isubr
         qrunh(i) = pdata(ibit+i)
      enddo
      do  i=isubr+1,isubr*2
         qevth(i-isubr) = pdata(ibit+i)
      enddo
      if ( qevth(11) .lt. 0.01 ) qevth(12) = 0.0

c - - - - - check on primary particle:
      if ( 0 .lt. qevth(3) .and. qevth(3) .lt. 200. ) then
         ip = int(qevth(3))
      elseif ( qevth(3) .lt. 6000. ) then
         ip = 200
         qpatext(200) = chemical(int(mod(qevth(3),100.)))//'       '
         if ( qevth(3).eq.201 ) qpatext(200) = ' Deuteron          '
         if ( qevth(3).eq.301 ) qpatext(200) = ' Tritium           '
      else
         ! write(*,*) '       invalid particle id ',qevth(3)
         ! should not but may occur by converted Aires simulations.
      endif
      chaspec = ' ' ! marker `*` for aires simulation run.

c - - - - - check models and date of simulation run - - - - -
      models = 0
      do  i=73,80
         models = models + 10**(80-i) * int(qevth(i))
      enddo
      imont = mod(qevth(45),1.e4)
      ijahr = qevth(45) * 1.00001e-4
      if ( ijahr < 30 ) ijahr = ijahr + 2000
      if ( ijahr < 100) ijahr = ijahr + 1900
      idate = 10000. * ijahr + imont
      if ( mod(idate,100) .gt. 31 ) idate = 31 + 100 * int(idate/100)
      if ( mod(idate,100) .eq.  0 ) idate =  1 + idate 
      if ( qevth(148) .eq. 0. ) qevth(148) = 1.
      czeile(1:11) = qpatext(ip)(1:11)
      if ( nbits(ifil) .eq. 64 ) czeile(9:11) = '_64'
      if ( qevth(152) .gt. 1. ) then
         qevth(152) = 0.01234 + log10(qevth(152)*1.e-2)
      else
         qevth(152) = 0.
      endif

c - - - - - set particle energy to lower limit of energy interval:
      if ( qevth(59) .ne. qevth(60) .and. nfsho(ifil) .gt. 1 )
     +   qevth(4) = qevth(59)
      if ( qevth(81) .ne. qevth(82) .and. nfsho(ifil) .gt. 1 ) then
         qevth(11) = qevth(81)*cpi180*1.001
         qevth(12) = qevth(82)*cpi180*1.001 
      endif

c - - - - - print out tabular quantities - - - - -
      if ( ip .gt. 0 ) then
        if ( qevth(159) .lt. 1.e-3 ) then
c - - - - - - - corsika simulation:
          if ( nfsho(ifil) .gt. 1 ) then
           if ( qevth(81) .ne. qevth(82) ) then
            write(*,'(a11,i5,a1,f7.2,i4,'':'',i2.2,'' 0:360'',i6,i8.6,
     +        f8.1,f9.2,f9.0,f6.1,f9.0,f6.1,f9.5,2i10,2f7.2,f7.2)')
     +        czeile(1:11),int(qevth(3)),chaspec,9.+log10(qevth(4)),
     +        int(c180pi*qevth(11)),int(c180pi*qevth(12)),nfsho(ifil),
     +        int(qevth(44)), fsize(ifil), 1.d-2*qevth(48),
     +        1.d-2*qevth(7), log10(qevth(148)+4.67735141e-34),
     +        qevth(151), qevth(152)+0.001, 1.e-5+qevth(46), models,
     +        idate, qevth(71), qevth(72)
cc   +        ,9.+log10(qevth(60)) ! do not print upper energy value.
           else
            write(*,'(a11,i5,a1,f7.2,f6.1,f7.1,i6,i8.6,f8.1,
     +        f9.2,f9.0,f6.1,f9.0,f6.1,f9.5,2i10,2f7.2,f7.2)')
     +        czeile(1:11),int(qevth(3)),chaspec,9.+log10(qevth(4)),
     +        c180pi*qevth(11), c180pi*qevth(12), nfsho(ifil),
     +        int(qevth(44)), fsize(ifil), 1.d-2*qevth(48),
     +        1.d-2*qevth(7), log10(qevth(148)+4.67735141e-34),
     +        qevth(151), qevth(152)+0.001, 1.e-5+qevth(46), models,
     +        idate, qevth(71), qevth(72)
cc   +        ,9.+log10(qevth(60))
           endif
          else         
            if ( qevth(59) .ne. qevth(60) ) then
              write(*,'(a11,i5,a1,f7.2,f6.1,f7.1,i6,i8.6,f8.1,
     +          f9.2,f9.0,f6.1,f9.0,f6.1,f9.5,2i10,2f7.2,f7.2)')
     +          czeile(1:11),int(qevth(3)),chaspec,9.+log10(qevth(4)),
     +          c180pi*qevth(11), c180pi*qevth(12), nfsho(ifil),
     +          int(qevth(44)), fsize(ifil), 1.d-2*qevth(48),
     +          1.d-2*qevth(7), log10(qevth(148)+4.67735141e-34),
     +          qevth(151), qevth(152)+0.001, 1.e-5+qevth(46), models,
     +          idate, qevth(71), qevth(72)
cc   +          ,9.+log10(qevth(60))
            else
              write(*,'(a11,i5,a1,f7.2,f6.1,f7.1,i6,i8.6,f8.1,
     +          f9.2,f9.0,f6.1,f9.0,f6.1,f9.5,2i10,2f7.2,f7.2)')
     +          czeile(1:11),int(qevth(3)),chaspec,9.+log10(qevth(4)),
     +          c180pi*qevth(11), c180pi*qevth(12), nfsho(ifil),
     +          int(qevth(44)), fsize(ifil), 1.d-2*qevth(48),
     +          1.d-2*qevth(7), log10(qevth(148)+4.67735141e-34),
     +          qevth(151), qevth(152)+0.001, 1.e-5+qevth(46), models,
     +          idate, qevth(71), qevth(72)
            endif
          endif
**   +,i3,3f7.0)') int(qevth(79)),qevth(160),qevth(161),qevth(162)
        else
c - - - - - - - converted aires simulation incl. Aires weight factor:
         write(*,'(a11,i5,a1,f7.2,f6.1,f8.2,i5,i8.6,f8.1,
     +   f9.2,f9.0,f6.1,f9.0,f6.1,f9.5,2i10,2f7.2,f11.6)')
     +   czeile(1:11),int(qevth(3)),chaspec,9.+log10(qevth(4)),
     +   c180pi*qevth(11), c180pi*qevth(12), nfsho(ifil),
     +   int(qevth(44)), fsize(ifil), 1.d-2*qevth(48),
     +   1.d-2*qevth(7), log10(qevth(148)+4.67735141e-34),
     +   qevth(151), qevth(152)+0.001, 1.e-5+qevth(46), models,
     +   idate, qevth(71), qevth(72)
     +   ,qevth(159)
        endif
      endif
      qpatext(200) = '                   '

c - - - - - end-of loop ifil=1,nfil.
  444 continue 
      close(unit=3)

c--print closing comment lines------------------------------------------
      if ( naires .le. 0 ) then
      write(*,'(/,9x,''primary   lg(E)  theta   phi   nsh'',
     +   3x,''runnr    size  obslvme  h1stme  thilev  wmax  thirad'',
     +   2x,''verspgm    models   rundate  Xmagn  Zmagn'',/)')
      else
      write(*,'(/,9x,''primary   lg(E)  theta   phi   nsh'',
     +   3x,''runnr    size  obslvme  h1stme  thilev  wmax  thirad'',
     +   2x,''verspgm    models   rundate  Xmagn  Zmagn'',
     +   2x,''Aireswfact'',/)')
      endif
  445 continue
      write(*,'(14x,''Total number of files:'',i7)') nstop

c--print explanation of model digits------------------------------------
      write(*,'(14x,''Appendix `_64`: simulation done by 64bit'',
     + '' executable, detected by 32bit processor, otherwise'',
     + '' no information possible.'')')
      write(*,'(14x,''Explanation of digits of `models`:'',
     + 2x,''(10^7): EGS flag;  (10^6): NKG flag;  (10^5):'',
     + '' lowEnergy flag, 1=gheisha, 2=urqmd,'')')
      write(*,'(14x,''3=fluka;  (10^4): highEnergy, 0=hdpm,'',
     + '' 1=venus, 2=sibyll, 3=qgsjet, 4=dpmjet, 5=nexus, 6=epos;'',
     + 2x,''(10^3): Cerenkov flag;'')')
      write(*,'(14x,''(10^2): Neutrino flag;   (10^1): Curved'',
     + '' flag, 0=standard, 1=Aires, 2=curved;   (10^0):'',
     + '' Computer, 3=unix, 4=macintosh.'')')
      write(*,'(14x,''thirad: radial thinning to 10^thirad meter.'',
     + '' An extra `*` at particle codes indicates a special primary'',
     + '' particle usage'')') 
      write(*,'(14x,''in the original Aires simulation; it has to be'',
     + '' converted to Corsika particle data structure (program '', 
     + '' `readciodemo6`)'')')
      write(*,'(14x,''before running `showsimulist`,'',
     + '' wherein the `models` reference number for Aires simulations'',
     + '' is 00130013.'',/)')

c--end of data----------------------------------------------------------
      stop
      end
