c=======================================================================
c
c  c o r s 2 i n p u t . f
c  =======================
c     creates a (possible) corsika steering (i.e. input) file using the
c     existing particle data output file, also already valid for 
c     corsika simulation data from 64 bit machines and parallel simul.
c     some steering quantities may be defined individually.
c     data files from `simprod` simulations are also correctly read.
c-----------------------------------------------------------------------
c compilation:
c     gfortran -fbounds-check cors2input.f -o cors2input
c     f77 -fbounds-check cors2input.f -o cors2input
c     ifort -C -check bounds cors2input.f -o cors2input
c execution:
c     ./cors2input
c          < enter particle data file name >
c-----------------------------------------------------------------------
c     "standard" Corsika:
c     output format for particle output (blocklength = 22932+8 fixed)
c     each block consists of 21 subblocks of 273 words (5733).
c     "thinning" Corsika:
c     output format for particle output (blocklength = 26208+8 fixed)
c     each block consists of 21 subblocks of 312 words (6552).
c-------------------------------------------
c     Transfer indices of thinning elements:
c         evth(149) = thin(1); evth(151) = thin(2); evth(152) = thin(3);
c         evth(148) = thin(1)/thinh(1); evth(150) = thin(2)/thinh(2);
c     Reading CER-file: bs = pdata(l) ! bunchsize
c     Reading DAT-file: if (pdata(l) .gt. 9.9e6) bs = mod(pdata(l),1.e5)
c-----------------------------------------------------------------------
c           RUNH = 211285.2812500000000000;
c           EVTH = 217433.0781250000000000;
c           LONG =  52815.2968750000000000;
c           EVTE =   3397.3918457031250000;
c           RUNE =   3301.3325195312500000;
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c=======================================================================

      program cors2input

      character chfile*240,chflag(0:1)*7,czeile*240,cinput*240,chuser*7
      character qpatext(0:200)*20,qnuclei(0:29)*2,csklst*240
      character cfmtint*5,crunh*4
      dimension pdata(5733),drunh(312),devth(312)
     +          ,nevent(16),ncount(16),nclast(16),nbit(0:2)
      double precision aatm(5),batm(5),catm(5),thickgr
      double precision aatmax(5),batmax(5),catmax(5),atmlay(5)
      real prunh 
      equivalence (crunh,prunh) 
      logical lexist, lstexi
      data crunh/'RUNH'/,nbit/32,32,64/
      common /atmos/aatm,batm,catm
      data aatm /-186.5562d0,  -94.919d0, 0.61289d0, 0.d0, .01128292d0/
      data batm /1222.6562d0,1144.9069d0, 1305.5948d0, 540.1778d0,0.d0/
      data catm /994186.38d0,878153.55d0,636143.04d0,772170.16d0,1.d-9/
      data aatmax
     + /-129.86987d0, -13.912578d0, 1.1379052d0, -4.5502d-4, 1.12829d-2/
      data batmax
     + /1170.07784d0, 1310.69613d0, 1490.6966d0, 503.613568d0, 1.0d0/
      data catmax
     + /971950.04d0, 682326.92d0, 615751.06d0, 795110.76d0, 1.d-9/
      data atmlay /10.7d5, 14.6d5, 36.6d5, 100.d5, 112.8292d5/
      data chflag/'      F','      T'/,c180pi/57.2957795/
      data chuser/'prakt  '/ ! chuser/'simprod'/, chuser/'maximo'/
      data nevent/
     +   1500, 1000, 400, 200, 100, 60, 30, 20, 11, 5, 2, 1, 1, 1, 1, 1/
      data ncount/
     +    18,   15,  22,  25,  28, 26, 29, 25, 25, 31,44,50,28,16, 9, 5/
      data nclast/
     +   1117,  811,  91,   0, 12, 21, 19,  0,  6,  3, 1, 0, 0, 0, 0, 0/
      data cfmtint/'(i10)'/ 
      data qnuclei/'  ',
     +             'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     +             'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     +             'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu'/
      data (qpatext(ip),ip=0,60)/'                   ',
     c' gamma             ',' positron          ',' electron          ',
     c' (stackin simulat.)',' muon +            ',' muon -            ',
     c' pion 0            ',' pion +            ',' pion -            ',
     c' Kaon 0 long       ',' Kaon +            ',' Kaon -            ',
     c' neutron           ',' proton            ',' anti proton       ',
     c' Kaon 0 short      ','                   ',' Lambda            ',
     c' Sigma +           ',' Sigma 0           ',' Sigma -           ',
     c' Xi 0              ',' Xi -              ',' Omega -           ',
     c' anti neutron      ',' anti Lambda       ',' anti Sigma -      ',
     c' anti Sigma 0      ',' anti Sigma +      ',' anti Xi 0         ',
     c' anti Xi +         ',' anti Omega +      ',16*'                ',
     c'                   ',' omega             ',' rho 0             ',
     c' rho +             ',' rho -             ',' Delta ++          ',
     c' Delta +           ',' Delta 0           ',' Delta -           ',
     c' anti Delta --     ',' anti Delta -      ',' anti Delta 0      '/
      data (qpatext(ip),ip=61,129)/
     c' anti Delta +      ',' Kaon * 0          ',' Kaon * +          ',
     c' Kaon * -          ',' anti Kaon * 0     ',' electron neutrino ',
     c' anti elec neutrino',' muon neutrino     ',' anti muon neutrino',
     c'                   ',' eta=>2*gamma      ',' eta=>3*pi0        ',
     c' eta=>pi+ pi- pi0  ',' eta=>pi+ pi- gamma',40*'                ',
     c'                   ',' D 0               ',' D +               ',
     c' anti D -          ',' anti D 0          ',' D s +             ',
     c' anti D s -        ',' eta c             ',' D * 0             ',
     c' D * +             ',' anti D * -        ',' anti D * 0        ',
     c' D * s +           ',' anti D * s -      ','                   '/
      data (qpatext(ip),ip=130,200)/
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

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - - read next name of a corsika data file:
   10 continue
      read(*,*,err=18,end=19) chfile
      i = 240 + 1
   11 continue
      i = i - 1
      if ( chfile(i:i) .eq. ' ' ) goto 11
      lenchf = i
      idat = index(chfile(1:lenchf),'DAT') - 1
      if ( idat .lt. 0 ) then
         if ( chfile(1:1).eq.'/' ) then
            i = lenchf
   12       continue
            i = i - 1
            if ( chfile(i:i) .ne. '/' ) goto 12
            idat = i ! position of last slash in the file name.
         endif
      endif
      itpa = index(chfile(1:lenchf),'.part')
      iusc = index(chfile(1:lenchf),'_')
      chuser = 'prakt  '
      if ( itpa .gt. 3 ) then
         cinput = chfile(idat+1:lenchf-5)//'.inp'
      elseif ( iusc .gt. 9 ) then ! simprod simulation found.
            chuser = 'simprod'
            write(cfmtint(3:4),'(i2)') lenchf-iusc
            read(chfile(iusc+1:lenchf),cfmtint) msimu
            if ( chfile(iusc-4:iusc-4).eq.'M' .or.
     +           chfile(iusc-4:iusc-4).eq.'m'
     +        .or. chfile(iusc-7:iusc-7).eq.'E' .or.
     +           chfile(iusc-7:iusc-7).eq.'e' )
     +         then
               cinput = chfile(iusc-9:lenchf)//'.inp'
            else
               cinput = chfile(1:lenchf)//'.inp'
            endif
      else ! idat>0 found.
         cinput = chfile(idat+1:lenchf)//'.inp'
      endif

c - - - - - - - - binary corsika data file:
      open(3,file=chfile(1:lenchf),status='old',form='unformatted')
      read(3,err=16) (pdata(i),i=1,5733)
      close(3)

c - - - - - check on 64 or 32 bit and standard or thinning corsika:
      if ( pdata(1).ge.211284. .and. pdata(1).le.211286. ) then 
         ibit = 0
         if ( pdata(313).ge.217432. .and. pdata(313).le.217434. ) then
            isubr = 312
         elseif (pdata(274).ge.217432..and.pdata(274).le.217434.) then
            isubr = 273
         endif
      elseif ( abs(pdata(1)).lt.1.e-6 ) then 
         ibit = 1
         if ( pdata(2).ge.211284. .and. pdata(2).le.211286. ) then
            if ( pdata(314).ge.217432. .and. pdata(314).le.217434.) then
               isubr = 312
            elseif (pdata(275).ge.217432..and.pdata(275).le.217434.)then
               isubr = 273
            endif
         endif
      elseif ( pdata(3).ge.2.0202 .and. pdata(3).lt.9.9999 ) then
         ibit = -1
         if ( pdata(312).ge.217432. .and. pdata(312).le.217434.) then
            isubr = 312
         elseif (pdata(273).ge.217432..and.pdata(273).le.217434.) then
            isubr = 273
         endif
      else
         write(*,*) '     pdata(1) =',pdata(1),' case should not occur!'
      endif
      ! copy run header to short array drunh:
      drunh(1) = prunh ! 211285.281
      do  i=2,isubr
         drunh(i) = pdata(ibit+i)
      enddo
      ! exclude invalid default dates, set date to 01.Jan.2011:
      if ( 22. .lt. drunh(3)/10000. .and. drunh(3)/10000. .lt. 99. )
     +   drunh(3) = 110101.
      ! copy event header to short array devth:
      do  i=isubr+1,isubr*2
         devth(i-isubr) = pdata(ibit+i)
      enddo
c - - - - - - use quantities of runh end evth subblocks:
      if ( devth(7) .lt. 0. ) devth(7) = -devth(7) ! curved version
      iday = mod(int(drunh(3)),100)
      imon = mod(int(drunh(3)),10000) / 100
      iyer = int(drunh(3)) / 10000 + 1900
      if ( iyer .lt. 1988 ) iyer = iyer + 100
      if ( drunh(4) .gt. 1.e6 ) drunh(4) = drunh(4) * 1.e-6
      write(*,'(/,4x,''runh'',9x,''runnumb'',6x,''date'',9x,'//
     +   '''versprog'',5x,''nobslev'',6x,''obslev1'',6x,''obslev2'','//
     +   '6x,''obslev3'')')
      write(*,'(1p,e13.5,0p,f11.0,2x,1p,6e13.5,/)') (drunh(l),l=1,8)
      write(*,'(4x,''obslev4'',6x,''obslev5'',6x,''obslev6'','//
     +   '6x,''obslev7'',6x,''obslev8'',6x,''obslev9'',5x,'//
     +   '''obsl0/nsh'',5x,''slope'')')
      write(*,'(1p,8e13.5,/)') (drunh(l),l=9,16)
      write(*,'(4x,''engymin'',6x,''engymax'',6x,''flagEGS4'','//
     +   '5x,''flagNKG'',6x,''ecut.hadr'',4x,''ecut.muon'','//
     +   '4x,''ecut.elec'',4x,''ecut.phot'')')
      write(*,'(1p,8e13.5,/)') (drunh(l),l=17,24)
      write(*,'(4x,''evth'',9x,''evtnumb'',6x,''particle'',5x,'//
     +   '''energy'',7x,''altit/gr'',5x,''nrtarget'',5x,'//
     +   '''height/cm'',4x,''moment.x'')')
      write(*,'(1p,8e13.5,/)') (devth(l),l=1,8)
      write(*,'(4x,''moment.y'',5x,''moment-z'',5x,''theta/deg'',4x,'//
     +   '''phi/deg'',6x,''nr.seeds'',5x,''seed1'',8x,''offset1a'','//
     +   '5x,''offset1b'')')
      write(*,'(1p,2e13.5,0p,f11.3,f13.3,2x,1p,4e13.5/)')
     +(devth(l),l=9,10),(devth(l)*57.29578,l=11,12),(devth(l),l=13,16)
      write(*,'(4x,''seed2'',8x,''offset2a'',5x,''offset2b'',5x,'//
     +   '''seed3'',8x,''offset3a'',5x,''offset3b'',5x,''seed4'','//
     +   '8x,''offset4a'')')
      write(*,'(1p,8e13.5,/)') (devth(l),l=17,24)
      if ( isubr .eq. 312 ) then
         write(*,'(i14,'' bit thinning simulation. lsubrec ='',i4)') 
     +      nbit(1+ibit),isubr
      elseif ( isubr .eq. 273 ) then
         write(*,'(i14,'' bit standard simulation. lsubrec ='',i4)') 
     +      nbit(1+ibit),isubr
      endif
      write(*,'(/,12x,a,/)') cinput(1:lenchf+4)

c - - - - - test on DATnnnnnn.lst file to get best steering infos:
      csklst = chfile(1:lenchf)//'.lst'
      inquire(file=csklst,exist=lstexi)
      iprm = int(devth(3))
      open(7,file=cinput,form='formatted',
     +     access='sequential',status='unknown')
**    if ( iprm .eq. 4 ) then

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if ( lstexi ) then
         open(unit=2,file=csklst,form='formatted',status='old')
   24    continue     
         read(2,'(a)',end=22,err=23) czeile(1:60)
         if ( index(czeile,'RUNNR ').le.0 ) goto 24
         ic = 60
   20    continue
         ic = ic - 1
         if ( czeile(ic:ic) .eq. ' ' ) goto 20
         write(7,'(a)') czeile(1:ic)
         do  iq=2,1234
            read(2,'(a)',end=22,err=23) czeile(1:60)
            ic = 60
   21       continue
            ic = ic - 1
            if ( czeile(ic:ic) .eq. ' ' ) goto 21
            write(7,'(a)') czeile(1:ic)
            if ( index(czeile,'EXIT').ge.1 ) then
               close(2)
               goto 23
            endif
         enddo
   22    continue
   23    continue
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      else
      write(7,'(''RUNNR  '',i10,4x,''version'',f9.6,'//
     +   '4x,''rundate'',i3.2,''.'',i2.2,''.'',i4)')
     +   int(drunh(2)),drunh(4),iday,imon,iyer
      if ( drunh(2) .lt. 10000. .and. devth(13) .ge. 6. ) then  
         icut = int(2.d0 * log10(drunh(17)*1.000001) - 14.d0) 
         drunh(80) = 4.5432d5 * 1.95583d0**icut
         write(7,'(''PARALLEL'',f10.0,f12.0,''   1   F'')')
     +      drunh(80)/1.d3,drunh(80)
      endif
      if ( chuser .ne. 'simprod' ) then
         write(7,'(''EVTNR  '',i10)') int(devth(2))
      else
         write(7,'(''EVTNR  '',i10,24x,a)')
     +      int(devth(2)),cinput(1:lenchf+4)
      endif
      is = devth(13)
      do  ii=1,is
         write(7,'(''SEED   '',3i10)') (int(devth(13+3*(ii-1)+i)),i=1,3)
      enddo
c - - - - - - - primary, energy, angles - - - - - - - - - - - - - - - -
      if ( ( 1 .le. iprm .and. iprm .le. 74 ) .or.
     +   ( 116 .le. iprm .and. iprm .le. 194 ) ) then
         write(7,'(''PRMPAR '',i10,23x,a)') iprm,qpatext(iprm)
      elseif ( 201 .le. iprm .and. iprm .le. 5629 ) then
         write(7,'(''PRMPAR '',i10,22x,i4,''-'',a)')
     +      iprm,int(iprm/100),qnuclei(mod(iprm,100))
      elseif ( ( mod(iprm,100) .le. 0 .or. 29 .lt. mod(iprm,100) ) .or.
     +   ( drunh(3) .lt. 100. ) ) then
         write(*,*) '  Warning: better using version cors2input64',
     +   ' for correct counting of elements. Stop! drunh(3)=',drunh(3)
         write(*,*) '  '
         goto 19
      endif
      write(7,'(''ESLOPE '',f13.2)') devth(58)
      write(7,'(''ERANGE '',1p,2e13.4)') devth(59),devth(60)
      write(7,'(''THETAP '',2f13.2)') devth(81),devth(82)
      write(7,'(''PHIP   '',2f13.2)') devth(83),devth(84)
      isho = 1
      iobs = devth(47)
      if ( devth(79) .eq. 1. .and. iobs .lt. 10 ) then
         if ( drunh(15) .gt. 0 ) isho = int(drunh(15))  
      elseif ( devth(79) .eq. 0. .or. devth(79) .eq. 2. ) then 
         if ( chuser .ne. 'simprod' ) then
            ! check lst file to get number of simulated showers:
            csklst = chfile(1:lenchf)//'.lst'
            inquire(file=csklst,exist=lstexi)
            if ( lstexi ) then 
               open(unit=2,file=csklst,form='formatted',status='old')
               do  iq=1,12345 
                  read(2,'(a)',end=13,err=13) czeile
                  if ( index(czeile,'NUMBER OF GENERATED').ge.1 ) then
                     close(2)
                     read(czeile(30:40),'(i11)') isho
                     goto 13
                  endif
               enddo
            endif
   13       continue
         else ! user 'simprod' special selection of shower numbers.
            is = 4. * log10(devth(4)) - 19.
            isho = nevent(is)
            if ( lenchf-iu .eq. 2 ) then
               if ( msimu .gt. ncount(is) ) isho = nclast(is)
            else 
               isho = 1
            endif
         endif
      endif
      write(7,'(''NSHOW  '',i10)') isho
      do  i=1,iobs
         gramms = thickgr(dble(devth(47+i)))
         if ( gramms .lt. 900. ) gramms = gramms + 0.03
         write(7,'(''OBSLEV '',f13.2,''E2'',f13.3,'' g/cm^2'')')
     +      devth(47+i)/100.,gramms
      enddo
c - - - - - - - cuts, prints, flags - - - - - - - - - - - - - - - - - -
      write(7,'(''ECUTS     '',2f11.3,1p,2e10.1)') (devth(i),i=61,64)
      write(7,'(''ECTMAP          9.87654E10'')')
      write(7,'(''RADNKG '',f13.2,''E2'')') devth(147)/100.
      write(7,'(''MAXPRT '',i10)') min(1,int(devth(2)))
      write(7,'(''HADFLG '',6i5)') (int(devth(i)),i=65,70)
c - - - - - - - model quantities - - - - - - - - - - - - - - - - - - - -
      write(7,'(''ELMFLG '',2(3x,a7))')
     +   chflag(int(devth(73))),chflag(int(devth(74)))
      if (devth(76).eq.0.) then
         write(7,'(''* HDPM            T         0'')')
      endif
      if (devth(76).eq.1.) then
         write(7,'(''VENUS           T         0'')')
         if (devth(145).eq.1.) then
            write(7,'(''VENSIG          T'')')
         elseif (devth(145).eq.2.) then
            write(7,'(''NEXSIG          T'')')
         else
            write(*,*) ' check elements `evth` 140,142,144 '
         endif
         write(7,'(''VENPAR ''''      ''''         0.'')')
      endif
      if (devth(76).eq.2.) then
         write(7,'(''SIBYLL          T         0'')')
         if (devth(140).ge.1.) then
            write(7,'(''SIBSIG          T'')')
         else
            write(*,*) ' check elements 142,144,145 '
         endif
      endif
      if (devth(76).eq.3.) then
         write(7,'(''QGSJET          T         0'')')
         if (devth(142).ge.1.) then
            write(7,'(''QGSSIG          T'')')
         else
            write(*,*) ' check elements 140,144,145 '
         endif
      endif
      if (devth(76).eq.4.) then
         write(7,'(''DPMJET          T         0'')')
         if (devth(144).ge.1.) then
            write(7,'(''DPJSIG          T'')')
         else
            write(*,*) ' check elements 140,142,145 '
         endif
      endif
      if (devth(76).eq.5.) then
         write(7,'(''NEXUS           T         0'')')
         if (devth(145).eq.1.) then
            write(7,'(''VENSIG          T'')')
         elseif (devth(145).eq.2.) then
            write(7,'(''NEXSIG          T'')')
            ! init input files for nexus use:
            write(7,'(''NEXPAR fname inics nexus/nexus.inics'')')
            write(7,'(''NEXPAR fname iniev nexus/nexus.iniev'')')
            write(7,'(''NEXPAR fname initl nexus/nexus.initl'')')
            write(7,'(''NEXPAR fname inirj nexus/nexus.inirj'')')
            ! dummy out files for epos (debug case):
            write(7,'(''NEXPAR fname check none'')')
            write(7,'(''NEXPAR fname histo none'')')
            write(7,'(''NEXPAR fname data  none'')')
            write(7,'(''NEXPAR fname copy  none'')')
         else
            write(*,*) ' check elements 140,142,144 '
         endif
      endif
      if (devth(76).eq.6.) then
         write(7,'(''EPOS            T         0'')')
         ! init input files for epos use:
         write(7,'(''EPOPAR input epos/epos.param      '')')
         write(7,'(''EPOPAR fname inics epos/epos.inics'')')
         write(7,'(''EPOPAR fname iniev epos/epos.iniev'')')
         write(7,'(''EPOPAR fname initl epos/epos.initl'')')
         write(7,'(''EPOPAR fname inirj epos/epos.inirj'')')
         write(7,'(''EPOPAR fname inihy epos/epos.ini1b'')')
         ! dummy out files for epos (debug case):
         write(7,'(''EPOPAR fname check none'')')
         write(7,'(''EPOPAR fname histo none'')')
         write(7,'(''EPOPAR fname data  none'')')
         write(7,'(''EPOPAR fname copy  none'')')
      endif 
c - - - - - - - logicals, hilow, longi, magnet - - - - - - - - - - - - -
      gramms = thickgr(dble(devth(7)))
      write(7,'(''MUMULT          T'')')
      write(7,'(''MUADDI    '',a7)') chflag(int(devth(94)))
      if ( devth(95) .lt. 1. ) devth(95) = 1.
      if ( devth(79) .eq. 1. .or.
     +   ( 2.02 .lt. drunh(4) .and. drunh(4) .lt. 2.12 ) ) devth(95)=1.
      write(7,'(''STEPFC '',f13.2)') devth(95)
      if ( devth(4) .lt. 3.2e6 .and. devth(13) .le. 3. .and.
     +     drunh(2) .lt. 100000. .and. chuser .ne. 'simprod' .and.
     +   ( devth(75) .eq. 1. .and. devth(76) .eq. 3 ) )
     +   write(7,'(''* PLOTSH          T'')') 
      if ( isho .eq. 1 ) then
         write(7,'(''* FIXHEI '',f13.2,''E2     0'',f13.3,'' g/cm^2'')')
     +      devth(7)/100.,gramms
         if ( devth(60)/devth(59) .gt. 1.0003 ) 
     +      write(7,'(''* ERANGE '',1p,2e13.4)') (devth(4),i=1,2)
         if ( devth(82)/devth(81) .gt. 1.0003 )
     +      write(7,'(''* THETAP '',2f13.2)') (devth(11)*c180pi,i=1,2)
         if ( devth(84)/devth(83) .gt. 1.0003 )
     +      write(7,'(''* PHIP   '',2f13.2)') (devth(12)*c180pi,i=1,2)
      endif 
      if ( devth(79) .eq. 1. .or.
     +   ( 2.02 .lt. drunh(4) .and. drunh(4) .lt. 2.12 ) )
     +   write(7,
     +      '(''* Aires simulation converted to corsika quantities.'')')
c - - - - - - - low energy model and hilow - - - - - - - - - - - - - - -
      if ( devth(75) .ge. 3. ) then
         hilow = 200.
         if ( devth(155) .ge. 80. ) hilow = devth(155)
         write(7,'(''* Low energy model Fluka used.'')')
      elseif ( devth(75) .eq. 2. ) then
         hilow = 80.
         if ( devth(155) .ge. 80. ) hilow = devth(155)
         write(7,'(''URQMD             T'')')
      elseif ( devth(75) .le. 1. ) then
         hilow = 80.
         if ( devth(155) .ge. 80. ) hilow = devth(155)
         write(7,'(''* Low energy model Gheisha used.'')')
      endif
      write(7,'(''HILOW'',f15.2)') hilow
c - - - - - - - check longitudinal step size - - - - - - - - - - - - - -
ccc   if ( drunh(2) .ge. 100000. ) then
         gramlong = 8.
         inquire(file=chfile(1:lenchf)//'.long',exist=lexist)
         if ( lexist ) then
            open(4,file=chfile(1:lenchf)//'.long',form='formatted',
     +           access='sequential',status='old')
            read(4,'(30x,i5,19x,f5.0)',err=15,end=15) longstep,gramlong
            goto 15
   14       continue
            gramlong = 12.34 
   15       continue
            close(4)
         endif
         write(7,'(''LONGI'',11x,''T'',f7.0,2(''     T''))') gramlong
ccc   endif
c - - - - - - - thinning quantities - - - - - - - - - - - - - - - - - -
      if (devth(150).gt.0.) then
         devth(149) = 1.00001d0 * devth(149)
         write(7,'(''THIN     '',1p,3e12.4)')
     +      devth(149),devth(151),devth(152)
         write(7,'(''THINH '',2f12.0)') 1.,100.
         !  devth(149)/devth(148),devth(151)/devth(150)
      endif
c - - - - - - - check quantities of magnetic field - - - - - - - - - - -   
      if ( devth(72) .gt. -17. .and. devth(72) .lt. -11. .and.
     +     devth(71) .gt.  18. .and. devth(71) .lt.  22. ) then
         write(7,'(''MAGNET  '',2f12.2,6x,''Auger'')')
     +      devth(71),devth(72)
      elseif ( devth(72) .gt. 41. .and. devth(72) .lt. 45. .and.
     +         devth(71) .gt. 18. .and. devth(71) .lt. 22. ) then
         write(7,'(''MAGNET  '',2f12.2,6x,''Karlsruhe'')')
     +      devth(71),devth(72)
      elseif ( abs(devth(71)) .lt. 1.e-2 .and.
     +         abs(devth(72)) .lt. 1.e-2 ) then
            write(7,'(''MAGNET  '',1p,2e13.3,''    NoMag'')')
     +         devth(71),devth(72)
      else
         write(7,'(''MAGNET  '',2f12.2)') devth(71),devth(72)
      endif
c - - - - - - - cherenkov quantities - - - - - - - - - - - - - - - - - -
      if (devth(77).gt.0.) then
         if (devth(96).gt.0.)
     +      write(7,'(''CWAVLG '',2f10.0)') devth(96),devth(97)
         if (devth(85).gt.0.)
     +      write(7,'(''CERSIZ '',f10.0)') devth(85)
         write(7,'(''CERFIL    '',a7)') chflag(int(devth(92)))
      endif
c - - - - - - - direct, host, user - - - - - - - - - - - - - - - - - - -
      if ( chuser.eq.'maximo' ) then
         write(7,'(''ATMOD   0'')')
         write(7,'(''ATMA  '',1p,4e15.7,e13.5)') (aatmax(i),i=1,5)
         write(7,'(''ATMB  '',1p,4e15.7,e13.5)') (batmax(i),i=1,5)
         write(7,'(''ATMC  '',1p,4e15.7,e13.5)') (catmax(i),i=1,5)
         write(7,'(''ATMLAY'',5(f9.3,''e5''))') (atmlay(i)*1.d-5,i=1,5)
         write(7,'(''DIRECT /data/corsdat5/maximo/'')')
         write(7,'(''USER   maximo'')')
      elseif ( chuser.eq.'simprod' ) then
         ! write(7,'(''DIRECT /lxdata/d1lx36/simprod/qgs_cont-2c/'')')
         write(7,'(''DIRECT /lxdata/d1lx36/simprod/'')')
         write(7,'(''USER   simprod'')')
      elseif ( devth(13) .ge. 6. ) then 
         ! write(7,'(''DIRECT /work/kit/ikp/jl5949/'')')
         write(7,'(''DIRECT csk'',i6.6,''/'')') int(drunh(2))
         write(7,'(''USER   joe'')')
      else 
         ! write(7,'(''DIRECT /fzk/cgwork/joe/aires/'')')
         write(7,'(''DIRECT /data/corsdat6/joe/corsika-run/'')')
         write(7,'(''USER   joe'')')
      endif
      if ( devth(13) .ge. 6. ) then
         write(7,'(''HOST   hc3n996'')')
      elseif ( devth(48) .ge. 1234.56e2 ) then
         write(7,'(''HOST   iwrcgvor2'')')
      else
         write(7,'(''HOST   iklx108'')')
      endif
      write(7,'(''EXIT'')')
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      endif
      close(7)
      goto 19 ! goto 10 ! only one file may be processed. 

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   16 continue
      write(*,'(''    err: file name: '',a,''!'')') chfile(1:lenchf)
      goto 19
   17 continue
      write(*,'(''    end: file name: '',a,''!'')') chfile(1:lenchf)
      goto 19
   18 continue
      write(*,'(''    err: file name: '',a,''!'')') chfile(1:lenchf)
   19 continue
      stop
      end
c=======================================================================

      double precision function heightcm( arg )

c-----------------------------------------------------------------------
c     height (cm) above sea level as function of gramms per cm^2
c-----------------------------------------------------------------------

      implicit double precision (a-h)
      double precision aatm(5),batm(5),catm(5),arg
      common /atmos/aatm,batm,catm

      if ( arg .gt. 631.1d0 ) then
        heightcm = catm(1) * log( batm(1) / (arg - aatm(1)) )
      elseif ( arg .gt. 271.7d0 ) then
        heightcm = catm(2) * log( batm(2) / (arg - aatm(2)) )
      elseif ( arg .gt. 3.0395d0 ) then
        heightcm = catm(3) * log( batm(3) / (arg - aatm(3)) )
      elseif ( arg .gt. 1.28292d-3 ) then
        heightcm = catm(4) * log( batm(4) / (arg - aatm(4)) )
      else
        heightcm = (aatm(5) - arg) / catm(5)
      endif

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
         thickgr = aatm(1) + batm(1) * exp( -arg / catm(1) )
      elseif ( arg .lt. 1.d6 ) then
         thickgr = aatm(2) + batm(2) * exp( -arg / catm(2) )
      elseif ( arg .lt. 4.d6 ) then
         thickgr = aatm(3) + batm(3) * exp( -arg / catm(3) )
      elseif ( arg .lt. 1.d7 ) then
         thickgr = aatm(4) + batm(4) * exp( -arg / catm(4) )
      else
         thickgr = aatm(5) - arg * catm(5)
      endif

      return
      end
