c=======================================================================
c
c  s h o w p a r a l l e l . f
c  ---------------------------
c     writing a tabular of simulation quantities to get an overview of
c     available parallel corsika simulation paths named csk00iiii/
c-----------------------------------------------------------------------
c compilation:
c     gfortran -fbounds-check showparallel.f -o showparallel
c     f77 -fbounds-check -m32 showparallel.f -o showparallel
c     ifort -C -check bounds showparallel.f -o showparallel
c execution:
c     # - - - - - create list of Job00*.out files:
c     ls -1 csk*/Job00*.out > showparallel.jobinfos
c     # - - - - - list pathes and Job info files:
c     ./showparallel < showparallel.jobinfos > showparallel.hc3-joblist
c-----------------------------------------------------------------------
c in case of csk002141/Job002141*.err does not contain number of files:
c ls -l csk002141/DAT00*.lst | wc | awk '{printf("%7d\n",$1)}' > csk002141/Job00*.err
c-----------------------------------------------------------------------
c     input-file: unit=*: showparallel.jobinfos
c                   csk001001/Job001001_390162.out
c                   csk001002/Job001002_390166.out
c-----------------------------------------------------------------------
c     Primary   lg(E)  theta    phi  runtsk  sizeGBy  procs
c           ectcut   ectmax  t(min)  files  RATIO  obslev  Xmagn  Zmagn
c               ecutha  ecutmu  ecutel  ecutga  thilev  wmax  lg(rad)
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
 
      program showparallel

      implicit double precision (a-h,o-z), integer (i-n)

      character chjobname(1000)*100, chzeile*100, chpfile*100
      character chtime*24, chtext*100, cfmtint*4, cfmtflt*7
      character cherrname(1000)*100, chstck*24

      real pdata(1000)
      dimension ecut(4), ecmnt(4), jcexp(4), jobnlen(1000)

      logical lexist

      data cfmtint/'(i9)'/, cfmtflt/'(f11.0)'/, mthi/0/, lthi/0/

c - - - - - - read all Job files and keep names in a character array - -
      do  job=1,1000
         read(*,'(a)',err=101,end=102) chjobname(job)
         jl = 100 + 1 
  100    jl = jl - 1
         if ( chjobname(job)(jl:jl) .eq. ' ' ) goto 100
         jobnlen(job) = jl       
      enddo
      goto 102
  101 continue
      write(*,*) ' e r r o r   reading chjobname ',chjobname(job)
  102 job = job - 1
      mcprev = 0
      enprev = 0.
      thprev = 0.
      phprev = 0.
      lthipr = 0

c - - - - - - read and check content of each Job file - - - - - - - - - 
      do  jj=1,job

      read(chjobname(jj)(4:9),'(i6)') mrun
      energy = 0.001
      ectcut = 1.e3
      ectmax = 1.e6
      sumgiga = 0.001
      do  j=1,4
         ecut(j) = 0.
      enddo   
      thinlev = 0.
      thiwmax = 0.
      thirad  = 0.
      obslev = 0. 
      bxmag = 1.e-5
      bzmag = 1.e-5
      theta = 0.
      phia = 0.
      phib = 0. 
      mprocs = 0 
      mcode = 0
      lthi = 0 
      inquire(file=chjobname(jj),exist=lexist)
      if ( .not. lexist ) goto 149 

      ! - - - get current number of .lst files from Job00*.err:
      jbl1 = index(chjobname(jj),'.out')
      cherrname(jj) = chjobname(jj)(1:jbl1)//'err'
      open(unit=1,file=cherrname(jj)(1:jobnlen(jj)),
     +            form='formatted',status='old')     
      ifiles = 50000
      read(1,*,end=103,err=103) ifiles
  103 continue
      close(unit=1)
      tminuts = 0.9876

      ! - - - read infos from Job00*.out: 
      open(unit=2,file=chjobname(jj)(1:jobnlen(jj)),
     +            form='formatted',status='old')
         
      do  lin=1,99 ! max number of lines in chjobname.

         read(2,'(a)',err=144,end=144) chzeile 

      ! - - - comment lines in steering file - - - - - - - - - - - - - -

         if ( chzeile(1:1) .eq. '*' .or. chzeile(1:2) .eq. 'c ' 
     +   .or. chzeile(1:2) .eq. 'C ' ) then
            jeng = index( chzeile, 'ERANGE')
            if ( jeng .ge. 2 ) then
               ! - - - - original primary energy (stck_in simulation)
               jbl1 = index( chzeile(jeng:100), ' ') + jeng - 1
  107          continue
               jbl1 = jbl1 + 1
               if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 107
               jbl2 = index( chzeile(jbl1:100), ' ') + jbl1 - 1
               write(cfmtflt(3:4),'(i2)') jbl2-jbl1
               read(chzeile(jbl1:jbl2-1),cfmtflt) energy
               deratio = 1.
               if ( ectmax .gt. 0. ) deratio = energy/ectmax
               energy = 9.+log10(energy)
            endif
            goto 143
         endif
         if ( chzeile(1:6) .eq. '      ' ) goto 143

      ! - - - Files: - - - - - - - - - - - - - - - - - - - - - - scc - -
         if ( index( chzeile, 'Real-time') .gt. 0 ) then
            read(chzeile(1:10),'(i10)') jtimbeg
            read(2,*,err=144,end=144) jtimend
            read(2,'(a)',err=144,end=144) chzeile
            jbl1 = index( chzeile, 'Files')
  108       continue
            jbl1 = jbl1 - 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 108
            write(cfmtflt(3:4),'(i2)') jbl1
            read(chzeile(1:jbl1),cfmtflt) qfiles
            read(2,*,err=144,end=144) tminuts
            mfiles = qfiles
            mprocs = qfiles
         endif

      ! - - - Exit and sum of file sizes - - - - - - - - - - - - - - - -
         if ( index( chzeile, 'EXIT') .ge. 1 ) then
            read(2,*,err=144,end=144) sumgiga
            goto 144
         endif

      ! - - - Command: - - - - - - - - - - - - - - - - - - - - - - - - - 
         if ( index( chzeile, 'Command:') .gt. 0 ) then
            if ( index( chzeile, ' -m') .le. 0 ) goto 149 
            ! - - - - number of processors - - - - - - - - - - - - - - - -
            jpcs = index( chzeile, ' -p') + 3 
            if ( jpcs .gt. 9 ) then
               do  j=0,100-jpcs
                  if ( chzeile(jpcs+j:jpcs+j) .ne. ' ' ) goto 109
               enddo
  109          continue  
               jpcs = jpcs + j
               jblk = index( chzeile(jpcs:100), ' ')
               write(cfmtint(3:3),'(i1)') jblk-1
               read(chzeile(jpcs:jpcs+jblk-2),cfmtint) mprocs
            endif
         endif

      ! - - - parallel parameters - - - - - - - - - - - - - - - - - - -
         jpar = index( chzeile, 'PARALLEL')
         if ( jpar .ge. 1 ) then
            jbl1 = index( chzeile, ' ')
  110       continue  
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 110 
            jbl2 = index( chzeile(jbl1:100), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1  
            read(chzeile(jbl1:jbl2-1),cfmtflt) ectcut
  111       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 111
            jbl3 = index( chzeile(jbl2:100), ' ') + jbl2 -1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) ectmax
         endif 

      ! - - - primary particle code - - - - - - - - - - - - - - - - - - 
         jcod = index( chzeile, 'PRMPAR')
         if ( jcod .eq. 1 ) then
            jbl1 = index( chzeile(jcod:100), ' ') ! + jcod
  112       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 112
            jbl2 = index( chzeile(jbl1:100), ' ') + jbl1 - 1
            write(cfmtint(3:3),'(i1)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtint) mcode
         endif

      ! - - - primary particle energy - - - - - - - - - - - - - - - - - 
         jeng = index( chzeile, 'ERANGE')
         if ( jeng .ge. 1 ) then
            jbl1 = index( chzeile(jeng:100), ' ')
  113       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 113
            jbl2 = index( chzeile(jbl1:100), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) energy
            deratio = 1.
            if ( ectmax .gt. 0. ) deratio = energy/ectmax
            energy = 9.+log10(energy)
         endif

      ! - - - theta angle in degrees - - - - - - - - - - - - - - - - - -
         jthe = index( chzeile, 'THETAP')
         if ( jthe .ge. 1 ) then
            jbl1 = index( chzeile, ' ')
  114       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 114
            jbl2 = index( chzeile(jbl1:100), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) theta
         endif

      ! - - - azimuth angle in degrees - - - - - - - - - - - - - - - - -
         jphi = index( chzeile, 'PHIP')
         if ( jphi .ge. 1 ) then
            jbl1 = index( chzeile, ' ')
  115       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 115
            jbl2 = index( chzeile(jbl1:100), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) phia
  116       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 116
            jbl3 = index( chzeile(jbl2:100), ' ') + jbl2 - 1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) phib
            if ( phia .lt. phib ) then
               phia = 123.45
               do  j=1,6
                  read(2,'(a)',err=144,end=144) chzeile
               enddo
               if ( index( chzeile, 'SEED') .ge. 1 ) then
                  jbl1 = index( chzeile, ' ')
  117             continue
                  jbl1 = jbl1 + 1
                  if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 117
                  jbl2 = index( chzeile(jbl1:100), ' ') + jbl1 - 1
                  write(cfmtint(3:3),'(i1)') jbl2-jbl1
                  read(chzeile(jbl1:jbl2-1),cfmtint) iseed
                  ! check particle data file for current phip.
                  chpfile = chjobname(jj)(1:19)//
     +               '-000000000-000000001'
                  chpfile(11:13) = 'DAT'
                  write(chpfile(21:29),'(i9.9)') iseed
                  open(unit=4,file=chpfile(1:39),
     +                 form='unformatted',status='old')
                  read(4,err=118,end=118) (pdata(j),j=1,546)
                  close(unit=4)
                  lblock = 273
                  if ( 217433.0 .lt. pdata(312+1) .and.
     +             pdata(312+1) .lt. 217433.2 ) lblock = 312 
                  phia = 57.2957795 * pdata(12+lblock) 
  118             continue          
               else
                  write(*,*) '   e r r o r   too few seeds.'
               endif
            endif
         endif

      ! - - - observation level - - - - - - - - - - - - - - - - - - - -
         jlev = index( chzeile, 'OBSLEV')
         if ( jlev .ge. 1 ) then
            jbl1 = index( chzeile, ' ')
  119       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 119
            jbl2 = index( chzeile(jbl1:100), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) obslev
            obslev = obslev / 100.
         endif

      ! - - - magnetic field components - - - - - - - - - - - - - - - -
         jmag = index( chzeile, 'MAGNET')
         if ( jmag .ge. 1 ) then
            jbl1 = index( chzeile, ' ')
  120       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 120
            jbl2 = index( chzeile(jbl1:100), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) bxmag
  121       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 121
            jbl3 = index( chzeile(jbl2:100), ' ') + jbl2 - 1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) bzmag
         endif 

      ! - - - energy cuts of four particle groups - - - - - - - - - - -
         ject = index( chzeile, 'ECUTS ')
         if ( ject .ge. 1 ) then
            jbl1 = index( chzeile, ' ')
  125       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 125
            jbl2 = index( chzeile(jbl1:100), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) ecut(1)
  126       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 126
            jbl3 = index( chzeile(jbl2:100), ' ') + jbl2 - 1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) ecut(2)
  127       continue
            jbl3 = jbl3 + 1
            if ( chzeile(jbl3:jbl3) .eq. ' ' ) goto 127
            jbl4 = index( chzeile(jbl3:100), ' ') + jbl3 - 1
            write(cfmtflt(3:4),'(i2)') jbl4-jbl3
            read(chzeile(jbl3:jbl4-1),cfmtflt) ecut(3)
  128       continue
            jbl4 = jbl4 + 1
            if ( chzeile(jbl4:jbl4) .eq. ' ' ) goto 128
            jbl5 = index( chzeile(jbl4:100), ' ') + jbl4 - 1
            write(cfmtflt(3:4),'(i2)') jbl5-jbl4
            read(chzeile(jbl4:jbl5-1),cfmtflt) ecut(4)
            do  j=1,4
               if ( ecut(j) .lt. 1. ) then
                  jcexp(j) = int(log10(ecut(j))*1.000003-0.999)
               else
                  jcexp(j) = int(log10(ecut(j))*1.000003)    
               endif  
               ecmnt(j) = ecut(j)/10.d0**jcexp(j)
            enddo  
         endif 

      ! - - - thinning specification - - - - - - - - - - - - - - - - - -
         jthi = index( chzeile, 'THIN ')
         if ( jthi .ge. 1 ) then
            jbl1 = index( chzeile, ' ')
  140       continue
            jbl1 = jbl1 + 1
            if ( chzeile(jbl1:jbl1) .eq. ' ' ) goto 140
            jbl2 = index( chzeile(jbl1:100), ' ') + jbl1 - 1
            write(cfmtflt(3:4),'(i2)') jbl2-jbl1
            read(chzeile(jbl1:jbl2-1),cfmtflt) thinlev
            thinlev = log10(thinlev+4.67735141e-34)
  141       continue
            jbl2 = jbl2 + 1
            if ( chzeile(jbl2:jbl2) .eq. ' ' ) goto 141
            jbl3 = index( chzeile(jbl2:100), ' ') + jbl2 - 1
            write(cfmtflt(3:4),'(i2)') jbl3-jbl2
            read(chzeile(jbl2:jbl3-1),cfmtflt) thiwmax
  142       continue
            jbl3 = jbl3 + 1
            if ( chzeile(jbl3:jbl3) .eq. ' ' ) goto 142
            jbl4 = index( chzeile(jbl3:100), ' ') + jbl3 - 1
            write(cfmtflt(3:4),'(i2)') jbl4-jbl3
            read(chzeile(jbl3:jbl4-1),cfmtflt) thirad
            thirad = log10(thirad+4.67735141e-34) - 2. ! meters
            mthi = mthi + 1 
            lthi = 1 
         endif       

      ! - - - end-of-loop lin=1,99 - - - - - - - - - - - - - - - - - - -
  143 continue
      enddo

  144 continue 
      close(unit=2)

      ! - - - if energy still 0.0: check DAT00nnnn.stck file - - - - - - 
      if ( energy .lt. 1.d0 ) then
         chstck = 'DAT000000.stck'
         write(chstck(4:9),'(i6.6)') mrun
         inquire(file=chstck,exist=lexist)
         if ( .not. lexist ) then 
            ! - - - - - - check stck file in subdirectory :
            chstck = 'csk000000/DAT000000.stck'
            write(chstck(4:9),'(i6.6)') mrun
            write(chstck(14:19),'(i6.6)') mrun
            inquire(file=chstck,exist=lexist)
         endif
         if ( lexist ) then
            open(unit=2,file=chstck,form='formatted',status='old')
            read(2,*) jstck, energy, jprim, h1stck 
            close(unit=2)
            deratio = 1.
            if ( ectmax .gt. 0. ) deratio = energy/ectmax
            energy = 9.+log10(energy)          
         endif 
      endif

c - - - - - - time analysis from file time.txt - - - - - - - hc3 - - - -
c - - - - - - or time.txt00nnnn or time00nnnn.txt  - - - - - - - - - - -
      if ( tminuts .lt. 1. ) then
      iest = 0
      chtime = chjobname(jj)(1:10)//'time.txt'
  145 continue
      inquire(file=chtime,exist=lexist)
      if ( lexist ) then
         open(unit=3,file=chtime,form='formatted',status='old')
         read(3,'(a)',err=148,end=148) chtext
         read(3,*,err=148,end=148) secst,secfi,tminuts
         read(3,'(a)',err=148,end=148) chtext
         read(3,'(a)',err=148,end=148) chtext
         jbl1 = index(chtext,'jobs =')
         jbl2 = jbl1+6
         ! count blanks before digits of mfiles:
  146    continue
         jbl2 = jbl2+1
         if ( chtext(jbl2:jbl2) .eq. ' ' ) goto 146       
         ! count valid digits of mfiles:
  147    continue
         jbl2 = jbl2+1            
         if ( chtext(jbl2:jbl2) .ne. ' ' ) goto 147
         write(cfmtint(3:3),'(i1)') jbl2-(jbl1+6)
         read(chtext(jbl1+6:jbl2-1),cfmtint) mfiles
  148    continue
         close(unit=3)
      else
         iest = iest + 1
         if ( iest .eq. 2 ) 
     +      chtime = chjobname(jj)(1:10)//'time.txt'//chjobname(jj)(4:9)
         if ( iest .eq. 3 ) chtime 
     +      = chjobname(jj)(1:10)//'time'//chjobname(jj)(4:9)//'.txt'
         if ( iest .le. 3 ) goto 145
      endif
      else
         ! tminuts already available from Job00*_scc.out (see above).
      endif

c - - - - - - - print tabular - - - - - - - - - - - - - - - - - - - - -
      if ( jj .eq. 1 ) then  
      ! - - - - - - print first title line - - - - - - - - - - - - - - -
      if ( lthi .ge. 1 ) then
         write(*,'(85x,''=E/ectmax'')')
         write(*,'(/,''primary  lg(E)  theta   phi  runtsk  sizeGBy'',
     +      2x,''procs   ectcut   ectmax  t(min)  files'',
     +      3x,''RATIO'',2x,''obslev'',3x,''Xmagn'',3x,''Zmagn'',
     +      2x,''ecutha'',2x,''ecutmu'',2x,''ecutel'',2x,''ecutga'', 
     +      2x,''thilev  wmax  lg(rad)'')')
      else 
         write(*,'(85x,''=E/ectmax'')')
         write(*,'(''primary  lg(E)  theta   phi  runtsk  sizeGBy'',
     +      2x,''procs   ectcut   ectmax  t(min)  files'',
     +      3x,''RATIO'',2x,''obslev'',3x,''Xmagn'',3x,''Zmagn'',
     +      2x,''ecutha'',2x,''ecutmu'',2x,''ecutel'',2x,''ecutga'')') 
      endif
      endif
      ! - - - - - - mcode=0 for stackin simulations - - - - - - - - - -
      if ( mcode .eq. 0 ) mcode = 4
      ! - - - - - - distinguish `thin` and `stnd` - - - - - - - - - - -
      if ( sumgiga .lt. 1.e-2 ) sumgiga = 1.e-2
      if ( mcode .ne. mcprev .or. energy .ne. enprev .or.
     +     theta .ne. thprev .or. phia   .ne. phprev .or.
     +      lthi .ne. lthipr ) then
           if ( mcprev .eq. 5626 .and. mcode .ne. 5626 ) write(*,'(1x)')
           if ( mcprev .ne. 703 .and. mcprev .ne. 5626 ) write(*,'(1x)')
           if ( mcprev .eq. 703 .and. mcode .ne. 703 ) write(*,'(1x)')
      endif

      if ( ifiles .lt. mfiles ) then ! almost all DAT files deleted:
        if ( lthi .ge. 1 ) then
          write(*,'(i6,''*'',f7.2,f6.1,f7.1,i8.6,f9.2,i6,1p,e10.1,e9.1,
     +      0p,f8.1,i7,2f8.1,2f8.2,4(ss,f5.1,''e'',sp,i2), 
     +      ss,f6.0,1p,e9.1,0p,f6.1)')
     +      mcode,energy,theta,phia,mrun,sumgiga,mprocs,ectcut,
     +      ectmax,tminuts,mfiles,deratio,obslev,bxmag,bzmag,
     +      (ecmnt(j),jcexp(j),j=1,4),thinlev,thiwmax,thirad
        else
          write(*,'(i6,''*'',f7.2,f6.1,f7.1,i8.6,f9.2,i6,1p,e10.1,e9.1,
     +      0p,f8.1,i7,2f8.1,2f8.2,4(ss,f5.1,''e'',sp,i2))')
     +      mcode,energy,theta,phia,mrun,sumgiga,mprocs,ectcut,
     +      ectmax,tminuts,mfiles,deratio,obslev,bxmag,bzmag,
     +      (ecmnt(j),jcexp(j),j=1,4)
        endif
      else ! all particle data files available:
        if ( lthi .ge. 1 ) then
          write(*,'(i6,f8.2,f6.1,f7.1,i8.6,f9.2,i6,1p,e10.1,e9.1,
     +      0p,f8.1,i7,2f8.1,2f8.2,4(ss,f5.1,''e'',sp,i2), 
     +      ss,f6.0,1p,e9.1,0p,f6.1)')
     +      mcode,energy,theta,phia,mrun,sumgiga,mprocs,ectcut,
     +      ectmax,tminuts,mfiles,deratio,obslev,bxmag,bzmag,
     +      (ecmnt(j),jcexp(j),j=1,4),thinlev,thiwmax,thirad
        else
          write(*,'(i6,f8.2,f6.1,f7.1,i8.6,f9.2,i6,1p,e10.1,e9.1,
     +      0p,f8.1,i7,2f8.1,2f8.2,4(ss,f5.1,''e'',sp,i2))')
     +      mcode,energy,theta,phia,mrun,sumgiga,mprocs,ectcut,
     +      ectmax,tminuts,mfiles,deratio,obslev,bxmag,bzmag,
     +      (ecmnt(j),jcexp(j),j=1,4)
        endif
      endif
      if ( 899 .le. mrun .and. mrun .le. 899 )
     +   write(*,'(51x,2(i8,''*''))') int(0.5+ectmax/ectcut)

c - - - - - - - copy new quantites to previous - - - - - - - - - - - - -
      mcprev = mcode 
      enprev = energy
      thprev = theta
      phprev = phia
      lthipr = lthi

  149 continue 

      enddo ! end-of-loop jj=1,job

c - - - - - - print last title line of tabular - - - - - - - - - - - - -
      if ( mthi .gt. 0 ) then
         write(*,'(/,''primary  lg(E)  theta   phi  runtsk  sizeGBy'',
     +      2x,''procs   ectcut   ectmax  t(min)  files'',
     +      3x,''RATIO'',2x,''obslev'',3x,''Xmagn'',3x,''Zmagn'',
     +      2x,''ecutha'',2x,''ecutmu'',2x,''ecutel'',2x,''ecutga'', 
     +      2x,''thilev  wmax  lg(rad)'',/,85x,''=E/ectmax'')')
      else
         write(*,'(/,''primary  lg(E)  theta   phi  runtsk  sizeGBy'',
     +      2x,''procs   ectcut   ectmax  t(min)  files'',
     +      3x,''RATIO'',2x,''obslev'',3x,''Xmagn'',3x,''Zmagn'',
     +      2x,''ecutha'',2x,''ecutmu'',2x,''ecutel'',2x,''ecutga'', 
     +      /,85x,''=E/ectmax'')')
      endif   

c - - - - - - end-of-program showparallel - - - - - - - - - - - - - -
      stop
      end
