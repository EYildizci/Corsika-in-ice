c=======================================================================
c
c  a c r e i n p h c 3 . f
c  -----------------------
c     Automatic creation of successive steering files to run corsika
c        simulations and of the corresponding shell script files.
c               (qgsjet, gheisha, single shower files)
c     ------------------------------------------------------------------
c     #!/bin/bash
c           f77 -fbounds-check -m32 acreinphc3.f -o acreinphc3
c           gfortran -fbounds-check acreinphc3.f -o acreinphc3
c           ifort -C -check bounds acreinphc3.f -o acreinphc3
c           ./acreinphc3
c           chmod +x jobhc3*
c     job_submit -p1 -cp -t1000 -m1000 showanalyhc3.sh001750
c           ls -1 DAT001750* | grep t -v > showanalyhc3.i001750
c           showanalyhc3 < showanalyhc3.i001750 > showanalyhc3.out001750
c           mv fort.9 showanalyhc3.fort001750
c     ------------------------------------------------------------------
c           mpi_corsika73723_stnd_qgsII4_gheisha_runner
c           mpi_corsika72495_thin_qgsII_urqmd_runner
c           mpi_corsika72499_stnd_qgsII3_urqmd_runner
c     ------------------------------------------------------------------
c     job_submit -cp -p16 -t60 -m1000 -oJob001517_%j.out
c                -eJob001517_%j.err jobhc3-001517 parallel-001517
c-----------------------------------------------------------------------
c                                   juergen.oehlschlaeger@kit.edu
c=======================================================================

      program acreinphc3
      implicit double precision (a-h,o-z)
      character cpinput*15,cplldir*10
      character cjobhc3*13,cquota*8,cquotb*8,crunmpi*26
      dimension mprim(8),mfanz(16),mshif(16)
      dimension qengy(0:16)

c - - - - - - set number of files, showers, values of energy - - - - - -
      data mprim/ 14, 402, 1206, 2814, 5626, 1, 1, 1/
      data mfanz/ 200, 160, 130, 100,  80,  20,  20,  20, 8*20/
      data mshif/   0, 200, 400, 600, 700, 800, 900, 950, 8*990/
      data qengy/1.0e7,1.78e7,3.16e7,5.62e7,1.00e8,1.78e8,3.16e8,5.62e8,
     +      1.00e9,1.78e9,3.16e9,5.62e9,1.00e10,1.78e10,3.16e10,5.62e10,
     +      1.00e11/
      data cquota/'"______ '/, cquotb/' ______"'/ 

c - - - - - - some initialisations for corsika files - - - - - - - - - -
      ecut1 = 0.1 
      ecut2 = 0.1
      ecut3 = 2.5e-4
      ecut4 = 2.5e-4 
      iauge = 1
      ifluk = 0
      eslop = -2.0
      themin = 8.0
      themax = 8.0
      phideg = 153.435
      dectcut = 500. ! 1.0e4
      dectmax = 5.e5 ! 1.0e7

c - - - - - - particle energy loop - - - - - - - - - - - - - - - - - - -
      do  29  iegy=0,0
      engya = 1.2345e+07 ! qengy(iegy)
      engyb = engya 
      write(*,'(''- - - - - - - - - - - - - - - - - - -'')')
      write(*,'(''   energ='',1p,e9.2,'' ...'',e9.2)')
     +    engya+6.,engyb+6.
      write(*,'(''- - - - - - - - - - - - - - - - - - -'')')
      ! isnmax = mfanz(iegy)
      isnmax = 1
      mstrt0 = 3456 ! + 5 * (iegy-1) + mshif(iegy)

c - - - - - - theta angle loop - - - - - - - -
      do  27  ithe=0,3
      themin = 8. + 10.*ithe
      themax = themin

c - - - - - - primary particle loop - - - - - - - -
      do  26  iprm=4,4,4
      mstrt = mstrt0  + ithe ! + 10*(iprm-1)
      iprim = mprim(iprm)
      write(*,'(10x,''prim='',i6)') mprim(iprm)

c - - - - - - (single) shower loop - - - - - - - -
      do  25  isnr=1,isnmax
         msrun = mstrt + isnr - 1
         if ( isnr .eq. 1 .or. isnr .eq. isnmax)
     +      write(*,*) '                     msrun=',msrun
         if ( isnr .eq. 2 .and. isnr .ne. isnmax)
     +      write(*,*) '                      . . . . . . . . .'
         irun = msrun
         im = mod(irun,1000000)

c - - - - - - some more quantities for shell scripts - - - - - - - -
         cpinput = 'parallel-000000'
         write(cpinput(10:15),'(I6.6)') im
         cplldir = 'csk000000/' 
         write(cplldir(4:9),'(I6.6)') im
         cjobhc3 = 'jobhc3-000000'
         write(cjobhc3(8:13),'(I6.6)') im
         crunmpi = 'runmpi000000.txt'
         write(crunmpi(7:12),'(I6.6)') im
         crunmpi = cplldir//crunmpi(1:16)

c - - - - - - - - create corsika steering file - - - - - - - - - - - -
         iseed = 1750*3 ! irun*3
         open(unit=9,FILE=cpinput,STATUS='UNKNOWN')
         write(9,'(''RUNNR '',I10)') irun
         write(9,'(''PARALLEL'',f9.0,f12.0,''  1  F'')') dectcut,dectmax
         write(9,'(''NSHOW '',I10)') 1
         write(9,'(''EVTNR '',I10)') 1
         write(9,'(''SEED  '',3I10)') iseed,0,0
         write(9,'(''SEED  '',3I10)') iseed+1,0,0
         write(9,'(''SEED  '',3I10)') iseed+2,0,0
         write(9,'(''SEED  '',3I10)') iseed+3,0,0
         write(9,'(''SEED  '',3I10)') iseed+4,0,0
         write(9,'(''SEED  '',3I10)') iseed+5,0,0
         write(9,'(''PRMPAR'',I10)') iprim
         if ( engya .lt. engyb ) write(9,'(''ESLOPE'',F10.2)') eslop
         write(9,'(''ERANGE'',6x,1P,2E12.4)') engya,engyb
         write(9,'(''THETAP'',2F12.2)') themin,themax
         write(9,'(''PHIP  '',2F12.2)') phideg,phideg ! 0.,360.
         if ( iauge .eq. 1 ) then
            write(9,'(''OBSLEV   1452.e2        870.000 g/cm^2'')')
            write(9,'(''MAGNET     19.71        -14.18     Auger'')')
         else
            write(9,'(''OBSLEV    110.e2       1022.647 g/cm^2'')')
            write(9,'(''MAGNET     20.40         43.20   Karlsruhe'')')
         endif
         write(9,'(''MAXPRT'',I10)') 1
         write(9,'(''ECTMAP     1.E11'')')
         write(9,'(''ECUTS    '',2F9.4,2f10.5)') ecut1,ecut2,ecut3,ecut4
         write(9,'(''RADNKG    200.E2'')')
         write(9,'(''HADFLG'',6I5)') (0,i=1,5),2
         write(9,'(''ELMFLG         T         T'')')
         ! write(9,'(''SIBYLL         T'',I10)') 0
         ! write(9,'(''SIBSIG         T'')')
         ! write(9,'(''QGSJET         T'',I10)') 0
         ! write(9,'(''QGSSIG         T'')')
         write(9,'(''MUADDI         T'')')
         write(9,'(''MUMULT         T'')')
         write(9,'(''STEPFC        1.'')')
         write(9,'(''LONGI          T     5.0     T     T'')')
         write(9,'(''HILOW       111.11'')')
         write(9,'(''DIRECT '',a)') cplldir
         write(9,'(''HOST   hc3.uni'')')
         write(9,'(''USER   joe'')')
         write(9,'(''EXIT  '')')
         close(unit=9)

c - - - - - - - - create shell script - - - - - - - - - - - - - - - - -
         itimec = 60
         if ( engya .gt. 3.3333e10 ) itimec = 68 * 60 ! max 4080 min.
         open(unit=9,FILE=cjobhc3,STATUS='UNKNOWN')
         write(9,'(''#!/bin/bash'')')
         write(9,'(''# script to run a parallel corsika simulation'','//
     +      '1x,''on `hc3.scc.kit.edu`'')')
         write(9,'(''if [ ! -e '',a,'' ] ; then'')') cplldir
         write(9,'(''   /bin/mkdir '',a)') cplldir
         write(9,'(''else'')')
         write(9,'(''   /bin/rm -f '',a,''*'')') cplldir
         write(9,'(''fi'')')
         write(9,'(''/bin/cat '',a)') cjobhc3 
         write(9,'(''/bin/cp '',a,1x,a)') cjobhc3,cplldir 
         write(9,'(''/bin/cp '',a,1x,a)') cpinput,cplldir
         write(9,'(''/bin/cp sum* '',a)') cplldir
         write(9,'(''/bin/cp totaltime* '',a)') cplldir
         write(9,'(''/bin/cp postprocess.sh '',a)') cplldir
         write(9,'(''/bin/cp readcsk* '',a)') cplldir
         write(9,'(''/bin/cp showanalyhc3* '',a)') cplldir
         write(9,'(''# '')')
         write(9,'(''job_submit -cd -p16 -t'',i2,'' -m2000 -oJob'',
     +      i6.6,''_%j.out -eJob'',i6.6,''_%j.err mpirun '',
     +      ''mpi_corsika73723_stnd_QGSII4_gheisha_runner '',a)')
     +      mod(itimec,100),im,im,cpinput
         write(9,'(''# '')')
         close(unit=9)
   25 continue
   26 continue
   27 continue
   28 continue
   29 continue
c - - - - - - end of all loops - - - - - - - - - - - - - - - - - - - - -

      stop
      end
