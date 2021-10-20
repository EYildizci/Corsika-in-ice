c=======================================================================
c
c  s h o w a n a l y h c 3 . f
c  ---------------------------
c     create first a text file containing printed histograms (see
c     following tabular) and isecond the bin contents of all histograms.
c compilation:
c     gfortran -fbounds-check showanalyhc3.f -o showanalyhc3
c     f77 -fbounds-check -m32 showanalyhc3.f -o showanalyhc3
c     ifort -C -check bounds showanalyhc3.f -o showanalyhc3
c execution:
c     # job_submit -p1 -cp -t2000 -m1000 showanalyhc3.sh000053
c     ls -1 DAT00* | grep t -v | n -v > showanalyhc3.i000053
c     ./showanalyhc3 < showanalyhc3.i000053 > showanalyhc3.out000053
c     mv fort.9 showanalyhc3.fort000053
c get number of histograms:
c     grep "  -99.8877" showanalyhc3.fort000053 | wc 
c-----------------------------------------------------------------------
c input example:
c           csk000053/DAT000053-882208370-000000564
c           csk000053/DAT000053-882426264-000000803
c           csk000053/DAT000053-882850130-000000136
c-----------------------------------------------------------------------
c           runh=211285.281   evth=217433.078
c           long=52815.2969   evte=3397.39185   rune=3301.33252
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
c
c       121 (1)   rec  density of gam   vs log10(r)  "always weight 1" 
c       122 (1)   rec  density of e+    vs log10(r)
c       123 (1)   rec  density of e-    vs log10(r)
c       124 (1)   rec  density of mu+   vs log10(r)
c       125 (1)   rec  density of mu-   vs log10(r)
c       126 (1)   rec  density of pi+   vs log10(r)
c       127 (1)   rec  density of pi-   vs log10(r)
c       128 (1)   rec  density of p     vs log10(r)
c       129 (1)   rec  density of n     vs log10(r)
c       130 (1)   rec  density of cerph vs log10(r)
c       131 (1)   rec  density of other vs log10(r)
c
c       181 (1)   N rec gam   vs log10(E)   "always weight 1"
c       182 (1)   N rec e+    vs log10(E)
c       183 (1)   N rec e-    vs log10(E)
c       184 (1)   N rec mu+   vs log10(E)
c       185 (1)   N rec mu-   vs log10(E)
c       186 (1)   N rec pi+   vs log10(E)
c       187 (1)   N rec pi-   vs log10(E)
c       188 (1)   N rec p     vs log10(E)
c       189 (1)   N rec n     vs log10(E)
c       190 (1)   N rec cerph vs log10(E)
c       191 (1)   N rec other vs log10(E)
c
c       221 (1)   N rec gam   vs log10(t)   "always weight 1"
c       222 (1)   N rec e+    vs log10(t)
c       223 (1)   N rec e-    vs log10(t)
c       224 (1)   N rec mu+   vs log10(t)
c       225 (1)   N rec mu-   vs log10(t)
c       226 (1)   N rec pi+   vs log10(t)
c       227 (1)   N rec pi-   vs log10(t)
c       228 (1)   N rec p     vs log10(t)
c       229 (1)   N rec n     vs log10(t)
c       230 (1)   N rec cerph vs log10(t)
c       231 (1)   N rec other vs log10(t)
c
c       299 (2)   particle codes vs log10(r)
c
c       341 (2)   N gam   vs log10(E) and log10(r)
c       342 (2)   N e+    vs log10(E) and log10(r)
c       343 (2)   N e-    vs log10(E) and log10(r)
c       344 (2)   N mu+   vs log10(E) and log10(r)
c       345 (2)   N mu-   vs log10(E) and log10(r)
c       346 (2)   N pi+   vs log10(E) and log10(r)
c       347 (2)   N pi-   vs log10(E) and log10(r)
c       348 (2)   N p     vs log10(E) and log10(r)
c       349 (2)   N n     vs log10(E) and log10(r)
c       350 (2)   N cerph vs log10(E) and log10(r)
c       351 (2)   N other vs log10(E) and log10(r)
c
c       361 (2)   N gam   vs log10(t) and log10(r)
c       362 (2)   N e+    vs log10(t) and log10(r)
c       363 (2)   N e-    vs log10(t) and log10(r)
c       364 (2)   N mu+   vs log10(t) and log10(r)
c       365 (2)   N mu-   vs log10(t) and log10(r)
c       366 (2)   N pi+   vs log10(t) and log10(r)
c       367 (2)   N pi-   vs log10(t) and log10(r)
c       368 (2)   N p     vs log10(t) and log10(r)
c       369 (2)   N n     vs log10(t) and log10(r)
c       370 (2)   N cerph vs log10(t) and log10(r)
c       371 (2)   N other vs log10(t) and log10(r)
c
c       381 (2)   N gam   vs log10(E) and log10(t)
c       382 (2)   N e+    vs log10(E) and log10(t)
c       383 (2)   N e-    vs log10(E) and log10(t)
c       384 (2)   N mu+   vs log10(E) and log10(t)
c       385 (2)   N mu-   vs log10(E) and log10(t)
c       386 (2)   N pi+   vs log10(E) and log10(t)
c       387 (2)   N pi-   vs log10(E) and log10(t)
c       388 (2)   N p     vs log10(E) and log10(t)
c       389 (2)   N n     vs log10(E) and log10(t)
c       390 (2)   N cerph vs log10(E) and log10(t)
c       391 (2)   N other vs log10(E) and log10(t)
c
c-----------------------------------------------------------------------

      program showanalyhc3

c-----------------------------------------------------------------------
c        hist vector header:
c     qhistos(1,ih) = ident number of the histogram
c     qhistos(2,ih) = number of bins (1st dim), nbin
c     qhistos(3,ih) = lower bound of 1st dim
c     qhistos(4,ih) = upper bound of 1st dim
c     qhistos(5,ih) = marker for all logarithmic axes
c     qhistos(6,ih) = number of bins (2nd dim), nbi2
c     qhistos(7,ih) = lower bound of 2nd dim
c     qhistos(8,ih) = upper bound of 2nd dim
c     qhistos(9,ih) = -99.8877
c     qhistos(10,ih) = not used and not written to unit 9
c     qhistos(11,ih) = first element of histogram contents
c        hist vector trailer 1-dim:                not written to unit 9
c     qhistos(10+nbin,ih) = last element of histogram contents
c     qhistos(10+nbin+3,ih) = min index of 1st dim with entries
c     qhistos(10+nbin+4,ih) = max index of 1st dim with entries
c     qhistos(10+nbin+5,ih) = minimum of number of entries
c     qhistos(10+nbin+6,ih) = maximum of number of entries
c        hist vector trailer 2-dim:                not written to unit 9
c     qhistos(10+nbin*nbi2,ih) = last element of histogram contents
c     qhistos(10+nbin*nbi2+3,ih) = min index of 1st dim with entries
c     qhistos(10+nbin*nbi2+4,ih) = max index of 1st dim with entries
c     qhistos(10+nbin*nbi2+5,ih) = minimum of number of entries
c     qhistos(10+nbin*nbi2+6,ih) = maximum of number of entries
c     qhistos(10+nbin*nbi2+7,ih) = min index of 2nd dim with entries    
c     qhistos(10+nbin*nbi2+8,ih) = max index of 2nd dim with entries
c-----------------------------------------------------------------------

      implicit double precision (a-h,o-z), integer (i-n)
      parameter (nblklen=312,nhigh=20)
      character chtitle(480)*40,chaziff(0:9)*1
      double precision qhistos(10020,480)
      common /histch/ chtitle
      common /histos/ qhistos,rpromin,rpromax
     +       ,eparmin,eparmax
     +       ,xmin348,xmax348,ymin348,ymax348
      common /buffer/ outbuf(nblklen,21)
      common /buffch/ crunh,crune,cevth,cevte,clong
      character chline(0:101),chaxis(0:101),chlin2(100),cheight(50)
      character*4 crunh,crune,cevth,cevte,clong
      real    ddata(6552),outbuf,chablk(5),curpar
      equivalence(chablk(1),crunh)
      data crunh/'RUNH'/, crune/'RUNE'/, clong/'LONG'/
      data cevth/'EVTH'/, cevte/'EVTE'/
      data chtitle/480*'  '/
      data chaziff/' ','1','2','3','4','5','6','7','8','9'/
      data cheight/'a','b','c','d','e','f','g','h','i','j',
     +             'k','l','m','n','o','p','q','r','s','t',
     +             'u','v','w','x','y',
     +             'A','B','C','D','E','F','G','H','I','J',
     +             'K','L','M','N','O','P','Q','R','S','T',
     +             'U','V','W','X','Y'/
c-----------------------------------------------------------------------
c  physical constants and run parameters used for simulation
      common / phys / obslev(10),enspec(3),ffegs,ffnkg,cut(4),
     *                consts(230),aatm(5),batm(5),catm(5),flag(4)
c-----------------------------------------------------------------------
c  additional run parameters in analysis job
c    ishow/ishowu: no of showers read/used
c    ishowf: no of showers read from current file
c    iblk  : number of blocks read
c    ifile : number of current data file
c    iret1/2: flags for error or runstop in subroutines particles/header
c    xoff, yoff  coordinate correction for inclined showers
      common /runpar/ heighp,thetap,phip,xoff,yoff,tcen,dintmod,dintcrs,
     *                nobslv,ishow,ishowu,ishowf,iblk,ifile,iret1,iret2,
     *                rmax,lfnkg,lfegs
      double precision heighp,thetap,phip,tcen,dintmod,dintcrs,rmax
      logical*1 lfnkg,lfegs
      double precision pama,pi
      common /pconst/pama(6000),pi
      data pi/3.1415926535979d0/
c-----------------------------------------------------------------------
c  curpar = current particle :
c  curpar (1)      =  ident * 1000  +  igen*10  +  iobslv
c    ident  = identification according to geant
c    igen   = generation
c    iobslv = observation level
c  curpar (2,3,4)  =  px,py,pz  in GeV/c
c  curpar (5,6,7)  =  x,y,t  in cm, ns
c  modified for cerenkov bunches
c     curpar(1)    =  99*10**5 + nphot*10 + 1
c     curpar(2,3)  =  x and y coordinate in cm
c     curpar(4,5)  =  direction cosini to x and y axis
c     curpar(6)    =  t in ns
c     curpar(7)    =  height of emission in cm
c     curpar(8)    =  weight of particle
c  epart,ppart,r   =  energy, momentum, distance to shower axis
c  ityp            =  particle group
c  tf              =  arrival time of shower front at point(x,y)
c  costhe, phi     =  cosine of polar angle, azimuth angle
      common /partic/ curpar(8),ident,igen,iobslv,ityp,
     *                epart,ppart,r,tf,costhe,phi
c-----------------------------------------------------------------------
c  shower parameters
c  factors needed to project particles into the plane perp. to shower axis
      common /shower/ profak1,profak2,profak3,
     *   ccla
      double precision profak1,profak2,profak3
      character*5 ccla(11)
      data ccla / 'gam  ','e+   ','e-   ','mu+  ','mu-  ',
     *            'pi+  ','pi-  ','p    ','n    ','cerph','other' /
      parameter (nmax=50000)
      common /steer/  lunmon,lundeb,lunin,lunhis,nfiles
      common /steerc/ chfile(nmax),chhist,chlong
      character*80    cdat,chfile,chhist,chlong
      logical lexist,first/.true./
      cheight(nhigh) = '@'
      iblklen = 273
      eparmin = 1.
      eparmax = 1.
      rpromin = 1.
      rpromax = 1.
      xmin348 = 1.
      xmax348 = 1.
      ymin348 = 1.
      ymax348 = 1.

c-----------------------------------------------------------------------

c  set defaults
      lunmon = 6
      lundeb = 6
      lunin  = 3
      lunhis = 6
      nobslv = 1
      nfiles = 0
      do ihist=1,480
      do i=1,10020
         qhistos(i,ihist) = 0.d0
      enddo
      enddo
      do i=1,nmax
        chfile(i) = ' '
      enddo
c  initialize particle masses
      call pamaf

c--read run parameters including file names-----------------------------
      do  ifile=1,nmax
         read(*,'(a)',end=427,err=427) cdat
         chfile(ifile) = cdat
      enddo
  427 continue
      nfiles = ifile - 1
      nfimod = 10
      if ( nfiles .gt. 10000 ) then
         nfimod = 200
      elseif ( nfiles .gt. 1000 ) then
         nfimod = 50
      elseif ( nfiles .gt. 300 ) then
         nfimod = 20
      endif

c--read first particle data file to test record length------------------
      open(lunin,file=chfile(1),status='old',
     *        form='unformatted',access='sequential',err=997)
      read(lunin) (ddata(i),i=1,5432)
      if ( 217433.0 .lt. ddata(273+1) .and.
     +                   ddata(273+1) .lt. 217433.2 ) then
         write(*,100)
  100 format(/,' ===========================================',/,
     *         ' standard  CORSIKA  particle  data  analysis',/,
     *         ' ===========================================')
         iblklen = 273
      elseif ( 217433.0 .lt. ddata(312+1) .and.
     +                       ddata(312+1) .lt. 217433.2 ) then
         write(*,101)
  101 format(/,' ===========================================',/,
     *         ' thinning  CORSIKA  particle  data  analysis',/,
     *         ' ===========================================')
         iblklen = 312
      endif
      close(lunin)

c--read particle data files---------------------------------------------
      do  ifile=1,nfiles
         if ( ifile .lt. 10 )
     +      write(*,'(i6,''. file: '',a)') ifile,chfile(ifile)
         if ( ifile .eq. 11 )
     +      write(*,'(''     ....................'')')
         inquire(file=chfile(ifile),exist=lexist)
         if ( .not.lexist ) then
            write(*,'(10x,''nr.'',i5,4x,a44,'' NOT found.'')')
     +         ifile,chfile(ifile)
            chfile(ifile) = 'dummy'   
            goto 431
         endif 
         open(lunin,file=chfile(ifile),status='old',
     *        form='unformatted',access='sequential',err=997)
         read(lunin,end=429,err=430) (ddata(i),i=1,iblklen*21)
         if ( ddata(1) .eq. 0. .and. ddata(iblklen+1) .eq. 0. ) then
            write(*,'(10x,''nr.'',i5,4x,a44,'' only ZEROES.'')')
     +         ifile,chfile(ifile)
            chfile(ifile) = 'dummy'
         endif 
         goto 431
  429    continue
         write(*,'(10x,''nr.'',i5,4x,a44,'' EOF detected.'')')
     +      ifile,chfile(ifile)
         chfile(ifile) = 'dummy'   
         goto 431 
  430    continue
         write(*,'(10x,''nr.'',i5,4x,a44,'' ERR detected.'')')
     +      ifile,chfile(ifile)
         chfile(ifile) = 'dummy'   
  431    continue
         close(lunin)
      enddo
      if ( nfiles .eq. 10 .or. ifile .gt. 10 )
     +   write(*,'(i6,''. file: '',a)') ifile-1,chfile(ifile-1)

c-----------------------------------------------------------------------
c  initialize some variables, histograms
      ifile  = 0
      ishow  = 0
      ishowu = 0
      iblk   = 0
      call hisini

c-----------------------------------------------------------------------
c  begin of loop over all files

   10 continue
      ifile = ifile + 1
      ishowf = 0

c  check on correctly read file and quantities
      if ( chfile(ifile)(1:5) .eq. 'dummy' ) then
         write(*,'(10x,''nr.'',i5,'' renamed to `dummy`.'')') ifile
         goto 10
      endif

c  open current file
      open(lunin, file=chfile(ifile), status='old',
     *     form='unformatted', access='sequential')
      if ( mod(ifile,nfimod) .eq. 0 )
     +   write(*,'(10x,''nr.'',i5,'' in progress: '',a)')
     +      ifile,chfile(ifile)

c  read run header (i.e. first data record)
      iret2 = 0
      read(lunin,end=11,err=22) (ddata(i),i=1,iblklen*21)
      iblk = iblk + 1
      do  isb=1,21
      do   ib=1,iblklen
         outbuf(ib,isb) = ddata(ib+iblklen*(isb-1))
      enddo
      enddo

c  not data tape start marker read
      if ( .not. ( 211285.2 .lt. outbuf(1,1) .and.
     +                           outbuf(1,1) .lt. 211285.4 ) ) then
        write(*,*) ' wrong header line of new run'
        iret2 = 1
        goto 88
      endif

c  define and print primary particle
      if ( ifile .eq. 1 ) then
         write(*,219) int(outbuf(2,1)),outbuf(3,2),outbuf(4,2),
     +      outbuf(7,2)/100.,outbuf(48,2)/100.,2.e7+outbuf(3,1),
     +      outbuf(4,1)
  219 format(/,' simulation run number   ',i15.6,'.',/,
     *         ' primary particle        ',f16.0,/,
     *         ' energy fixed [GeV]',1p,e22.4,0p,/,
     *         ' height of 1st int. [m]  ',f16.0,/,
     *         ' observation level [m]   ',f16.0,/,
     *         ' date of simulation      ',f16.0,/,
     *         ' version                 ',f16.6)
         write(*,224) (outbuf(i,1),i=21,24)
  224 format(' cuts [GeV] hadron, muon, electr, photon: ',4f9.5)
         obslev(1) = outbuf(48,2)
         lfnkg = ( ffnkg .eq. 1. )
         lfegs = ( ffegs .eq. 1. )
         nflran = flag(1)
         nfldif = flag(2)
         nflpi0 = mod( int(flag(3)),100 )
         nflpif = int(flag(3)) / 100
         nflche = mod( int(flag(4)),100 )
         nfragm = int(flag(4)) / 100
      endif
      goto 88

   11 continue
      write(*,*) ' end of file at reading run header '
      iret2 = 1
      goto 88
   22 continue
      write(*,*) ' error at reading run header '
      iret2 = 1

c = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   88 continue
      if ( iret2 .ne. 0 ) then
        write(*,*) ' error or eof in reading header '
        goto 99
      endif

c-----------------------------------------------------------------------
c  read data blocks from current file and call analysis routines
      iret1 = 0
      call particles(iblklen)
      if ( iret2 .ne. 0 ) then
        write(*,*) ' error or eof in header '
        goto 99
      endif

c  close current file
      close(lunin)
      if ( iret1 .ne. 0 ) then
          
         if ( iret1 .eq. 4 ) write(*,*) ' unexpected end-of-file at ',
     +      ifile,'th file: ',chfile(ifile) 

c     write(*,*) 'eof detected when reading lunin '
c     iret1 = 5 
c     return
c     write(*,*) ' error at reading lunin '
c     iret1 = 3
c     return

        write(*,1001) ifile,nfiles
 1001   format(' error or end in particles at ',
     *         i3,'. file of',i5,' input files ')
      endif

c  does another input file exist?
      if ( nfiles .gt. ifile ) goto 10

c-----------------------------------------------------------------------
c  end of event loop
c-----------------------------------------------------------------------

   99 continue
      write(*,'(/,5x,''total number of records'',i10,/)') iblk

c-----------------------------------------------------------------------
c  printer plots of histograms at the end of the run
c     histext
c     nhist, nbin, xlow, xupr, vlog, [nbi2, ylow, yupr, dummy]
c     plots / table of contents
c  histogram quantities will be written to fort.9 (without statistics)
c-----------------------------------------------------------------------

      do  ihist=121,389
      if ( chtitle(ihist)(1:2) .ne. '  ' ) then

c - - - - - - print title and optionally contents of the histogram:
          nbin = qhistos(2,ihist)
          nbi2 = max(1,int(qhistos(6,ihist)))
          write(*,'(55('' =''))')
          write(*,'(1x,a40)') chtitle(ihist)
          write(9,'(1x,a40)') chtitle(ihist)
          if ( qhistos(6,ihist) .eq. 0. ) then ! 1-dim.
            write(*,'(2(0p,2f8.0,1p,2g12.4))') 
     +         (qhistos(i,ihist),i=1,5)
            write(9,'(2(0p,2f8.0,1p,2g12.4),0p,f14.4)')
     +         (qhistos(i,ihist),i=1,5), 0.,-99.,-99., -99.8877
            write(9,'(1p,10e12.4)') (qhistos(i,ihist),i=11,10+nbin)
          else ! 2-dim.
            write(*,'(2(0p,2f8.0,1p,2g12.4))')
     +         (qhistos(i,ihist),i=1,8)
            write(9,'(2(0p,2f8.0,1p,2g12.4),0p,f14.4)')
     +         (qhistos(i,ihist),i=1,8), -99.8877
            write(9,'(1p,10e12.4)') (qhistos(i,ihist),i=11,10+nbin*nbi2)
          endif
c - - - - - - calculate minimum and maximum of the histogram:        
          qhistos(10+nbin*nbi2+3,ihist) = -1.
          qhistos(10+nbin*nbi2+4,ihist) = -1.
          qhmin = qhistos(11,ihist)
          qhmax = qhistos(11,ihist)
          do  ib=11,10+nbin*nbi2
             if ( qhistos(ib,ihist) .gt. qhmax )
     +          qhmax = qhistos(ib,ihist)
             if ( qhistos(ib,ihist) .lt. qhmin )
     +          qhmin = qhistos(ib,ihist)
             if ( qhistos(10+nbin*nbi2+3,ihist) .eq. -1. .and.
     +            qhistos(ib,ihist) .gt. 0. ) then
                qhistos(10+nbin*nbi2+3,ihist) = -10 + ib
             endif
             if ( qhistos(10+nbin*nbi2+4,ihist) .eq. -1. .and.
     +            qhistos(21+nbin*nbi2-ib,ihist) .gt. 0. ) then
                qhistos(10+nbin*nbi2+4,ihist) = 11. + nbin*nbi2 - ib
             endif
          enddo
          qhistos(10+nbin*nbi2+5,ihist) = qhmin
          qhistos(10+nbin*nbi2+6,ihist) = qhmax

          if ( qhmin .eq. qhmax ) goto 109

c-----------------------------------------------------------------------
          if ( ihist .lt. 299 ) then

c - - - - - - - calculate stretching factor (5 or 2):
             qhlog = log10(qhmax)
             if ( qhlog .gt. 0. ) then
                ihlog = int(qhlog)
                qhlog = qhlog - ihlog
                qscal = 1. - ihlog
             else
                ihlog = 1. + int(-qhlog)
                qhlog = qhlog + ihlog
                qscal = 1.*ihlog
             endif
             if ( qhlog .gt. 0.69897d0 ) then
                qfact = 1.
             elseif ( qhlog .gt. 0.30103d0 ) then
                qfact = 2.
             else
                qfact = 5.
             endif
             qscal = 10.d0**qscal
             write(*,'(42x,''hmin='',1p,g12.5,''   hmax='',g12.5,
     +          9x,''(content *'',g9.2,'')'')') qhmin,qhmax,1.d0/qscal
c - - - - - - - plot title axis and scale quantities:
             do  i=1,100
                chaxis(i) = '-'
             enddo
             if ( qfact .eq. 5. ) then
                do  i=10,90,10
                   chaxis(i) = 'o'
                enddo
                chaxis(25) = '+'
                chaxis(50) = '|'
                chaxis(75) = '+'
                chaxis(100) = '|'
             elseif ( qfact .eq. 2. ) then
                do  i=10,90,20
                   chaxis(i) = '+'
                   chaxis(i+10) = '|'
                enddo
             else
                do  i=10,100,10
                   chaxis(i) = '|'
                enddo
             endif
             if ( qfact .eq. 5. ) then
                write(*,'(i8,4i25)') 0,(i,i=5,20,5)
             elseif ( qfact .eq. 2. ) then
                write(*,'(i8,5i20)') 0,(10*i,i=1,5)
             else
                write(*,'(i8,10i10)') 0,(10*i,i=1,10)
             endif
             write(*,'(6x,'' !'',100a1)') (chaxis(i),i=1,100)
c - - - - - - - plot histogram contents (to the right):
            do  ib=11,10+nbin
               do  is=1,100
                  chline(is) = ' '
               enddo
               if ( qfact .eq. 5. ) then
                  do  is=25,100,25
                     chline(is) = '|'
                  enddo
               elseif ( qfact .eq. 2. ) then
                  do  is=20,100,20
                     chline(is) = '|'
                  enddo
               else
                  do  is=10,100,10
                     chline(is) = '|'
                  enddo
               endif
               is = qscal * qfact * qhistos(ib,ihist)
               if ( is .eq. 0 .and. qhistos(ib,ihist) .gt. 0. ) is = 1
               chline(is) = '*'
               write(*,'(f6.2,'' |'',100a1)')
     +            qhistos(3,ihist)+(qhistos(4,ihist)-qhistos(3,ihist))/
     +            qhistos(2,ihist)*(ib-11),(chline(i),i=1,is)
            enddo 
c - - - - - - - plot closing axis and scale quantities:
            write(*,'(6x,'' !'',100a1)') (chaxis(i),i=1,100)
            if ( qfact .eq. 5. ) then
               write(*,'(i8,4i25)') 0,(i,i=5,20,5)  
            elseif ( qfact .eq. 2. ) then
               write(*,'(i8,5i20)') 0,(10*i,i=1,5)  
            else
               write(*,'(i8,10i10)') 0,(10*i,i=1,10) 
            endif

c-----------------------------------------------------------------------
          elseif ( ihist .eq. 299 ) then

c - - - - - - - hist.299: plot title axis and scale quantities:
             write(*,'(42x,''hmin='',1p,g12.5,''   hmax='',g12.5)')
     +         qhmin,qhmax
             ! - - - - - linear scale for contents:
             qfact = 0.99999d0 * nhigh / qhmax
          !write(*,'(59x,''hmax/'',i2,''='',1p,g12.5)')nhigh,qhmax/nhigh
             ! - - - - - logarithmic scale for contents:
             qfact = (qhmax+1.d-2)**(1.d0/nhigh)
             write(*,'(62x,''hmax^(1/'',i2,'')='',1p,g12.5)')
     +         nhigh,qfact
             if ( qhmin .eq. qhmax ) goto 105
             ! - - - - - print only for qhmax > qhmin:
             do  i=1,100
                chaxis(i) = '-'
             enddo
             do  i=20,100,20
                chaxis(i) = '|'
                chaxis(i-10) = '+'
             enddo
             write(*,'(f9.1,10f10.1)') (qhistos(7,ihist)+
     +          (qhistos(8,ihist)-qhistos(7,ihist))/10.*i,i=0,10)
             write(*,'(6x,'' !'',100a1)') (chaxis(i),i=1,100)
c - - - - - - - plot 2-dim-histogram contents using cheight steps:
             do  ib=0,nbin-1
                do  iy=1,nbi2
                   chlin2(iy) = ' '
                   if ( ib .gt. 0 ) then 
                   if ( qhistos(10+ib+(iy-1)*nbin,ihist) .gt. 0. ) then
                      ! - - - - - - - - - - linear steps:
                    ! ih = 1. + qhistos(10+ib+(iy-1)*nbin,ihist) * qfact
                    ! if ( ih .eq. 0 ) ih = 1  
                    ! chlin2(iy) = cheight(ih)
                    ! if ( qhistos(10+ib+(iy-1)*nbin,ihist).le.9. ) then
                    !    ih = qhistos(10+ib+(iy-1)*nbin,ihist)
                    !    chlin2(iy) = chaziff(ih)
                    ! endif
                      ! - - - - - - - - - - logarithmic steps:
                      qsqua = (qhmax+1.d-2)
                      ih = 21
  104                 continue
                      qsqua = qsqua / qfact
                      ih = ih - 1
                      if ( qhistos(10+ib+(iy-1)*nbin,ihist) .le. qsqua )
     +                   goto 104                 
                      if ( ih .eq. 0 ) ih = 1
                      chlin2(iy) = cheight(ih)      
                   endif
                   endif
                enddo
                write(*,'(i6,'' |'',100a1,''|'')')
     +             ib,(chlin2(iy),iy=1,nbi2)
             enddo
c - - - - - - - plot closing axis and scale quantities:
             write(*,'(6x,'' !'',100a1)') (chaxis(i),i=1,100)
             write(*,'(f9.1,10f10.1)') (qhistos(7,ihist)+
     +          (qhistos(8,ihist)-qhistos(7,ihist))/10.*i,i=0,10)
  105        continue

c-----------------------------------------------------------------------
          else ! ihist > 299

             if ( ihist .eq. 348 ) then
                write(*,'(42x,''xmin='',1p,g12.5,''   xmax='',g12.5)')
     +             xmin348,xmax348
                write(*,'(42x,''ymin='',1p,g12.5,''   ymax='',g12.5)')
     +             ymin348,ymax348
             endif 

c - - - - - - - print histogram 348 (2) N p   vs log10(E) and log10(r)
c - - - - - - - - - - - - - - - 368 (2) N p   vs log10(t) and log10(r)
c - - - - - - - - - - - - - - - 388 (2) N p   vs log10(E) and log10(t)
             write(*,'(42x,''hmin='',1p,g12.5,''   hmax='',g12.5)')
     +          qhmin,qhmax
             write(*,'(59x,''hmax/'',i2,''='',1p,g12.5)')
     +          nhigh,qhmax/nhigh
             qfact = 0.99999d0 * nhigh / qhmax
             do  i=1,100
                chaxis(i) = '-'
             enddo
             do  i=20,100,20
                chaxis(i) = '|'
                chaxis(i-10) = '+'
             enddo
             write(*,'(f9.1,10f10.1)') (qhistos(7,ihist)+
     +          (qhistos(8,ihist)-qhistos(7,ihist))/10.*i,i=0,10)
             write(*,'(6x,'' !'',100a1)') (chaxis(i),i=1,100)
             do  ib=1,nbin
                do  iy=1,nbi2
                   chlin2(iy) = ' '
                   if ( qhistos(10+ib+(iy-1)*nbin,ihist) .gt. 0. ) then
                      ih = 1. + qhistos(10+ib+(iy-1)*nbin,ihist) * qfact
                      chlin2(iy) = cheight(ih)
                   endif
                enddo
                write(*,'(f6.1,'' |'',100a1,''|'')')
     +             qhistos(3,ihist)+(qhistos(4,ihist)-qhistos(3,ihist))/
     +             qhistos(2,ihist)*(ib-1),(chlin2(iy),iy=1,nbi2)
             enddo
             write(*,'(6x,'' !'',100a1)') (chaxis(i),i=1,100)
             write(*,'(f9.1,10f10.1)') (qhistos(7,ihist)+
     +          (qhistos(8,ihist)-qhistos(7,ihist))/10.*i,i=0,10)

          endif ! end-of ihist > 299

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - end-of individual histograms.

  109   continue
        endif 
      enddo

      goto 999

c-----------------------------------------------------------------------

  997 write(*,'(a)') ' problems with open data file ==> run is stopped'
      write(*,'(a,4x,a)') ' chfile =', chfile(ifile)

  999 continue
      stop
      end
c=======================================================================

      subroutine hisini

c-----------------------------------------------------------------------
c  book histograms
c-----------------------------------------------------------------------
c  definitions
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (nblklen=312)
      character chtitle(480)*40
      double precision qhistos(10020,480)
      common /histch/ chtitle
      common /histos/ qhistos,rpromin,rpromax
     +       ,eparmin,eparmax
     +       ,xmin348,xmax348,ymin348,ymax348
      common /buffer/ outbuf(nblklen,21)
      common /buffch/ crunh,crune,cevth,cevte,clong
      character*4 crunh,crune,cevth,cevte,clong
      real      outbuf,chablk(5),curpar
      equivalence(chablk(1),crunh)
      double precision pama,pi
      common /pconst/pama(6000),pi
      common /partic/ curpar(8),ident,igen,iobslv,ityp,
     *                epart,ppart,r,tf,costhe,phi
      common / phys / obslev(10),enspec(3),ffegs,ffnkg,cut(4),
     *                consts(230),aatm(5),batm(5),catm(5),flag(4)
      common /runpar/ heighp,thetap,phip,xoff,yoff,tcen,dintmod,dintcrs,
     *                nobslv,ishow,ishowu,ishowf,iblk,ifile,iret1,iret2,
     *                rmax,lfnkg,lfegs
      double precision heighp,thetap,phip,tcen,dintmod,dintcrs,rmax
      logical*1 lfnkg,lfegs
      common /shower/ profak1,profak2,profak3,
     *   ccla
      double precision profak1,profak2,profak3
      character*5 ccla(11)
      parameter (nmax=50000)
      common /steer/  lunmon,lundeb,lunin,lunhis,nfiles
      common /steerc/ chfile(nmax),chhist,chlong
      character*80    chfile,chhist,chlong,htit
      double precision tl, ul, tl2, ul2, ulog, area
      integer nb,nb2

c  1-dim histograms
c  histograms as function of log10(r)
      nb = 100
      tl = 0.
      ul = 5.
c  record density
      htit = 'rec  density of       vs log10(r)'
      nn = 120
      ulog = 1.
      do i=1,09
        htit(17:21) = ccla(i)
        call hbook1(nn+i,htit,nb,tl,ul,ulog)
      enddo

c  histograms as function of log10(E)
      nb = 100
      tl = -4.
      ul =  6.
c  energy distribution for recs
      htit = 'N rec       vs log10(E)'
      nn = 180
      ulog = 1.
      do i=1,09
        htit(7:11) = ccla(i)
        call hbook1(nn+i,htit,nb,tl,ul,ulog)
      enddo

c  histograms as function of log10(t)
      nb = 100
      tl = -5.
      ul =  5.
c  time distribution for recs
      htit = 'N rec       vs log10(t)'
      nn = 220
      ulog = 1.
      do i=1,09
        htit(7:11) = ccla(i)
        call hbook1(nn+i,htit,nb,tl,ul,ulog)
      enddo

c  2-dim histograms
c  particle numbers as function of particle code and r
      nb = 30
      tl =  0.
      ul = 30.
      nb2 = 100
      tl2 = -4.
      ul2 =  6.
      ulog = 2.
      htit = 'particle codes vs log10(r)'
      nn = 299
      call hbook2(nn,htit,nb,tl,ul,nb2,tl2,ul2,ulog)

c  particle numbers as function of energy and r
      nb = 80 ! 100
      tl = -4. ! -4.
      ul =  4. ! 6.
      nb2 = 100
      tl2 = -4.
      ul2 =  6.
      nn = 340
      ulog = 3.
      htit = 'N       vs log10(E) and log10(r)'
      do i=1,09
        htit(3:7) = ccla(i)
        call hbook2(nn+i,htit,nb,tl,ul,nb2,tl2,ul2,ulog)
      enddo

c  particle numbers as function of time and r
      nb = 100
      tl = -5.
      ul =  5.
      nb2 = 100
      tl2 = -4.
      ul2 = 6.
      htit = 'N       vs log10(t) and log10(r)'
      nn = 360
      ulog = 3.
      do i=1,09
        htit(3:7) = ccla(i)
        call hbook2(nn+i,htit,nb,tl,ul,nb2,tl2,ul2,ulog)
      enddo

c  particle numbers as function of energy and t
      nb = 80 ! 100
      tl = -4. ! -4.
      ul =  4. ! 6.
      nb2 = 100
      tl2 = -5.
      ul2 =  5.
      htit = 'N       vs log10(E) and log10(t)'
      nn = 380
      ulog = 3.
      do i=1,09
        htit(3:7) = ccla(i)
        call hbook2(nn+i,htit,nb,tl,ul,nb2,tl2,ul2,ulog)
      enddo

      return
      end
c=======================================================================

      subroutine particles(iblklen)

c-----------------------------------------------------------------------
c  reads all data blocks from unit 'lunin' for current input file
c  and calls analyses routines for showers to be used
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z), integer (i-n)
      parameter (nblklen=312)
      character chtitle(480)*40
      double precision qhistos(10020,480)
      common /histch/ chtitle
      common /histos/ qhistos,rpromin,rpromax
     +       ,eparmin,eparmax
     +       ,xmin348,xmax348,ymin348,ymax348
      common /buffer/ outbuf(nblklen,21)
      common /buffch/ crunh,crune,cevth,cevte,clong
      character*4 crunh,crune,cevth,cevte,clong
      real     ddata(6552),outbuf,chablk(5),curpar
      equivalence(chablk(1),crunh)
c-----------------------------------------------------------------------
      double precision pama,pi
      common /pconst/pama(6000),pi
c-----------------------------------------------------------------------
c  variables are the same as used in corsika :
c  czx1   =  lateral distribution in x - direction for 1. level
c   ...
c  czxy2  =  lateral distribution in xy- direction for 2. level etc.
c   ...
c  nenkg1(i) =  Number of electrons in steps of 100 g/cm**2
c  age1  (i) =  age of shower in steps of 100 g/cm**2
c  dist  (i) =  distance bins for NKG lateral dist
c  lage1 (i) =  local age level 1
c  nenkg2(i) =  nenkg1(i);  dummy
c  age2  (i) =  age1(i);  dummy
c  tlev  (i) =  levels for NKG output in g/cm^2
c  cmlev (i) =  levels for NKG output in cm
c  dummy     =  unused words
      parameter (npd=10)
      common /nkg/czx1 (-npd:npd),  czy1 (-npd:npd),
     *            czxy1(-npd:npd),  czyx1(-npd:npd),
     *            czx2 (-npd:npd),  czy2 (-npd:npd),
     *            czxy2(-npd:npd),  czyx2(-npd:npd),
     *            nenkg1(npd), age1(npd), dist(npd), lage1(npd),
     *            tlev(npd),cmlev(npd),
     *            dummy(2*npd)
      real*4 nkgpar(248),nenkg1,age1,dist,lage1,tlev,cmlev
      equivalence (nkgpar(1),czx1(-npd))
      common /partic/ curpar(8),ident,igen,iobslv,ityp,
     *                epart,ppart,r,tf,costhe,phi
      common / phys / obslev(10),enspec(3),ffegs,ffnkg,cut(4),
     *                consts(230),aatm(5),batm(5),catm(5),flag(4)
      common /runpar/ heighp,thetap,phip,xoff,yoff,tcen,dintmod,dintcrs,
     *                nobslv,ishow,ishowu,ishowf,iblk,ifile,iret1,iret2,
     *                rmax,lfnkg,lfegs
      double precision heighp,thetap,phip,tcen,dintmod,dintcrs,rmax
      logical*1 lfnkg,lfegs
      common /shower/ profak1,profak2,profak3,
     *   ccla
      double precision profak1,profak2,profak3
      character*5 ccla(11)
      parameter (nmax=50000)
      common /steer/  lunmon,lundeb,lunin,lunhis,nfiles
      common /steerc/ chfile(nmax),chhist,chlong
      character*80    chfile,chhist,chlong
      character*8 use(2) /'used    ','ignored '/
      logical*1 first /.true./

c-----------------------------------------------------------------------
c  beginning of new shower
c-----------------------------------------------------------------------

      isb  = 2
    1 continue
      ishow  = ishow  + 1
      ishowf = ishowf + 1

c  decide whether shower is used or ignored
      if ( ishowf .lt. 0 ) then
        iu = 2
      else
        iu = 1
        ishowu=ishowu + 1
      endif

      if ( iu .eq. 1 ) then

c  calculate offsets of shower fronts for the observation level
        heighp = outbuf( 7,isb)
        thetap = outbuf(11,isb)
        phip   = outbuf(12,isb)
        xoff   = -(heighp-obslev(1)) * tan(thetap) * cos(phip)
        yoff   = -(heighp-obslev(1)) * tan(thetap) * sin(phip)

c  get maximal radius for radial thinning (in m)
        rmax = 1.d-2 * outbuf(152,isb)

c  arrival time of shower core at observation level
        tcen = sqrt((heighp-obslev(1))**2+xoff**2+yoff**2)/29.9792425

c  factors for projection to the plane perpendicular to the shower plane
        profak1 = cos(phip)**2*cos(thetap)+sin(phip)**2
        profak2 = cos(phip)*sin(phip)*(cos(thetap)-1.)
        profak3 = sin(phip)**2*cos(thetap)+cos(phip)**2

c  get interaction model and cross sections
c  venus  = 1, qgsjet = 2, sibyll = 3, dpmjet = 4, hdpm = 5
        if ( outbuf(76,isb) .eq. 1. ) then
          dintmod = 1.
        elseif ( outbuf(139,isb) .eq. 1. ) then
          dintmod = 3.
        elseif ( outbuf(141,isb) .eq. 1. ) then
          dintmod = 2.
        elseif ( outbuf(143,isb) .eq. 1. ) then
          dintmod = 4.
        else
          dintmod = 5.
        endif
        if ( outbuf(145,isb) .eq. 1. ) then
          dintcrs = 1.
        elseif ( outbuf(140,isb) .eq. 1. ) then
          dintcrs = 3.
        elseif ( outbuf(142,isb) .eq. 1. ) then
          dintcrs = 2.
        elseif ( outbuf(144,isb) .eq. 1. ) then
          dintcrs = 4.
        else
          dintcrs = 5.
        endif

      endif

c-----------------------------------------------------------------------
c  get next subblock

    2 continue
      isb = isb + 1
      if ( isb .gt. 21 ) then
        read(lunin,end=112,err=114) (ddata(i),i=1,iblklen*21)
        iblk = iblk + 1
        do  isb=1,21
        do   ib=1,iblklen
           outbuf(ib,isb) = ddata(ib+iblklen*(isb-1))
        enddo
        enddo
        isb = 1
      endif

c  shower finished ? 
      if ( 3397.3 .lt. outbuf(1,isb) .and.
     +                 outbuf(1,isb) .lt. 3397.5 ) goto 5
      lh = 0
      if ( iu .eq. 2 ) lh = 38

c-----------------------------------------------------------------------
c  next particle. transfer particle of datablock to current particle
    3 continue
      iblpart = 7
      if ( iblklen .eq. 312 ) iblpart = 8
      do i=1,iblpart
        curpar(i) = outbuf(iblpart*lh+i,isb)
      enddo

c  shower finished if 0 in the word where particle code should be
      if ( curpar(1) .eq. 0. ) then
        goto 4
      endif

c-----------------------------------------------------------------------
c  call standard analysis routines and
c  optional user analysis routines for each particles
c  if the shower is to be used for the analysis

      if ( iu .eq. 1 ) then
        call analys(iblklen)
      endif

      if ( lh .eq. 38 ) goto 2
      lh = lh + 1
      goto 3

c-----------------------------------------------------------------------
c  shower finished
c-----------------------------------------------------------------------
    4 continue
      isb = isb + 1
      if ( isb .gt. 21 ) then
        read(lunin,end=112,err=114) (ddata(i),i=1,iblklen*21)
        iblk = iblk + 1
        do  isb=1,21
        do   ib=1,iblklen
           outbuf(ib,isb) = ddata(ib+iblklen*(isb-1))
        enddo
        enddo
        isb = 1
      endif

c  next subblock a LONG block?
      if ( 52815.2 .lt. outbuf(1,isb) .and.
     +                  outbuf(1,isb) .lt. 52815.4 ) then
        write(*,*) 'LONG block detected, ignore for the moment'
        goto 4
      endif
      if ( .not. ( 3397.3 .lt. outbuf(1,isb) .and.
     +                         outbuf(1,isb) .lt. 3397.5 ) ) then
        iret1 = 7
        write(*,*) 'no EVTE = evt end after particle code 0 '
        return
      endif

    5 continue

c  print statistics for shower (from end of event block)
      if ( ifile .eq. 1 ) write(*,'(''  '')')
      if ( ifile .eq. 1 .and. nfiles .eq. 1 ) then
         write(*,'(i6,''. gam,el,had,mu,sum:'',1p,5e14.6)')
     +      ifile,(outbuf(i,isb),i=3,7)
      else 
      if ( ifile .lt. 10 ) then
         if ( nfiles .gt. 1 ) then
            if ( outbuf(7,isb) .gt. 0. )
     +         write(*,'(i6,''. gam,el,had,mu,sum:'',1p,5e14.6)')
     +            ifile,(outbuf(i,isb),i=3,7)
         endif
      elseif ( nfiles .ge. 10 .and. ifile .eq. nfiles ) then
         write(*,'(i6,''. gam,el,had,mu,sum:'',1p,5e14.6)')
     +      ifile,(outbuf(i,isb),i=3,7)
      endif
      endif

c-----------------------------------------------------------------------

c  define additional shower parameter from nkg-routine
      if ( lfnkg ) then
        do i=1,248
          nkgpar(i) = outbuf(i+7,isb)
        enddo
      endif

c  next shower to be read ?
      if ( 0 .eq. ishowf ) then
        iret1 = 6
        return
      endif

c  get next subblock
      isb = isb + 1
      if ( isb .gt. 21 ) then
        read(lunin,end=113,err=114) (ddata(i),i=1,iblklen*21)
        iblk = iblk + 1
        do  isb=1,21
        do   ib=1,iblklen
           outbuf(ib,isb) = ddata(ib+iblklen*(isb-1))
        enddo
        enddo
        isb = 1
      endif

      if ( 217433.0 .lt. outbuf(1,isb) .and.
     +                   outbuf(1,isb) .lt. 217433.2 ) then
c  new event header
        goto 1
      elseif ( 3301.3 .lt. outbuf(1,isb) .and.
     +                     outbuf(1,isb) .lt. 3301.5 ) then
c  run end found 
        iret1 = 0
        return
      endif

c-----------------------------------------------------------------------
c  read marks: error end
  111 continue
      ! write(*,*) 'eof detected when reading '
      iret1 = 5 
      return
  112 continue
      ! write(*,*) ' unexpected end-of-file '
      ishowu = ishowu - 1
      iret1 = 4
      return
  113 continue
      ! write(*,*) ' eof at reading lunin after "evte" '
      return
  114 continue
      ! write(*,*) ' error at reading '
      ishowu = ishowu - 1
      iret1 = 3
      return

      end
c=======================================================================

      subroutine analys(iblklen)

c-----------------------------------------------------------------------
c  analysis routine:  called for each particle from subroutine next
c                     if shower is used
c-----------------------------------------------------------------------

      implicit double precision (a-h,o-z), integer (i-n)
      parameter (nblklen=312)
      character chtitle(480)*40
      double precision qhistos(10020,480)
      common /histch/ chtitle
      common /histos/ qhistos,rpromin,rpromax
     +       ,eparmin,eparmax
     +       ,xmin348,xmax348,ymin348,ymax348
      common /buffer/ outbuf(nblklen,21)
      common /buffch/ crunh,crune,cevth,cevte,clong
      character*4 crunh,crune,cevth,cevte,clong
      real      outbuf,chablk(5),curpar
      equivalence(chablk(1),crunh)
c-----------------------------------------------------------------------
      double precision pama,pi
      common /pconst/pama(6000),pi
      common /partic/ curpar(8),ident,igen,iobslv,ityp,
     *                epart,ppart,r,tf,costhe,phi
      common / phys / obslev(10),enspec(3),
     *                ffegs,ffnkg,cut(4),
     *                consts(230),aatm(5),batm(5),catm(5),flag(4)
      common /runpar/ heighp,thetap,phip,xoff,yoff,tcen,dintmod,dintcrs,
     *                nobslv,ishow,ishowu,ishowf,iblk,ifile,iret1,iret2,
     *                rmax,lfnkg,lfegs
      double precision heighp,thetap,phip,tcen,dintmod,dintcrs,rmax
      logical*1 lfnkg,lfegs
      common /shower/ profak1,profak2,profak3,
     *   ccla
      double precision profak1,profak2,profak3
      character*5 ccla(11)
      parameter (nmax=50000)
      common /steer/  lunmon,lundeb,lunin,lunhis,nfiles
      common /steerc/ chfile(nmax),chhist,chlong
      character*80    chfile,chhist,chlong
      double precision dident,rrplog,timlog,timlogw
      integer icount /0/

c-----------------------------------------------------------------------

c do not analyze the muon additional information
      if ( int(curpar(1)/1000) .eq. 75 .or.
     *     int(curpar(1)/1000) .eq. 76 ) return
      if ( int(curpar(1)/1000) .eq. 85 .or.
     *     int(curpar(1)/1000) .eq. 86 ) return
      if ( int(curpar(1)/1000) .eq. 95 .or.
     *     int(curpar(1)/1000) .eq. 96 ) return
      if ( curpar(1) .lt. 1001. ) return

c  is particle from the correct observation level ?
      iobslv= 1 ! mod( int(curpar(1)) , 10 )

c  unpack information for normal particles
      if ( curpar(1) .lt. 1.e6 ) then
        igen  = mod( int(curpar(1)), 1000 ) / 10
        ident = curpar(1) / 1000
        px = curpar(2)
        py = curpar(3)
        pz = curpar(4)
        x  = curpar(5)
        y  = curpar(6)
        t  = curpar(7)
        wt = 1.
        if ( iblklen .eq. 312 ) wt = curpar(8)
c  calculate momentum and energy in GeV
        psq   = px**2 + py**2 + pz**2
        ppart = sqrt( psq )
        epart = sqrt( pama(ident)**2 + psq )
        eplog = log10(epart)
c  calculate cosine of theta and phi
        costhe = pz / ppart
        if ( px .eq. 0.  .and.  py .eq. 0. ) then
          phi = 0.
        else
          phi = atan2( py , px )
        endif

c  unpack information for cherenkov bunches
c  no momentum, just direction available
      elseif ( curpar(1) .ge. 1.e6 ) then
        igen = 0
        ident = int(curpar(1)/1.e3)
        xnpho = int((curpar(1)-ident*1.e3)/10.)
        www = curpar(4)**2 + curpar(5)**2
        costhe = -sqrt(1. - www)
        theta = acos(costhe)
        if ( www .le. 0. ) then
          phi = 0.
        else
          phi = -atan2(curpar(5),curpar(4))
        endif
        x  = curpar(2)
        y  = curpar(3)
        cx = curpar(4) ! direction cosine to x-axis.
        cy = curpar(5) ! direction cosine to y-axis.
        t  = curpar(6)
        wt = 1.
        if ( iblklen .eq. 312 ) wt = (curpar(1)-99.e5)/10.
        curpar(8) = wt
      endif

c  calculate arrival time of spherical shower front at point(x,y)
      tf = sqrt(  ( heighp-obslev(1) ) **2  +
     *             ( x-xoff ) **2  + ( y-yoff ) **2  ) / 29.9792425
c  get the arrival time of particle with respect to the spherical shower front
c  this a good quantity to get a local time distribution at one detector
      trel = t - tf
c  however, with trel values one cannot reconstruct a shower front
c  since tf is subtracted from the times and tf is dependent on the
c  position in x and y.
c  To get the correct arrival times at different detectors subtract from
c  the simple arrival time t the time tcen when the shower axis hits the ground.
c  (see evaluation of tcen in routine next)
c  this difference is negative for about half the particles in inclined showers
      tabs = t - tcen
      if (trel .gt. 0 ) then
        timlog = log10(trel)
      else
        timlog = -1.e30
      endif

c  calculate distance from shower axis in meters
      xx = x * 1.e-2
      yy = y * 1.e-2
      rr = sqrt(xx**2 + yy**2)
c  calculate the radius projected on to the plane perpendicular to the shower axis
      xxproj = xx * profak1 + yy * profak2
      yyproj = xx * profak2 + yy * profak3
      rrproj = sqrt(xxproj**2 + yyproj**2)
c  use for plotting the radius in the plane perpendicular to the shower axis
      rrplog = log10(rrproj)

c  now filling of histograms starts with overall 2dim-histogram
      if ( 1 .le. ident .and. ident .lt. 30 ) then
        dident = ident
        call hfill2(299,dident,rrplog,wt)
      endif

c  define particle type and count particles
c  ...... EGS photons
      idi = 11
      if     ( ident .eq. 1 ) then
        idi = 1
c  ...... positrons
      elseif ( ident .eq. 2 ) then
        idi = 2
c  ...... electrons
      elseif ( ident .eq. 3 ) then
        idi = 3
c  ...... muons+
      elseif ( ident .eq. 5 ) then
        idi = 4
c  ...... muons-
      elseif ( ident .eq. 6 ) then
        idi = 5
c  ...... pions+
      elseif ( ident .eq. 8 ) then
        idi = 6
c  ...... pions-
      elseif ( ident .eq. 9 ) then
        idi = 7
c  ...... proton/antiproton
      elseif ( ident .eq. 14 .or. ident .eq. 15 ) then
        idi = 8
        if ( rrplog .lt. ymin348 ) ymin348 = rrplog
        if ( rrplog .gt. ymax348 ) ymax348 = rrplog 
        if ( epart .lt. eparmin ) eparmin = epart
        if ( epart .gt. eparmax ) eparmax = epart
        if ( eplog .lt. xmin348 ) xmin348 = eplog
        if ( eplog .gt. xmax348 ) xmax348 = eplog
c  ...... neutron/antineutron
      elseif ( ident .eq. 13 .or. ident .eq. 25 ) then
        idi = 9
c  ...... cerenkov photons
      elseif ( ident .ge. 9900 ) then
c       idi = 10
      else
c  ...... other particles
c       idi = 11
      endif

c  fill histograms
c  (attention: there is a slight inconsistency in that corsika
c  applies radial thinning for particles within rmax from the
c  shower core in the observation level.
c  here however we plot radial distributions where the radius
c  is calculated in the plane perpendicular to the shower axis.)
c  make the radius cut only for those histograms that do not show
c  dependence on radius
        call hfill1(120+idi,rrplog,0.d0,wt)
        call hfill2(340+idi,eplog,rrplog,wt)
        call hfill2(360+idi,timlog,rrplog,wt)
***     if ( rrproj .gt. rmax ) then
          call hfill1(180+idi,eplog,0.d0,wt)
          call hfill1(220+idi,timlog,0.d0,wt)
          call hfill2(380+idi,eplog,timlog,wt)
***     endif

      return
      end
c=======================================================================

      subroutine pamaf

c-----------------------------------------------------------------------
c  pa(rticle) ma(ss) f(illing)
c
c  fills particle mass for particle ip in array pama
c  resonances and strange baryons included
c  particle masses according to geant table,
c  taken from the periodic table
c  or calculated with the mass formula of weizsaecker
c  this subroutine is called from start
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z), integer (i-n)

      double precision pama,pi
      common /pconst/pama(6000),pi

      double precision masses(75),charge(75),amus(59,14)
      data masses /
     * 0.0        , 0.000511   , 0.000511   , 0.0        , 0.105658   ,
     * 0.105658   , 0.134973   , 0.139568   , 0.139568   , 0.497671   ,
     * 0.493646   , 0.493646   , 0.939566   , 0.938272   , 0.938272   ,
     * 0.497671   , 0.5488     , 1.11563    , 1.18937    , 1.19255    ,
     * 1.19743    , 1.3149     , 1.32132    , 1.67243    , 0.939566   ,
     * 1.11563    , 1.18937    , 1.19255    , 1.19743    , 1.3149     ,
     * 1.32132    , 1.67243    , 1.7841     , 1.7841     , 1.8693     ,
     * 1.8693     , 1.8645     , 1.8645     , 1.9693     , 1.9693     ,
     * 2.2852     , 80.6       , 80.6       , 91.161     , 1.877      ,
     * 2.817      , 3.755      , 0.0        , 0.0        , 0.0        ,
     * 0.7669     , 0.7681     , 0.7681     , 1.2309     , 1.2323     ,
     * 1.2336     , 1.2349     , 1.2309     , 1.2323     , 1.2336     ,
     * 1.2349     , 0.89624    , 0.89209    , 0.89209    , 0.89624    ,
     * 0.0        , 0.0        , 0.0        , 0.0        , 0.0        ,
     * 0.5488     , 0.5488     , 0.5488     , 0.5488     , 0.0        /

c  isotope masses calculated from: atomic data and nucl.data tables 39
c  (1988) 289, (wapstra's values, corrected for electron masses)
      data ((amus(i,l),i=1,59),l=1,7) /
     * 1.8756  ,  2.8089  ,                                    57*0.  ,
     * 2.8083  ,  3.7273  ,  4.6678  ,  5.6054  ,  6.5454  ,   54*0.  ,
     * 2*0.  ,  5.6014  ,  6.5337  ,  7.4712  ,  8.4067  ,
     *              9.3471  , 10.2856  ,                       51*0.  ,
     * 2*0.  ,  6.5341  ,  7.4547  ,  8.3926  ,  9.3253  ,
     *             10.2644  , 11.2008  ,                       51*0.  ,
     * 2*0.  ,  7.4722  ,  8.3932  ,  9.3243  , 10.2524  ,
     *             11.1886  , 12.1232  , 13.0618  , 13.9986  , 49*0.  ,
     * 2*0.  ,  8.4091  ,  9.3274  , 10.2538  , 11.1747  , 12.1093  ,
     *             13.0406  , 13.9790  , 14.9143  , 15.8531  , 48*0.  ,
     * 4*0.  , 11.1915  , 12.1110  , 13.0400  , 13.9687  , 14.9057  ,
     *             15.8394  , 16.7761  , 17.7104  ,            47*0.  /
      data ((amus(i,l),i=1,59),l=8,14) /
     * 4*0.  , 12.1282  , 13.0446  , 13.9709  , 14.8948  , 15.8302  ,
     *             16.7617  , 17.6973  , 18.6293  , 19.5650  , 46*0.  ,
     * 7*0.  , 15.8325  , 16.7629  , 17.6920  , 18.6429  , 19.5564  ,
     *             20.4907  , 21.4227  , 22.3587  ,            44*0.  ,
     * 6*0.  , 15.8464  , 16.7668  , 17.6947  , 18.6174  , 19.5502  ,
     *  20.4794  , 21.4137  , 22.3444  , 23.2839  , 24.2138  , 43*0.  ,
     * 8*0.  , 18.6308  , 19.5532  , 20.4817  , 21.4088  , 22.3414  ,
     *  23.2720  , 24.2059  , 25.1387  , 26.0746  , 27.0099  ,
     *  27.9469  , 28.8820  , 29.8173  , 30.7546  , 31.6913  , 36*0.  ,
     * 7*0.  , 18.6410  , 19.5658  , 20.4860  , 21.4124  , 22.3354  ,
     *  23.2676  , 24.1961  , 25.1292  , 26.0602  , 26.9961  ,
     *  27.9291  , 28.8660  , 29.7994  , 30.7376  ,            38*0.  ,
     * 9*0.  , 21.4241  , 22.3488  , 23.2714  , 24.1996  , 25.1261  ,
     *  26.0579  , 26.9880  , 27.9218  , 28.8541  , 29.7894  ,
     *  30.7233  , 31.6599  , 32.5944  , 33.5316  ,            36*0.  ,
     * 9*0.  , 22.3591  , 23.2836  , 24.2041  , 25.1304  , 26.0527  ,
     *  26.9838  , 27.9128  , 28.8457  , 29.7761  , 30.7111  ,
     *  31.6431  , 32.5803  , 33.5128  , 34.4505  , 35.3837  , 35*0.  /

c-----------------------------------------------------------------------

c  geant particles  including rho, k*, and delta
      do ip = 1,75
        pama(ip) = masses(ip)
      enddo
      do ip = 76,6000
        pama(ip) = 0.
      enddo

      do ia = 1,59
      do ic = 1,ia
        in = ia - ic
        ip = ia * 100 + ic
        if ( ic .le. 14 ) then
c  masses from mass table for isotopes
          if ( in .eq. 0 ) then
            pama(ip) = ic * pama(14)
          else
            pama(ip) = amus(in,ic)
          endif
c  simple sum of proton and neutron masses
          if ( pama(ip) .eq. 0.   )
     *               pama(ip) = pama(14) * ic + in * pama(13)
        else
c  weizsaeckers mass formula gives binding energy in mev
          b1 =  14.1  * ia
          b2 = -13.   * ia**(2./3.)
          b3 = -0.595 * ic**2 / ia**(1./3.)
          b4 = -19.   * (ic-in)**2 / ia
          b5 =  33.5  / ia**0.75
          if ( mod(ic,2) .eq. 0 .and. mod(in,2) .eq. 0 ) then
            ss =  1.
          elseif ( mod(ic,2) .eq. 1 .and. mod(in,2) .eq. 1 ) then
            ss = -1.
          else
            ss =  0.
          endif
          bind = (b1 + b2 + b3 + b4 + ss*b5) * 1.e-3
          pama(ip) = masses(13) * in + masses(14) * ic - bind
        endif
      enddo
      enddo

c  masses of multineutron clusters
      do in = 1,59
        ip = 100 * in
        pama(ip) = in * pama(13)
      enddo

      return
      end
c=======================================================================

      subroutine hbook1(ihist,htit,nbin,xlow,xupr,ulog)

      implicit double precision (a-h,o-z), integer (i-n)
      character chtitle(480)*40,htit*(*)
      double precision qhistos(10020,480)
      common /histch/ chtitle
      common /histos/ qhistos,rpromin,rpromax
     +       ,eparmin,eparmax
     +       ,xmin348,xmax348,ymin348,ymax348

      chtitle(ihist) = '                                        '
      chtitle(ihist) = htit
      qhistos(1,ihist) = 1.d0*ihist
      qhistos(2,ihist) = 1.d0*nbin
      qhistos(3,ihist) = xlow
      qhistos(4,ihist) = xupr 
      qhistos(5,ihist) = ulog
      qhistos(6,ihist) = 0.d0
      qhistos(7,ihist) = 0.d0
      qhistos(8,ihist) = 0.d0
      qhistos(9,ihist) = -99.8877d0
      do  i=10,nbin+10
         qhistos(i,ihist) = 0.
      enddo

      return
      end
c=======================================================================

      subroutine hbook2(ihist,htit,nbin,xlow,xupr,nbi2,xlo2,xup2,ulog)

      implicit double precision (a-h,o-z), integer (i-n)
      character chtitle(480)*40,htit*(*)
      double precision qhistos(10020,480)
      common /histch/ chtitle
      common /histos/ qhistos,rpromin,rpromax
     +       ,eparmin,eparmax
     +       ,xmin348,xmax348,ymin348,ymax348

      chtitle(ihist) = '                                        '
      chtitle(ihist) = htit
      qhistos(1,ihist) = 1.d0*ihist
      qhistos(2,ihist) = 1.d0*nbin
      qhistos(3,ihist) = xlow
      qhistos(4,ihist) = xupr 
      qhistos(5,ihist) = ulog
      qhistos(6,ihist) = 1.d0*nbi2
      qhistos(7,ihist) = xlo2
      qhistos(8,ihist) = xup2 
      qhistos(9,ihist) = -99.8877d0
      do  i=10,nbin*nbi2+10
         qhistos(i,ihist) = 0.
      enddo

      return
      end
c=======================================================================

      subroutine hfill1(ihist,xval,yval,wval)

      implicit double precision (a-h,o-z), integer (i-n)
      character chtitle(480)*40
      double precision qhistos(10020,480)
      common /histch/ chtitle
      common /histos/ qhistos,rpromin,rpromax
     +       ,eparmin,eparmax
     +       ,xmin348,xmax348,ymin348,ymax348

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      nbin = int(qhistos(2,ihist)) 
      if ( nbin .gt. 0 ) then
      if ( 121 .le. ihist .and. ihist .le. 231 ) then
c 128 (1)   rec  density of p     vs log10(r)
c 188 (1)   N rec p     vs log10(energy)
c 228 (1)   N rec p     vs log10(time)
      if ( xval .le. qhistos(4,ihist) ) then   
         if ( xval .gt. qhistos(3,ihist) ) then ! regular interval.
            ix = (xval - qhistos(3,ihist)) * qhistos(2,ihist) /
     /         (qhistos(4,ihist) - qhistos(3,ihist))
            if ( ix .ge. 0 ) then
               qhistos(11+ix,ihist) = qhistos(11+ix,ihist) + wval 
            else
               qhistos(10,ihist) = qhistos(10,ihist) + wval
            endif
         else ! underflow.
            qhistos(10,ihist) = qhistos(10,ihist) + wval
         endif
      else ! overflow.
         write(*,*) '    ih=',ihist,xval,wval,'  overflow'
         qhistos(10+nbin+1,ihist) = qhistos(10+nbin+1,ihist) + wval  
      endif
      endif
      endif

      return
      end
c=======================================================================

      subroutine hfill2(ihist,xval,yval,wval)

      implicit double precision (a-h,o-z), integer (i-n)
      character chtitle(480)*40
      double precision qhistos(10020,480)
      common /histch/ chtitle
      common /histos/ qhistos,rpromin,rpromax
     +       ,eparmin,eparmax
     +       ,xmin348,xmax348,ymin348,ymax348

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if ( ihist .eq. 299 ) then
c particle codes (0,...,29) vs log10(r)              
c    299.     30.  0.0000       30.00          2.    100.  -4.000       6.000
         if ( yval .lt. qhistos(8,ihist) ) then
         if ( yval .gt. qhistos(7,ihist) ) then
            ix = int(xval)
            iy = 1. + (yval - qhistos(7,ihist)) * qhistos(6,ihist) /
     /         (qhistos(8,ihist) - qhistos(7,ihist))
            iq = ix + (iy-1) * int(qhistos(2,ihist))
            qhistos(10+iq,ihist) = qhistos(10+iq,ihist) + wval
         else ! y-underflow
            iq = 10
            qhistos(iq,ihist) = qhistos(iq,ihist) + wval
         endif    
         else ! y-overflow
            iq = 10+int(qhistos(2,ihist))*int(qhistos(6,ihist))+1
            qhistos(iq,ihist) = qhistos(iq,ihist) + wval
         endif
      endif 

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if ( ( 341 .le. ihist .and. ihist .le. 349 ) .or.
     +     ( 361 .le. ihist .and. ihist .le. 369 ) .or.
     +     ( 381 .le. ihist .and. ihist .le. 389 ) ) then
c N p     vs log10(E) and log10(r)        
c    348.    80.  -4.000       4.000          3.    100.  -4.000      6.000
c N p     vs log10(t) and log10(r)        
c    368.    100.  -4.000       6.000          3.    100.   0.000       5.000
c N p     vs log10(E) and log10(t)        
c    388.    100.  -4.000       6.000          3.    100.   0.000       5.000
         if ( yval .lt. qhistos(8,ihist) ) then
            if ( yval .gt. qhistos(7,ihist) ) then
               if (  xval .lt. qhistos(4,ihist) ) then
                  if ( xval .ge. qhistos(3,ihist) ) then
                  ix = 1. + (xval-qhistos(3,ihist)) * qhistos(2,ihist) /
     /                 (qhistos(4,ihist)-qhistos(3,ihist))
                  iy = 1. + (yval-qhistos(7,ihist)) * qhistos(6,ihist) /
     /                 (qhistos(8,ihist)-qhistos(7,ihist))
                     iq = 10 + ix + (iy-1) * int(qhistos(2,ihist))
                     qhistos(iq,ihist) = qhistos(iq,ihist) + wval
                  else ! x-underflow
                     iq = 10
                     qhistos(iq,ihist) = qhistos(iq,ihist) + wval
                  endif 
               else ! x-overflow
                  iq = 10+int(qhistos(2,ihist))*int(qhistos(6,ihist))+1
                  qhistos(iq,ihist) = qhistos(iq,ihist) + wval
               endif
            else ! y-underflow
               iq = 10
               qhistos(iq,ihist) = qhistos(iq,ihist) + wval
            endif
         else ! y-overflow
            iq = 10+int(qhistos(2,ihist))*int(qhistos(6,ihist))+1
            qhistos(iq,ihist) = qhistos(iq,ihist) + wval
         endif
      endif 

      return
      end
c=======================================================================

      subroutine tabhistos

c-----------------------------------------------------------------------
c  Full tabular of histograms including ident number and dimension
c-----------------------------------------------------------------------
c
c         1 (1)   longi dist of gammas
c         2 (1)   longi dist of e+
c         3 (1)   longi dist of e-
c         4 (1)   longi dist of mu+
c         5 (1)   longi dist of mu-
c         6 (1)   longi dist of hadrons
c         7 (1)   longi dist of charged
c         8 (1)   longi dist of nuclei
c         9 (1)   longi dist of cerenk
c
c        11 (1)   longi edep of gammas
c        12 (1)   longi edep of em ioniz
c        13 (1)   longi edep of em cut
c        14 (1)   longi edep of mu ioniz
c        15 (1)   longi edep of mu cut
c        16 (1)   longi edep of had ioniz
c        17 (1)   longi edep of had cut
c        18 (1)   longi edep of neutrino
c        19 (1)   longi edep of sum
c
c        21 (1)   longi dist of nucleon
c        22 (1)   longi dist of p
c        23 (1)   longi dist of n
c        24 (1)   longi dist of pi+/-
c        25 (1)   longi dist of K+/-
c        26 (1)   longi dist of Kl
c        27 (1)   longi dist of Ks
c        28 (1)   longi dist of Kl/s
c  
c        99 (2)   particle codes  vs  log10(r)    (here on 299)
c
c       101 (1)   part density of gam    vs  log10(r)
c       102 (1)   part density of e+     vs  log10(r)
c       103 (1)   part density of e-     vs  log10(r)
c       104 (1)   part density of mu+    vs  log10(r)
c       105 (1)   part density of mu-    vs  log10(r)
c       106 (1)   part density of pi+    vs  log10(r)
c       107 (1)   part density of pi-    vs  log10(r)
c       108 (1)   part density of p      vs  log10(r)
c       109 (1)   part density of n      vs  log10(r)
c       110 (1)   part density of cerph  vs  log10(r)
c       111 (1)   part density of other  vs  log10(r)
c
c       121 (1)   rec  density of gam    vs  log10(r)  "always weight 1" 
c       122 (1)   rec  density of e+     vs  log10(r)
c       123 (1)   rec  density of e-     vs  log10(r)
c       124 (1)   rec  density of mu+    vs  log10(r)
c       125 (1)   rec  density of mu-    vs  log10(r)
c       126 (1)   rec  density of pi+    vs  log10(r)
c       127 (1)   rec  density of pi-    vs  log10(r)
c       128 (1)   rec  density of p      vs  log10(r)
c       129 (1)   rec  density of n      vs  log10(r)
c       130 (1)   rec  density of cerph  vs  log10(r)
c       131 (1)   rec  density of other  vs  log10(r)
c
c       141 (1)   ener density of gam    vs  log10(r)
c       142 (1)   ener density of e+     vs  log10(r)
c       143 (1)   ener density of e-     vs  log10(r)
c       144 (1)   ener density of mu+    vs  log10(r)
c       145 (1)   ener density of mu-    vs  log10(r)
c       146 (1)   ener density of pi+    vs  log10(r)
c       147 (1)   ener density of pi-    vs  log10(r)
c       148 (1)   ener density of p      vs  log10(r)
c       149 (1)   ener density of n      vs  log10(r)
c       150 (1)   ener density of cerph  vs  log10(r)
c       151 (1)   ener density of other  vs  log10(r)
c
c       161 (1)   N gam      vs  log10(energy)
c       162 (1)   N e+       vs  log10(energy)
c       163 (1)   N e-       vs  log10(energy)
c       164 (1)   N mu+      vs  log10(energy)
c       165 (1)   N mu-      vs  log10(energy)
c       166 (1)   N pi+      vs  log10(energy)
c       167 (1)   N pi-      vs  log10(energy)
c       168 (1)   N p        vs  log10(energy)
c       169 (1)   N n        vs  log10(energy)
c       170 (1)   N cerph    vs  log10(energy)
c       171 (1)   N other    vs  log10(energy)
c
c       181 (1)   N rec gam    vs  log10(energy)   "always weight 1"
c       182 (1)   N rec e+     vs  log10(energy)
c       183 (1)   N rec e-     vs  log10(energy)
c       184 (1)   N rec mu+    vs  log10(energy)
c       185 (1)   N rec mu-    vs  log10(energy)
c       186 (1)   N rec pi+    vs  log10(energy)
c       187 (1)   N rec pi-    vs  log10(energy)
c       188 (1)   N rec p      vs  log10(energy)
c       189 (1)   N rec n      vs  log10(energy)
c       190 (1)   N rec cerph  vs  log10(energy)
c       191 (1)   N rec other  vs  log10(energy)
c
c       201 (1)   N gam    vs  log10(time)
c       202 (1)   N e+     vs  log10(time)
c       203 (1)   N e-     vs  log10(time)
c       204 (1)   N mu+    vs  log10(time)
c       205 (1)   N mu-    vs  log10(time)
c       206 (1)   N pi+    vs  log10(time)
c       207 (1)   N pi-    vs  log10(time)
c       208 (1)   N p      vs  log10(time)
c       209 (1)   N n      vs  log10(time)
c       210 (1)   N cerph  vs  log10(time)
c       211 (1)   N other  vs  log10(time)
c
c       221 (1)   N rec gam    vs  log10(time)   "always weight 1"
c       222 (1)   N rec e+     vs  log10(time)
c       223 (1)   N rec e-     vs  log10(time)
c       224 (1)   N rec mu+    vs  log10(time)
c       225 (1)   N rec mu-    vs  log10(time)
c       226 (1)   N rec pi+    vs  log10(time)
c       227 (1)   N rec pi-    vs  log10(time)
c       228 (1)   N rec p      vs  log10(time)
c       229 (1)   N rec n      vs  log10(time)
c       230 (1)   N rec cerph  vs  log10(time)
c       231 (1)   N rec other  vs  log10(time)
c
c       241 (1)   N gam    vs  theta
c       242 (1)   N e+     vs  theta
c       243 (1)   N e-     vs  theta
c       244 (1)   N mu+    vs  theta
c       245 (1)   N mu-    vs  theta
c       246 (1)   N pi+    vs  theta
c       247 (1)   N pi-    vs  theta
c       248 (1)   N p      vs  theta
c       249 (1)   N n      vs  theta
c       250 (1)   N cerph  vs  theta
c       251 (1)   N other  vs  theta
c
c       261 (1)   N rec gam    vs  theta   "always weight 1"
c       262 (1)   N rec e+     vs  theta
c       263 (1)   N rec e-     vs  theta
c       264 (1)   N rec mu+    vs  theta
c       265 (1)   N rec mu-    vs  theta
c       266 (1)   N rec pi+    vs  theta
c       267 (1)   N rec pi-    vs  theta
c       268 (1)   N rec p      vs  theta
c       269 (1)   N rec n      vs  theta
c       270 (1)   N rec cerph  vs  theta
c       271 (1)   N rec other  vs  theta
c
c       281 (1)   N gam    vs  phi
c       282 (1)   N e+     vs  phi
c       283 (1)   N e-     vs  phi
c       284 (1)   N mu+    vs  phi
c       285 (1)   N mu-    vs  phi
c       286 (1)   N pi+    vs  phi
c       287 (1)   N pi-    vs  phi
c       288 (1)   N p      vs  phi
c       289 (1)   N n      vs  phi
c       290 (1)   N cerph  vs  phi
c       291 (1)   N other  vs  phi
c
c     ( 295 (1)   ring areas for particle densities )
c
c       299 (2)   particle codes  vs  log10(r)
c
c       301 (1)   N rec gam    vs  phi   "always weight 1"
c       302 (1)   N rec e+     vs  phi
c       303 (1)   N rec e-     vs  phi
c       304 (1)   N rec mu+    vs  phi
c       305 (1)   N rec mu-    vs  phi
c       306 (1)   N rec pi+    vs  phi
c       307 (1)   N rec pi-    vs  phi
c       308 (1)   N rec p      vs  phi
c       309 (1)   N rec n      vs  phi
c       310 (1)   N rec cerph  vs  phi
c       311 (1)   N rec other  vs  phi
c
c       321 (1)   N gam    vs  log10(weight)
c       322 (1)   N e+     vs  log10(weight)
c       323 (1)   N e-     vs  log10(weight)
c       324 (1)   N mu+    vs  log10(weight)
c       325 (1)   N mu-    vs  log10(weight)
c       326 (1)   N pi+    vs  log10(weight)
c       327 (1)   N pi-    vs  log10(weight)
c       328 (1)   N p      vs  log10(weight)
c       329 (1)   N n      vs  log10(weight)
c       330 (1)   N cerph  vs  log10(weight)
c       331 (1)   N other  vs  log10(weight)
c
c       341 (2)   N gam    vs  log10(E) and log10(r)
c       342 (2)   N e+     vs  log10(E) and log10(r)
c       343 (2)   N e-     vs  log10(E) and log10(r)
c       344 (2)   N mu+    vs  log10(E) and log10(r)
c       345 (2)   N mu-    vs  log10(E) and log10(r)
c       346 (2)   N pi+    vs  log10(E) and log10(r)
c       347 (2)   N pi-    vs  log10(E) and log10(r)
c       348 (2)   N p      vs  log10(E) and log10(r)
c       349 (2)   N n      vs  log10(E) and log10(r)
c       350 (2)   N cerph  vs  log10(E) and log10(r)
c       351 (2)   N other  vs  log10(E) and log10(r)
c
c       361 (2)   N gam    vs  log10(t) and log10(r)
c       362 (2)   N e+     vs  log10(t) and log10(r)
c       363 (2)   N e-     vs  log10(t) and log10(r)
c       364 (2)   N mu+    vs  log10(t) and log10(r)
c       365 (2)   N mu-    vs  log10(t) and log10(r)
c       366 (2)   N pi+    vs  log10(t) and log10(r)
c       367 (2)   N pi-    vs  log10(t) and log10(r)
c       368 (2)   N p      vs  log10(t) and log10(r)
c       369 (2)   N n      vs  log10(t) and log10(r)
c       370 (2)   N cerph  vs  log10(t) and log10(r)
c       371 (2)   N other  vs  log10(t) and log10(r)
c
c       381 (2)   N gam    vs  log10(E) and log10(t)
c       382 (2)   N e+     vs  log10(E) and log10(t)
c       383 (2)   N e-     vs  log10(E) and log10(t)
c       384 (2)   N mu+    vs  log10(E) and log10(t)
c       385 (2)   N mu-    vs  log10(E) and log10(t)
c       386 (2)   N pi+    vs  log10(E) and log10(t)
c       387 (2)   N pi-    vs  log10(E) and log10(t)
c       388 (2)   N p      vs  log10(E) and log10(t)
c       389 (2)   N n      vs  log10(E) and log10(t)
c       390 (2)   N cerph  vs  log10(E) and log10(t)
c       391 (2)   N other  vs  log10(E) and log10(t)
c
c-----------------------------------------------------------------------

      end
