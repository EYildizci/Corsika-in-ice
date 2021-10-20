c=======================================================================
c
c  s u m l o n g i f i l e s . f
c  -----------------------------
c     sum up content of a number of `*.long` files,
c     sum up content of a lot of `*.long` files for one single shower
c     in the case of a parallelized corsika simulation,
c     and print the sum of distribution and energy quantities.
c-----------------------------------------------------------------------
c compilation:
c     f77 -fbounds-check sumlongifiles.f -o sumlongifiles
c     gfortran -fbounds-check sumlongifiles.f -o sumlongifiles
c     ifort -C -check bounds sumlongifiles.f -o sumlongifiles
c execution:
c     ls -1 DAT002345-*.long > sumlongifiles.002345
c     ./sumlongifiles < sumlongifiles.i002345 > sumlongifiles.out002345
c     more sumlongifiles.sum002345
c-----------------------------------------------------------------------
cLONGITUDINAL DISTRIBUTION IN   ivs VERTICAL STEPS OF  fgr. G/CM**2 FOR SHOWER
cDEPTH     GAMMAS   POSITRONS   ELECTRONS         MU+         MU-     HADRONS 
c   5. 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00 0.00000e+00 1.00000e+00 
c  10. 1.19640e+04 1.03913e+02 1.57065e+02 0.00000e+00 0.00000e+00 1.00000e+00
c  15. 1.89921e+04 1.02062e+02 6.02637e+02 0.00000e+00 0.00000e+00 1.00000e+00
c  20. 2.21814e+04 6.25936e+02 7.22297e+02 0.00000e+00 0.00000e+00 1.00000e+00
c  25. 2.33861e+04 2.67479e+02 6.90617e+02 0.00000e+00 0.00000e+00 1.00000e+00
c-----------------------------------------------------------------------
c                                     juergen.oehlschlaeger@kit.edu
c-----------------------------------------------------------------------
 
      program sumlongifiles

      implicit double precision (a-h,o-z), integer (i-n)
      parameter (ldim=2800)
      character clongname(50000)*120, czeilong(4)*120, czeile*120
      character chlongsum*23, chpartic*29, chsteer*15, cformat*7
      dimension qdistrb(0:9,ldim), qenergy(0:9,ldim), qzeile(0:9)
      dimension qxdistr(ldim), qxenerg(ldim)
      dimension longlen(50000)
      real pdata(400)
      
c - - - - - - read all long files and keep names in a character array - -
      do  long=1,12345
         read(*,'(a)',err=101,end=102) clongname(long)
         jl = 120 + 1 
  100    continue 
         jl = jl - 1
         if ( clongname(long)(jl:jl) .eq. ' ' ) goto 100
         longlen(long) = jl       
      enddo
      goto 102
  101 continue
      write(*,*) ' e r r o r   reading clongname ',clongname(long)
  102 continue  
      long = long - 1
      read(clongname(1)(4:9),'(i6)') lrunnr      
      chlongsum = 'sumlongifiles.sum000000'
      write(chlongsum(18:23),'(i6.6)') lrunnr

c - - - - - - read quantities from parallel steering file  - - - - - - -
      parengy = 0.
      if ( long .ge. 2 ) then
         write(chsteer, '(''parallel-'',i6.6)') lrunnr 
         open(unit=2,file=chsteer,form='unformatted',status='old')
         do  is=1,lrunnr   
            read(2,'(a)',end=103,err=103) czeile
            if ( index(czeile,'ERANGE') .gt. 0 ) then
               ia = index(czeile(7:120),' ')
               ib = index(czeile(ia+6:120),'.')
               if ( ib .le. 0 ) ib = max(
     +          index(czeile(ia+6:120),'e'),index(czeile(ia+6:120),'E'),
     +          index(czeile(ia+6:120),'d'),index(czeile(ia+6:120),'D'))
               ic = index(czeile(ia+6:120),' ') - 1
               write(cformat,'(''(f'',i2,''.2)'')') ic-ia+1
               read(czeile(ia+6:ic),cformat) parengy 
               goto 103
            endif
         enddo
  103    continue
         close(unit=2)
      endif

c - - - - - - read first long file - - - - - - - - - - - - - - - - - - -
      open(unit=1,file=clongname(1)(1:longlen(1)),
     +            form='formatted',status='old')
      ! - - - - - clear all arrays:
      do  is=1,ldim
      do  il=0,9
         qdistrb(il,is) = 0.d0
         qenergy(il,is) = 0.d0
      enddo
      qxdistr(is) = 0.d0
      qxenerg(is) = 0.d0
      enddo
      ! - - - - - read first long table (check format):
      read(1,'(a)',end=119,err=118) czeilong(1)
      read(1,'(a)',end=119,err=118) czeilong(2)
      write(czeilong(1)(81:86),'(i6.6)') lrunnr
      ! - - - - - get number of steps from title line and check: 
      read(czeilong(1)(31:35),'(i5)') lsteps
      if ( lsteps .le. 0 .or. lsteps .ge. ldim ) lsteps = 3
      ! - - - - - read first 3 lines from long table and check grstep:
      do  is=1,3
         read(1,*) (qzeile(il),il=0,9)
         do  il=1,9
            qdistrb(il,is) = qzeile(il)
         enddo 
         qxdistr(is) = qzeile(0)
      enddo
      qgrams = qxdistr(3) - qxdistr(2)
      grstep = 0.
      if ( index(czeilong(1)(55:59),'.') .gt. 0 ) then
         read(czeilong(1)(55:59),'(f5.0)') grstep
      else 
         read(czeilong(1)(55:59),'(i5)') jl
         grstep = jl
      endif
      if ( grstep .le. 1. .or. grstep .ge. 345. ) grstep = qgrams
      ! - - - - - read other lines from long table:
      do  is=4,lsteps
         read(1,*,end=110,err=110) (qzeile(il),il=0,9)
         do  il=1,9
            qdistrb(il,is) = qzeile(il)
         enddo 
         qxdistr(is) = qzeile(0)
      enddo
      ! - - - - - read title lines of second long table:
      read(1,'(a)',end=119,err=118) czeilong(3)
      read(1,'(a)',end=119,err=118) czeilong(4)
      goto 111
      ! - - - - - error exit, lsteps does not match defined stepsize:
  110 continue    
      lsteps = is - 1
      lsteps = 2
      if ( lsteps .le. 2 ) then
         write(*,'(a)') czeilong(1)(1:86)
         write(*,'(a)') czeilong(2)(1:54)
         do  jl=1,3 
         write(*,'(f6.1,1p,3e12.5)') qxdistr(jl),(qdistrb(i,jl),i=1,3)
         enddo
         write(*,'('' ERROR:''
     +      '' mismatching number of steps and stepsize.'')')
         goto 119
      endif 
      ! - - - - - insert correct step quantities:
      write(czeilong(1)(31:35),'(i5)') lsteps
      write(czeilong(1)(54:59),'(f6.1)') grstep
      ! - - - - - create 3rd title line of second long table:
      czeilong(3)( 1:14) = czeilong(1)(1:14)
      czeilong(3)(15:28) = 'ENERGY DEPOSIT'
      czeilong(3)(29:88) = czeilong(1)(27:86)
      ! - - - - - read 4th title line of second long table:
      read(1,'(a)',end=119,err=118) czeilong(4)
  111 continue
      ! - - - - - read numbers of the second long table:
      do  is=1,lsteps
         read(1,*) (qzeile(i),i=0,9) 
         do  il=1,9
            qenergy(il,is) = qzeile(il)
         enddo 
         qxenerg(is) = qzeile(0)
      enddo
      close(unit=1)
      write(*,'(4x,''sumlongifiles.out'',i6.6,4x,''lsteps ='',i4,4x,
     +   ''grstep ='',f6.1)') lrunnr,lsteps,grstep

c - - - - - - read following long files and sum content - - - - - - - -
      if ( long .ge. 2 ) then
      do  ifile=2,long
         if ( mod(ifile,100) .eq. 0 ) write(*,*) '       ifile =',ifile
         open(unit=1,file=clongname(ifile)(1:longlen(ifile)),
     +               form='formatted',status='old')
         read(1,'(a)',end=119,err=118) czeile
         read(1,'(a)',end=119,err=118) czeile
         do  is=1,lsteps
            read(1,*,end=114,err=114) (qzeile(il),il=0,9)
            do  il=1,9
               qdistrb(il,is) = qdistrb(il,is) + qzeile(il)
            enddo
         enddo
         read(1,'(a)',end=119,err=118) czeile
         read(1,'(a)',end=119,err=118) czeile
         do  is=1,lsteps
            read(1,*,end=114,err=114) (qzeile(i),i=0,9)
            do  il=1,9
               qenergy(il,is) = qenergy(il,is) + qzeile(il)
            enddo
         enddo
  114    continue
         close(unit=1)
      enddo
      if ( mod(ifile-1,100) .ne. 0 ) write(*,*) '       ifile =',ifile-1
      endif

c - - - - - - write summed tables  - - - - - - - - - - - - - - - - - - -
      open(unit=8,file=chlongsum,form='formatted',status='unknown')
      write(8,'(a)') czeilong(1)
      write(8,'(a)') czeilong(2)
      do  is=1,lsteps
         write(8,'(f6.1,1p,9e12.5)') qxdistr(is),(qdistrb(il,is),il=1,9)
      enddo
      write(8,'(a)') czeilong(3)
      write(8,'(a)') czeilong(4)
      engysum = 0.d0
      do  is=1,lsteps
         engysum = engysum + qenergy(9,is) 
         write(8,'(f6.1,1p,9e12.5)') qxenerg(is),(qenergy(il,is),il=1,9)
      enddo
      close(unit=8)

c - - - - - - calculate energy sum of all levels - - - - - - - - - - - -
      if ( long .ge. 2 ) then
         if ( parengy .gt. 1. ) then
            write(*,*) '     parengy =',parengy,' GeV'
            write(*,*) '     engysum =',engysum,' * ',parengy/engysum
         endif
      endif

c - - - - - - end-of program sumlongifiles.
  118 continue
  119 continue
      stop
      end
