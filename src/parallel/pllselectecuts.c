
/* pllselectecuts.c:
**************************************************************************
            gcc -fbounds-check pllselectecuts.c -o pllselectecuts -lm
            ./pllselectecuts  Energy  [nproc]
*************************************************************************/
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <math.h>

void energyratio( double Engy, double Rati, double *Erat);

/* main program to display parallel corsika submit parameters:
*************************************************************************/

int main(int argc, char *argv[])
{
   double Energy=7.50e7, Eratio, Ratio, ratio, elog10, efiles;
   double qselect[22][8]={ 
          {  0.00,     0.,    0.,   0.,    0.,   0.,   0.,  0.},
          { 16.00,    16.,   60.,  24.,   45.,   8., 120.,  0.},
          { 16.25,    16.,   90.,   0.,    0.,   0.,   0.,  0.},
          { 16.50,    16.,  120.,  24.,   80.,  32.,  60.,  0.},
          { 16.75,    24.,  120.,   0.,    0.,   0.,   0.,  0.},
          { 17.00,    24.,  150.,  32.,  110.,  48.,  80.,  0.},
          { 17.25,    32.,  150.,   0.,    0.,   0.,   0.,  0.},
          { 17.50,    32.,  360.,  64.,  200., 120., 120.,  0.},
          { 17.75,    56.,  360.,   0.,    0.,   0.,   0.,  0.},
          { 18.00,    56.,  400.,  80.,  250., 120., 160.,  0.},
          { 18.25,    80.,  500.,   0.,    0.,   0.,   0.,  0.},
          { 18.50,   120.,  800., 160.,  700., 240., 400.,  0.},
          { 18.75,   144., 1400.,   0.,    0.,   0.,   0.,  0.},
          { 19.00,   160., 1700., 240., 1300., 320., 900.,  0.},
          { 19.25,   240., 2000.,   0.,    0.,   0.,   0.,  0.},
          { 19.50,   400., 2700., 360., 3000.,   0.,   0.,  0.},
          { 19.75,   480., 3200.,   0.,    0.,   0.,   0.,  0.},
          { 20.00,   800., 3500.,   0.,    0.,   0.,   0.,  0.},
          { 20.25,   880., 4000.,   0.,    0.,   0.,   0.,  0.},
          { 20.50,   920., 4300.,   0.,    0.,   0.,   0.,  0.},
          { 20.75,   960., 4310.,   0.,    0.,   0.,   0.,  0.},
          { 21.00,   992., 4320.,   0.,    0.,   0.,   0.,  0.} };
   int ieng, icnt, iprc, irat, ifiles, mfiles, mprocs=16;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

   if ( 1 <= argc && argc <= 3 ) {
   if ( argc >= 2 ) Energy = atof(argv[1]);
   if ( argc == 3 ) mprocs = atoi(argv[2]);

   // for( ieng=2; ieng<12; ieng++ ) { Energy = 15.376 + 0.5 * ieng;

   if ( Energy <  23.456789 ) Energy = pow(10.,Energy); // value to eV.
   if ( Energy >= 1.2345e12 ) Energy = 1.e-9 * Energy; // value to GeV.
   elog10 = 9.+log10(Energy); // log10 of GeV.
   // - - - - energy < 1.00e7 GeV:
   if ( Energy < 1.0e7 ) {
      printf("  \n");
      printf("          Energy = %10.4e GeV   log10(Energy/eV) = %7.4f \n",
         Energy,elog10);
      printf("               ==> run shower simulation as single"
             " sequential job: \n");
      printf("                   ./corsika73756_..... < inp310023"
             " > DAT310023.lst \n \n");
   }
   // - - - - energy >= 1.00e7 GeV:
   else {
      printf("  \n");
      printf("  log10(E)    Energy    nproc    ectcut      ectmax"
             "    cpumin  nfiles   Ratio \n \n");
      // - - - - - - one parameter given (energy): 
      icnt = (0.001234 + elog10 - 16. + 2.*0.25) / 0.25;
      if ( elog10 > 18.4987 ) {
         /* * * * * calculate submit quantitiy of energy[icnt-2] */
         iprc = qselect[icnt-2][1];
         Eratio = 2.4567/4.5678 * Energy / iprc;
         ratio = Energy / Eratio;
         energyratio( Energy, ratio, &Eratio );
         if ( elog10 > 19.5 ) Eratio = Eratio - 1.4e7; 
         if ( 2.4567*Energy/Eratio < 7000. ) 
           printf(
             " %8.4f   %10.4e %5d. %9.1f  %11.3e  %6.0f. %6.0f. %7.1f \n",
               elog10, Energy, iprc, Eratio*1.e-3, Eratio,
               qselect[icnt][2], 2.4567*Energy/Eratio, Energy/Eratio);
         /* * * * * calculate submit quantities of energy[icnt-1] */
         iprc = qselect[icnt-1][1];
         ifiles = 2.4567*Energy/Eratio; 
         Eratio = 2.4567/4.5678 * Energy / iprc;
         ratio = Energy / Eratio;
         energyratio( Energy, ratio, &Eratio );
         if ( elog10 > 19.5 ) Eratio = Eratio - 1.4e7; 
         mfiles = 2.4567*Energy/Eratio;
         if ( 2.4567*Energy/Eratio < 7000. && mfiles != ifiles ) 
           printf(
             " %8.4f   %10.4e %5d. %9.1f  %11.3e  %6.0f. %6.0f. %7.1f \n",
               elog10, Energy, iprc, Eratio*1.e-3, Eratio,
               qselect[icnt][2], 2.4567*Energy/Eratio, Energy/Eratio);
      }
      /* * * * * calculate submit quantities of original energy */
      iprc = qselect[icnt][1];
      for( irat=-1; irat<=1; irat++ ) {  
         Eratio = 2.4567/4.5678 * Energy / iprc;
         ratio = Energy / Eratio * (1.4567+0.4567*irat);
         energyratio( Energy, ratio, &Eratio );
         if ( elog10 > 19.5 ) Eratio = Eratio - 1.4e7; 
         ifiles =  2.4567*Energy/Eratio; 
         if ( 2.4567*Energy/Eratio < 7000. && mfiles != ifiles ) 
           printf(
             " %8.4f   %10.4e %5d. %9.1f  %11.3e  %6.0f. %6.0f. %7.1f \n",
               elog10, Energy, iprc, Eratio*1.e-3, Eratio,
               qselect[icnt][2], 2.4567*Energy/Eratio, Energy/Eratio);
         mfiles = ifiles;
      }
   printf("\n");
   printf("  log10(E)    Energy    nproc    ectcut      ectmax"
          "    cpumin  nfiles   Ratio \n \n");
   }

   // } // end-of optional for loop for testing selections.

   }

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

   return 0;
}

/* energyratio.c: function to calculate ratio=energy/ectmax
**************************************************************************
      double qvalue[13]={1.0,1.5,2.0,2.5,3.0,3.5,4.,5.,6.,7.,8.,9.,10.};*/

void energyratio( double Energy, double ratio, double *Eratio) {

   double qlog10[13]={0.0, 0.17609125905568124, 0.30102999566398120,     
      0.39794000867203760, 0.47712125471966244, 0.54406804435027567,     
      0.60205999132796240, 0.69897000433601886, 0.77815125038364363,     
      0.84509804001425681, 0.90308998699194354, 0.95424250943932487}; 
   double eratlog, fratlog, dlog;

   int i, iratlog, ilog;

      eratlog = log10( Energy / ratio );
      iratlog = eratlog;
      fratlog = eratlog - iratlog;
      dlog = 0.5 * qlog10[1];
      for( i=1; i<=12; i++ ) {
         if ( fabs(qlog10[i]-fratlog) < dlog ) {
            dlog = fabs(qlog10[i]-fratlog);
            ilog = i;
         } 
      }
      eratlog = qlog10[ilog] + iratlog;

   *Eratio = pow(10.,eratlog);

}

/****************************************************************************
  log10(E)    Energy    nproc    ectcut      ectmax    cpumin  nfiles   Ratio
 
  16.0000   1.0000e+07    16.     350.0    3.500e+05      90.     70.    28.6
  16.0000   1.0000e+07    16.     250.0    2.500e+05      90.     98.    40.0
  16.0000   1.0000e+07    16.     200.0    2.000e+05      90.    123.    50.0
 
  16.5000   3.1623e+07    24.     700.0    7.000e+05     120.    111.    45.2
  16.5000   3.1623e+07    24.     500.0    5.000e+05     120.    155.    63.2
  16.5000   3.1623e+07    24.     350.0    3.500e+05     120.    222.    90.4
 
  17.0000   1.0000e+08    32.    1500.0    1.500e+06     150.    164.    66.7
  17.0000   1.0000e+08    32.    1000.0    1.000e+06     150.    246.   100.0
  17.0000   1.0000e+08    32.     900.0    9.000e+05     150.    273.   111.1
 
  17.5000   3.1623e+08    56.    3000.0    3.000e+06     360.    259.   105.4
  17.5000   3.1623e+08    56.    2000.0    2.000e+06     360.    388.   158.1
  17.5000   3.1623e+08    56.    1500.0    1.500e+06     360.    518.   210.8
 
  18.0000   1.0000e+09    80.    7000.0    7.000e+06     500.    351.   142.9
  18.0000   1.0000e+09    80.    5000.0    5.000e+06     500.    491.   200.0
  18.0000   1.0000e+09    80.    3500.0    3.500e+06     500.    702.   285.7

  18.5000   3.1623e+09    80.   20000.0    2.000e+07    1400.    388.   158.1
  18.5000   3.1623e+09   120.   15000.0    1.500e+07    1400.    518.   210.8
  18.5000   3.1623e+09   144.   10000.0    1.000e+07    1400.    777.   316.2

  19.0000   1.0000e+10   144.   35000.0    3.500e+07    2000.    702.   285.7
  19.0000   1.0000e+10   240.   25000.0    2.500e+07    2000.    983.   400.0
  19.0000   1.0000e+10   240.   15000.0    1.500e+07    2000.   1638.   666.7

  19.5000   3.1623e+10   240.   70000.0    7.000e+07    3200.   1110.   451.8
  19.5000   3.1623e+10   400.   40000.0    4.000e+07    3200.   1942.   790.6
  19.5000   3.1623e+10   480.   35000.0    3.500e+07    3200.   2220.   903.5

  log10(E)    Energy    nproc    ectcut      ectmax    cpumin  nfiles   Ratio
****************************************************************************/
