//========================================================================
// 
//  r c o r s i k a . c p p
//  =======================
//    This program is able to read corsika particle data files of 
//    "standard" corsika or "thinning" corsika simulations and also
//    of (translated) aires simulations; valid for 32bit compilation
//    and for 64bit compilation; (here without Root/C++ graphics;) 
//    running the program writes a protocol file to standard out
//    containing 
//          primary, energy, runnumb, simnrsh, #hadrs, #muons, #gammas,
//          #elecs, #nkgelecs, obslev, theta, phi, h1km, h1gr. 
// - - - - - - - - - - - - - - CompLinkRun - - - - - - - - - - - - - - - -
// CompLink:
//          g++ -fbounds-check rcorsika.cpp -o rcorsika -lm
//          g++ -fbounds-check -m64 rcorsika.cpp -o rcorsika64 -lm
// RunProgr:
//          ./rcorsika < rcorsika.i [ > rcorsika.out ]
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    input-files:
//    unit=1: current particle data file. 
//    unit=*: steering input file (rcorsika.i):
//                    1        '_total_number_of_showers_'
//                    1        '_total_number_of_files___'
//         '/data/corsdata/joe/corsika-6990/run/DAT370375'
//                    1
//    unit=*: steering input file (rcorsika.l):
//                    1        '_total_number_of_showers_'
//                    1        '_total_number_of_files___'
//         '/data/corsdata/joe/corsika-6990/run/DAT370378'
//                    1
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//         runh=211285.281   evth=217433.078
//         long=52815.2969   evte=3397.39185   rune=3301.33252
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//         standard simulation 22932 == 3.2134576383896705E-41
//         thinning simulation 26208 == 3.6725230153024806E-41
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// (1)=black, (2)=red, (3)=green, (4)=blue, (5)=yellow, (6)=pink, (7)=cyan
// (8)=green, (9)=blue, (10)=white, (12)=d'grey, (14)=m'grey, (17)=l'grey.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//   naming conventions of corsika particles: 
//    1   gamma           24   Omega           64   K* -
//    2   positron        25   anti neutron    65   anti K* 0
//    3   electron        26   anti Lambda     66   electron neutrino
//    4   neutrino        27   anti Sigma -    67   electron anti neutrino
//    5   muon +          28   anti Sigma 0    68   muon neutrino
//    6   muon -          29   anti Sigma +    69   muon anti neutrino
//    7   pion 0          30   anti Xi 0       71   eta-> 2*gam        
//    8   pion +          31   anti Xi +       72   eta-> 3*pi0
//    9   pion -          32   anti Omega      73   eta-> pi+ + pi- + pi0
//   10   Kaon 0 long     50   omega           74   eta-> pi+ + pi- + gam
//   11   Kaon            51   rho 0           201   Deuteron  
//   12   Kaon -          52   rho +           301   Tritium
//   13   neutron         53   rho -           402   alpha
//   14   proton          54   Delta ++       1206   Carbon
//   15   anti proton     55   Delta +        1407   Nitrogen
//   16   Kaon 0 short    56   Delta 0        1608   Oxygen  
//   17   eta (71..74)    57   Delta -        2713   Aluminium  
//   18   Lambda          58   anti Delta --  2814   Silicon 
//   19   Sigma +         59   anti Delta -   3216   Sulfur
//   20   Sigma 0         60   anti Delta 0   5626   Iron    
//   21   Sigma -         61   anti Delta +   9900   Cherenkov photons   
//   22   Xi 0            62   K* 0       
//   23   Xi -            63   K* +     
//  116   D 0            131   tau +           150   anti Xi c -
//  117   D +            132   tau -           151   anti Xi c 0
//  118   anti D -       133   tau neutrino    152   anti Sigma c --
//  119   anti D 0       134   anti tau neutr  153   anti Sigma c -
//  120   D s +          137   Lambda c +      154   anti Sigma c 0
//  121   anti D s -     138   Xi c +          155   anti Xi c prime -
//  122   eta c          139   Xi c 0          156   anti Xi c prime 0
//  123   D*0            140   Sigma c ++      157   anti Omega c 0
//  124   D*+            141   Sigma c +       161   Sigma c * ++
//  125   anti D*-       142   Sigma c 0       162   Sigma c * +
//  126   anti D*0       143   Xi c prime +    163   Sigma c * 0
//  127   D* s +         144   Xi c prime 0    171   anti Sigma c * --
//  128   anti D* s -    145   Omega c 0       172   anti Sigma c * -
//  130   J/psi          149   anti Lambda c-  173   anti Sigma c * 0
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <math.h>
using namespace std;

// = = = = = = = = = = = =
// = = = = = = = = = = = = function definitions = = = = = = = = = = = = = =

void analyzeEvte(); // analyze subblock EVTE.

void analyzeEvth(); // analyze subblock EVTH.

void analyzeLong(); // analyze subblock LONG.

void analyzePart(); // analyze particle data subblock.

void analyzeRune(); // analyze subblock RUNE.

void analyzeRunh(); // analyze subblock RUNH.
 
double heightcm(double g); // calculate height (cm) for a given thickness.

void pamafill(); // fill particle masses up to 59_Ni.

void printable(int m); // formatted output table of shower quantities.

void printpartic(int m); // print particles of the first subblock.

double thickgr(double h); // thickness of atmosphere depending on height.

// = = = = = = = = = = = =
// = = = = = = = = = = = = = = global constants = = = = = = = = = = = = = =

   const int nrecstd = 22940;           // "standard corsika" record size,
                                        //        pdata[0] = 6.08270e-311.
                                        //    or  pdata[0] = 3.213456e-41.
   const int nrecthi = 26216;           // "thinning corsika" record size,
                                        //        pdata[0] = 6.95166e-311.
                                        //    or  pdata[0] = 3.672523e-41.
   const int ndim = 20000000;
   const int nsubblo = 21;
   const int nsblstd = 273;
   const int nsblthi = 312;
   const int nmaxrec = 2180000;              // valid up to 50 GBytes.
   const int numbstd = nrecstd / 4;          // =  5735.
   const int numbthi = nrecthi / 4;          // =  6554.
   int ndatpar, nreclen, nsblock, lbit,
       nxpix=750, nypix=650, nprimry, lobslev, jobslev, ifil, iruntyp=1;

   float pdata[numbthi]; // to read a single corsika record thinned data.
   float rdata[numbthi]; // to keep first corsika record thinned data.
   float sdata[numbstd]; // to read a single corsika record standard data.
   float hdata[819]; // to read end of first corsika record thinned data.
   float zdata[2]; // to read record length information on 64bit machines. 

   double qrunh[275], qevth[275], qevte[275], qdata[101];

   const double velight = 29.9792458; // velocity of light in cm/nsec.

   const double AATM[6] = { 0.e0,       // param for US standard atmosph.
                   -186.5562e0, -94.919e0, 0.61289e0, 0.e0, 1.128292e-2 };
   const double BATM[6] = { 0.e0,       // param for US standard atmosph.
                1222.6562e0, 1144.9069e0, 1305.5948e0, 540.1778e0, 0.e0 };
   const double CATM[6] = { 0.e0,       // param for US standard atmosph.
              994186.38e0, 878153.55e0, 636143.04e0, 772170.16e0, 1.e-9 };
 
   const double parmas[202] = { 0.0e0,               // particle masses.
         0.e0     ,5.10998902e-4,5.10998902e-4, 0.e0     ,0.105658357e0,
         0.105658357e0,0.1349766e0,0.13957018e0,0.13957018e0,0.497672e0,
         0.493677e0, 0.493677e0, 0.93956533e0,0.93827200e0,0.93827200e0,
         0.497672e0 , 0.54730e0  , 1.115683e0 , 1.18937e0 , 1.192642e0 ,
         1.197449e0 , 1.31483e0  , 1.32131e0  , 1.67245e0 ,0.93956533e0,
         1.115683e0 , 1.18937e0  , 1.192642e0 , 1.197449e0, 1.31483e0  ,
         1.32131e0  , 1.67245e0  , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.78257e0  ,
         0.7690e0   , 0.7665e0   , 0.7665e0   , 1.2305e0  , 1.2318e0   ,
         1.2331e0   , 1.2344e0   , 1.2309e0   , 1.2323e0  , 1.2336e0   ,
         1.2349e0   , 0.89610e0  , 0.89166e0  , 0.89166e0 , 0.89610e0  ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.54730e0  , 0.54730e0  , 0.54730e0  , 0.54730e0 , 0.105658e0 ,
         0.105658e0 , 0.0e0      , 0.0e0      , 0.0e0     , 0.0e0      ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         1864.5e0   , 1869.3e0   , 1869.3e0   , 1864.5e0  , 1968.6e0   ,
         1968.5e0   , 2979.7e0   , 2006.7e0   , 2010.0e0  , 2010.0e0   ,
         2006.7e0   , 2112.4e0   , 2112.4e0   , 3510.51e0 , 3096.87e0  ,
         1776.99e0  , 1776.99e0  , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 2284.9e0   , 2466.3e0   , 2471.8e0  , 2452.6e0   ,
         2451.3e0   , 2452.2e0   , 2574.1e0   , 2578.8e0  , 2697.5e0   ,
         0.e0       , 0.e0       , 0.e0       , 2284.9e0  , 2466.3e0   ,
         2471.8e0   , 2452.6e0   , 2451.3e0   , 2452.2e0  , 2574.1e0   ,
         2578.8e0   , 2697.5e0   , 0.e0       , 0.e0      , 0.e0       ,
         2519.4e0   , 2515.9e0   , 2517.5e0   , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         2519.4e0   , 2515.9e0   , 2517.5e0   , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0, 0.e0 };

   const double chargs[202] = { 0.,                 // particle charges.
           0.  ,+1.e0,-1.e0, 0.  ,+1.e0,-1.e0, 0.  ,+1.e0,-1.e0, 0.  ,
          +1.e0,-1.e0, 0.  ,+1.e0,-1.e0, 0.  , 0.  , 0.  ,+1.e0, 0.  ,
          -1.e0, 0.  ,-1.e0,-1.e0, 0.  , 0.  ,-1.e0, 0.  ,+1.e0, 0.  ,
          +1.e0,+1.e0, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  ,+1.e0,-1.e0,+2.e0,+1.e0, 0.  ,-1.e0,-2.e0,-1.e0, 0.  ,
          +1.e0, 0.  ,+1.e0,-1.e0, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,+1.e0,-1.e0, 0.  ,+1.e0,
          -1.e0, 0.  , 0.  ,+1.e0,-1.e0, 0.  ,+1.e0,-1.e0, 0.  , 0.  ,
          +1.e0,-1.e0, 0.  , 0.  , 0.  , 0.  ,+1.e0,+1.e0, 0.  ,+2.e0,
          +1.e0, 0.  ,+1.e0, 0.  , 0.  , 0.  , 0.  , 0.  ,-1.e0,-1.e0,
           0.  ,-2.e0,-1.e0, 0.  ,-1.e0, 0.  , 0.  , 0.  , 0.  , 0.  ,
          +2.e0,+1.e0, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
          -2.e0,-1.e0, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  };

// = = = = = = = = = = = =
// = = = = = = = = = = = = = = global quantities = = = = = = = = = = = = =
  
   char corsfile[21][120];
  
   int mprimary, ibytes, mpinull, mpiplus, mpiminus;
   int nshowf[21], nfils, nshow, nshof, ishow, ishof, isho, irec, itext;
   int idpa, iobs, igen, icod, mutr, ii, la, lb, lz, lp, lthi=1;
  
   double epart, eplog, equad; // extra memory for quantities.
  
   double pmass[6001], psign[6001], prest[6001], rx, ry;
  
   double slopest, engylow, enghigh, thetlow, thehigh, obslev1;
 
// = = = = = = = = = = = =
// = = = = = = = = = = = = = = = = = main= = = = = = = = = = = = = = = = =
//          
//=== I.a - read total number of showers (nshow)
//=== I.b - read total number of files (nfils)
//=== II. - read file names and corresponding numbers of showers - - - - -
//=== III.- loop over all file names read in
//==== III.a.- - test next corsika file
//==== III.b.- - loop over all data records
//==== III.b.1.- - - read a single record from corsika file
//==== III.b.2.- - - analyze subblocks of the record
//=== IV. - after all showers and particles
//
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int main(int argc, char *argv[])
{
   char chtext[120];
   time_t tseconds;
   tm *datetime;
   pamafill();  // fill particle masses up to 59_Ni.
   if ( argc > 1 ) iruntyp = atoi(argv[1]);

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//=== I.a - read total number of showers (nshow) - - - - - - - - - - - - -
   cin >> nshow >> chtext;
   itext = strlen(chtext);
   chtext[0] = ' ';
   chtext[itext-1] = '\0';
   cout << "   " << endl;
   cout << chtext << setw(6) << " number of showers:" <<  nshow << endl;

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//=== I.b - read total number of files (nfils) - - - - - - - - - - - - - -
   cin >> nfils >> chtext;
   itext = strlen(chtext);
   chtext[0] = ' ';
   chtext[itext-1] = '\0';
   cout << chtext << setw(6) << " Number of files:" << nfils << endl;
   if (nfils > 20) {
      nfils = 20;
      cout << " Number of corsika files, nfils, reduced to 20,"
           << " but may be increased carefully. " << endl;
   }
   tseconds = time(&tseconds);
   datetime = localtime(&tseconds);
   printf("        _date_of_analysis_ ");
   printf(" %.2d. %.2d. %.4d.  %.2d.%.2d%.2d hms  (%ld sec) \n",
      datetime->tm_mday, datetime->tm_mon+1, datetime->tm_year+1900,
      datetime->tm_hour, datetime->tm_min, int(tseconds)%60, tseconds); 

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//=== II. - read file names and corresponding numbers of showers - - - - -
   cout << "   " << endl;
   for( ifil=1; ifil<nfils+1; ifil++ ) {
      cin >> chtext >> isho;
         itext = strlen(chtext);
         chtext[0] = ' ';
         chtext[itext-1] = '\0';
         strcpy(corsfile[ifil],&chtext[1]);
      if ( iruntyp > 1 )
         cout << " corsika.file:   " << corsfile[ifil]
              << "          " << isho << " sh. " << endl;
      if ( isho > nshow ) isho=nshow;
      nshowf[ifil] = isho;
   }

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//=== III. - loop over all file names read in- - - - - - - - - - - - - - -
   nshof = 0;
   ishow = 0;
   for( ifil=1; ifil<=nfils; ifil++ )
   {
   nshof = nshof + nshowf[ifil]; // total showers till end of current file
   ishof = 0;                    // counting showers of current file.

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//==== III.b.1.- - read/check first record from current corsika file - - -
   ifstream inputcors(corsfile[ifil],ios::in); // allocate file.
   if ( inputcors.peek() == EOF ) {
      printf("         file ==> does not exist at this path! Stop!\n");
      inputcors.clear();
      inputcors.close();
      break;
   }
   else {
      inputcors.read((char*)&sdata,sizeof(sdata));
      if ( 217433.0 < sdata[nsblstd+1] && sdata[nsblstd+1] < 217433.2 ) {
         lbit = 32;
         lthi = 0;
         for( ii=numbstd-1; ii<=numbthi; ii++ ) rdata[ii] = 0.;
         for( ii=0; ii<=numbstd; ii++ ) rdata[ii] = sdata[ii];
      }
      else if ( 217433.0 < sdata[nsblstd+2] && sdata[nsblstd+2] < 217433.2 ) {
         lbit = 64;
         lthi = 0;
         inputcors.read((char*)&zdata,sizeof(zdata));
         for( ii=numbstd-1; ii<=numbthi; ii++ ) rdata[ii] = 0.;
         for( ii=0; ii<=numbstd; ii++ ) rdata[ii] = sdata[ii+1];
      }
      else if ( 217433.0 < sdata[nsblthi+1] && sdata[nsblthi+1] < 217433.2 ) {
         lbit = 32;
         lthi = 1;
         inputcors.read((char*)&hdata,sizeof(hdata));
         for( ii=0; ii<numbstd; ii++ ) rdata[ii] = sdata[ii];
         for( ii=numbstd; ii<numbthi; ii++ ) rdata[ii] = hdata[ii-numbstd];
      }
      else if ( 217433.0 < sdata[nsblthi+2] && sdata[nsblthi+2] < 217433.2 ) {
         lbit = 64;
         lthi = 1;
         inputcors.read((char*)&hdata,sizeof(hdata));
         inputcors.read((char*)&zdata,sizeof(zdata));
         for( ii=0; ii<numbstd; ii++ ) rdata[ii] = sdata[ii+1];
         for( ii=numbstd; ii<numbthi; ii++ ) rdata[ii] = hdata[ii-numbstd+1];
      }
      lbit = (lbit-32) / 32;
      nsblock = nsblstd + lthi*39;
      ndatpar = nsblock / 39;
   }

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//==== III.b - - loop over all data records of current file  - - - - - - -
   for( irec=1; irec<nmaxrec; irec++ ) {

// _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//==== III.b.1.- - - get record of corsika particle data file: - - - - - -
   if ( irec > 1 ) {
   // - - - - - - - - read data irec>1 of known corsika record length:
      if ( lthi == 1 ) { // - - - - - - - long record (thinning).
         inputcors.read((char*)&pdata,sizeof(pdata));
         if ( pdata[0] < 1.21345763838967e-41 ) { // EOF found.
            break;
         }
         // - - - - - - - not yet at end of file:
         if ( lbit == 1 ) inputcors.read((char*)&zdata,sizeof(zdata));
         pdata[0] = 0.21345763838967e-41;
         if ( lbit == 1 ) {
            for( ii=1; ii<numbthi; ii++ ) pdata[ii] = pdata[ii+lbit];
            pdata[numbthi] = zdata[0];
         }
         for( ii=1; ii<=numbstd; ii++ ) sdata[ii] = pdata[ii];
      }
      else if ( lthi == 0 ) { // - - - - - - - short record (standard).
         inputcors.read((char*)&sdata,sizeof(sdata));
         if ( sdata[0] < 1.21345763838967e-41 ) { // EOF found.
            break;
         }
         // - - - - - - - not yet at end of file:
         if ( lbit == 1 ) inputcors.read((char*)&zdata,sizeof(zdata));
         sdata[0] = 0.21345763838967e-41;
         for( ii=numbstd; ii<=numbthi; ii++ ) pdata[ii] = 0.;
         for( ii=1; ii<=numbstd; ii++ ) pdata[ii] = sdata[ii+lbit];
      }
   }
   // - - - - - - - - take data of first record from read/check:
   else if ( irec == 1 ) {
      if ( iruntyp > 1 ) {
         printf("         nsblock = %d ",nsblock);
         if ( lbit == 1 ) printf("     64bit file");
         printf(" \n");
      }
      for( ii=0; ii<numbthi ; ii++ ) pdata[ii] = rdata[ii];
      for( lz=1; lz<=nsblstd; lz++ ) qrunh[lz] = rdata[lz];
      for( lz=1; lz<=nsblstd; lz++ ) qevth[lz] = rdata[lz+nsblock];
      // - - - - - - - - optional test print to verify contents:
      if ( lbit > 2 ) {
         printf("        rdata[0-1]   %18.15e  %18.15e \n",rdata[0],rdata[1]);
         printf("        rdata[2-6]   %f  %f  %f  %f  %f \n",rdata[2],rdata[3],
             rdata[4],rdata[5],rdata[6]);
         printf("       rdata[7-11]   %f  %f  %f  %f  %f \n",rdata[7],rdata[8],
             rdata[9],rdata[10],rdata[11]);
         if ( nsblock == 273 ) {
            printf("  rdata[5725-5729]   %f  %f  %f  %f  %f \n",rdata[5725],
             rdata[5726],rdata[5727],rdata[5728],rdata[5729]);
            printf("  rdata[5730-5734]   %f  %f  %f  %f  %f \n",rdata[5730],
             rdata[5731],rdata[5732],rdata[5733],rdata[5734]);
            printf("  rdata[5735-5736]   %18.15e  %18.15e \n",
             rdata[5735],rdata[5736]);
         }
      }
      // - - - - - - - - optional print out of quantities of first record:
      if ( lbit < 2 ) {
         if ( iruntyp > 1 ) 
         for( lp=1; lp<nsblock*4; lp+=nsblock )
         printf(" subl%.2d:%14.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e"//
            " %13.6e\n",int((lp+nsblock)/nsblock),rdata[lp],rdata[lp+1],
            rdata[lp+2],rdata[lp+3],rdata[lp+4],rdata[lp+5],rdata[lp+6],
            rdata[lp+7]);
      }
   }

// _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//==== III.b.2.- - - analyze subblocks of the record - - - - - - - - - - -
      for( lp=1; lp<nsblock*21; lp+=nsblock )
      {

//========== RUNH subblock found =========================================
      if ( ( 211285.1 < pdata[lp] && pdata[lp] < 211285.4 ) ||
         ( 1111110. < pdata[lp] && pdata[lp] < 1111112. ) ) {
         if ( ifil < 0 ) cout << "        RUNH subblock found. " << endl;
         analyzeRunh(); 
      }

//========== EVTH subblock found =========================================
      else if ( ( 217433.0 < pdata[lp] && pdata[lp] < 217433.2 ) ||
                ( 3333332. < pdata[lp] && pdata[lp] < 3333334. ) ) {
         if ( ifil < 0 ) cout << "        EVTH subblock found. " << endl;      
         analyzeEvth();
      }

//========== LONG subblock found =========================================
      else if ( ( 52815.2 < pdata[lp] && pdata[lp] < 52815.3 ) ||
               ( 5555554. < pdata[lp] && pdata[lp] < 5555556. ) ) {
         if ( ifil < 0 ) cout << "        LONG subblock ignored." << endl;
         analyzeLong();
      }

//========== particle data subblock ======================================
      else if ( ( 1000.0 < pdata[lp] && pdata[lp] < 3301.3 ) ||
                ( 3301.4 < pdata[lp] && pdata[lp] < 3397.3 ) ||
                ( 3397.4 < pdata[lp] && pdata[lp] < 52815.2 ) ||
               ( 52815.3 < pdata[lp] && pdata[lp] < 77000.0 ) ||
              ( 116000.0 < pdata[lp] && pdata[lp] < 174000.0 )) {
         analyzePart();
      }

//========== EVTE subblock found =========================================
      else if ( ( 3397.3 < pdata[lp] && pdata[lp] < 3397.5 ) ||
              ( 7777776. < pdata[lp] && pdata[lp] < 7777778. ) ) {
         if ( ifil < 0 ) cout << "        EVTE subblock found. " << endl;
         analyzeEvte();
         if ( ishow == nshow ) break; // all showers analyzed. 
      }  

//========== RUNE subblock found =========================================
      else if ( ( 3301.2 < pdata[lp] && pdata[lp] < 3301.5 ) ||
              ( 9999990. < pdata[lp] && pdata[lp] < 1.e7 ) ) {
         if (ifil < 0 ) cout << "        RUNE subblock found. " << endl;
         analyzeRune();
      }

      } // end_of_loop of analyzing subblocks - - - - - - - - - - - - - - -

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//==== III.c.- - - continue with next data record or close file. 
 
      if ( ishow == nshow ) break; // all showers analyzed.

    } // endfor_loop over all data records - - - - - - - - - - - - - - - -

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

    inputcors.close(); // close corsika file also after break at EOF;
 
    if ( ishow == nshow ) break; // all showers analyzed.

   } // endfor_loop over all file names read in  - - - - - - - - - - - - -

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//=== IV. - after all showers and particles- - - - - - - - - - - - - - - -
  
   cout << "   " << endl;
   cout << " End of particle data file:"; 
   if ( ishow == 1 )
      cout << " 1 shower processed.";                
   else if ( ishow > 1 )
      cout << " " << ishow << " showers processed.";
   else 
      cout << " WARNING, newer " << (lbit+1)*32
           << "bit simulation, may be better using `rcorsika"
           << (lbit+1)*32 <<".cpp`."; 
   cout << endl << "  " << endl;
   
   return 0;
}


// = = = = = = = = = = = =
// = = = = = = = = = = = = = = = analyzeEvte = = = = = = = = = = = = = = =
//        analyze or check corsika data subblock event end
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void analyzeEvte()
{
   ishow++;
   if ( ishow <= nshow ) {
      // - - test-print of some quantities of evte subblock - -
      cout.precision(5);
      cout.setf(ios::scientific,ios::fixed);
      if ( iruntyp > 1 )
      cout << "                                 EVTE   hadr= "
           << pdata[lp+4] << "   allp= " << pdata[lp+6]
           << "   irec= " << irec << endl;
   }
   if (pdata[lp+4]+pdata[lp+5] > 0.) {
      // - - fill electron numbers also at NKG simulations - -
      if ( pdata[lp+184] < 1. ) {
         la = 1;
         pdata[lp+184] = pdata[lp+184-la];
         while( pdata[lp+184] < 1. )
            { la++; pdata[lp+184] = pdata[lp+184-la]; } 
      }
      if ( pdata[lp+3] < 1. ) pdata[lp+3] = pdata[lp+184];
      if ( pdata[lp+194] < 0.25 ) {
         la = 1;
         pdata[lp+194] = pdata[lp+194-la];
         while( pdata[lp+194] < 0.25)
            { la++; pdata[lp+194] = pdata[lp+194-la]; } 
      }
   }
//- - - - - - - - - - number of hadrons at obslev - - - - - - - - - - - -
   qdata[5] = pdata[lp+4]; 
//- - - - - - - - - - number of muons at obslev - - - - - - - - - - - - -
   qdata[6] = pdata[lp+5]; 
//- - - - - - - - - - number of photons at obslev - - - - - - - - - - - -
   qdata[7] = pdata[lp+2]; 
//- - - - - - - - - - number of electrons at obslev - - - - - - - - - - -
   qdata[8] = pdata[lp+3]; 
//- - - - - - - - - - NKG number of electrons at obslev - - - - - - - - -
   qdata[9] = pdata[lp+184]; 
//- - - - - - - - - - shower age at obslev- - - - - - - - - - - - - - - -
   qdata[10] = pdata[lp+194]; 
   printable(ishow);           // print shower information to file.
   return;
}


// = = = = = = = = = = = =
// = = = = = = = = = = = = = = = analyzeEvth = = = = = = = = = = = = = = = =
//        analyze or check corsika data subblock event header
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void analyzeEvth()
{
//- - - - - - - - - - number of shower in current file - - - - - - - - - -
   ishof++;
//- - - - - - - - - - set counters to zero - - - - - - - - - - - - - - - -
   mutr = 0; 
   mpinull = 0;
   mpiplus = 0;
   mpiminus = 0;
//- - - - - - - - - - number of observation level - - - - - - - - - - - -
   lobslev = (int) pdata[lp+46];
//- - - - - - - - - - primary particle code - - - - - - - - - - - - - - -
   qdata[1] = pdata[lp+2];
//- - - - - - - - - - primary particle energy in GeV  - - - - - - - - - -
   qdata[2] = pdata[lp+3];
//- - - - - - - - - - - run number  - - - - - - - - - - - - - - - - - - -
   qdata[3] = pdata[lp+43];
//- - - - - - - - - - simulated shower number - - - - - - - - - - - - - -
   qdata[4] = pdata[lp+1];
//- - - - - - - - - - copy some parameters and all observation levels - -
   for( la=44; la<61; la++ ) { qdata[la] = pdata[lp+la-1]; }
//- - - - - - - - - - - observation level (meter) - - - - - - - - - - - -
   qdata[11] = 1.e-2 * pdata[lp+47];
//- - - - - - - - - - - theta in degrees- - - - - - - - - - - - - - - - -
   qdata[12] = pdata[lp+10] * 57.295779513e0;
//- - - - - - - - - - - phi in degrees- - - - - - - - - - - - - - - - - -
   qdata[13] = pdata[lp+11] * 57.295779513e0;
//- - - - - - - - - - height of first interaction (km)  - - - - - - - - -
   if ( pdata[lp+6] < 0. ) pdata[lp+6] = fabs(pdata[lp+6]);
   qdata[14] = 1.e-5 * pdata[lp+6];
//- - - - - - - - - - height of first interaction (gramm/cm2) - - - - - -
   qdata[15] = thickgr(pdata[lp+6]);
//- - - - - - - - - - runtime of light from first intact to obslev- - - -
   qdata[16] = ( pdata[lp+6] - pdata[lp+46+lobslev] )
             / ( velight * cos(pdata[lp+10]) );
   cout.precision(5);
   cout.setf(ios::scientific,ios::fixed);
   if ( iruntyp > 1 )
   cout << " ish= " << ishow+1 << "   run= " << int(pdata[lp+43])
        << "   evt= " << pdata[lp+1] << "   EVTH"
        << "   prim= " << pdata[lp+2] << "   engy= " << pdata[lp+3]
        << "   thet= " << pdata[lp+10] * 57.295779513e0
        << "   phi= " << pdata[lp+11] * 57.295779513e0 << endl;
   return;
}


// = = = = = = = = = = = =
// = = = = = = = = = = = = = = = analyzeLong = = = = = = = = = = = = = = = =
//        analyze or check corsika data subblock of longitudinal tables
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void analyzeLong()
{
   // analyze longitudinal tables of 
   // vertical depth of step n, number of gammas, positrons,
   // electrons, muons(+), muons(-), hadrons, all charged,
   // nuclei, and cherenkov photons at step n, test fit parameters.    
   // see separate program code `rcorsik1.cpp`.
   // cout << "        LONG subblock ignored. " << endl;
   return;
}


// = = = = = = = = = = = =
// = = = = = = = = = = = = = = = analyzePart = = = = = = = = = = = = = = = =
//        analyze or check corsika particle data subblock
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void analyzePart()
{
   for( lz=lp; lz<lp+nsblock; lz+=ndatpar ) {
 
      // - - - - - - some print out of first particles:
      if ( irec == 1 && ishow == 0 && ifil == 1 ) {
         if ( lz < 3*nsblock ) {
            // if ( lz == lp ) cerr << "   " << endl;
            // if ( pdata[lz] < 77000. ) printpartic(lz);
            // if ( lz == lp+nsblock-ndatpar ) cerr << "   " << endl;
         }
      } 
 
      // - - - - - - valid particle id found: - - - - - - - - - - - - - -
      if ( pdata[lz] > 0. ) {
         icod = (int) pdata[lz]; // - integer copy of `full` particle id.
         iobs = icod % 10;    
         if ( iobs == 0 ) iobs = 10; 
         igen = (icod % 1000) / 10; // - generation of interaction.
         idpa = icod / 1000; // - particle id in corsika notation.

         // - - - - - - - - calculate total energy also for nuclei:
         if ( idpa < 200 ) { 
            equad = parmas[idpa]*parmas[idpa] + pdata[lz+1]*pdata[lz+1]
                    + pdata[lz+2]*pdata[lz+2] + pdata[lz+3]*pdata[lz+3];
         }
         else if ( idpa == 9900 ) {
            equad = parmas[1]*parmas[1] + pdata[lz+1]*pdata[lz+1]
                    + pdata[lz+2]*pdata[lz+2] + pdata[lz+3]*pdata[lz+3];
         }
         else if ( idpa < 6000 ) {
            la = idpa%100; // number of protons in nucleus.
            lb = int(double(idpa/100)); // number of neutrons in nucleus.
            equad = parmas[14]*double(la) + parmas[13]*double(lb-la);
            equad = equad*equad + pdata[lz+1]*pdata[lz+1]
                  + pdata[lz+2]*pdata[lz+2] + pdata[lz+3]*pdata[lz+3];
         }
         epart = sqrt(equad);
         eplog = log10(epart);

         // - - - - - - - - switch coordinates to meters instead of cm. 
         rx = 1.e-2 * pdata[lz+4];
         ry = 1.e-2 * pdata[lz+5];
 
         // =========== muons =========== 
         if ( idpa == 5 || idpa == 6 ) {
            // - - - - - - - - - - - count truncated muons:
            if (1600. <= rx*rx+ry*ry && rx*rx+ry*ry <= 40000.) mutr++;
         } 

         // ========== hadrons ========== 
         if ( 7 <= idpa && idpa <= 65 ) { 
  
         // ========== pion(0) ========== 
            if ( idpa == 7 ) {
                mpinull++;
            } 
  
         // ========== pion(+) ========== 
            if ( idpa == 8 )  {
               mpiplus++;
            } 
  
         // ========== pion(-) ========== 
            if ( idpa == 9 ) {
               mpiminus++;
            } 
  
         } // end_of_hadrons. 
 
      }
      // - - - - - - invalid particle id found: - - - - - - - - - - - - -
      if ( pdata[lz] < 0. ) {
         // cerr << "   " << "invalid particle id found!" << "   " << endl; 
      }
 
   }
   
   return;
}


// = = = = = = = = = = = =
// = = = = = = = = = = = = = = = analyzeRune = = = = = = = = = = = = = = = =
//        analyze or check corsika data subblock run end
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void analyzeRune()
{
   // cout << "      RUNE " << endl;
   return;
}


// = = = = = = = = = = = =
// = = = = = = = = = = = = = = = analyzeRunh = = = = = = = = = = = = = = = =
//        analyze or check corsika data subblock run header
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void analyzeRunh()
{
   // cout << "      RUNH " << endl;
   return;
}


// = = = = = = = = = = = =
// = = = = = = = = = = = = = = = heightcm  = = = = = = = = = = = = = = = =
//     calculate height (cm) a.s.l. for a given thickness (gramms/cm**2)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double heightcm(double g)
{
   double height;
   if (g > 631.1e0)
      { height = CATM[1] * log( BATM[1] / (g-AATM[1]) ); }
   else if (g > 271.7e0)
      { height = CATM[2] * log( BATM[2] / (g-AATM[2]) ); }
   else if (g > 3.0395e0)
      { height = CATM[3] * log( BATM[3] / (g-AATM[3]) ); }
   else if (g > 1.28292e-3)
      { height = CATM[4] * log( BATM[4] / (g-AATM[4]) ); }
   else
      { height = (AATM[5] - g) / CATM[5]; }
   return height;
}


// = = = = = = = = = = = =
// = = = = = = = = = = = = = = = pamafill= = = = = = = = = = = = = = = = =
//                    fill particle masses up to 59_Ni.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void pamafill()
{
   double az,eb;
   int ia,in,ip,is;

// - - - - - - initialize particle masses, charges, names:
   const int ipartyp[202] = {
                 0,1,2,3,0,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
           21,22,23,24,25,26,27,28,29,30,31,32,17*0,50,51,52,53,54,55,
         56,57,58,59,60,61,62,63,64,65,66,67,68,69,0,71,72,73,74,41*0,
        116,117,118,119,120,121,122,123,124,125,126,127,128,0,130,131,
         132,133,134,2*0,137,138,139,140,141,142,143,144,145,3*0,149,
           150,151,152,153,154,155,156,157,3*0,161,162,163,7*0,171,
                 172,173,25*0,199,0,0 };
   const double parmas[202] = { 0.0e0,               // particle masses.
         0.e0     ,5.10998902e-4,5.10998902e-4, 0.e0     ,0.105658357e0,
         0.105658357e0,0.1349766e0,0.13957018e0,0.13957018e0,0.497672e0,
         0.493677e0, 0.493677e0, 0.93956533e0,0.93827200e0,0.93827200e0,
         0.497672e0 , 0.54730e0  , 1.115683e0 , 1.18937e0 , 1.192642e0 ,
         1.197449e0 , 1.31483e0  , 1.32131e0  , 1.67245e0 ,0.93956533e0,
         1.115683e0 , 1.18937e0  , 1.192642e0 , 1.197449e0, 1.31483e0  ,
         1.32131e0  , 1.67245e0  , 0.e0       ,  16*0.e0  , 0.78257e0  ,
         0.7690e0   , 0.7665e0   , 0.7665e0   , 1.2305e0  , 1.2318e0   ,
         1.2331e0   , 1.2344e0   , 1.2309e0   , 1.2323e0  , 1.2336e0   ,
         1.2349e0   , 0.89610e0  , 0.89166e0  , 0.89166e0 , 0.89610e0  ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.54730e0  , 0.54730e0  , 0.54730e0  , 0.54730e0 , 0.105658e0 ,
         0.105658e0 , 0.0e0      , 0.0e0      , 0.0e0     ,  36*0.0e0  ,
         1864.5e0   , 1869.3e0   , 1869.3e0   , 1864.5e0  , 1968.6e0   ,
         1968.5e0   , 2979.7e0   , 2006.7e0   , 2010.0e0  , 2010.0e0   ,
         2006.7e0   , 2112.4e0   , 2112.4e0   , 3510.51e0 , 3096.87e0  ,
         1776.99e0  , 1776.99e0  , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 2284.9e0   , 2466.3e0   , 2471.8e0  , 2452.6e0   ,
         2451.3e0   , 2452.2e0   , 2574.1e0   , 2578.8e0  , 2697.5e0   ,
         0.e0       , 0.e0       , 0.e0       , 2284.9e0  , 2466.3e0   ,
         2471.8e0   , 2452.6e0   , 2451.3e0   , 2452.2e0  , 2574.1e0   ,
         2578.8e0   , 2697.5e0   , 0.e0       , 0.e0      ,     0.e0   ,
         2519.4e0   , 2515.9e0   , 2517.5e0   , 0.e0      ,   6*0.e0   ,
         2519.4e0   , 2515.9e0   , 2517.5e0   , 0.e0      ,  27*0.e0  }; 
   const double chargs[202] = { 0.,                 // particle charges.
           0.  ,+1.e0,-1.e0, 0.  ,+1.e0,-1.e0, 0.  ,+1.e0,-1.e0, 0.  ,
          +1.e0,-1.e0, 0.  ,+1.e0,-1.e0, 0.  , 0.  , 0.  ,+1.e0, 0.  ,
          -1.e0, 0.  ,-1.e0,-1.e0, 0.  , 0.  ,-1.e0, 0.  ,+1.e0, 0.  ,
          +1.e0,+1.e0, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,11*0.,
           0.  ,+1.e0,-1.e0,+2.e0,+1.e0, 0.  ,-1.e0,-2.e0,-1.e0, 0.  ,
          +1.e0, 0.  ,+1.e0,-1.e0, 0.  , 0.  , 0.  , 0.  , 0.  ,41*0.,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,+1.e0,-1.e0, 0.  ,+1.e0,
          -1.e0, 0.  , 0.  ,+1.e0,-1.e0, 0.  ,+1.e0,-1.e0, 0.  , 0.  ,
          +1.e0,-1.e0, 0.  , 0.  , 0.  , 0.  ,+1.e0,+1.e0, 0.  ,+2.e0,
          +1.e0, 0.  ,+1.e0, 0.  , 0.  , 0.  , 0.  , 0.  ,-1.e0,-1.e0,
           0.  ,-2.e0,-1.e0, 0.  ,-1.e0, 0.  , 0.  , 0.  , 0.  , 0.  ,
          +2.e0,+1.e0, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
          -2.e0,-1.e0, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,22*0.}; 
   char chptext[202][20] = {"                   ", 
      " gamma             "," positron          "," electron          ",
      "                   "," muon +            "," muon -            ",
      " pion 0            "," pion +            "," pion -            ",
      " Kaon 0 long       "," Kaon +            "," Kaon -            ",
      " neutron           "," proton            "," anti proton       ",
      " Kaon 0 short      ","                   "," Lambda            ",
      " Sigma +           "," Sigma 0           "," Sigma -           ",
      " Xi 0              ","  Xi -             "," Omega -           ",
      " anti neutron      "," anti Lambda       "," anti Sigma -      ",
      " anti Sigma 0      "," anti Sigma +      "," anti Xi 0         ",
                            // 0 to 30.
      " anti Xi +         "," anti Omega +      ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   "," omega             "," rho 0             ",
      " rho +             "," rho -             "," Delta ++          ",
      " Delta +           "," Delta 0           "," Delta -           ",
      " anti Delta --     "," anti Delta -      "," anti Delta 0      ",
                            // 31 to 60.
      " anti Delta +      "," Kaon * 0          "," Kaon * +          ",
      " Kaon * -          "," anti Kaon * 0     "," electron neutrino ",
      " anti elec neutrino"," muon neutrino     "," anti muon neutrino",
      "                   "," eta=>2*gamma      "," eta=>3*pi0        ",
      " eta=>pi+ pi- pi0  "," eta=>pi+ pi- gamma","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
                            // 61 to 90.
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   "," D 0               "," D +               ",
      " anti D -          "," anti D 0          "," D s +             ",
                            // 91 to 120.
      " anti D s -        "," eta c             "," D * 0             ",
      " D * +             "," anti D * -        "," anti D * 0        ",
      " D * s +           "," anti D * s -      ","                   ",
      " J/psi             "," tau +             "," tau -             ",
      " tau neutrino      "," anti tau neutrino ","                   ",
      "                   "," Lambda c +        "," Xi c +            ",
      " Xi c 0            "," Sigma c ++        "," Sigma c +         ",
      " Sigma c 0         "," Xi c prime +      "," Xi c prime 0      ",
      " Omega c 0         ","                   ","                   ",
      "                   "," anti Lambda c -   "," anti Xi c -       ",
                            // 121 to 150.
      " anti Xi c 0       "," anti Sigma c --   "," anti Sigma c -    ",
      " anti Sigma c 0    "," anti Xi c prime - "," anti Xi c prime 0 ",
      " anti Omega c 0    ","                   ","                   ",
      "                   "," Sigma c * ++      "," Sigma c * +       ",
      " Sigma c * 0       ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   "," anti Sigma c * -- ",
      " anti Sigma c * -  "," anti Sigma c * 0  ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
                            // 151 to 180.
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      " Cherenkov photon  ","                   ","                   "};

// - - - - - - calculate particle (nuclei) masses:   
   for( ia=0; ia<=76; ia++ )
      { pmass[ia] = parmas[ia]; psign[ia] = chargs[ia]; }
   for( ia=77; ia<=6000; ia++ )
      { pmass[ia] = 0.; psign[ia] = 0.; }
   for( ia=1; ia<60; ia++ ) {
      for( ip=1; ip<=ia; ip++ ) {
         in = ia - ip;
         is = ia * 100 + ip;
         // - - without binding energy effects.
         pmass[is] = parmas[13] * in + parmas[14] * ip;
         // - - nuclei are assumed to be fully ionized.
         psign[is] = ip;
         // - - applying binding energy of nuclei.
         az = ia;
         eb = 14.1 * az - 0.595 * ip*ip * pow(az,-1./3.)
            - 13. * pow(az,2./3.) - 19. * (ip-in)*(ip-in) / az;
         if ( ip%2 == 0 && in%2 == 0 )
            { eb = eb + 33.5 * pow(az,-0.75); }
         else
            if ( ip%2 == 1 && in%2 == 1 )
            { eb = eb - 33.5 * pow(az,-0.75); }
         if ( eb > 0. ) eb=1.e-3*eb;
         else eb=0.; 
         pmass[is] = pmass[is]-eb;
      }
   }
   // masses of multi-neutron clusters.
   for( in=1; in<60; in++ ) {
      is = in * 100;
      pmass[is] = parmas[13] * in;
      psign[is] = 0.;
   }
}


// = = = = = = = = = = = =
// = = = = = = = = = = = = = = = printable = = = = = = = = = = = = = = = =
//        formatted output table of shower quantities (to txt-file):
//    primary, energy, runnumb, simnrsh, #hadrs, #muons, #gammas, #elecs,
//            #nkgelecs, obslev, theta, phi, h1km, h1gr. 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void printable(int ish)
{
   if ( ish == 1 ) cerr << " primary     energy     runnumb  simnrsh    "
      << "   #hadrs       #muons      #gammas       #elecs    #nkgelecs  " 
      << "  obslev   theta    phi     h1km    h1gr" << endl; 
   cerr.precision(1);
   cerr.setf(ios::fixed,ios::floatfield);
   cerr << setw(8) << qdata[1];
   cerr.precision(4);
   cerr.setf(ios::scientific,ios::floatfield);
   cerr << setw(13) << qdata[2];
   cerr.precision(1);
   cerr.setf(ios::fixed,ios::floatfield);
   cerr << setw(10) << qdata[3];
   cerr << setw(9)  << qdata[4];
   cerr.precision(5);
   cerr.setf(ios::scientific,ios::floatfield);
   cerr << setw(13) << qdata[5];
   cerr << setw(13) << qdata[6];
   cerr << setw(13) << qdata[7];
   cerr << setw(13) << qdata[8];
   cerr << setw(13) << qdata[9];
   cerr.precision(1);
   cerr.setf(ios::fixed,ios::floatfield);
   cerr << setw(10) << qdata[11];
   cerr.precision(2);
   cerr << setw(8) << qdata[12];
   cerr << setw(8) << qdata[13];
   cerr << setw(8) << qdata[14];
   cerr << setw(8) << qdata[15];
   cerr << endl; 
   return;
}


// = = = = = = = = = = = =
// = = = = = = = = = = = = = = = printpartic = = = = = = = = = = = = = = =
//         print particle quantities of the first particle data subblock
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void printpartic(int m)
{
   cerr << setw(5) << (int) (pdata[m]/1000.) << "." << ((int)pdata[m])%10;
   cerr.precision(5);
   for( int i=1; i<ndatpar; i++ ) {
      cerr.setf(ios::scientific,ios::floatfield);
      cerr << setw(13) << pdata[i+m];
   }
   cerr.precision();
   cerr << endl;
   return;
}


// = = = = = = = = = = = =
// = = = = = = = = = = = = = = = = thickgr = = = = = = = = = = = = = = = =
//       thickgr (gramms/cm**2) of atmosphere depending on height (cm)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double thickgr(double h)
{
   double thickn;
   if (h < 4.e5)
      { thickn = AATM[1] + BATM[1] * exp(-h/CATM[1]); }
   else if (h < 1.e6)
      { thickn = AATM[2] + BATM[2] * exp(-h/CATM[2]); }
   else if (h < 4.e6)
      { thickn = AATM[3] + BATM[3] * exp(-h/CATM[3]); }
   else if (h < 1.e7)
      { thickn = AATM[4] + BATM[4] * exp(-h/CATM[4]); }
   else
      { thickn = AATM[5] - h*CATM[5]; }
   return thickn;
}

