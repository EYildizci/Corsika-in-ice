//************************************************************************
// Original code by Juergen Oeschlaeger to write DAT file using C code
// 14/11/2012
// Bug fixes by Konrad Bernloehr
// 26/02/2014
//************************************************************************

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if HAVE_CONFIG_H
#include "config.h"
#endif

FILE *fmpatap;
#if __CERENKOV__
FILE *fmcetap;
#endif

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __REDHAT__
#if __THIN__
 void fwritempatap__( int *maxbuf, int *nsubbl, float outvect[6552] ) {
#else
 void fwritempatap__( int *maxbuf, int *nsubbl, float outvect[5733] ) {
#endif
#else
#if __THIN__
 void fwritempatap_( int *maxbuf, int *nsubbl, float outvect[6552] ) {
#else
 void fwritempatap_( int *maxbuf, int *nsubbl, float outvect[5733] ) {
#endif
#endif

  //  printf("input %d, %d",*maxbuf,*nsubbl);

   union fltint {
         float funi;
         int iuni;
   } intflt4;

   intflt4.iuni = *maxbuf * *nsubbl * 4; // record length 22932 wihtout thin or 26208 with thin.

   //   printf("outvect: %d, %f",intflt4.iuni,outvect[0]);

   if ( fmpatap != NULL ) /* Only try to write if we have an output file */
   {
      fwrite(&intflt4.funi, sizeof(float), 1, fmpatap);
      fwrite(outvect, intflt4.iuni, 1, fmpatap);
      fwrite(&intflt4.funi, sizeof(float), 1, fmpatap);
      /* Does it make sense to continue if there is an error (e.g. disk full) ? */
      /* if ( ferror(fmaptap) ) ... */
   }

}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __REDHAT__
 void fopenmpatap__( const char *corsika_name, int *ibl ) {
#else
 void fopenmpatap_( const char *corsika_name, int *ibl ) {
#endif

   unsigned long length= (unsigned long) (*ibl);
   const char *s = strchr(corsika_name,' ');
   if ( s != NULL && (s-corsika_name) < length )
      length = (unsigned long) (s-corsika_name);
   char *file_name=malloc(length+1);

   strncpy(file_name,corsika_name,length);
   file_name[length]='\0';
   if ( length != (unsigned long) (*ibl) )
   {
      fprintf(stderr,"\n\nFile name in fopenmpatap() function reported to have length %d\n", *ibl);
      fprintf(stderr,"but actually is only of length %lu: '%s'\n\n", length, file_name);
   }

   //     printf("file name: %s end\n",file_name) ;

   /* If the null device was given, it is more efficient to not open anything */
   fmpatap = NULL;
   if ( strcmp(file_name,"/dev/null") != 0 )
   {
      fmpatap = fopen( file_name, "w");
      if ( fmpatap == NULL )
         perror(file_name);
   }

}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __REDHAT__
void fclosempatap__() {
#else
void fclosempatap_() {
#endif

   if ( fmpatap != NULL )
      fclose(fmpatap); // close test version of particle data file.
   fmpatap = NULL;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __CERENKOV__
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __REDHAT__
#if __THIN__
 void fwritemcetap__( int *maxbf2, int *nsubbl, float outvect[6552] ) {
#else
 void fwritemcetap__( int *maxbf2, int *nsubbl, float outvect[5733] ) {
#endif
#else
#if __THIN__
 void fwritemcetap_( int *maxbf2, int *nsubbl, float outvect[6552] ) {
#else
 void fwritemcetap_( int *maxbf2, int *nsubbl, float outvect[5733] ) {
#endif
#endif

  //  printf("input %d, %d",*maxbf2,*nsubbl);

   union fltint {
         float funi;
         int iuni;
   } intflt4;

   intflt4.iuni = *maxbf2 * *nsubbl * 4; // record length 22932 wihtout thin or 26208 with thin.

   //   printf("outvect: %d, %f",intflt4.iuni,outvect[0]);

   if ( fmcetap != NULL ) /* Only try to write if we have an output file */
   {
      fwrite(&intflt4.funi, sizeof(float), 1, fmcetap);
      fwrite(outvect, intflt4.iuni, 1, fmcetap);
      fwrite(&intflt4.funi, sizeof(float), 1, fmcetap);
      /* Does it make sense to continue if there is an error (e.g. disk full) ? */
      /* if ( ferror(fmcetap) ) ... */
   }

}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __REDHAT__
 void fopenmcetap__( const char *corsika_name, int *ibl ) {
#else
 void fopenmcetap_( const char *corsika_name, int *ibl ) {
#endif

   unsigned long length= *ibl;
   const char *s = strchr(corsika_name,' ');
   if ( s != NULL && (s-corsika_name) < length )
      length = (unsigned long) (s-corsika_name);
   char *file_name=malloc(length+1);

   strncpy(file_name,corsika_name,length);
   file_name[length]='\0';
   if ( length != (unsigned long) (*ibl) )
   {
      fprintf(stderr,"\n\nFile name in fopenmcetap() function reported to have length %d\n", *ibl);
      fprintf(stderr,"but actually is only of length %lu: '%s'\n\n", length, file_name);
   }

   //     printf("file name: %s end\n",file_name) ;

   /* If the null device was given, it is more efficient to not open anything */
   fmcetap = NULL;
   if ( strcmp(file_name,"/dev/null") != 0 )
   {
      fmcetap = fopen( file_name, "w");
      if ( fmpatap == NULL )
         perror(file_name);
   }
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __REDHAT__
void fclosemcetap__() {
#else
void fclosemcetap_() {
#endif

   if ( fmcetap != NULL )
      fclose(fmcetap); // close test version of particle data file.
   fmcetap = NULL;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#endif
