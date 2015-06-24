/* ############################################################

     FILE:     cooling.h

     PURPOSE:  contains shared definitions with scope
               limited to the cooling module ONLY

   ############################################################ */

/* ##############################################################

                   P R O T O T Y P I N G

   ############################################################## */

void   CompEquil(double n, double T, double *v);
double GetMaxRate (double *, double *, double);
double MeanMolecularWeight  (double *);
void   Radiat (double *, double *);
void   NormalizeIons (double *);

/* ############################################################

         Fractions and Atomic Wts.

 ############################################################## */

/* Proto-Solar Mass Fractions for Hydrogen and Helium  
 * (Lodders, ApJ 591, 2003 )                                       */

#define H_MASS_FRAC       0.7110
#define He_MASS_FRAC      0.2741

#define NIONS    3
#define X_HI     (NFLX)
#define X_H2     (NFLX+1)
#define X_HII    (NFLX+2)
