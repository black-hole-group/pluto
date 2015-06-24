/* ############################################################
      
     FILE:     cooling.h

     PURPOSE:  contains common definitions for the 
               whole CODE

   ############################################################ */

#define NIONS  1
#define X_HI    NFLX   

double GetMaxRate (double *, double *, double);
double MeanMolecularWeight  (double *);
double H_MassFrac (void);
double CompEquil (double, double, double *);
void Radiat (double *, double *);










