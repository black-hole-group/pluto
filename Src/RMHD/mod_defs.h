/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Module header file for relativistic MHD (RMHD).

  Set label, indexes and basic prototyping for the relativistic 
  MHD module.


  \authors A. Mignone (mignone@ph.unito.it)\n
           C. Zanni   (zanni@oato.inaf.it)\n
  \date    Oct 5, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#ifndef RESISTIVE_RMHD
 #define RESISTIVE_RMHD   NO
#endif

/* Set variable name labels.
  We make extra vector components, when not needed, point 
  to the last element (255) of the array stored by startup.c.  */

enum {

 #if COMPONENTS == 1

  RHO, MX1, BX1, ENG, PRS = ENG,
  #if MHD_FORMULATION == DIV_CLEANING
   PSI_GLM,
  #endif
  MX2 = 255, BX2 = 255, MX3 = 255, BX3 = 255,

 #elif COMPONENTS == 2

  RHO, MX1, MX2, BX1, BX2, ENG, PRS = ENG,
  #if MHD_FORMULATION == DIV_CLEANING
   PSI_GLM,
  #endif
  MX3 = 255, BX3 = 255,

 #elif COMPONENTS == 3

  RHO, MX1, MX2, MX3, BX1, BX2, BX3, ENG, PRS = ENG,
  #if MHD_FORMULATION == DIV_CLEANING
   PSI_GLM,
  #endif

 #endif

 VX1 = MX1, VX2 = MX2, VX3 = MX3, 

/* -- backward compatibility -- */

 DN = RHO, PR = PRS, EN = ENG,
 VX = VX1, VY = VX2, VZ = VX3,
 MX = MX1, MY = MX2, MZ = MX3,
 BX = BX1, BY = BX2, BZ = BX3

};


/* **************************************
     add the PSI_GLM label if necessary 
   ************************************** */
/*
#if MHD_FORMULATION == DIV_CLEANING 
 #define PSI_GLM (2 + 2*COMPONENTS)
 #define NFLX (2 + 2*COMPONENTS + 1)
#else
 #define NFLX (2 + 2*COMPONENTS)
#endif
*/

#define NFLX (2 + 2*COMPONENTS + (MHD_FORMULATION==DIV_CLEANING))

/* *********************************************************
    Label the different waves in increasing order 
    following the number of vector components.

    IMPORTANT: the KPSI_GLMM & KPSI_GLMP modes are 
               present only in the MHD-GLM formulation.
               We keep them at the END of the enumeration
               so we can skip them in unnecessary loops.
               Please do NOT change them !
   ********************************************************* */

enum KWAVES {
 KFASTM, KFASTP, KENTRP

 #if MHD_FORMULATION != DIV_CLEANING
  , KDIVB
 #endif

 #if COMPONENTS >= 2
  , KSLOWM, KSLOWP
  #if COMPONENTS == 3
   , KALFVM, KALFVP
  #endif
 #endif

 #if MHD_FORMULATION == DIV_CLEANING  
  , KPSI_GLMM, KPSI_GLMP 
 #endif
};


#define SUBTRACT_DENSITY   YES  /**< By turning SUBTRACT_DENSITY to YES, 
                                     we let PLUTO evolve the total energy 
                                     minus the mass density contribution. */
                                     
/* ********************************************************************* */
/*! The Map_param structure is used to pass input/output arguments 
    during the conversion from conservative to primitive variables 
    operated by the ConsToPrim() function in the relativistic modules
    (RHD and RMHD).
    The output parameter, rho, W, lor and p, must be set at the end
    of every root-finder routine (EnergySolve(), EntropySolve() and
    PressureFix()).
    Additionally, some of the input parameters must be re-computed in
    EntropySolve() and PressureFix().
   ********************************************************************* */
typedef struct MAP_PARAM{
 double D;       /**< Lab density       (input). */
 double sigma_c; /**< Conserved entropy (input). */
 double E;       /**< Total energy      (input). */
 double m2;      /**< Square of total momentum (input). */
 double S;       /**< m<dot>B                  (input). */
 double S2;      /**< Square of S              (input). */
 double B2;      /**< Square of magnetic field (input). */

 double rho;     /**< proper density     (output)  */
 double W;       /**< D*h*lor            (output). */
 double lor;     /**< Lorentz factor     (output). */
 double prs;     /**< Thermal pressure   (output). */
} Map_param;

/* ******************************************************
     Vector potential: these labels are and MUST only
     be used in the STARTUP / INIT  functions;
     they're convenient in obtaining a discretization 
     that preserve divB since the beginning.   
   ****************************************************** */

#define   AX1  (NVAR + 1)
#define   AX2  (NVAR + 2)
#define   AX3  (NVAR + 3)

#define AX  AX1  /* backward compatibility */
#define AY  AX2
#define AZ  AX3

/* *************************************************
     Now define more convenient and user-friendly 
     pointer labels for geometry setting      
   ************************************************* */

#if GEOMETRY == CYLINDRICAL
                               
 #define iVR    VX1
 #define iVZ    VX2
 #define iVPHI  VX3
                               
 #define iMR    MX1
 #define iMZ    MX2
 #define iMPHI  MX3
                               
 #define iBR    BX1
 #define iBZ    BX2
 #define iBPHI  BX3
                               
#endif

#if GEOMETRY == POLAR
                               
 #define iVR    VX1
 #define iVPHI  VX2
 #define iVZ    VX3
                               
 #define iMR    MX1
 #define iMPHI  MX2
 #define iMZ    MX3
                               
 #define iBR    BX1
 #define iBPHI  BX2
 #define iBZ    BX3
                               
#endif

#if GEOMETRY == SPHERICAL
                    
 #define iVR     VX1
 #define iVTH    VX2
 #define iVPHI   VX3
                             
 #define iMR    MX1
 #define iMTH   MX2
 #define iMPHI  MX3
                               
 #define iBR    BX1
 #define iBTH   BX2
 #define iBPHI  BX3
 
#endif

/*      Prototyping goes here          */

int  ConsToPrim   (double **, double **, int, int, unsigned char *);
void PRIM_EIGENVECTORS (double *, double, double, double *, double **, double **);

/*
int  EnergySolve  (double *, double *);
int  EntropySolve (double *, double *);
int  PressureFix  (double *, double *);
*/

int  EntropySolve (Map_param *);
int  EnergySolve  (Map_param *);
int  PressureFix  (Map_param *);
 
void Flux      (double **, double **, double *, double **, double *, int, int);
void HLL_Speed (double **, double **, double *, double *, double *, double *,
                double *, double *, int, int);
int MaxSignalSpeed (double **, double *, double *, double *, double *, int, int);

void PrimToCons   (double **, double **, int, int);
void VelocityLimiter (double *, double *, double *);

int  Magnetosonic (double *vp, double cs2, double h, double *lambda);
int  QuarticSolve (double, double, double, double, double *);
int  CubicSolve   (double, double, double, double *);

Riemann_Solver LF_Solver, HLL_Solver, HLLC_Solver, HLLD_Solver, 
               MUSTA_Solver, RusanovDW_Solver;

#if MHD_FORMULATION == EIGHT_WAVES

 void POWELL_DIVB_SOURCE(const State_1D *, int, int, Grid *);
 void HLL_DIVB_SOURCE (const State_1D *, double **, int, int, Grid *);

#elif MHD_FORMULATION == DIV_CLEANING

 #include "MHD/GLM/glm.h"

#elif MHD_FORMULATION == CONSTRAINED_TRANSPORT

  #include "MHD/CT/ct.h"

#endif
